#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

#define solve_greedy_l0 deletion_solve_greedy_l0
#include "greedy_l0.cpp"
#undef solve_greedy_l0

using namespace Rcpp;

static std::vector<double> signed_deletion_rank_scores(
    const std::vector<double>& beta,
    const std::vector<double>& evidence) {
  std::vector<double> scores(beta.size(), 0.0);
  std::vector<int> selected;
  selected.reserve(beta.size());
  for (int i = 0; i < static_cast<int>(beta.size()); ++i) {
    if (R_finite(beta[i]) && beta[i] != 0.0 &&
        R_finite(evidence[i]) && evidence[i] >= 0.0) {
      selected.push_back(i);
    }
  }
  std::sort(selected.begin(), selected.end(), [&](int left, int right) {
    if (evidence[left] != evidence[right]) {
      return evidence[left] > evidence[right];
    }
    return left < right;
  });
  int groups = 0;
  for (int i = 0; i < static_cast<int>(selected.size()); ) {
    const double anchor = evidence[selected[i]];
    int j = i + 1;
    while (j < static_cast<int>(selected.size()) &&
           evidence[selected[j]] == anchor) {
      ++j;
    }
    ++groups;
    i = j;
  }
  int group_index = 0;
  for (int i = 0; i < static_cast<int>(selected.size()); ) {
    const double anchor = evidence[selected[i]];
    int j = i + 1;
    while (j < static_cast<int>(selected.size()) &&
           evidence[selected[j]] == anchor) {
      ++j;
    }
    const double magnitude = 1.0 -
      (static_cast<double>(group_index) + 0.5) /
      static_cast<double>(groups);
    for (int k = i; k < j; ++k) {
      const int index = selected[k];
      scores[index] = beta[index] < 0.0 ? -magnitude : magnitude;
    }
    ++group_index;
    i = j;
  }
  return scores;
}

// [[Rcpp::export]]
DataFrame infer_network(
    NumericMatrix expression,
    CharacterVector gene_names,
    NumericMatrix pseudotime,
    List params
) {
  const int p = expression.nrow();
  const int n = expression.ncol();
  if (p == 0 || n < 2) {
    return DataFrame::create(
      _["regulator"] = CharacterVector(),
      _["target"] = CharacterVector(),
      _["standardized_beta"] = NumericVector(),
      _["deletion_delta_bic"] = NumericVector(),
      _["weight"] = NumericVector(),
      _["stringsAsFactors"] = false
    );
  }

  int max_support = static_cast<int>(list_numeric_or_default(params, "max_support_size", 0.0));
  if (max_support <= 0) {
    max_support = p - 1;
  } else {
    max_support = std::min(max_support, p - 1);
  }
  max_support = std::max(0, max_support);
  const double min_improvement = std::max(0.0, list_numeric_or_default(params, "min_improvement", 1e-10));
  const double lag_fraction = list_numeric_or_default(params, "pseudotime_lag_fraction", 0.05);
  const int requested_lag_steps = list_int_or_default(params, "pseudotime_lag_steps", 0);
  const int cores = std::max(1, list_int_or_default(params, "cores", 1));

  if (pseudotime.ncol() > 0 && pseudotime.nrow() != n) {
    stop("`pseudotime` must contain one row per cell");
  }
  std::vector<LagBranch> lag_branches;
  lag_branches.reserve(pseudotime.ncol());
  for (int branch = 0; branch < pseudotime.ncol(); ++branch) {
    NumericVector branch_time = pseudotime(_, branch);
    std::vector<std::vector<int> > branch_states = pseudotime_states(branch_time);
    const int branch_steps = effective_lag_steps(
      static_cast<int>(branch_states.size()),
      lag_fraction,
      requested_lag_steps
    );
    if (branch_steps <= 0) {
      continue;
    }
    if (static_cast<int>(branch_states.size()) - branch_steps < 2) {
      continue;
    }
    LagBranch lag_branch;
    lag_branch.states.swap(branch_states);
    lag_branch.steps = branch_steps;
    lag_branches.push_back(lag_branch);
  }
  std::sort(
    lag_branches.begin(),
    lag_branches.end(),
    [](const LagBranch& left, const LagBranch& right) {
      if (left.states != right.states) {
        return std::lexicographical_compare(
          left.states.begin(), left.states.end(),
          right.states.begin(), right.states.end()
        );
      }
      return left.steps < right.steps;
    }
  );
  lag_branches.erase(
    std::unique(
      lag_branches.begin(),
      lag_branches.end(),
      [](const LagBranch& left, const LagBranch& right) {
        return left.steps == right.steps && left.states == right.states;
      }
    ),
    lag_branches.end()
  );
  std::map<std::vector<int>, int> transition_multiplicity;
  for (const LagBranch& branch : lag_branches) {
    const int lag_n = static_cast<int>(branch.states.size()) - branch.steps;
    for (int row = 0; row < lag_n; ++row) {
      ++transition_multiplicity[lag_transition_key(branch, row)];
    }
  }
  for (LagBranch& branch : lag_branches) {
    const int lag_n = static_cast<int>(branch.states.size()) - branch.steps;
    branch.transition_weight.assign(lag_n, 1.0);
    for (int row = 0; row < lag_n; ++row) {
      const int multiplicity = transition_multiplicity[lag_transition_key(branch, row)];
      branch.transition_weight[row] = 1.0 / static_cast<double>(multiplicity);
    }
  }
  const bool use_lag = !lag_branches.empty();

  const std::vector<char> regulator_mask = gene_mask_from_param(params, "regulators", gene_names, p);
  const std::vector<char> target_mask = gene_mask_from_param(params, "targets", gene_names, p);

  std::vector<std::string> regulators;
  std::vector<std::string> targets;
  std::vector<double> standardized_betas;
  std::vector<double> deletion_delta_bics;
  std::vector<double> weights;
  regulators.reserve(static_cast<size_t>(p) * std::max(max_support, 0));
  targets.reserve(static_cast<size_t>(p) * std::max(max_support, 0));
  standardized_betas.reserve(static_cast<size_t>(p) * std::max(max_support, 0));
  deletion_delta_bics.reserve(static_cast<size_t>(p) * std::max(max_support, 0));
  weights.reserve(static_cast<size_t>(p) * std::max(max_support, 0));

  std::vector<double> core_raw(static_cast<size_t>(p) * p, 0.0);
  std::vector<double> core_evidence(static_cast<size_t>(p) * p, 0.0);

  std::vector<double> static_gram;
  std::vector<char> static_valid(p, 0);
  if (!use_lag) {
    std::vector<double> standardized(
      static_cast<size_t>(p) * static_cast<size_t>(n),
      0.0
    );
    #pragma omp parallel for num_threads(cores)
    for (int gene = 0; gene < p; ++gene) {
      double sum = 0.0;
      int finite_count = 0;
      for (int cell = 0; cell < n; ++cell) {
        const double value = expression(gene, cell);
        if (R_finite(value)) {
          sum += value;
          ++finite_count;
        }
      }
      if (finite_count < 2) {
        continue;
      }
      const double mean = sum / static_cast<double>(finite_count);
      double ss = 0.0;
      for (int cell = 0; cell < n; ++cell) {
        const double value = expression(gene, cell);
        const double centered_value = R_finite(value) ? value - mean : 0.0;
        standardized[static_cast<size_t>(gene) * n + cell] = centered_value;
        ss += centered_value * centered_value;
      }
      if (ss <= 0.0) {
        continue;
      }
      const double scale = std::sqrt(ss / static_cast<double>(finite_count - 1));
      static_valid[gene] = 1;
      for (int cell = 0; cell < n; ++cell) {
        standardized[static_cast<size_t>(gene) * n + cell] /= scale;
      }
    }

    symmetric_crossproduct(standardized, n, p, static_gram);
  }

  const int lag_n_total = static_cast<int>(transition_multiplicity.size());
  std::vector<double> lag_predictor_gram;
  std::vector<double> lag_cross;
  std::vector<double> lag_y_ss_vec;
  std::vector<char> lag_predictor_valid(p, 0);
  std::vector<char> lag_response_valid(p, 0);
  if (use_lag) {
    lag_predictor_gram.assign(static_cast<size_t>(p) * p, 0.0);
    lag_cross.assign(static_cast<size_t>(p) * p, 0.0);
    lag_y_ss_vec.assign(p, 0.0);

    for (const LagBranch& branch : lag_branches) {
      const int lag_n = static_cast<int>(branch.states.size()) - branch.steps;

      std::vector<double> lag_x(static_cast<size_t>(p) * lag_n, 0.0);
      std::vector<double> lag_y(static_cast<size_t>(p) * lag_n, 0.0);
      #pragma omp parallel for num_threads(cores)
      for (int gene = 0; gene < p; ++gene) {
        double x_sum = 0.0;
        double y_sum = 0.0;
        double x_weight_sum = 0.0;
        double y_weight_sum = 0.0;
        std::vector<char> x_row_valid(lag_n, 0);
        std::vector<char> y_row_valid(lag_n, 0);
        for (int row = 0; row < lag_n; ++row) {
          const std::vector<int>& x_cells = branch.states[row];
          const std::vector<int>& y_cells = branch.states[row + branch.steps];
          double x_state_sum = 0.0;
          double y_state_sum = 0.0;
          int x_state_n = 0;
          int y_state_n = 0;
          for (int cell : x_cells) {
            const double value = expression(gene, cell);
            if (R_finite(value)) {
              x_state_sum += value;
              ++x_state_n;
            }
          }
          for (int cell : y_cells) {
            const double value = expression(gene, cell);
            if (R_finite(value)) {
              y_state_sum += value;
              ++y_state_n;
            }
          }
          if (x_state_n > 0) {
            const double x_value = x_state_sum / static_cast<double>(x_state_n);
            lag_x[static_cast<size_t>(gene) * lag_n + row] = x_value;
            x_row_valid[row] = 1;
            x_sum += branch.transition_weight[row] * x_value;
            x_weight_sum += branch.transition_weight[row];
          }
          if (y_state_n > 0) {
            const double y_value = y_state_sum / static_cast<double>(y_state_n);
            lag_y[static_cast<size_t>(gene) * lag_n + row] = y_value;
            y_row_valid[row] = 1;
            y_sum += branch.transition_weight[row] * y_value;
            y_weight_sum += branch.transition_weight[row];
          }
        }
        const double x_mean = x_weight_sum > 0.0
          ? x_sum / x_weight_sum
          : 0.0;
        const double y_mean = y_weight_sum > 0.0
          ? y_sum / y_weight_sum
          : 0.0;
        double x_ss = 0.0;
        double y_ss = 0.0;
        for (int row = 0; row < lag_n; ++row) {
          const double x_value = lag_x[static_cast<size_t>(gene) * lag_n + row];
          const double y_value = lag_y[static_cast<size_t>(gene) * lag_n + row];
          if (x_row_valid[row]) {
            const double centered_x = x_value - x_mean;
            lag_x[static_cast<size_t>(gene) * lag_n + row] = centered_x;
            x_ss += branch.transition_weight[row] * centered_x * centered_x;
          }
          if (y_row_valid[row]) {
            const double centered_y = y_value - y_mean;
            lag_y[static_cast<size_t>(gene) * lag_n + row] = centered_y;
            y_ss += branch.transition_weight[row] * centered_y * centered_y;
          }
        }
        if (x_weight_sum > 0.0 && x_ss > 0.0) {
          const double x_scale = std::sqrt(x_ss / x_weight_sum);
          lag_predictor_valid[gene] = 1;
          for (int row = 0; row < lag_n; ++row) {
            lag_x[static_cast<size_t>(gene) * lag_n + row] /= x_scale;
          }
        }
        if (y_weight_sum > 0.0 && y_ss > 0.0) {
          const double y_scale = std::sqrt(y_ss / y_weight_sum);
          lag_response_valid[gene] = 1;
          for (int row = 0; row < lag_n; ++row) {
            lag_y[static_cast<size_t>(gene) * lag_n + row] /= y_scale;
          }
        }
      }

      accumulate_weighted_crossproducts(
        lag_x,
        lag_y,
        branch.transition_weight,
        lag_n,
        p,
        lag_predictor_gram,
        lag_cross,
        lag_y_ss_vec
      );
    }
  }

  const std::vector<double>& active_gram = use_lag
    ? lag_predictor_gram
    : static_gram;
  const std::vector<double>& active_cross = use_lag
    ? lag_cross
    : static_gram;
  const std::vector<char>& active_predictor_valid = use_lag
    ? lag_predictor_valid
    : static_valid;
  const std::vector<char>& active_response_valid = use_lag
    ? lag_response_valid
    : static_valid;
  const int active_n = use_lag ? lag_n_total : n;
  const int active_groups = use_lag
    ? static_cast<int>(lag_branches.size())
    : 1;
  const int active_max_support = std::max(
    0,
    std::min(max_support, active_n - active_groups - 1)
  );

  #pragma omp parallel for num_threads(cores) schedule(dynamic)
  for (int target = 0; target < p; ++target) {
    if (!target_mask[target] || !active_response_valid[target]) {
      continue;
    }
    std::vector<double> xty(p, 0.0);
    for (int regulator = 0; regulator < p; ++regulator) {
      xty[regulator] = active_cross[static_cast<size_t>(regulator) * p + target];
    }
    const double y_ss = use_lag
      ? lag_y_ss_vec[target]
      : static_gram[static_cast<size_t>(target) * p + target];
    if (y_ss <= 0.0) {
      continue;
    }
    std::vector<char> allowed(p, 1);
    allowed[target] = 0;
    for (int regulator = 0; regulator < p; ++regulator) {
      if (!regulator_mask[regulator] || !active_predictor_valid[regulator]) {
        allowed[regulator] = 0;
      }
    }

    const TargetFit fit = fit_target_greedy_from_gram(
      active_gram,
      xty,
      allowed,
      target,
      p,
      active_n,
      active_max_support,
      min_improvement,
      y_ss
    );
    if (!fit.support.empty()) {
      std::vector<double> support_inverse;
      if (!invert_subset_gram(active_gram, fit.support, p, support_inverse)) {
        stop("Selected support Gram matrix is not identifiable.");
      }
      const int selected_count = static_cast<int>(fit.support.size());
      for (int selected = 0; selected < selected_count; ++selected) {
        const int regulator = fit.support[selected];
        const double beta = fit.coefficient[regulator];
        const double inverse_diagonal = support_inverse[
          static_cast<size_t>(selected) * selected_count + selected
        ];
        if (!R_finite(beta) || beta == 0.0 ||
            !R_finite(inverse_diagonal) || inverse_diagonal <= 0.0) {
          stop("Invalid selected coefficient or inverse Gram diagonal.");
        }
        const double removed_rss = fit.rss + beta * beta / inverse_diagonal;
        double delta_bic = subset_bic(
          removed_rss, selected_count - 1, active_n
        ) - fit.bic;
        const double evidence_tolerance =
          1e-8 * (1.0 + std::fabs(fit.bic));
        if (delta_bic < -evidence_tolerance) {
          stop("Selected support is not deletion-local-optimal.");
        }
        const size_t edge_index =
          static_cast<size_t>(regulator) * p + target;
        core_raw[edge_index] = beta;
        core_evidence[edge_index] = std::max(0.0, delta_bic);
      }
    }
  }

  const std::vector<double> network_scores =
    signed_deletion_rank_scores(core_raw, core_evidence);
  std::vector<NetworkEdge> edges;
  edges.reserve(static_cast<size_t>(p) * std::max(max_support, 0));
  for (int target = 0; target < p; ++target) {
    if (!target_mask[target]) {
      continue;
    }
    for (int regulator = 0; regulator < p; ++regulator) {
      const double weight = network_scores[static_cast<size_t>(regulator) * p + target];
      if (regulator != target && weight != 0.0) {
        edges.push_back({regulator, target, weight});
      }
    }
  }
  std::sort(edges.begin(), edges.end(), [](const NetworkEdge& left, const NetworkEdge& right) {
    const double left_abs = std::fabs(left.weight);
    const double right_abs = std::fabs(right.weight);
    if (left_abs != right_abs) {
      return left_abs > right_abs;
    }
    if (left.regulator != right.regulator) {
      return left.regulator < right.regulator;
    }
    return left.target < right.target;
  });
  for (const NetworkEdge& edge : edges) {
    regulators.push_back(as<std::string>(gene_names[edge.regulator]));
    targets.push_back(as<std::string>(gene_names[edge.target]));
    const size_t edge_index =
      static_cast<size_t>(edge.regulator) * p + edge.target;
    standardized_betas.push_back(core_raw[edge_index]);
    deletion_delta_bics.push_back(core_evidence[edge_index]);
    weights.push_back(edge.weight);
  }

  return DataFrame::create(
    _["regulator"] = regulators,
    _["target"] = targets,
    _["standardized_beta"] = standardized_betas,
    _["deletion_delta_bic"] = deletion_delta_bics,
    _["weight"] = weights,
    _["stringsAsFactors"] = false
  );
}

// [[Rcpp::export]]
List solve_greedy_l0_batch(
    NumericMatrix gram_matrix,
    NumericMatrix xty_matrix,
    NumericVector y_ss,
    List candidate_predictors,
    int n_obs,
    int max_support_size,
    double min_improvement
) {
  const int p = gram_matrix.nrow();
  const int targets = xty_matrix.ncol();
  if (p < 1 || gram_matrix.ncol() != p || xty_matrix.nrow() != p ||
      y_ss.size() != targets || candidate_predictors.size() != targets ||
      n_obs < 3) {
    stop("Invalid batch dimensions.");
  }

  std::vector<int> edge_predictor;
  std::vector<int> edge_target;
  std::vector<double> edge_beta;
  std::vector<double> edge_delta_bic;
  IntegerVector candidate_size(targets);
  IntegerVector support_size(targets);
  NumericVector rss(targets, NA_REAL);
  NumericVector bic(targets, NA_REAL);
  CharacterVector status(targets, "unprocessed");

  for (int target = 0; target < targets; ++target) {
    IntegerVector supplied = candidate_predictors[target];
    std::vector<int> candidate;
    candidate.reserve(supplied.size());
    std::vector<char> seen(p, 0);
    for (int index = 0; index < supplied.size(); ++index) {
      const int predictor = supplied[index] - 1;
      if (predictor < 0 || predictor >= p || seen[predictor]) {
        stop("Candidate predictor indices must be unique and one-based.");
      }
      seen[predictor] = 1;
      const double diagonal = gram_matrix(predictor, predictor);
      const double cross = xty_matrix(predictor, target);
      if (!R_finite(diagonal) || !R_finite(cross)) {
        stop("Batch Gram and cross-products must be finite.");
      }
      if (diagonal > 0.0) candidate.push_back(predictor);
    }
    candidate_size[target] = static_cast<int>(candidate.size());
    if (candidate.empty()) {
      status[target] = "no_variable_candidate";
      continue;
    }
    if (!R_finite(y_ss[target]) || y_ss[target] <= 0.0) {
      status[target] = "constant_response";
      continue;
    }

    const int k = static_cast<int>(candidate.size());
    std::vector<double> local_gram(static_cast<size_t>(k) * k, 0.0);
    std::vector<double> local_xty(k, 0.0);
    std::vector<char> allowed(k, 1);
    for (int left = 0; left < k; ++left) {
      local_xty[left] = xty_matrix(candidate[left], target);
      for (int right = 0; right < k; ++right) {
        const double value = gram_matrix(candidate[left], candidate[right]);
        if (!R_finite(value)) stop("Batch Gram matrix must be finite.");
        local_gram[static_cast<size_t>(left) * k + right] = value;
      }
    }
    const int requested_support = max_support_size <= 0
      ? k
      : std::min(max_support_size, k);
    const int support_cap = std::max(0, std::min(requested_support, n_obs - 2));
    const TargetFit fit = fit_target_greedy_from_gram(
      local_gram, local_xty, allowed, -1, k, n_obs, support_cap,
      std::max(0.0, min_improvement), y_ss[target]
    );
    rss[target] = fit.rss;
    bic[target] = fit.bic;
    support_size[target] = static_cast<int>(fit.support.size());
    status[target] = fit.support.empty() ? "empty_support" : "fit";
    if (fit.support.empty()) continue;

    std::vector<double> support_inverse;
    if (!invert_subset_gram(local_gram, fit.support, k, support_inverse)) {
      stop("Selected support Gram matrix is not identifiable.");
    }
    const int selected_count = static_cast<int>(fit.support.size());
    for (int selected = 0; selected < selected_count; ++selected) {
      const int local_predictor = fit.support[selected];
      const double beta = fit.coefficient[local_predictor];
      const double inverse_diagonal = support_inverse[
        static_cast<size_t>(selected) * selected_count + selected
      ];
      if (!R_finite(beta) || beta == 0.0 ||
          !R_finite(inverse_diagonal) || inverse_diagonal <= 0.0) {
        stop("Invalid selected coefficient or inverse Gram diagonal.");
      }
      const double removed_rss = fit.rss + beta * beta / inverse_diagonal;
      double delta_bic = subset_bic(
        removed_rss, selected_count - 1, n_obs
      ) - fit.bic;
      const double tolerance = 1e-8 * (1.0 + std::fabs(fit.bic));
      if (delta_bic < -tolerance) {
        stop("Selected support is not deletion-local-optimal.");
      }
      edge_predictor.push_back(candidate[local_predictor] + 1);
      edge_target.push_back(target + 1);
      edge_beta.push_back(beta);
      edge_delta_bic.push_back(std::max(0.0, delta_bic));
    }
  }

  return List::create(
    _["predictor_index"] = wrap(edge_predictor),
    _["target_index"] = wrap(edge_target),
    _["standardized_beta"] = wrap(edge_beta),
    _["deletion_delta_bic"] = wrap(edge_delta_bic),
    _["candidate_size"] = candidate_size,
    _["support_size"] = support_size,
    _["rss"] = rss,
    _["bic"] = bic,
    _["status"] = status,
    _["n_obs"] = n_obs
  );
}
