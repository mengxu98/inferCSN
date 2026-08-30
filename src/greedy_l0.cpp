#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

struct NetworkEdge {
  int regulator;
  int target;
  double weight;
};

struct TargetFit {
  std::vector<double> coefficient;
  double r2;
  double rss;
  double bic;
  std::vector<int> support;

  TargetFit() : r2(0.0), rss(0.0), bic(std::numeric_limits<double>::infinity()) {}
};

struct LagBranch {
  std::vector<std::vector<int> > states;
  std::vector<double> transition_weight;
  int steps;
};

[[maybe_unused]] static std::vector<int> lag_transition_key(
    const LagBranch& branch,
    int row
) {
  const std::vector<int>& source = branch.states[row];
  const std::vector<int>& destination = branch.states[row + branch.steps];
  std::vector<int> key;
  key.reserve(source.size() + destination.size() + 1);
  key.insert(key.end(), source.begin(), source.end());
  key.push_back(-1);
  key.insert(key.end(), destination.begin(), destination.end());
  return key;
}

static double list_numeric_or_default(const List& params, const char* name, double value) {
  if (!params.containsElementNamed(name)) {
    return value;
  }
  SEXP item = params[name];
  if (Rf_isNull(item)) {
    return value;
  }
  NumericVector x(item);
  if (x.size() == 0 || !R_finite(x[0])) {
    return value;
  }
  return x[0];
}

[[maybe_unused]] static int list_int_or_default(const List& params, const char* name, int value) {
  return static_cast<int>(list_numeric_or_default(params, name, static_cast<double>(value)));
}

[[maybe_unused]] static std::vector<char> gene_mask_from_param(
    const List& params,
    const char* name,
    const CharacterVector& gene_names,
    int p
) {
  std::vector<char> mask(p, 1);
  if (!params.containsElementNamed(name)) {
    return mask;
  }
  SEXP item = params[name];
  if (Rf_isNull(item)) {
    return mask;
  }
  CharacterVector requested(item);
  if (requested.size() == 0) {
    return mask;
  }

  std::map<std::string, int> index;
  for (int i = 0; i < p; ++i) {
    if (!CharacterVector::is_na(gene_names[i])) {
      index[as<std::string>(gene_names[i])] = i;
    }
  }

  std::vector<char> selected(p, 0);
  int selected_count = 0;
  for (int i = 0; i < requested.size(); ++i) {
    if (CharacterVector::is_na(requested[i])) {
      continue;
    }
    const std::string gene = as<std::string>(requested[i]);
    std::map<std::string, int>::const_iterator found = index.find(gene);
    if (found == index.end()) {
      continue;
    }
    if (!selected[found->second]) {
      selected[found->second] = 1;
      ++selected_count;
    }
  }
  return selected_count > 0 ? selected : mask;
}

[[maybe_unused]] static int effective_lag_steps(int state_count, double lag_fraction, int lag_steps) {
  if (state_count < 2) {
    return 0;
  }
  if (lag_steps > 0) {
    return std::min(lag_steps, state_count - 1);
  }
  int fraction_steps = 0;
  if (lag_fraction > 0.0 && R_finite(lag_fraction)) {
    fraction_steps = static_cast<int>(
      std::ceil(static_cast<double>(state_count) * lag_fraction)
    );
  }
  return std::min(std::max(1, fraction_steps), state_count - 1);
}

[[maybe_unused]] static std::vector<std::vector<int> > pseudotime_states(
    const NumericVector& pseudotime
) {
  std::vector<std::pair<double, int> > ordered;
  ordered.reserve(pseudotime.size());
  for (int i = 0; i < pseudotime.size(); ++i) {
    if (R_finite(pseudotime[i])) {
      ordered.push_back(std::make_pair(pseudotime[i], i));
    }
  }
  std::sort(
    ordered.begin(),
    ordered.end(),
    [](const std::pair<double, int>& left,
       const std::pair<double, int>& right) {
      if (left.first == right.first) {
        return left.second < right.second;
      }
      return left.first < right.first;
    }
  );

  std::vector<std::vector<int> > states;
  int i = 0;
  while (i < static_cast<int>(ordered.size())) {
    int j = i + 1;
    while (j < static_cast<int>(ordered.size()) &&
           ordered[j].first == ordered[i].first) {
      ++j;
    }
    std::vector<int> cells;
    cells.reserve(j - i);
    for (int pos = i; pos < j; ++pos) {
      cells.push_back(ordered[pos].second);
    }
    states.push_back(cells);
    i = j;
  }
  return states;
}

static bool solve_linear_system(
    std::vector<double> a,
    std::vector<double> b,
    int k,
    std::vector<double>& solution
) {
  solution.assign(k, 0.0);
  if (k <= 0) {
    return true;
  }
  double matrix_scale = 0.0;
  for (int i = 0; i < k; ++i) {
    matrix_scale = std::max(
      matrix_scale,
      std::fabs(a[static_cast<size_t>(i) * k + i])
    );
  }
  const double pivot_tolerance = 1e-10 * std::max(1.0, matrix_scale);
  for (int col = 0; col < k; ++col) {
    int pivot = col;
    double pivot_abs = std::fabs(a[static_cast<size_t>(col) * k + col]);
    for (int row = col + 1; row < k; ++row) {
      const double candidate = std::fabs(a[static_cast<size_t>(row) * k + col]);
      if (candidate > pivot_abs) {
        pivot = row;
        pivot_abs = candidate;
      }
    }
    if (pivot_abs <= pivot_tolerance) {
      return false;
    }
    if (pivot != col) {
      for (int j = col; j < k; ++j) {
        std::swap(a[static_cast<size_t>(col) * k + j], a[static_cast<size_t>(pivot) * k + j]);
      }
      std::swap(b[col], b[pivot]);
    }
    const double diag = a[static_cast<size_t>(col) * k + col];
    for (int row = col + 1; row < k; ++row) {
      const double factor = a[static_cast<size_t>(row) * k + col] / diag;
      if (factor == 0.0) {
        continue;
      }
      a[static_cast<size_t>(row) * k + col] = 0.0;
      for (int j = col + 1; j < k; ++j) {
        a[static_cast<size_t>(row) * k + j] -= factor * a[static_cast<size_t>(col) * k + j];
      }
      b[row] -= factor * b[col];
    }
  }
  for (int row = k - 1; row >= 0; --row) {
    double rhs = b[row];
    for (int col = row + 1; col < k; ++col) {
      rhs -= a[static_cast<size_t>(row) * k + col] * solution[col];
    }
    const double diag = a[static_cast<size_t>(row) * k + row];
    if (std::fabs(diag) <= pivot_tolerance) {
      return false;
    }
    solution[row] = rhs / diag;
  }
  return true;
}

static bool fit_subset_coefficients(
    const std::vector<double>& gram,
    const std::vector<double>& xty,
    const std::vector<int>& subset,
    int p,
    std::vector<double>& beta
) {
  const int k = static_cast<int>(subset.size());
  std::vector<double> lhs(static_cast<size_t>(k) * k, 0.0);
  std::vector<double> rhs(k, 0.0);
  for (int i = 0; i < k; ++i) {
    rhs[i] = xty[subset[i]];
    for (int j = 0; j < k; ++j) {
      lhs[static_cast<size_t>(i) * k + j] = gram[static_cast<size_t>(subset[i]) * p + subset[j]];
    }
  }
  return solve_linear_system(lhs, rhs, k, beta);
}

static bool invert_subset_gram(
    const std::vector<double>& gram,
    const std::vector<int>& subset,
    int p,
    std::vector<double>& inverse
) {
  const int k = static_cast<int>(subset.size());
  inverse.assign(static_cast<size_t>(k) * k, 0.0);
  if (k == 0) {
    return true;
  }

  std::vector<double> work(static_cast<size_t>(k) * k, 0.0);
  double matrix_scale = 0.0;
  for (int row = 0; row < k; ++row) {
    for (int column = 0; column < k; ++column) {
      work[static_cast<size_t>(row) * k + column] =
        gram[static_cast<size_t>(subset[row]) * p + subset[column]];
    }
    inverse[static_cast<size_t>(row) * k + row] = 1.0;
    matrix_scale = std::max(
      matrix_scale,
      std::fabs(work[static_cast<size_t>(row) * k + row])
    );
  }

  const double pivot_tolerance = 1e-10 * std::max(1.0, matrix_scale);
  for (int column = 0; column < k; ++column) {
    int pivot = column;
    double pivot_abs = std::fabs(
      work[static_cast<size_t>(column) * k + column]
    );
    for (int row = column + 1; row < k; ++row) {
      const double candidate = std::fabs(
        work[static_cast<size_t>(row) * k + column]
      );
      if (candidate > pivot_abs) {
        pivot = row;
        pivot_abs = candidate;
      }
    }
    if (pivot_abs <= pivot_tolerance) {
      return false;
    }
    if (pivot != column) {
      for (int item = 0; item < k; ++item) {
        std::swap(
          work[static_cast<size_t>(column) * k + item],
          work[static_cast<size_t>(pivot) * k + item]
        );
        std::swap(
          inverse[static_cast<size_t>(column) * k + item],
          inverse[static_cast<size_t>(pivot) * k + item]
        );
      }
    }

    const double diagonal = work[static_cast<size_t>(column) * k + column];
    for (int item = 0; item < k; ++item) {
      work[static_cast<size_t>(column) * k + item] /= diagonal;
      inverse[static_cast<size_t>(column) * k + item] /= diagonal;
    }
    for (int row = 0; row < k; ++row) {
      if (row == column) {
        continue;
      }
      const double factor = work[static_cast<size_t>(row) * k + column];
      if (factor == 0.0) {
        continue;
      }
      for (int item = 0; item < k; ++item) {
        work[static_cast<size_t>(row) * k + item] -=
          factor * work[static_cast<size_t>(column) * k + item];
        inverse[static_cast<size_t>(row) * k + item] -=
          factor * inverse[static_cast<size_t>(column) * k + item];
      }
    }
  }
  return true;
}

static void residual_candidate_statistics(
    const std::vector<double>& gram,
    const std::vector<double>& xty,
    const std::vector<int>& subset,
    const std::vector<double>& beta,
    const std::vector<double>& inverse,
    const std::vector<char>& allowed,
    const std::vector<char>& used,
    int p,
    std::vector<double>& projection,
    std::vector<double>& residual_variance,
    std::vector<double>& residual_cross,
    std::vector<char>& valid
) {
  const int k = static_cast<int>(subset.size());
  projection.assign(static_cast<size_t>(p) * k, 0.0);
  residual_variance.assign(p, 0.0);
  residual_cross.assign(p, 0.0);
  valid.assign(p, 0);

  if (k > 0) {
    std::vector<double> inverse_column_major(static_cast<size_t>(k) * k, 0.0);
    std::vector<double> selected_cross(static_cast<size_t>(k) * p, 0.0);
    for (int column = 0; column < k; ++column) {
      for (int row = 0; row < k; ++row) {
        inverse_column_major[static_cast<size_t>(column) * k + row] =
          inverse[static_cast<size_t>(row) * k + column];
      }
    }
    for (int regulator = 0; regulator < p; ++regulator) {
      for (int row = 0; row < k; ++row) {
        selected_cross[static_cast<size_t>(regulator) * k + row] =
          gram[static_cast<size_t>(subset[row]) * p + regulator];
      }
    }
    const char no_transpose = 'N';
    const double one = 1.0;
    const double zero = 0.0;
    F77_CALL(dgemm)(
      &no_transpose,
      &no_transpose,
      &k,
      &p,
      &k,
      &one,
      inverse_column_major.data(),
      &k,
      selected_cross.data(),
      &k,
      &zero,
      projection.data(),
      &k FCONE FCONE
    );
  }

  for (int regulator = 0; regulator < p; ++regulator) {
    if (!allowed[regulator] || used[regulator]) {
      continue;
    }
    double this_cross = xty[regulator];
    double this_variance =
      gram[static_cast<size_t>(regulator) * p + regulator];
    for (int index = 0; index < k; ++index) {
      const double gram_value =
        gram[static_cast<size_t>(regulator) * p + subset[index]];
      this_cross -= gram_value * beta[index];
      this_variance -= gram_value *
        projection[static_cast<size_t>(regulator) * k + index];
    }
    const double variance_tolerance = 1e-10 * std::max(
      1.0,
      std::fabs(gram[static_cast<size_t>(regulator) * p + regulator])
    );
    residual_variance[regulator] = this_variance;
    residual_cross[regulator] = this_cross;
    if (this_variance <= variance_tolerance || !R_finite(this_variance)) {
      continue;
    }
    valid[regulator] = 1;
  }
}

[[maybe_unused]] static void symmetric_crossproduct(
    const std::vector<double>& matrix,
    int rows,
    int columns,
    std::vector<double>& result
) {
  result.assign(static_cast<size_t>(columns) * columns, 0.0);
  if (rows <= 0 || columns <= 0) {
    return;
  }
  const char transpose = 'T';
  const char no_transpose = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dgemm)(
    &transpose,
    &no_transpose,
    &columns,
    &columns,
    &rows,
    &one,
    matrix.data(),
    &rows,
    matrix.data(),
    &rows,
    &zero,
    result.data(),
    &columns FCONE FCONE
  );
}

[[maybe_unused]] static void accumulate_weighted_crossproducts(
    std::vector<double>& predictors,
    std::vector<double>& responses,
    const std::vector<double>& weights,
    int rows,
    int columns,
    std::vector<double>& predictor_gram,
    std::vector<double>& predictor_response,
    std::vector<double>& response_ss
) {
  for (int row = 0; row < rows; ++row) {
    const double root_weight = std::sqrt(std::max(0.0, weights[row]));
    for (int column = 0; column < columns; ++column) {
      predictors[static_cast<size_t>(column) * rows + row] *= root_weight;
      responses[static_cast<size_t>(column) * rows + row] *= root_weight;
    }
  }

  const char transpose = 'T';
  const char no_transpose = 'N';
  const double one = 1.0;
  F77_CALL(dgemm)(
    &transpose,
    &no_transpose,
    &columns,
    &columns,
    &rows,
    &one,
    predictors.data(),
    &rows,
    predictors.data(),
    &rows,
    &one,
    predictor_gram.data(),
    &columns FCONE FCONE
  );

  F77_CALL(dgemm)(
    &transpose,
    &no_transpose,
    &columns,
    &columns,
    &rows,
    &one,
    responses.data(),
    &rows,
    predictors.data(),
    &rows,
    &one,
    predictor_response.data(),
    &columns FCONE FCONE
  );

  for (int target = 0; target < columns; ++target) {
    const double* response = &responses[static_cast<size_t>(target) * rows];
    double sum_squares = 0.0;
    for (int row = 0; row < rows; ++row) {
      sum_squares += response[row] * response[row];
    }
    response_ss[target] += sum_squares;
  }
}

static double subset_sse(
    const std::vector<double>& gram,
    const std::vector<double>& xty,
    const std::vector<int>& subset,
    const std::vector<double>& beta,
    double y_ss,
    int p
) {
  double explained = 0.0;
  for (int i = 0; i < static_cast<int>(subset.size()); ++i) {
    explained += 2.0 * beta[i] * xty[subset[i]];
    for (int j = 0; j < static_cast<int>(subset.size()); ++j) {
      explained -= beta[i] * beta[j] * gram[static_cast<size_t>(subset[i]) * p + subset[j]];
    }
  }
  const double sse = y_ss - explained;
  return std::max(0.0, sse);
}

static double subset_bic(double sse, int support_size, int n_obs) {
  const double n = static_cast<double>(std::max(1, n_obs));
  return n * std::log(std::max(sse / n, 1e-12)) +
    static_cast<double>(support_size) * std::log(n);
}

static TargetFit fit_target_greedy_from_gram(
    const std::vector<double>& gram,
    const std::vector<double>& xty,
    const std::vector<char>& allowed,
    int target,
    int p,
    int n_obs,
    int max_support,
    double min_improvement,
    double y_ss
) {
  TargetFit out;
  out.coefficient.assign(p, 0.0);
  const double baseline_bic = subset_bic(y_ss, 0, n_obs);
  out.rss = y_ss;
  out.bic = baseline_bic;
  if (y_ss <= 0.0 || max_support <= 0 || n_obs < 2) {
    return out;
  }

  std::vector<int> selected;
  std::vector<char> used(p, 0);
  if (target >= 0 && target < p) {
    used[target] = 1;
  }
  double path_sse = y_ss;
  double path_bic = baseline_bic;
  double best_path_bic = baseline_bic;
  std::vector<int> best_path_selected;
  std::vector<double> path_beta;
  const double forward_bic_tolerance = 1e-10;

  for (int step = 0; step < max_support; ++step) {
    std::vector<double> path_inverse;
    if (!invert_subset_gram(gram, selected, p, path_inverse)) {
      break;
    }
    std::vector<double> path_projection;
    std::vector<double> path_variance;
    std::vector<double> path_cross;
    std::vector<char> path_valid;
    residual_candidate_statistics(
      gram,
      xty,
      selected,
      path_beta,
      path_inverse,
      allowed,
      used,
      p,
      path_projection,
      path_variance,
      path_cross,
      path_valid
    );
    int best_regulator = -1;
    double best_sse = path_sse;
    for (int regulator = 0; regulator < p; ++regulator) {
      if (!path_valid[regulator]) {
        continue;
      }
      const double sse = std::max(
        0.0,
        path_sse - path_cross[regulator] * path_cross[regulator] /
          path_variance[regulator]
      );
      if (sse >= path_sse - min_improvement) {
        continue;
      }
      const double tie_tolerance = 1e-12 * (1.0 + std::fabs(best_sse));
      if (best_regulator < 0 || sse < best_sse - tie_tolerance ||
          (std::fabs(sse - best_sse) <= tie_tolerance && regulator < best_regulator)) {
        best_sse = sse;
        best_regulator = regulator;
      }
    }
    if (best_regulator < 0) {
      break;
    }
    std::vector<int> candidate_selected = selected;
    candidate_selected.push_back(best_regulator);
    std::vector<double> candidate_beta;
    if (!fit_subset_coefficients(
      gram,
      xty,
      candidate_selected,
      p,
      candidate_beta
    )) {
      break;
    }
    best_sse = subset_sse(
      gram,
      xty,
      candidate_selected,
      candidate_beta,
      y_ss,
      p
    );
    if (best_sse >= path_sse - min_improvement) {
      break;
    }
    const double candidate_bic = subset_bic(
      best_sse,
      static_cast<int>(candidate_selected.size()),
      n_obs
    );
    const double current_tolerance = forward_bic_tolerance *
      (1.0 + std::fabs(path_bic));
    if (candidate_bic >= path_bic - current_tolerance) {
      break;
    }
    selected.swap(candidate_selected);
    path_beta.swap(candidate_beta);
    used[best_regulator] = 1;
    path_sse = best_sse;
    path_bic = candidate_bic;

    if (path_bic < best_path_bic) {
      best_path_bic = path_bic;
      best_path_selected = selected;
    }
  }

  if (best_path_selected.empty()) {
    return out;
  }

  std::vector<int> current_support = best_path_selected;
  std::sort(current_support.begin(), current_support.end());
  std::vector<double> current_beta;
  if (!fit_subset_coefficients(gram, xty, current_support, p, current_beta)) {
    return out;
  }
  double current_sse = subset_sse(
    gram,
    xty,
    current_support,
    current_beta,
    y_ss,
    p
  );
  double current_bic = subset_bic(
    current_sse,
    static_cast<int>(current_support.size()),
    n_obs
  );
  const double bic_tolerance = 1e-10;

  while (true) {
    std::vector<int> best_support = current_support;
    std::vector<double> best_beta = current_beta;
    double best_sse = current_sse;
    double best_bic = current_bic;
    bool best_beta_valid = true;

    const int support_size = static_cast<int>(current_support.size());
    std::vector<char> current_used(p, 0);
    for (int regulator : current_support) {
      current_used[regulator] = 1;
    }
    std::vector<double> current_inverse;
    if (!invert_subset_gram(
      gram,
      current_support,
      p,
      current_inverse
    )) {
      break;
    }
    std::vector<double> outside_projection;
    std::vector<double> outside_variance;
    std::vector<double> outside_cross;
    std::vector<char> outside_valid;
    residual_candidate_statistics(
      gram,
      xty,
      current_support,
      current_beta,
      current_inverse,
      allowed,
      current_used,
      p,
      outside_projection,
      outside_variance,
      outside_cross,
      outside_valid
    );

    const auto consider_scored = [
      &best_support,
      &best_sse,
      &best_bic,
      &best_beta_valid,
      current_bic,
      bic_tolerance,
      n_obs
    ](std::vector<int> trial_support, double trial_sse) {
      std::sort(trial_support.begin(), trial_support.end());
      const double trial_bic = subset_bic(
        trial_sse,
        static_cast<int>(trial_support.size()),
        n_obs
      );
      const double current_tol = bic_tolerance * (1.0 + std::fabs(current_bic));
      if (trial_bic >= current_bic - current_tol) {
        return;
      }
      const double best_tol = bic_tolerance * (1.0 + std::fabs(best_bic));
      if (trial_bic < best_bic - best_tol ||
          (std::fabs(trial_bic - best_bic) <= best_tol &&
           trial_support < best_support)) {
        best_support.swap(trial_support);
        best_sse = trial_sse;
        best_bic = trial_bic;
        best_beta_valid = false;
      }
    };

    const auto score_is_competitive = [
      &best_bic,
      current_bic,
      bic_tolerance
    ](double trial_bic) {
      const double current_tol = bic_tolerance * (1.0 + std::fabs(current_bic));
      if (trial_bic >= current_bic - current_tol) {
        return false;
      }
      const double best_tol = bic_tolerance * (1.0 + std::fabs(best_bic));
      return trial_bic <= best_bic + best_tol;
    };

    std::vector<double> removal_sse(support_size, current_sse);
    std::vector<double> removal_variance(support_size, 0.0);
    for (int remove = 0; remove < support_size; ++remove) {
      const double inverse_diagonal =
        current_inverse[static_cast<size_t>(remove) * support_size + remove];
      if (inverse_diagonal <= 0.0 || !R_finite(inverse_diagonal)) {
        continue;
      }
      const double residual_variance = 1.0 / inverse_diagonal;
      removal_variance[remove] = residual_variance;
      removal_sse[remove] = std::max(
        0.0,
        current_sse + current_beta[remove] * current_beta[remove] *
          residual_variance
      );
      const double trial_bic = subset_bic(
        removal_sse[remove],
        support_size - 1,
        n_obs
      );
      if (!score_is_competitive(trial_bic)) {
        continue;
      }
      std::vector<int> trial = current_support;
      trial.erase(trial.begin() + remove);
      consider_scored(trial, removal_sse[remove]);
    }

    if (support_size < max_support) {
      for (int regulator = 0; regulator < p; ++regulator) {
        if (!outside_valid[regulator]) {
          continue;
        }
        const double trial_sse = std::max(
          0.0,
          current_sse - outside_cross[regulator] * outside_cross[regulator] /
            outside_variance[regulator]
        );
        const double trial_bic = subset_bic(
          trial_sse,
          support_size + 1,
          n_obs
        );
        if (!score_is_competitive(trial_bic)) {
          continue;
        }
        std::vector<int> trial = current_support;
        trial.push_back(regulator);
        consider_scored(trial, trial_sse);
      }
    }

    if (support_size > 0) {
      for (int remove = 0; remove < support_size; ++remove) {
        const double removed_variance = removal_variance[remove];
        if (removed_variance <= 0.0) {
          continue;
        }
        for (int regulator = 0; regulator < p; ++regulator) {
          if (!allowed[regulator] || current_used[regulator]) {
            continue;
          }
          const double removed_projection = outside_projection[
            static_cast<size_t>(regulator) * support_size + remove
          ];
          const double trial_variance = outside_variance[regulator] +
            removed_projection * removed_projection * removed_variance;
          const double trial_variance_tolerance = 1e-10 * std::max(
            1.0,
            std::fabs(gram[static_cast<size_t>(regulator) * p + regulator])
          );
          if (trial_variance <= trial_variance_tolerance ||
              !R_finite(trial_variance)) {
            continue;
          }
          const double trial_cross = outside_cross[regulator] +
            current_beta[remove] * removed_projection * removed_variance;
          const double trial_sse = std::max(
            0.0,
            removal_sse[remove] - trial_cross * trial_cross / trial_variance
          );
          if (trial_sse >= current_sse) {
            continue;
          }
          const double trial_bic = subset_bic(
            trial_sse,
            support_size,
            n_obs
          );
          if (!score_is_competitive(trial_bic)) {
            continue;
          }
          std::vector<int> trial = current_support;
          trial[remove] = regulator;
          consider_scored(trial, trial_sse);
        }
      }
    }

    if (best_support == current_support) {
      break;
    }
    if (!best_beta_valid) {
      if (!fit_subset_coefficients(
        gram,
        xty,
        best_support,
        p,
        best_beta
      )) {
        break;
      }
      best_sse = subset_sse(
        gram,
        xty,
        best_support,
        best_beta,
        y_ss,
        p
      );
      best_bic = subset_bic(
        best_sse,
        static_cast<int>(best_support.size()),
        n_obs
      );
      const double current_tol = bic_tolerance *
        (1.0 + std::fabs(current_bic));
      if (best_bic >= current_bic - current_tol) {
        break;
      }
    }
    current_support.swap(best_support);
    current_beta.swap(best_beta);
    current_sse = best_sse;
    current_bic = best_bic;
  }

  for (int idx = 0; idx < static_cast<int>(current_support.size()); ++idx) {
    out.coefficient[current_support[idx]] = current_beta[idx];
  }
  out.r2 = std::max(0.0, 1.0 - current_sse / y_ss);
  out.rss = current_sse;
  out.bic = current_bic;
  out.support = current_support;
  return out;
}

// [[Rcpp::export]]
List solve_greedy_l0(
    NumericMatrix x,
    NumericVector y,
    int max_support_size,
    double min_improvement
) {
  const int n_total = x.nrow();
  const int p = x.ncol();
  NumericVector coefficient(p, 0.0);
  if (p == 0 || n_total < 3 || y.size() != n_total) {
    return List::create(
      _["coefficient"] = coefficient,
      _["support"] = IntegerVector(),
      _["r_squared"] = 0.0,
      _["rss"] = NA_REAL,
      _["bic"] = NA_REAL,
      _["n_obs"] = 0
    );
  }

  std::vector<int> rows;
  rows.reserve(n_total);
  for (int row = 0; row < n_total; ++row) {
    if (R_finite(y[row])) {
      rows.push_back(row);
    }
  }
  const int n_obs = static_cast<int>(rows.size());
  if (n_obs < 3) {
    return List::create(
      _["coefficient"] = coefficient,
      _["support"] = IntegerVector(),
      _["r_squared"] = 0.0,
      _["rss"] = NA_REAL,
      _["bic"] = NA_REAL,
      _["n_obs"] = n_obs
    );
  }

  double y_mean = 0.0;
  for (int row : rows) {
    y_mean += y[row];
  }
  y_mean /= static_cast<double>(n_obs);

  std::vector<double> standardized_y(n_obs, 0.0);
  double y_ss = 0.0;
  for (int i = 0; i < n_obs; ++i) {
    const double value = y[rows[i]] - y_mean;
    standardized_y[i] = value;
    y_ss += value * value;
  }
  if (y_ss > 0.0) {
    const double y_scale = std::sqrt(y_ss / static_cast<double>(n_obs - 1));
    y_ss = 0.0;
    for (int i = 0; i < n_obs; ++i) {
      standardized_y[i] /= y_scale;
      y_ss += standardized_y[i] * standardized_y[i];
    }
  }

  std::vector<double> standardized_x(static_cast<size_t>(p) * n_obs, 0.0);
  for (int column = 0; column < p; ++column) {
    double sum = 0.0;
    int finite_count = 0;
    for (int row : rows) {
      const double value = x(row, column);
      if (R_finite(value)) {
        sum += value;
        ++finite_count;
      }
    }
    const double mean = finite_count > 0 ? sum / static_cast<double>(finite_count) : 0.0;
    double ss = 0.0;
    for (int i = 0; i < n_obs; ++i) {
      const double value = x(rows[i], column);
      const double centered_value = R_finite(value) ? value - mean : 0.0;
      standardized_x[static_cast<size_t>(column) * n_obs + i] = centered_value;
      ss += centered_value * centered_value;
    }
    if (finite_count > 1 && ss > 0.0) {
      const double scale = std::sqrt(ss / static_cast<double>(finite_count - 1));
      for (int i = 0; i < n_obs; ++i) {
        standardized_x[static_cast<size_t>(column) * n_obs + i] /= scale;
      }
    }
  }

  std::vector<double> gram(static_cast<size_t>(p) * p, 0.0);
  std::vector<double> xty(p, 0.0);
  for (int left = 0; left < p; ++left) {
    for (int i = 0; i < n_obs; ++i) {
      xty[left] += standardized_x[static_cast<size_t>(left) * n_obs + i] * standardized_y[i];
    }
    for (int right = 0; right <= left; ++right) {
      double dot = 0.0;
      for (int i = 0; i < n_obs; ++i) {
        dot += standardized_x[static_cast<size_t>(left) * n_obs + i] *
          standardized_x[static_cast<size_t>(right) * n_obs + i];
      }
      gram[static_cast<size_t>(left) * p + right] = dot;
      gram[static_cast<size_t>(right) * p + left] = dot;
    }
  }

  std::vector<char> allowed(p, 0);
  for (int column = 0; column < p; ++column) {
    allowed[column] = gram[static_cast<size_t>(column) * p + column] > 0.0;
  }
  const int requested_support = max_support_size <= 0
    ? p
    : std::min(max_support_size, p);
  const int support_cap = std::max(
    0,
    std::min(requested_support, n_obs - 2)
  );
  const TargetFit fit = fit_target_greedy_from_gram(
    gram,
    xty,
    allowed,
    -1,
    p,
    n_obs,
    support_cap,
    std::max(0.0, min_improvement),
    y_ss
  );
  for (int column = 0; column < p; ++column) {
    coefficient[column] = fit.coefficient[column];
  }
  IntegerVector support(fit.support.size());
  for (int i = 0; i < static_cast<int>(fit.support.size()); ++i) {
    support[i] = fit.support[i] + 1;
  }
  return List::create(
    _["coefficient"] = coefficient,
    _["support"] = support,
    _["r_squared"] = fit.r2,
    _["rss"] = fit.rss,
    _["bic"] = fit.bic,
    _["n_obs"] = n_obs
  );
}
