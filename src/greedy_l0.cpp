#include "greedy_l0.h"

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
  NumericVector deletion_delta_bic(p, NA_REAL);
  if (p == 0 || n_total < 3 || y.size() != n_total) {
    return List::create(
      _["coefficient"] = coefficient,
      _["deletion_delta_bic"] = deletion_delta_bic,
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
      _["deletion_delta_bic"] = deletion_delta_bic,
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
  if (!fit.support.empty()) {
    std::vector<double> support_inverse;
    if (!invert_subset_gram(gram, fit.support, p, support_inverse)) {
      stop("Selected support Gram matrix is not identifiable.");
    }
    const int selected_count = static_cast<int>(fit.support.size());
    for (int selected = 0; selected < selected_count; ++selected) {
      const int column = fit.support[selected];
      deletion_delta_bic[column] = selected_deletion_delta_bic(
        fit.coefficient[column],
        support_inverse[static_cast<size_t>(selected) * selected_count + selected],
        fit, selected_count, n_obs
      );
    }
  }
  IntegerVector support(fit.support.size());
  for (int i = 0; i < static_cast<int>(fit.support.size()); ++i) {
    support[i] = fit.support[i] + 1;
  }
  return List::create(
    _["coefficient"] = coefficient,
      _["deletion_delta_bic"] = deletion_delta_bic,
    _["support"] = support,
    _["r_squared"] = fit.r2,
    _["rss"] = fit.rss,
    _["bic"] = fit.bic,
    _["n_obs"] = n_obs
  );
}
