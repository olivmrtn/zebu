#include "functions.h"
#include <Rcpp.h>
using namespace Rcpp;

//' @title Local Association Measures
//'
//' @description Subroutines called by \code{\link[zebu]{lassie}} to compute
//' local and global association measures from a list of probabilities.
//'
//' @param x list of probabilities as outputted by \code{\link[zebu]{estimate_prob}}.
//' @param measure name of measure to be used:
//' \itemize{
//' \item'chisq': Chi-squared residuals.
//' \item'd': Lewontin's D.
//' \item'z': Ducher's 'z'.
//' \item'pmi': Pointwise mutual information (in bits).
//' \item'npmi': Normalized pointwise mutual information (Bouma).
//' \item'npmi2': Normalized pointwise mutual information (Multivariate).
//' }
//' @param nr number of rows/samples. Only used to estimate chi-squared residuals.
//'
//' @return List containing the following values:
//' \itemize{
//' \item local: local association array (may contain NA, NaN and Inf values).
//' \item global: global association numeric value.
//' }
//'
//'
//' @seealso \code{\link[zebu]{lassie}}
//'
//' @example /inst/examples/lassie.R
//'
//' @name local_association
//'
//' @export
// [[Rcpp::export]]
List local_association(const List &x,
                       const std::string &measure = "chisq",
                       const double &nr = 1) {

  List out(2);

  if (measure == "z") {
    out = duchers_z(x);
  } else if (measure == "d") {
    out = lewontin_d(x);
  } else if (measure == "pmi") {
    out = pmi(x, 0);
  } else if (measure == "npmi") {
    out = pmi(x, 1);
  } else if (measure == "npmi2") {
    out = pmi(x, 2);
  } else if (measure == "chisq") {
    out = chisq(x, nr);
  } else {
    stop("Invalid 'measure' argument.");
  }

  return out;
}


//' @rdname local_association
//' @export
// [[Rcpp::export]]
List lewontin_d(const List &x) {

  // Retrieve probabilities
  NumericVector observed = x["observed"];
  NumericVector expected = x["expected"];
  List margins = x["margins"];

  // Dimensions
  int number_of_dims = margins.length();
  IntegerVector number_of_levels = observed.attr("dim");
  int number_of_cells = 1;
  for (int i = 0; i < number_of_dims; i++) {
    number_of_cells *= number_of_levels[i];
  }

  // Compute local association
  NumericVector local(number_of_cells);
  local.attr("dim") = number_of_levels;
  local.attr("dimnames") = observed.attr("dimnames");
  for (int i = 0; i < number_of_cells; ++i) {
    local[i] = observed[i] - expected[i];
  }

  // Compute global association
  double global = 0.0;
  for (int i = 0; i != number_of_cells; ++i) {
    global += local[i] * observed[i];
  }

  List out = Rcpp::List::create(
    Rcpp::Named("global")  = global,
    Rcpp::Named("local") = local);
  return out;
}

//' @rdname local_association
//' @export
// [[Rcpp::export]]
List duchers_z(const List &x) {

  // Retrieve R functions
  Environment pkg = Environment::namespace_env("zebu");
  Function compute_theoretical_min_prob = pkg[".compute_theoretical_min_prob"];
  Function compute_theoretical_max_prob = pkg[".compute_theoretical_max_prob"];

  // Retrieve probabilities
  NumericVector observed = x["observed"];
  NumericVector expected = x["expected"];
  List margins = x["margins"];

  // Dimensions
  int number_of_dims = margins.length();
  IntegerVector number_of_levels = observed.attr("dim");
  int number_of_cells = 1;
  for (int i = 0; i < number_of_dims; i++) {
    number_of_cells *= number_of_levels[i];
  }

  // Compute minimum marginal probability
  NumericVector mini = compute_theoretical_min_prob(margins);
  NumericVector maxi = compute_theoretical_max_prob(margins);

  // Compute local association
  NumericVector local(number_of_cells);
  local.attr("dim") = number_of_levels;
  local.attr("dimnames") = observed.attr("dimnames");
  for (int i = 0; i < number_of_cells; ++i) {
    local[i] = observed[i] - expected[i];
    if (local[i] > 0) {
      local[i] /= (maxi[i] - expected[i]);
    } else if (local[i] < 0) {
      local[i] /= (expected[i] - mini[i]);
    }
  }

  // Compute global association
  double global = 0.0;
  for (int i = 0; i != number_of_cells; ++i) {
    global += local[i] * observed[i];
  }

  List out = Rcpp::List::create(
    Rcpp::Named("global")  = global,
    Rcpp::Named("local") = local);
  return out;
}

//' @rdname local_association
//' @param normalize 0 for pmi, 1 for npmi, 2 for npmi2
//' @export
// [[Rcpp::export]]
List pmi(const List &x, const int &normalize) {

  // Retrieve probabilities
  NumericVector observed = x["observed"];
  NumericVector expected = x["expected"];
  List margins = x["margins"];

  // Dimensions
  int number_of_dims = margins.length();
  IntegerVector number_of_levels = observed.attr("dim");
  int number_of_cells = 1;
  for (int i = 0; i < number_of_dims; i++) {
    number_of_cells *= number_of_levels[i];
  }

  // Compute unnormalized local association
  NumericVector local(number_of_cells);
  local.attr("dim") = number_of_levels;
  local.attr("dimnames") = observed.attr("dimnames");
  for (int i = 0; i < number_of_cells; ++i) {
    local[i] = log2(observed[i]) - log2(expected[i]);
  }

  LogicalVector isinf = is_infinite(local);


  // Normalize
  if (normalize == 1) {
    for (int i = 0; i < number_of_cells; ++i) {
      if (isinf[i]) {
        local[i] = -1.0;
      } else {
        local[i] /= -log2(observed[i]);
      }
    }
  } else if (normalize == 2) {

    Environment pkg = Environment::namespace_env("zebu");
    Function compute_theoretical_max_prob = pkg[".compute_theoretical_max_prob"];
    NumericVector maxi = compute_theoretical_max_prob(margins);

    for (int i = 0; i < number_of_cells; ++i) {
      if (isinf[i]) {
        local[i] = -1.0;
      } else if (local[i] > 0) {
        local[i] /= (log2(maxi[i]) - log2(expected[i]));
      } else if (local[i] < 0) {
        local[i] /= -log2(observed[i]);
      }
    }
  }

  // Compute global association
  double global = 0.0;
  for (int i = 0; i != number_of_cells; ++i) {
    if (! isinf[i]) {
      global += local[i] * observed[i];
    }
  }

  List out = Rcpp::List::create(
    Rcpp::Named("global")  = global,
    Rcpp::Named("local") = local);
  return out;
}

//' @rdname local_association
//' @export
// [[Rcpp::export]]
List chisq(const List &x, const double &nr) {

  // Retrieve probabilities
  NumericVector observed = x["observed"];
  NumericVector expected = x["expected"];
  List margins = x["margins"];

  // Dimensions
  int number_of_dims = margins.length();
  IntegerVector number_of_levels = observed.attr("dim");
  int number_of_cells = 1;
  for (int i = 0; i < number_of_dims; i++) {
    number_of_cells *= number_of_levels[i];
  }

  // Sqrt of number of rows
  double sqrt_nr = sqrt(nr);

  // Compute local association
  NumericVector local(number_of_cells);
  local.attr("dim") = number_of_levels;
  local.attr("dimnames") = observed.attr("dimnames");
  for (int i = 0; i != number_of_cells; ++i) {
    local[i] = sqrt_nr * (observed[i] - expected[i]) / sqrt(expected[i]);
  }

  // Compute global association
  double global = 0.0;
  for (int i = 0; i != number_of_cells; ++i) {
    global += pow(local[i], 2);
  }

  List out = Rcpp::List::create(
    Rcpp::Named("global")  = global,
    Rcpp::Named("local") = local);
  return out;
}
