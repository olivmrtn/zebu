#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_observed_probability(const DataFrame &x,
                                           const int &number_of_samples,
                                           const int &number_of_dims,
                                           const int &number_of_cells,
                                           const IntegerVector &number_of_levels) {


  // Index offset for computing observed probabilities
  IntegerVector offset(number_of_dims, 1);
  for (int i = 1; i != number_of_dims; ++i) {
    for (int j = i; j != number_of_dims; ++j) {
      offset[j] *= number_of_levels[i-1];
    }
  }

  // Compute observed probabilities
  int index;
  double unit = 1.0 / number_of_samples;
  NumericVector observed(number_of_cells);
  observed.attr("dim") = Dimension(number_of_levels);
  for (int i = 0; i != number_of_samples; ++i) {
    index = 0;
    for (int j = 0; j != number_of_dims; ++j) {
      index += (as<IntegerVector>(x[j])[i] - 1) * offset[j];
    }
    observed[index] += unit;
  }

  return observed;
}

// [[Rcpp::export]]
NumericVector compute_expected_probability(const List &margins) {
  Environment pkg = Environment::namespace_env("zebu");
  Function f = pkg[".compute_expected_probability"];
  return f(margins);
}

// [[Rcpp::export]]
List compute_marginal_probability(const DataFrame &x,
                                  const int &number_of_samples,
                                  const int &number_of_dims,
                                  const int &number_of_cells,
                                  const IntegerVector &number_of_levels) {

  int index;
  double unit = 1.0 / number_of_samples;

  IntegerVector col(number_of_samples);
  List margins(number_of_dims);

  for (int j = 0; j != number_of_dims; ++j) {
    NumericVector margin(number_of_levels[j], 0.0);

    col = x[j];
    for (int i = 0; i != number_of_samples; ++i) {
      index = col[i] - 1;
      margin[index] += unit;
    }

    margins[j] = margin;
  }

  return margins;

}

//' Estimate marginal and multivariate probabilities
//'
//' Maximum-likelihood estimation of marginal and multivariate observed and expected independence probabilities. Marginal probability refers to probability of each factor per individual column. Multivariate probability refer to cross-classifying factors for all columns.
//'
//' @param x data.frame or matrix.
//'
//' @return List containing the following values:
//'   \itemize{
//'   \item margins: a list of marginal probabilities. Names correspond to colnames(x).
//'   \item observed: observed multivariate probability array.
//'   \item expected: expected multivariate probability array
//'   }
//'
//' @example /inst/examples/lassie.R
//'
//' @export
// [[Rcpp::export]]
List estimate_prob(const DataFrame &x) {

  // Retrieve dimensions
  int number_of_samples = x.nrow();
  int number_of_dims    = x.ncol();

  // Check if there are at least two variables
  if (number_of_dims < 2) {
    stop("Invalid 'x' argument: needs to have at least two columns.");
  }

  // Number of levels per variable and in total
  IntegerVector number_of_levels(number_of_dims);
  int number_of_cells = 1;
  for (int j = 0; j != number_of_dims; ++j) {
    number_of_levels[j] = as<CharacterVector>((as<IntegerVector>(x[j]).attr("levels"))).length();
    number_of_cells *= number_of_levels[j];
  }

  // Retrieve dimnames
  List dim_names(number_of_dims);
  dim_names.names() = x.names();
  for (int j = 0; j != number_of_dims; ++j) {
    dim_names[j] = as<IntegerVector>(x[j]).attr("levels");
  }

  // Compute marginal probabilities
  List margins(number_of_dims);
  margins = compute_marginal_probability(x,
                                         number_of_samples,
                                         number_of_dims,
                                         number_of_cells,
                                         number_of_levels);
  margins.names() = dim_names.names();
  for (int j = 0; j != number_of_dims; ++j) {
    as<NumericVector>(margins[j]).names() = dim_names[j];
  }

  // Compute observed probabilities
  NumericVector observed(number_of_cells);
  observed.attr("dim") = Dimension(number_of_levels);
  observed = compute_observed_probability(x,
                                          number_of_samples,
                                          number_of_dims,
                                          number_of_cells,
                                          number_of_levels);
  observed.attr("dimnames") = dim_names;

  // Compute expected probabilities
  NumericVector expected(number_of_cells);
  expected.attr("dim") = Dimension(number_of_levels);
  expected = compute_expected_probability(margins);
  expected.attr("dimnames") = dim_names;


  // Return results
  List out = Rcpp::List::create(
    Rcpp::Named("margins")  = margins,
    Rcpp::Named("observed") = observed,
    Rcpp::Named("expected") = expected);

  return out;
}
