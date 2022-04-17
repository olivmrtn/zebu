#include "functions.h"
#include <Rcpp.h>
#include <algorithm>
#include <random>
#include <chrono>

using namespace Rcpp;


// [[Rcpp::export]]
DataFrame permute_dataframe(DataFrame &x, const List &group) {

  // Find seeds for each column group
  CharacterVector gr;
  CharacterVector col_names = x.names();
  std::vector<unsigned int> seeds(x.ncol());
  unsigned int seed;

  for (int i = 0; i != group.length(); ++i) {
    gr = group[i];
    seed = std::chrono::system_clock::now()
      .time_since_epoch()
      .count();

    for (int j = 0; j != x.ncol(); ++j) {
      for (int k = 0; k != gr.length(); ++k) {
        if (col_names[j] == gr[k]) {
          seeds[j] = seed;
        }
      }
    }
  }

  // Permute columns
  IntegerVector perm(x.nrow());
  for (int j = 0; j != x.ncol(); ++j) {
    perm = x[j];
    std::shuffle(perm.begin(), perm.end(), std::default_random_engine(seeds[j]));
    x[j] = perm;
  }

  return x;
}

// [[Rcpp::export]]
List permtest_rcpp(const List &x, const int &nb, const List &group) {

  // Retrieve data
  DataFrame df = clone(as<DataFrame>(as<List>(x["data"])["pp"]));
  std::string measure = as<List>(x["lassie_params"])["measure"];
  double nr = df.nrow();

  // Permutations
  List prob(3);
  List results(nb);
  for (int i = 0; i != nb; ++i) {
    df = permute_dataframe(df, group);
    prob = estimate_prob(df);
    results[i] = local_association(prob, measure, nr);
  }

  // List out;
  return results;
}
