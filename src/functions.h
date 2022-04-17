#include <Rcpp.h>
using namespace Rcpp;

List estimate_prob(const DataFrame &x);

List local_association(const List &x,
                       const std::string &measure,
                       const double &nr);

List lewontin_d(const List &prob);

List duchers_z(const List &prob);

List pmi(const List &prob, const int &normalize);

List chisq(const List &prob, const double &nr);
