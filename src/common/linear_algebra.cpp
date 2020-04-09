#include "common/linear_algebra.hpp"


namespace robotcgmres {
namespace linearalgebra {

double InnerProduct(const int dim, const double *vec1, const double *vec2) {
  double ans = 0;
  for (int i=0; i<dim; ++i) {
    ans += vec1[i] * vec2[i];
  }
  return ans;
}

double SquaredNorm(const int dim, const double *vec) {
  double ans = 0;
  for (int i=0; i<dim; ++i) {
    ans += vec[i] * vec[i];
  }
  return ans;
}

} // namespace linearalgebra
} // namespace robotcgmres
