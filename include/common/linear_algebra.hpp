#ifndef LINEAR_ALGEBRA_HPP_
#define LINEAR_ALGEBRA_HPP_

#include <cmath>


namespace robotcgmres {
namespace linearalgebra {

// Returns inner product of vec_1 and vec_2.
double InnerProduct(const int dim, const double *vec1, const double *vec2);

// Returns squared norm of vec.
double SquaredNorm(const int dim, const double *vec);

} // namespace linearalgebra
} // namespace robotcgmres


#endif // LINEAR_ALGEBRA_HPP_