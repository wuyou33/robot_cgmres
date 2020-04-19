#ifndef ROBOTCGMRES_FDGMRES_HPP_
#define ROBOTCGMRES_FDGMRES_HPP_

#include <iostream>
#include <cmath>
#include <limits>

#include "common/memory_manager.hpp"
#include "common/linear_algebra.hpp"


namespace robotcgmres {

template <class NewtonTypeSystem>
class FDGMRES {
public:

  // Allocates robot_ptrarrays used in FDGMRES iteration.
  // Argments:
  //   dim_linear_problem: The dimension of the solution of the Newton system's 
  //      equation.
  //   kmax: The maximum number of GMRES iteration.
  FDGMRES(const unsigned int dim_linear_problem, const unsigned int kmax)
    : dim_linear_problem_(dim_linear_problem), 
      kmax_(kmax), 
      hessenberg_mat_(memorymanager::NewMatrix(kmax+1, kmax+1)), 
      basis_mat_(memorymanager::NewMatrix(kmax+1, dim_linear_problem)), 
      r_(memorymanager::NewVector(dim_linear_problem)), 
      givens_c_(memorymanager::NewVector(kmax+1)), 
      givens_s_(memorymanager::NewVector(kmax+1)), 
      g_(memorymanager::NewVector(kmax+1)) {
    if (kmax > dim_linear_problem) {
      kmax_ = dim_linear_problem;
    }
  }

  // Frees arrays used in FDGMRES iteration.
  ~FDGMRES() {
    memorymanager::DeleteMatrix(hessenberg_mat_);
    memorymanager::DeleteMatrix(basis_mat_);
    memorymanager::DeleteVector(r_);
    memorymanager::DeleteVector(givens_c_);
    memorymanager::DeleteVector(givens_s_);
    memorymanager::DeleteVector(g_);
  }

  // Computes the update quantity of the solution of newton_type_system under 
  // current t, q, v, solution. The result is stored in solution_update and 
  // it must have intial guess value before call this function.
  // Argments: 
  //   t: Time.
  //   q: Generalized configuration. The size must be dimq.
  //   v: Generalized velocity. The size must be dimv.
  //   solution: Current iterate solution. The size must be dim_linear_problem.
  //   solution_update: Update of the solution. The size must be dim_linear_problem.
  void computeSolutionUpdate(NewtonTypeSystem& newton_type_system,
                             const double t, const double* q, const double* v,
                             const double* solution, double* solution_update) {
    // Initializes vectors for QR factrization by Givens rotation.
    // Set givens_c_, givens_s_, g_ as zero.
    for (int i=0; i<kmax_+1; ++i) givens_c_[i] = 0;
    for (int i=0; i<kmax_+1; ++i) givens_s_[i] = 0;
    for (int i=0; i<kmax_+1; ++i) g_[i] = 0;
    // Generates the initial basis of the Krylov subspace.
    newton_type_system.computeInitialResidual(t, q, v, solution, 
                                              solution_update, r_);
    g_[0] = std::sqrt(linearalgebra::SquaredNorm(dim_linear_problem_, r_));
    // basis_mat_[0] = r_ / g_[0]
    for (int i=0; i<dim_linear_problem_; ++i) {
      basis_mat_[0][i] = r_[i] / g_[0];
    }
    // k : the dimension of the Krylov subspace at the current iteration.
    int k;
    for (k=0; k<kmax_; ++k) {
      newton_type_system.computeAx(t, q, v, solution_update, basis_mat_[k], 
                                    basis_mat_[k+1]);
      for (int j=0; j<=k; ++j) {
        hessenberg_mat_[k][j] = linearalgebra::InnerProduct(
            dim_linear_problem_, basis_mat_[k+1], basis_mat_[j]);
        // basis_mat_[k+1] -= hessenberg_mat_[k][j] * basis_mat_[j];
        for (int i=0; i<dim_linear_problem_; ++i) {
          basis_mat_[k+1][i] -= hessenberg_mat_[k][j] * basis_mat_[j][i];
        }
      }
      hessenberg_mat_[k][k+1] = std::sqrt(linearalgebra::SquaredNorm(
          dim_linear_problem_, basis_mat_[k+1]));
      if (std::abs(hessenberg_mat_[k][k+1]) 
          < std::numeric_limits<double>::epsilon()) {
#ifdef RUNTIME_CHECKS
        std::cout << "The modified Gram-Schmidt breakdown at k = " << k 
                  << std::endl;
#endif // RUNTIME_CHECKS
        break;
      }
      // basis_mat_[k+1] = basis_mat_[k+1] / hessenberg_mat_[k][k+1];
      for (int i=0; i<dim_linear_problem_; ++i) {
        basis_mat_[k+1][i] /= hessenberg_mat_[k][k+1];
      }
      // Givens Rotation for QR factrization of the least squares problem.
      for (int j=0; j<k; ++j) {
        applyGivensRotation(hessenberg_mat_[k], j);
      }
      double nu = std::sqrt(hessenberg_mat_[k][k]*hessenberg_mat_[k][k]
                            +hessenberg_mat_[k][k+1]*hessenberg_mat_[k][k+1]);
      givens_c_[k] = hessenberg_mat_[k][k] / nu;
      givens_s_[k] = - hessenberg_mat_[k][k+1] / nu;
      hessenberg_mat_[k][k] = givens_c_[k] * hessenberg_mat_[k][k] 
                              - givens_s_[k] * hessenberg_mat_[k][k+1];
      hessenberg_mat_[k][k+1] = 0;
      applyGivensRotation(g_, k);
    }
    // Computes solution_vec by solving hessenberg_mat_ * y = g_.
    for (int i=k-1; i>=0; --i) {
      double tmp = g_[i];
      for (int j=i+1; j<k; ++j) {
        tmp -= hessenberg_mat_[j][i] * givens_c_[j];
      }
      givens_c_[i] = tmp / hessenberg_mat_[i][i];
    }
    for (int i=0; i<dim_linear_problem_; ++i) {
      double tmp = 0;
      for (int j=0; j<k; ++j) { 
        tmp += basis_mat_[j][i] * givens_c_[j];
      }
      solution_update[i] += tmp;
    }
  }

  // Prohibits copy constructor.
  FDGMRES(const FDGMRES&) = delete;

  // Prohibits copy operator.
  FDGMRES& operator=(const FDGMRES&) = delete;

private:
  int dim_linear_problem_, kmax_;
  double **hessenberg_mat_, **basis_mat_, *r_, *givens_c_, *givens_s_,*g_;

  // Applies the Givens rotation for row_index element and 
  // row_index+1 element of column_vec, which is a column vector of 
  // a matrix. 
  // Argments: 
  //   column_vec: A column vector of a matrix that is factrized by 
  //      Givens rotation.
  //   row_index: Index of the row where Givens rotation applied.
  inline void applyGivensRotation(double* column_vec, const int row_index) {
    double tmp1 = givens_c_[row_index] * column_vec[row_index] 
                      - givens_s_[row_index] * column_vec[row_index+1];
    double tmp2 = givens_s_[row_index] * column_vec[row_index] 
                      + givens_c_[row_index] * column_vec[row_index+1];
    column_vec[row_index] = tmp1;
    column_vec[row_index+1] = tmp2;
  }

};

} // namespace robotcgmres


#endif // ROBOTCGMRES_FDGMRES_HPP_ 