#include <gtest/gtest.h>
#include "Eigen/Core"
#include "Eigen/LU"
#include <cmath>
#include <cstdlib>
#include <limits>

#include "ocp/fdgmres.hpp"
#include "common/memory_manager.hpp"
#include "common/linear_algebra.hpp"


namespace robotcgmres {

class LinearProblemSystem {
public:
  LinearProblemSystem(const int dim, const Eigen::MatrixXd& A_mat)
    : dim_(dim),
      A_mat_(A_mat) {
    if (A_mat_.rows() != dim) {
      std::abort();
    }
    if (A_mat_.cols() != dim) {
      std::abort();
    }
  }

  ~LinearProblemSystem() {
  }

  void computeInitialResidual(const double t, const double* q, const double* v, 
                              const double* solution, 
                              const double* initial_direction, 
                              double* initial_residual) {
    Eigen::Map<Eigen::VectorXd>(initial_residual, dim_) 
        = Eigen::Map<const Eigen::VectorXd>(solution, dim_) 
          - A_mat_ * Eigen::Map<const Eigen::VectorXd>(initial_direction, dim_);
  }

  void computeAx(const double t, const double* q, const double* v, 
                 const double* solution, const double* direction, double* Ax) {
    Eigen::Map<Eigen::VectorXd>(Ax, dim_) 
        = A_mat_ * Eigen::Map<const Eigen::VectorXd>(direction, dim_);
  }

private:
  int dim_;
  Eigen::MatrixXd A_mat_;

};

class FDGMRESTest: public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    dim_ = 50 + (int)std::abs(Eigen::VectorXd::Random(2)[1]) % 50;
    solution_ = memorymanager::NewVector(dim_);
    Eigen::Map<Eigen::VectorXd>(solution_, dim_) = Eigen::VectorXd::Random(dim_);
    solution_update_ = memorymanager::NewVector(dim_);
    Eigen::Map<Eigen::VectorXd>(solution_update_, dim_) = Eigen::VectorXd::Random(dim_);
    A_mat_ = Eigen::MatrixXd::Random(dim_, dim_);
    // A_mat_ must be nonsingular.
    while (A_mat_.determinant() == 0) {
      A_mat_ = Eigen::MatrixXd::Random(dim_, dim_);
    }
  }

  virtual void TearDown() {
    memorymanager::DeleteVector(solution_);
    memorymanager::DeleteVector(solution_update_);
  }

  int dim_;
  double *solution_, *solution_update_;
  Eigen::MatrixXd A_mat_;
};


TEST_F(FDGMRESTest, computeSolutionUpdate) {
  LinearProblemSystem linear_problem(dim_, A_mat_);
  double residual = (Eigen::Map<Eigen::VectorXd>(solution_, dim_) 
      - A_mat_ * Eigen::Map<Eigen::VectorXd>(solution_update_, dim_)).norm();
  double old_residual = residual;
  double *empty_ptr = nullptr;
  for (int kmax=1; kmax<=dim_; ++kmax) {
    FDGMRES<LinearProblemSystem> gmres(dim_, kmax);
    gmres.computeSolutionUpdate(linear_problem, 0, empty_ptr, empty_ptr, 
                                solution_, solution_update_);
    old_residual = residual;
    residual = (Eigen::Map<Eigen::VectorXd>(solution_, dim_) 
      - A_mat_ * Eigen::Map<Eigen::VectorXd>(solution_update_, dim_)).norm();
    EXPECT_TRUE(residual <= old_residual);
  }
  EXPECT_TRUE((A_mat_.inverse()*Eigen::Map<Eigen::VectorXd>(solution_, dim_))
              .isApprox(Eigen::Map<Eigen::VectorXd>(solution_update_, dim_), 
                        std::sqrt(std::numeric_limits<double>::epsilon())));
}

} // namespace robotcgmres


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}