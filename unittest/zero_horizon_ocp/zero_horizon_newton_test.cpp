#include <string>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robot/robot.hpp"
#include "zero_horizon_ocp/zero_horizon_newton.hpp"
#include "zero_horizon_ocp/zero_horizon_ocp.hpp"
#include "solver/fdgmres.hpp"
#include "solver/newton_system_for_ocp.hpp"
#include "common/memory_manager.hpp"
#include "common/linear_algebra.hpp"

#include "iiwa14/cost_function_iiwa14.hpp"


namespace robotcgmres {

class ZeroHorizonNewtonTest: public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    kmax_ = 5;
    fixed_urdf_file_name_ = "../../../examples/iiwa14/iiwa14.urdf";
    pinocchio::urdf::buildModel(fixed_urdf_file_name_, fixed_model_);
    q_ = memorymanager::NewVector(fixed_model_.nq);
    v_ = memorymanager::NewVector(fixed_model_.nv);
    initial_solution_ = memorymanager::NewVector(fixed_model_.nv);
    Eigen::Map<Eigen::VectorXd>(q_, fixed_model_.nq) 
        = pinocchio::randomConfiguration(fixed_model_, 
                                         -Eigen::VectorXd::Ones(fixed_model_.nq), 
                                         Eigen::VectorXd::Ones(fixed_model_.nq));
    Eigen::Map<Eigen::VectorXd>(v_, fixed_model_.nv) 
        = Eigen::VectorXd::Random(fixed_model_.nv);
    Eigen::Map<Eigen::VectorXd>(initial_solution_, fixed_model_.nv) 
        = 10*Eigen::VectorXd::Ones(fixed_model_.nv);
    t_ = Eigen::VectorXd::Random(2)[1];
    baumgarte_alpha_ = Eigen::VectorXd::Random(2)[1];
    baumgarte_beta_ = Eigen::VectorXd::Random(2)[1];
    finite_difference_ = 1.0e-08;
  }

  virtual void TearDown() {
    memorymanager::DeleteVector(q_);
    memorymanager::DeleteVector(v_);
    memorymanager::DeleteVector(initial_solution_);
  }

  std::string fixed_urdf_file_name_;
  pinocchio::Model fixed_model_;
  int kmax_;
  double baumgarte_alpha_, baumgarte_beta_, t_, finite_difference_;
  double *q_, *v_, *initial_solution_;
};


TEST_F(ZeroHorizonNewtonTest, dim) {
  Robot fixed_robot(fixed_urdf_file_name_);
  CostFunction_iiwa14 cost_function(&fixed_robot);
  ZeroHorizonOCP ocp(&fixed_robot, &cost_function);
  NewtonSystemForOCP<ZeroHorizonOCP> newton(&fixed_robot, &cost_function, 
                                            finite_difference_);
  FDGMRES<NewtonSystemForOCP<ZeroHorizonOCP>> fdgmres(newton.dim_solution(), 
                                                      kmax_);
  ZeroHorizonNewton zero_horizon_newton(&fixed_robot, &cost_function, 
                                        finite_difference_, kmax_);
  EXPECT_EQ(zero_horizon_newton.dimq(), newton.dimq());
  EXPECT_EQ(zero_horizon_newton.dimv(), newton.dimv());
  EXPECT_EQ(zero_horizon_newton.dim_solution(), newton.dim_solution());
}


TEST_F(ZeroHorizonNewtonTest, solveZeroHorizonOCP) {
  Robot fixed_robot(fixed_urdf_file_name_);
  CostFunction_iiwa14 cost_function(&fixed_robot);
  ZeroHorizonOCP ocp(&fixed_robot, &cost_function);
  NewtonSystemForOCP<ZeroHorizonOCP> newton(&fixed_robot, &cost_function, 
                                            finite_difference_);
  FDGMRES<NewtonSystemForOCP<ZeroHorizonOCP>> fdgmres(newton.dim_solution(), 
                                                      kmax_);
  ZeroHorizonNewton zero_horizon_newton(&fixed_robot, &cost_function, 
                                        finite_difference_, kmax_);
  double *solution = memorymanager::NewVector(zero_horizon_newton.dim_solution());
  zero_horizon_newton.solveZeroHorizonOCP(t_, q_, v_, solution);
  double *solution_ref = memorymanager::NewVector(newton.dim_solution());
  for (int i=0; i<newton.dim_solution(); ++i) {
    solution_ref[i] = 1.0;
  }
  double *solution_update = memorymanager::NewVector(newton.dim_solution());
  for (int i=0; i<newton.dim_solution(); ++i) {
    solution_update[i] = 1.0;
  }
  double error = newton.residualNorm(t_, q_, v_, solution_ref);
  int num_newton = 0;
  while (num_newton < 50 && error > 1.0e-03) {
    fdgmres.computeSolutionUpdate(newton, t_, q_, v_, solution_ref, 
                                  solution_update);
    newton.integrateSolution(solution_update, 1.0, solution_ref);
    std::cout << "error[" << num_newton << "] = " << error << std::endl;
    error = newton.residualNorm(t_, q_, v_, solution_ref);
    ++num_newton;
  }
  std::cout << "error[" << num_newton << "] = " << error << std::endl;
  EXPECT_TRUE(Eigen::Map<Eigen::VectorXd>(solution, zero_horizon_newton.dim_solution()).
              isApprox(Eigen::Map<Eigen::VectorXd>(solution_ref, newton.dim_solution())));
  memorymanager::DeleteVector(solution);
  memorymanager::DeleteVector(solution_ref);
  memorymanager::DeleteVector(solution_update);
}


TEST_F(ZeroHorizonNewtonTest, solveZeroHorizonOCPWithSettingParameters) {
  Robot fixed_robot(fixed_urdf_file_name_);
  CostFunction_iiwa14 cost_function(&fixed_robot);
  ZeroHorizonOCP ocp(&fixed_robot, &cost_function);
  NewtonSystemForOCP<ZeroHorizonOCP> newton(&fixed_robot, &cost_function, 
                                            finite_difference_);
  FDGMRES<NewtonSystemForOCP<ZeroHorizonOCP>> fdgmres(newton.dim_solution(), 
                                                      kmax_);
  ZeroHorizonNewton zero_horizon_newton(&fixed_robot, &cost_function, 
                                        finite_difference_, kmax_);
  double *solution = memorymanager::NewVector(zero_horizon_newton.dim_solution());
  zero_horizon_newton.setInitialGuessSolution(initial_solution_);
  zero_horizon_newton.setCriteriaForTermination(30, 1.0e-05);
  zero_horizon_newton.solveZeroHorizonOCP(t_, q_, v_, solution);
  double *solution_ref = memorymanager::NewVector(newton.dim_solution());
  for (int i=0; i<newton.dim_solution(); ++i) {
    solution_ref[i] = initial_solution_[i];
  }
  double *solution_update = memorymanager::NewVector(newton.dim_solution());
  for (int i=0; i<newton.dim_solution(); ++i) {
    solution_update[i] = initial_solution_[i];
  }
  double error = newton.residualNorm(t_, q_, v_, solution_ref);
  int num_newton = 0;
  while (num_newton < 30 && error > 1.0e-05) {
    fdgmres.computeSolutionUpdate(newton, t_, q_, v_, solution_ref, 
                                  solution_update);
    newton.integrateSolution(solution_update, 1.0, solution_ref);
    std::cout << "error[" << num_newton << "] = " << error << std::endl;
    error = newton.residualNorm(t_, q_, v_, solution_ref);
    ++num_newton;
  }
  std::cout << "error[" << num_newton << "] = " << error << std::endl;
  EXPECT_TRUE(Eigen::Map<Eigen::VectorXd>(solution, zero_horizon_newton.dim_solution()).
              isApprox(Eigen::Map<Eigen::VectorXd>(solution_ref, newton.dim_solution())));
  memorymanager::DeleteVector(solution);
  memorymanager::DeleteVector(solution_ref);
  memorymanager::DeleteVector(solution_update);
}


} // namespace robotcgmres


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}