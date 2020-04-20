#include <string>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robot/robot.hpp"
#include "solver/newton_system_for_ocp.hpp"
#include "zero_horizon_ocp/zero_horizon_ocp.hpp"
#include "common/memory_manager.hpp"
#include "common/linear_algebra.hpp"

#include "iiwa14/cost_function_iiwa14.hpp"


namespace robotcgmres {

class ZeroHorizonNewtonSystemTest: public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_urdf_file_name_ = "../../../examples/iiwa14/iiwa14.urdf";
    pinocchio::urdf::buildModel(fixed_urdf_file_name_, fixed_model_);
    q_ = memorymanager::NewVector(fixed_model_.nq);
    v_ = memorymanager::NewVector(fixed_model_.nv);
    a_ = memorymanager::NewVector(fixed_model_.nv);
    u_ = memorymanager::NewVector(fixed_model_.nv);
    Eigen::Map<Eigen::VectorXd>(q_, fixed_model_.nq) 
        = pinocchio::randomConfiguration(fixed_model_, 
                                         -Eigen::VectorXd::Ones(fixed_model_.nq), 
                                         Eigen::VectorXd::Ones(fixed_model_.nq));
    Eigen::Map<Eigen::VectorXd>(v_, fixed_model_.nv) 
        = Eigen::VectorXd::Random(fixed_model_.nv);
    Eigen::Map<Eigen::VectorXd>(a_, fixed_model_.nv) 
        = Eigen::VectorXd::Random(fixed_model_.nv);
    t_ = Eigen::VectorXd::Random(2)[1];
    baumgarte_alpha_ = Eigen::VectorXd::Random(2)[1];
    baumgarte_beta_ = Eigen::VectorXd::Random(2)[1];
    finite_difference_ = 1.0e-06;
  }

  virtual void TearDown() {
    memorymanager::DeleteVector(q_);
    memorymanager::DeleteVector(v_);
    memorymanager::DeleteVector(a_);
    memorymanager::DeleteVector(u_);
  }

  std::string fixed_urdf_file_name_;
  pinocchio::Model fixed_model_;
  double baumgarte_alpha_, baumgarte_beta_, t_, finite_difference_;
  double *q_, *v_, *a_, *u_, *beta_, *gamma_, *dRNEA_da_dot_beta_, 
         *dRNEA_dfext_dot_beta_, *dBaum_dq_dot_alpha_, *dBaum_dv_dot_alpha_, 
         *dBaum_da_dot_alpha_;
};


TEST_F(ZeroHorizonNewtonSystemTest, dim_solution) {
  Robot fixed_robot(fixed_urdf_file_name_);
  CostFunction_iiwa14 cost_function(&fixed_robot);
  ZeroHorizonOCP ocp(&fixed_robot, &cost_function);
  NewtonSystemForOCP<ZeroHorizonOCP> newton(&fixed_robot, &cost_function, 
                                            finite_difference_);
  EXPECT_EQ(newton.dimq(), ocp.dimq());
  EXPECT_EQ(newton.dimv(), ocp.dimv());
  EXPECT_EQ(newton.dim_solution(), ocp.dim_solution());
}


TEST_F(ZeroHorizonNewtonSystemTest, computeInitialResidual) {
  Robot fixed_robot(fixed_urdf_file_name_);
  CostFunction_iiwa14 cost_function(&fixed_robot);
  ZeroHorizonOCP ocp(&fixed_robot, &cost_function);
  NewtonSystemForOCP<ZeroHorizonOCP> newton(&fixed_robot, &cost_function, 
                                            finite_difference_);
  const int dim_solution = ocp.dim_solution();
  double *solution = memorymanager::NewVector(dim_solution);
  double *direction = memorymanager::NewVector(dim_solution);
  double *residual = memorymanager::NewVector(dim_solution);
  Eigen::Map<Eigen::VectorXd>(solution, dim_solution) = Eigen::VectorXd::Random(dim_solution);
  Eigen::Map<Eigen::VectorXd>(direction, dim_solution) = Eigen::VectorXd::Random(dim_solution);
  newton.computeInitialResidual(t_, q_, v_, solution, direction, residual);
  double *residual_ref = memorymanager::NewVector(dim_solution);
  double *incremented_solution = memorymanager::NewVector(dim_solution);
  double *opt_res = memorymanager::NewVector(dim_solution);
  double *incremented_opt_res = memorymanager::NewVector(dim_solution);
  for (int i=0; i<ocp.dim_solution(); ++i) {
    incremented_solution[i] = solution[i] + finite_difference_ * direction[i];
  }
  ocp.computeOptimalityResidual(t_, q_, v_, solution, opt_res);
  ocp.computeOptimalityResidual(t_, q_, v_, incremented_solution, incremented_opt_res);
  for (int i=0; i<ocp.dim_solution(); ++i) {
    residual_ref[i] = - opt_res[i] 
        - (incremented_opt_res[i]-opt_res[i]) / finite_difference_;
  }
  EXPECT_TRUE(Eigen::Map<Eigen::VectorXd>(residual, newton.dim_solution()).
              isApprox(Eigen::Map<Eigen::VectorXd>(residual_ref, ocp.dim_solution())));
  memorymanager::DeleteVector(solution);
  memorymanager::DeleteVector(direction);
  memorymanager::DeleteVector(residual);
  memorymanager::DeleteVector(residual_ref);
  memorymanager::DeleteVector(incremented_solution);
  memorymanager::DeleteVector(opt_res);
  memorymanager::DeleteVector(incremented_opt_res);
}


TEST_F(ZeroHorizonNewtonSystemTest, computeAx) {
  Robot fixed_robot(fixed_urdf_file_name_);
  CostFunction_iiwa14 cost_function(&fixed_robot);
  ZeroHorizonOCP ocp(&fixed_robot, &cost_function);
  NewtonSystemForOCP<ZeroHorizonOCP> newton(&fixed_robot, &cost_function, 
                                            finite_difference_);
  const int dim_solution = ocp.dim_solution();
  double *solution = memorymanager::NewVector(dim_solution);
  double *direction = memorymanager::NewVector(dim_solution);
  double *ax = memorymanager::NewVector(dim_solution);
  Eigen::Map<Eigen::VectorXd>(solution, dim_solution) = Eigen::VectorXd::Random(dim_solution);
  Eigen::Map<Eigen::VectorXd>(direction, dim_solution) = Eigen::VectorXd::Random(dim_solution);
  newton.computeInitialResidual(t_, q_, v_, solution, direction, ax);
  newton.computeAx(t_, q_, v_, solution, direction, ax);
  double *ax_ref = memorymanager::NewVector(dim_solution);
  double *incremented_solution = memorymanager::NewVector(dim_solution);
  double *opt_res = memorymanager::NewVector(dim_solution);
  double *incremented_opt_res = memorymanager::NewVector(dim_solution);
  for (int i=0; i<ocp.dim_solution(); ++i) {
    incremented_solution[i] = solution[i] + finite_difference_ * direction[i];
  }
  ocp.computeOptimalityResidual(t_, q_, v_, solution, opt_res);
  ocp.computeOptimalityResidual(t_, q_, v_, incremented_solution, incremented_opt_res);
  for (int i=0; i<ocp.dim_solution(); ++i) {
    ax_ref[i] = (incremented_opt_res[i]-opt_res[i]) / finite_difference_;
  }
  EXPECT_TRUE(Eigen::Map<Eigen::VectorXd>(ax, newton.dim_solution()).
              isApprox(Eigen::Map<Eigen::VectorXd>(ax_ref, ocp.dim_solution())));
  memorymanager::DeleteVector(solution);
  memorymanager::DeleteVector(direction);
  memorymanager::DeleteVector(ax);
  memorymanager::DeleteVector(ax_ref);
  memorymanager::DeleteVector(incremented_solution);
  memorymanager::DeleteVector(opt_res);
  memorymanager::DeleteVector(incremented_opt_res);
}


TEST_F(ZeroHorizonNewtonSystemTest, residualNorm) {
  Robot fixed_robot(fixed_urdf_file_name_);
  CostFunction_iiwa14 cost_function(&fixed_robot);
  ZeroHorizonOCP ocp(&fixed_robot, &cost_function);
  NewtonSystemForOCP<ZeroHorizonOCP> newton(&fixed_robot, &cost_function, 
                                            finite_difference_);
  const int dim_solution = ocp.dim_solution();
  double *solution = memorymanager::NewVector(dim_solution);
  Eigen::Map<Eigen::VectorXd>(solution, dim_solution) = Eigen::VectorXd::Random(dim_solution);
  const double residual_norm = newton.residualNorm(t_, q_, v_, solution);
  double *opt_res = memorymanager::NewVector(dim_solution);
  ocp.computeOptimalityResidual(t_, q_, v_, solution, opt_res);
  const double residual_norm_ref 
      = std::sqrt(linearalgebra::SquaredNorm(ocp.dim_solution(), opt_res));
  EXPECT_DOUBLE_EQ(residual_norm, residual_norm_ref);
  memorymanager::DeleteVector(solution);
  memorymanager::DeleteVector(opt_res);
}


TEST_F(ZeroHorizonNewtonSystemTest, integrateSolution) {
  Robot fixed_robot(fixed_urdf_file_name_);
  CostFunction_iiwa14 cost_function(&fixed_robot);
  ZeroHorizonOCP ocp(&fixed_robot, &cost_function);
  NewtonSystemForOCP<ZeroHorizonOCP> newton(&fixed_robot, &cost_function, 
                                            finite_difference_);
  const int dim_solution = ocp.dim_solution();
  double *solution = memorymanager::NewVector(dim_solution);
  double *solution_ref = memorymanager::NewVector(dim_solution);
  double *solution_update = memorymanager::NewVector(dim_solution);
  const double integration_length =  Eigen::VectorXd::Random(2)[1];
  Eigen::Map<Eigen::VectorXd>(solution, dim_solution) = Eigen::VectorXd::Random(dim_solution);
  Eigen::Map<Eigen::VectorXd>(solution_ref, dim_solution) = Eigen::Map<Eigen::VectorXd>(solution, dim_solution);
  Eigen::Map<Eigen::VectorXd>(solution_update, dim_solution) = Eigen::VectorXd::Random(dim_solution);
  newton.integrateSolution(solution_update, integration_length, solution);
  ocp.integrateSolution(solution_update, integration_length, solution_ref);
  EXPECT_TRUE(Eigen::Map<Eigen::VectorXd>(solution, newton.dim_solution()).
              isApprox(Eigen::Map<Eigen::VectorXd>(solution_ref, ocp.dim_solution())));
  memorymanager::DeleteVector(solution);
  memorymanager::DeleteVector(solution_ref);
  memorymanager::DeleteVector(solution_update);
}


TEST_F(ZeroHorizonNewtonSystemTest, getControlInput) {
  Robot fixed_robot(fixed_urdf_file_name_);
  CostFunction_iiwa14 cost_function(&fixed_robot);
  ZeroHorizonOCP ocp(&fixed_robot, &cost_function);
  NewtonSystemForOCP<ZeroHorizonOCP> newton(&fixed_robot, &cost_function, 
                                            finite_difference_);
  const int dim_solution = ocp.dim_solution();
  double *solution = memorymanager::NewVector(dim_solution);
  double *control_input = memorymanager::NewVector(ocp.dimv());
  double *control_input_ref = memorymanager::NewVector(ocp.dimv());
  const double integration_length =  Eigen::VectorXd::Random(2)[1];
  Eigen::Map<Eigen::VectorXd>(solution, dim_solution) = Eigen::VectorXd::Random(dim_solution);

  newton.getControlInputFromSolution(q_, v_, solution, control_input);
  ocp.getControlInputFromSolution(q_, v_, solution, control_input_ref);
  EXPECT_TRUE(Eigen::Map<Eigen::VectorXd>(control_input, ocp.dimv()).
              isApprox(Eigen::Map<Eigen::VectorXd>(control_input_ref, ocp.dimv())));
  memorymanager::DeleteVector(solution);
  memorymanager::DeleteVector(control_input);
  memorymanager::DeleteVector(control_input_ref);
}

} // namespace robotcgmres


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}