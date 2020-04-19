#include <string>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "Eigen/Core"

#include "robot/robot.hpp"
#include "ocp/zero_horizon_ocp.hpp"
#include "cost_function_iiwa14.hpp"
#include "constraints_iiwa14.hpp"
#include "common/memory_manager.hpp"
#include "common/linear_algebra.hpp"



namespace robotcgmres {

class ZeroHorizonOCPTest: public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    fixed_urdf_file_name_ = "../../../examples/iiwa14/iiwa14.urdf";
    pinocchio::urdf::buildModel(fixed_urdf_file_name_, fixed_model_);
    q_ = memorymanager::NewVector(fixed_model_.nq);
    v_ = memorymanager::NewVector(fixed_model_.nv);
    a_ = memorymanager::NewVector(fixed_model_.nv);
    u_ = memorymanager::NewVector(fixed_model_.nv);
    beta_ = memorymanager::NewVector(fixed_model_.nv);
    gamma_ = memorymanager::NewVector(fixed_model_.nv);
    dRNEA_da_dot_beta_ = memorymanager::NewVector(fixed_model_.nv);
    dRNEA_dfext_dot_beta_ = memorymanager::NewVector(fixed_model_.nv);
    dBaum_dq_dot_alpha_ = memorymanager::NewVector(fixed_model_.nv);
    dBaum_dv_dot_alpha_ = memorymanager::NewVector(fixed_model_.nv);
    dBaum_da_dot_alpha_ = memorymanager::NewVector(fixed_model_.nv);
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
  }

  virtual void TearDown() {
    memorymanager::DeleteVector(q_);
    memorymanager::DeleteVector(v_);
    memorymanager::DeleteVector(a_);
    memorymanager::DeleteVector(u_);
    memorymanager::DeleteVector(beta_);
    memorymanager::DeleteVector(gamma_);
    memorymanager::DeleteVector(dRNEA_da_dot_beta_);
    memorymanager::DeleteVector(dRNEA_dfext_dot_beta_);
    memorymanager::DeleteVector(dBaum_dq_dot_alpha_);
    memorymanager::DeleteVector(dBaum_dv_dot_alpha_);
    memorymanager::DeleteVector(dBaum_da_dot_alpha_);
  }

  std::string fixed_urdf_file_name_;
  pinocchio::Model fixed_model_;
  double baumgarte_alpha_, baumgarte_beta_, t_;
  double *q_, *v_, *a_, *u_, *beta_, *gamma_, *dRNEA_da_dot_beta_, 
         *dRNEA_dfext_dot_beta_, *dBaum_dq_dot_alpha_, *dBaum_dv_dot_alpha_, 
         *dBaum_da_dot_alpha_;
};


TEST_F(ZeroHorizonOCPTest, dims) {
  Robot fixed_robot(fixed_urdf_file_name_);
  CostFunction_iiwa14 cost_function(&fixed_robot);
  Constraints_iiwa14 constraints(&fixed_robot);
  ZeroHorizonOCP ocp_(&fixed_robot, &cost_function, &constraints);
  EXPECT_EQ(ocp_.dimq(), fixed_robot.dimq());
  EXPECT_EQ(ocp_.dimv(), fixed_robot.dimv());
  EXPECT_EQ(ocp_.dimf(), fixed_robot.dimf());
  EXPECT_EQ(ocp_.dim_passive(), fixed_robot.dim_passive());
  EXPECT_EQ(ocp_.dim_constraints(), constraints.dim_constraints());
  const int dim_solution = fixed_robot.dimv()+2*fixed_robot.dimf()+fixed_robot.dim_passive()+constraints.dim_constraints();
  EXPECT_EQ(ocp_.dim_solution(), dim_solution);
}


TEST_F(ZeroHorizonOCPTest, computeOptimalityResidualWithoutFext) {
  Robot fixed_robot(fixed_urdf_file_name_);
  CostFunction_iiwa14 cost_function(&fixed_robot);
  Constraints_iiwa14 constraints(&fixed_robot);
  ZeroHorizonOCP ocp(&fixed_robot, &cost_function, &constraints);
  const int dim_solution = ocp.dim_solution();
  EXPECT_EQ(ocp.dimf(), 0);
  EXPECT_EQ(dim_solution, ocp.dimv()+ocp.dim_passive()+ocp.dim_constraints());
  double *solution = memorymanager::NewVector(dim_solution);
  Eigen::Map<Eigen::VectorXd>(solution, dim_solution) = Eigen::VectorXd::Random(dim_solution);
  double* residual_ref = memorymanager::NewVector(ocp.dim_solution());
  double* residual = memorymanager::NewVector(ocp.dim_solution());
  ocp.computeOptimalityResidual(t_, q_, v_, solution, residual);

  double *a = solution;
  double *nu = &(solution[fixed_robot.dimv()]);
  double *mu = &(solution[fixed_robot.dimv()+fixed_robot.dim_passive()]);
  // Compute the generalized torque for fully actuated robot u_
  fixed_robot.RNEA(q_, v_, a, u_);
  // Compute hu and store it in beta_
  cost_function.lu(&fixed_robot, t_, q_, v_, a, u_, beta_);
  fixed_robot.addVecToPassiveIndices(nu, beta_);
  constraints.addCuDotVec(&fixed_robot, t_, q_, v_, a, u_, mu, beta_);
  // Residuals with respect to the generalized acceleration
  fixed_robot.RNEADerivativesTransDotVec(q_, beta_, dRNEA_da_dot_beta_);
  cost_function.la(&fixed_robot, t_, q_, v_, a, u_, residual_ref);
  constraints.addCaDotVec(&fixed_robot, t_, q_, v_, a, u_, mu, residual_ref);
  cost_function.phiv(&fixed_robot, t_, q_, v_, gamma_);
  for (int i=0; i<fixed_robot.dimv(); ++i) {
    residual_ref[i] += gamma_[i] + dRNEA_da_dot_beta_[i];
  }
  // Compute the residual of the passive joints constraints.
  fixed_robot.passiveTorqueViolation(u_, &(residual_ref[fixed_robot.dimv()]));
  // Compute the residual of the equality constraints.
  constraints.residual(&fixed_robot, t_, q_, v_, a, u_, 
                       &(residual_ref[fixed_robot.dimv()+fixed_robot.dim_passive()]));
  EXPECT_TRUE(Eigen::Map<Eigen::VectorXd>(residual, ocp.dim_solution()).
              isApprox(Eigen::Map<Eigen::VectorXd>(residual_ref, ocp.dim_solution())));

  memorymanager::DeleteVector(solution);
  memorymanager::DeleteVector(residual_ref);
  memorymanager::DeleteVector(residual);
}


TEST_F(ZeroHorizonOCPTest, computeOptimalityResidualWithFext) {
  Robot fixed_robot(fixed_urdf_file_name_);
  fixed_robot.addPointContact(18, baumgarte_alpha_, baumgarte_beta_);
  CostFunction_iiwa14 cost_function(&fixed_robot);
  Constraints_iiwa14 constraints(&fixed_robot);
  ZeroHorizonOCP ocp(&fixed_robot, &cost_function, &constraints);
  const int dim_solution = ocp.dim_solution();
  EXPECT_EQ(ocp.dimf(), 3);
  EXPECT_EQ(dim_solution, 
            ocp.dimv()+2*ocp.dimf()+ocp.dim_passive()+ocp.dim_constraints());
  double *solution = memorymanager::NewVector(dim_solution);
  Eigen::Map<Eigen::VectorXd>(solution, dim_solution) = Eigen::VectorXd::Random(dim_solution);
  double* residual = memorymanager::NewVector(ocp.dim_solution());
  ocp.computeOptimalityResidual(t_, q_, v_, solution, residual);

  // Compute ref residual
  double* residual_ref = memorymanager::NewVector(ocp.dim_solution());
  double* a = solution;
  double* fext = &(solution[fixed_robot.dimv()]);
  double* alpha = &(solution[fixed_robot.dimv()+fixed_robot.dimf()]); 
  double* nu = &(solution[fixed_robot.dimv()+2*fixed_robot.dimf()]); 
  double* mu = &(solution[fixed_robot.dimv()+2*fixed_robot.dimf()+fixed_robot.dim_passive()]); 
  // Compute the generalized torque u_
  fixed_robot.setFext(fext);
  fixed_robot.RNEA(q_, v_, a, fext, u_);
  // Compute hu and store in beta_
  cost_function.lu(&fixed_robot, t_, q_, v_, a, u_, fext, beta_);
  fixed_robot.addVecToPassiveIndices(nu, beta_);
  constraints.addCuDotVec(&fixed_robot, t_, q_, v_, a, u_, fext, mu, beta_);
  // Residuals of acceleration, fext, and contact constraints.
  fixed_robot.updateKinematics(q_, v_, a);
  fixed_robot.RNEADerivativesTransDotVec(q_, beta_, dRNEA_da_dot_beta_, 
                                     dRNEA_dfext_dot_beta_);
  cost_function.la(&fixed_robot, t_, q_, v_, a, u_, fext, residual_ref);
  constraints.addCaDotVec(&fixed_robot, t_, q_, v_, a, u_, fext, mu, 
                           residual_ref);
  fixed_robot.baumgarteDerivativesDotVec(alpha, dBaum_dq_dot_alpha_, 
                                         dBaum_dv_dot_alpha_, 
                                         dBaum_da_dot_alpha_);
  cost_function.phiv(&fixed_robot, t_, q_, v_, gamma_);
  for (int i=0; i<fixed_robot.dimv(); ++i) {
    residual_ref[i] += gamma_[i] + dRNEA_da_dot_beta_[i] 
                              + dBaum_da_dot_alpha_[i];
  }
  cost_function.lf(&fixed_robot, t_, q_, v_, a, u_, fext, 
                     &(residual_ref[fixed_robot.dimv()]));
  constraints.addCfDotVec(&fixed_robot, t_, q_, v_, a, u_, fext, mu,
                            &(residual_ref[fixed_robot.dimv()]));
  for (int i=0; i<fixed_robot.dimf(); ++i) {
    residual_ref[fixed_robot.dimv()+i] += dRNEA_dfext_dot_beta_[i];
  }
  // Compute the residual of the contact.
  fixed_robot.baumgarteResidual(&(residual_ref[fixed_robot.dimv()+fixed_robot.dimf()]));
  // Compute the residual of the passive joints constraints.
  fixed_robot.passiveTorqueViolation(u_, 
                                     &(residual_ref[fixed_robot.dimv()+2*fixed_robot.dimf()]));
  // Compute the residual of the equality constraints.
  constraints.residual(&fixed_robot, t_, q_, v_, a, u_, fext, 
                        &(residual_ref[fixed_robot.dimv()+2*fixed_robot.dimf()+fixed_robot.dim_passive()]));
  EXPECT_TRUE(Eigen::Map<Eigen::VectorXd>(residual, ocp.dim_solution()).
              isApprox(Eigen::Map<Eigen::VectorXd>(residual_ref, ocp.dim_solution())));

  memorymanager::DeleteVector(solution);
  memorymanager::DeleteVector(residual);
  memorymanager::DeleteVector(residual_ref);
}


TEST_F(ZeroHorizonOCPTest, integrateSolution) {
  Robot fixed_robot(fixed_urdf_file_name_);
  CostFunction_iiwa14 cost_function(&fixed_robot);
  Constraints_iiwa14 constraints(&fixed_robot);
  ZeroHorizonOCP ocp(&fixed_robot, &cost_function, &constraints);
  const int dim_solution = ocp.dim_solution();
  EXPECT_EQ(ocp.dimf(), 0);
  EXPECT_EQ(dim_solution, ocp.dimv()+ocp.dim_passive()+ocp.dim_constraints());
  double *solution = memorymanager::NewVector(dim_solution);
  double *solution_ref = memorymanager::NewVector(dim_solution);
  double *solution_update = memorymanager::NewVector(dim_solution);
  const double integration_length = Eigen::VectorXd::Random(2)[1];
  Eigen::Map<Eigen::VectorXd>(solution, dim_solution) = Eigen::VectorXd::Random(dim_solution);
  Eigen::Map<Eigen::VectorXd>(solution_ref, dim_solution) = Eigen::Map<Eigen::VectorXd>(solution, dim_solution);
  Eigen::Map<Eigen::VectorXd>(solution_update, dim_solution) = Eigen::VectorXd::Random(dim_solution);

  ocp.integrateSolution(solution_update, integration_length, solution);
  for (int i=0; i<ocp.dim_solution(); ++i) {
    solution_ref[i] += integration_length * solution_update[i];
  }

  EXPECT_TRUE(Eigen::Map<Eigen::VectorXd>(solution, ocp.dim_solution()).
              isApprox(Eigen::Map<Eigen::VectorXd>(solution_ref, ocp.dim_solution())));
  memorymanager::DeleteVector(solution);
  memorymanager::DeleteVector(solution_ref);
  memorymanager::DeleteVector(solution_update);
}


TEST_F(ZeroHorizonOCPTest, getControlInputWithoutFext) {
  Robot fixed_robot(fixed_urdf_file_name_);
  CostFunction_iiwa14 cost_function(&fixed_robot);
  Constraints_iiwa14 constraints(&fixed_robot);
  ZeroHorizonOCP ocp(&fixed_robot, &cost_function, &constraints);
  const int dim_solution = ocp.dim_solution();
  EXPECT_EQ(ocp.dimf(), 0);
  EXPECT_EQ(dim_solution, ocp.dimv()+ocp.dim_passive()+ocp.dim_constraints());
  double *solution = memorymanager::NewVector(dim_solution);
  Eigen::Map<Eigen::VectorXd>(solution, dim_solution) = Eigen::VectorXd::Random(dim_solution);
  double *control_input = memorymanager::NewVector(ocp.dimv());
  double *control_input_ref = memorymanager::NewVector(ocp.dimv());
  ocp.getControlInputFromSolution(q_, v_, solution, control_input);
  fixed_robot.RNEA(q_, v_, solution, control_input_ref);
  fixed_robot.setPassiveTorques(control_input_ref);
  EXPECT_TRUE(Eigen::Map<Eigen::VectorXd>(control_input, ocp.dimv()).
              isApprox(Eigen::Map<Eigen::VectorXd>(control_input_ref, ocp.dimv())));
  memorymanager::DeleteVector(solution);
  memorymanager::DeleteVector(control_input);
  memorymanager::DeleteVector(control_input_ref);
}


TEST_F(ZeroHorizonOCPTest, getControlInputWithFext) {
  Robot fixed_robot(fixed_urdf_file_name_);
  fixed_robot.addPointContact(18, baumgarte_alpha_, baumgarte_beta_);
  CostFunction_iiwa14 cost_function(&fixed_robot);
  Constraints_iiwa14 constraints(&fixed_robot);
  ZeroHorizonOCP ocp(&fixed_robot, &cost_function, &constraints);
  const int dim_solution = ocp.dim_solution();
  EXPECT_EQ(ocp.dimf(), 3);
  double *solution = memorymanager::NewVector(dim_solution);
  Eigen::Map<Eigen::VectorXd>(solution, dim_solution) = Eigen::VectorXd::Random(dim_solution);
  double *control_input = memorymanager::NewVector(ocp.dimv());
  double *control_input_ref = memorymanager::NewVector(ocp.dimv());
  ocp.getControlInputFromSolution(q_, v_, solution, control_input);
  fixed_robot.RNEA(q_, v_, solution, &(solution[fixed_robot.dimv()]), control_input_ref);
  fixed_robot.setPassiveTorques(control_input_ref);
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