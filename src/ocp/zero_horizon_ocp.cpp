#include "ocp/zero_horizon_ocp.hpp"


namespace robotcgmres {

ZeroHorizonOCP::ZeroHorizonOCP(const Robot* robot_ptr) 
  : robot_(robot_ptr),
    passive_joint_constraint_(robot_ptr),
    cost_function_(robot_ptr),
    constraints_(robot_ptr),
    dimq_(robot_ptr->dimq()),
    dimv_(robot_ptr->dimv()),
    dimf_(robot_ptr->dimf()),
    dim_passive_(passive_joint_constraint_.dim_passive()),
    dim_constraints_(constraints_.dim_constraints()),
    dim_solution_(robot_ptr->dimv()+2*robot_ptr->dimf()
                  +assive_joint_constraint_.dim_passive()
                  +constraints_.dim_constraints()),
    u_(memorymanager::NewVector(robot_ptr->dimv())),
    beta_(memorymanager::NewVector(robot_ptr->dimv())),
    gamma_(memorymanager::NewVector(robot_ptr->dimv())),
    dRNEA_da_dot_beta_(memorymanager::NewVector(robot_ptr->dimv())),
    dRNEA_dfext_dot_beta_(memorymanager::NewVector(robot_ptr->dimf())),
    dBaum_dq_dot_alpha_(memorymanager::NewVector(robot_ptr->dimv())),
    dBaum_dv_dot_alpha_(memorymanager::NewVector(robot_ptr->dimv())),
    dBaum_da_dot_alpha_(memorymanager::NewVector(robot_ptr->dimv())) {
}

ZeroHorizonOCP::~ZeroHorizonOCP() {
  memorymanager::DeleteVector(u_);
  memorymanager::DeleteVector(beta_);
  memorymanager::DeleteVector(gamma_);
  memorymanager::DeleteVector(dRNEA_da_dot_beta_);
  memorymanager::DeleteVector(dRNEA_dfext_dot_beta_);
  memorymanager::DeleteVector(dBaum_dq_dot_alpha_);
  memorymanager::DeleteVector(dBaum_dv_dot_alpha_);
  memorymanager::DeleteVector(dBaum_da_dot_alpha_);
}

void ZeroHorizonOCP::computeOptimalityResidual(const double t, const double* q, 
                                               const double* v, 
                                               const double* solution,
                                               double* optimality_residual) {
  if (dimf_ > 0) {
    robot_->setFext(&(solution[dimv_]));
    robot_->RNEA(q, v, solution, &(solution[dimv_]), u_);
  } 
  else {
    robot_->RNEA(q, v, solution, u_);
  }

  cost_function_.lu(t, q, v, solution, u_, &(solution[dimv_]), beta_);
  constraints_.addCu(t, q, v, solution, u_, &(solution[dimv_]),
                     &(solution[dimv_+2*dimf_]), beta_);

  if (dimf_ > 0) {
    robot_->updateKinematics(q, v, solution);
    robot_->RNEADerivativesTransDotVector(q, beta_, dRNEA_da_dot_beta_, 
                                          dRNEA_dfext_dot_beta_);
    robot_->baumgarteResidual(q, v, solution, 
                              &(optimality_residual[dimv_+dimf_]));
    robot_->baumgarteDerivativesDotVec(q, v, solution, &(solution[dimv_+dimf_]), 
                                       dBaum_dq_dot_alpha_, dBaum_dv_dot_alpha_,
                                       dBaum_da_dot_alpha_);
    cost_function_.la(t, q, v, solution, u_, &(solution[dimv_]), 
                      optimality_residual);
    constraints_.addCa(t, q, v, solution, u_, &(solution[dimv_]), 
                      &(solution[dimv_+2*dimf_]), optimality_residual);
    cost_function_.phiv(t, q, v, gamma_);
    for (int i=0; i<dimv_; ++i) {
      optimality_residual[i] += gamma_[i] + dRNEA_da_dot_beta_[i] 
                                + dBaum_da_dot_alpha_[i];
  } 
  else {
    robot_->RNEADerivativeTransDotVector(q, beta_, dRNEA_da_dot_beta_);
    cost_function_.la(t, q, v, solution, u_, &(solution[dimv_]), 
                      optimality_residual);
    constraints_.addCa(t, q, v, solution, u_, &(solution[dimv_]), 
                      &(solution[dimv_+2*dimf_]), optimality_residual);
    cost_function_.phiv(t, q, v, gamma_);
    for (int i=0; i<dimv_; ++i) {
      optimality_residual[i] += gamma_[i] + dRNEA_da_dot_beta_[i];
    }
  }

  computeLfCf(time, state_vec, &(state_vec[dim_configuration_]), solution_vec, 
              &(solution_vec[dim_velocity_]), u_vec_,
              &(solution_vec[dim_velocity_+2*dim_fext_]), 
              &(optimality_residual[dim_velocity_]));
  for (int i=0; i<dim_fext_; ++i) {
    optimality_residual[dim_velocity_+i] += dRNEA_dfext_dot_beta_[i];
  }

  constraints_.C(t, q, v, solution, u_, &(solution[dimv_]), 
                 &(solution[dimv_+2*dimf_]));
}

} // namespace robotcgmres