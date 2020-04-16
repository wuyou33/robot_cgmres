#include "ocp/zero_horizon_ocp.hpp"


namespace robotcgmres {

ZeroHorizonOCP::ZeroHorizonOCP(const Robot* robot_ptr) 
  : robot_(robot_ptr),
    passive_joints_(robot_ptr),
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
                                               const double* v, const double* a, 
                                               const double* mutipliers, 
                                               double* optimality_residual) {
  // Compute the generalized torque for fully actuated robot u_
  robot_->RNEA(q, v, a, u_);
  // Compute hu and store it in beta_
  cost_function_.lu(t, q, v, a, u_, beta_);
  passive_joints_.addPassiveConstraintsDerivativeDotVec(mutipliers, beta_);
  constraints_.addCuDotVec(t, q, v, a, u_, &(mutipliers[dim_passive_]), beta_);
  // Residuals with respect to the generalized acceleration
  robot_->RNEADerivativesTransDotVec(q, beta_, dRNEA_da_dot_beta_);
  cost_function_.la(t, q, v, a, u_, optimality_residual);
  constraints_.addCaDotVec(t, q, v, a, u_, &(mutipliers[dim_passive_]), 
                           optimality_residual);
  cost_function_.phiv(t, q, v, gamma_);
  for (int i=0; i<dimv_; ++i) {
    optimality_residual[i] += gamma_[i] + dRNEA_da_dot_beta_[i];
  }
  // Compute the residual of the passive joints constraints.
  passive_joints_.passiveConstraintsResidual(u_, &(optimality_residual[dimv_]));
  // Compute the residual of the equality constraints.
  constraints_.C(t, q, v, a, u_, &(optimality_residual[dimv_+dim_passive_]));
}

void ZeroHorizonOCP::computeOptimalityResidual(const double t, const double* q, 
                                               const double* v, const double* a, 
                                               const double* fext, 
                                               const double* mutipliers, 
                                               double* optimality_residual) {
  // Compute the generalized torque u_
  robot_->setFext(fext);
  robot_->RNEA(q, v, a, fext, u_);
  // Compute hu and store in beta_
  cost_function_.lu(t, q, v, a, u_, fext, beta_);
  passive_joints_.addPassiveConstraintsDerivativeDotVec(mutipliers, beta_);
  constraints_.addCuDotVec(t, q, v, a, u_, fext, 
                           &(mutipliers[dim_passive_+dimf_]), beta_);
  // Residuals of acceleration, fext, and contact constraints.
  robot_->updateKinematics(q, v, a);
  robot_->RNEADerivativesTransDotVec(q, beta_, dRNEA_da_dot_beta_, 
                                     dRNEA_dfext_dot_beta_);
  cost_function_.la(t, q, v, a, u_, fext, optimality_residual);
  constraints_.addCaDotVec(t, q, v, a, u_, fext, 
                           &(mutipliers[dim_passive_+dimf_]), 
                           optimality_residual);
  robot_->baumgarteDerivativesDotVec(&(mutipliers[dim_passive_]),
                                     dBaum_dq_dot_alpha_, dBaum_dv_dot_alpha_,
                                     dBaum_da_dot_alpha_);
  cost_function_.phiv(t, q, v, gamma_);
  for (int i=0; i<dimv_; ++i) {
    optimality_residual[i] += gamma_[i] + dRNEA_da_dot_beta_[i] 
                              + dBaum_da_dot_alpha_[i];
  }
  cost_function_.lf(t, q, v, a, u_, fext, &(optimality_residual[dimv_]));
  constraints_.addCfDotVec(t, q, v, a, u_, fext,  
                           &(mutipliers[dim_passive_+dimf_])
                           &(optimality_residual[dimv_]));
  for (int i=0; i<dimv_; ++i) {
    optimality_residual[dimv_+i] += dRNEA_dfext_dot_beta_[i];
  }
  // Compute the residual of the passive joints constraints.
  passive_joints_.passiveConstraintsResidual(
      u_, &(optimality_residual[dimv_+dimf_]));
  // Compute the residual of the contact.
  robot_->baumgarteResidual(&(optimality_residual[dimv_+dimf_+dim_passive_]));
  // Compute the residual of the equality constraints.
  constraints_.C(t, q, v, a, u_, fext, 
                 &(optimality_residual[dimv_+2*dimf_+dim_passive_]));
}

void ZeroHorizonOCP::integrateSolution(const double* solution_update, 
                                       const double integration_length, 
                                       double* solution) {
  for (int i=0; i<dim_solution_; ++i) {
    solution[i] += integration_length * solution_update[i];
  }
}

int ZeroHorizonOCP::dimq() const {
  return dimq_;
}

int ZeroHorizonOCP::dimv() const {
  return dimv_;
}

int ZeroHorizonOCP::dimf() const {
  return dimf_;
}

int ZeroHorizonOCP::dim_passive() const {
  return dim_passive_;
}

int ZeroHorizonOCP::dim_constraints() const {
  return dim_constraints_;
}

int ZeroHorizonOCP::dim_solution() const {
  return dim_solution_;
}

} // namespace robotcgmres