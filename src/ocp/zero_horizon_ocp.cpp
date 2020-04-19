#include "ocp/zero_horizon_ocp.hpp"


namespace robotcgmres {

ZeroHorizonOCP::ZeroHorizonOCP(const Robot* robot_ptr,  
                               const CostFunctionInterface* cost_function,
                               const ConstraintsInterface* constraints) 
  : robot_(const_cast<Robot*>(robot_ptr)),
    cost_function_(const_cast<CostFunctionInterface*>(cost_function)),
    constraints_(const_cast<ConstraintsInterface*>(constraints)),
    dimq_(robot_ptr->dimq()),
    dimv_(robot_ptr->dimv()),
    dimf_(robot_ptr->dimf()),
    dim_passive_(robot_ptr->dim_passive()),
    dim_constraints_(constraints_->dim_constraints()),
    dim_solution_(robot_ptr->dimv()+2*robot_ptr->dimf()
                  +robot_ptr->dim_passive()+constraints_->dim_constraints()),
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
  // If there are contacts, the OCP with contacts is called.
  if (dimf_ > 0) {
    computeOptimalityResidualWithContacts(t, q, v, solution, 
                                          &(solution[dimv_]), 
                                          &(solution[dimv_+dimf_]), 
                                          &(solution[dimv_+2*dimf_]), 
                                          &(solution[dimv_+2*dimf_+dim_passive_]), 
                                          optimality_residual);
  } 
  // If there are no contacts, the OCP without contacts is called.
  else {
    computeOptimalityResidualWithoutContacts(t, q, v, solution, 
                                             &(solution[dimv_]), 
                                             &(solution[dimv_+dim_passive_]), 
                                             optimality_residual);
  }
}

void ZeroHorizonOCP::integrateSolution(const double* solution_update, 
                                       const double integration_length, 
                                       double* solution) {
  for (int i=0; i<dim_solution_; ++i) {
    solution[i] += integration_length * solution_update[i];
  }
}

void ZeroHorizonOCP::getControlInputFromSolution(
    const double* q, const double* v, const double* solution, 
    double* control_input_torques) {
  // If there are contacts, RNEA with contacts is called.
  if (dimf_ > 0) {
    robot_->RNEA(q, v, solution, &(solution[dimv_]), control_input_torques);
  } 
  // If there are no contacts, RNEA without contacts is called.
  else {
    robot_->RNEA(q, v, solution, control_input_torques);
  }
  // Sets the torques corresnponding to passive joints zero.
  robot_->setPassiveTorques(control_input_torques);
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

inline void ZeroHorizonOCP::computeOptimalityResidualWithoutContacts(
    const double t, const double* q, const double* v, const double* a, 
    const double* nu, const double* mu, double* optimality_residual) {
  // Compute the generalized torque for fully actuated robot u_
  robot_->RNEA(q, v, a, u_);
  // Compute hu and store it in beta_
  cost_function_->lu(robot_, t, q, v, a, u_, beta_);
  robot_->addVecToPassiveIndices(nu, beta_);
  constraints_->addCuDotVec(robot_, t, q, v, a, u_, mu, beta_);
  // Residuals with respect to the generalized acceleration
  robot_->RNEADerivativesTransDotVec(q, beta_, dRNEA_da_dot_beta_);
  cost_function_->la(robot_, t, q, v, a, u_, optimality_residual);
  constraints_->addCaDotVec(robot_, t, q, v, a, u_, mu, optimality_residual);
  cost_function_->phiv(robot_, t, q, v, gamma_);
  for (int i=0; i<dimv_; ++i) {
    optimality_residual[i] += gamma_[i] + dRNEA_da_dot_beta_[i];
  }
  // Compute the residual of the passive joints constraints.
  robot_->passiveTorqueViolation(u_, &(optimality_residual[dimv_]));
  // Compute the residual of the equality constraints.
  constraints_->residual(robot_, t, q, v, a, u_, 
                         &(optimality_residual[dimv_+dim_passive_]));
}

inline void ZeroHorizonOCP::computeOptimalityResidualWithContacts(
    const double t, const double* q, const double* v, const double* a, 
    const double* fext, const double* alpha, const double* nu, 
    const double* mu, double* optimality_residual) {
  // Compute the generalized torque u_
  robot_->setFext(fext);
  robot_->RNEA(q, v, a, fext, u_);
  // Compute hu and store in beta_
  cost_function_->lu(robot_, t, q, v, a, u_, fext, beta_);
  robot_->addVecToPassiveIndices(nu, beta_);
  constraints_->addCuDotVec(robot_, t, q, v, a, u_, fext, mu, beta_);
  // Residuals of acceleration, fext, and contact constraints.
  robot_->updateKinematics(q, v, a);
  robot_->RNEADerivativesTransDotVec(q, beta_, dRNEA_da_dot_beta_, 
                                     dRNEA_dfext_dot_beta_);
  cost_function_->la(robot_, t, q, v, a, u_, fext, optimality_residual);
  constraints_->addCaDotVec(robot_, t, q, v, a, u_, fext, mu, 
                            optimality_residual);
  robot_->baumgarteDerivativesDotVec(alpha, dBaum_dq_dot_alpha_, 
                                     dBaum_dv_dot_alpha_, dBaum_da_dot_alpha_);
  cost_function_->phiv(robot_, t, q, v, gamma_);
  for (int i=0; i<dimv_; ++i) {
    optimality_residual[i] += gamma_[i] + dRNEA_da_dot_beta_[i] 
                              + dBaum_da_dot_alpha_[i];
  }
  cost_function_->lf(robot_, t, q, v, a, u_, fext, 
                     &(optimality_residual[dimv_]));
  constraints_->addCfDotVec(robot_, t, q, v, a, u_, fext, mu,
                            &(optimality_residual[dimv_]));
  for (int i=0; i<dimf_; ++i) {
    optimality_residual[dimv_+i] += dRNEA_dfext_dot_beta_[i];
  }
  // Compute the residual of the contact.
  robot_->baumgarteResidual(&(optimality_residual[dimv_+dimf_]));
  // Compute the residual of the passive joints constraints.
  robot_->passiveTorqueViolation(u_, &(optimality_residual[dimv_+2*dimf_]));
  // Compute the residual of the equality constraints.
  constraints_->residual(robot_, t, q, v, a, u_, fext, 
                         &(optimality_residual[dimv_+2*dimf_+dim_passive_]));
}

} // namespace robotcgmres