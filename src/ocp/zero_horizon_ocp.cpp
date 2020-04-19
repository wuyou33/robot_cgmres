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
    hu_(memorymanager::NewVector(robot_ptr->dimv())),
    phiv_(memorymanager::NewVector(robot_ptr->dimv())),
    dRNEA_da_dot_hu_(memorymanager::NewVector(robot_ptr->dimv())),
    dRNEA_dfext_dot_hu_(memorymanager::NewVector(robot_ptr->dimf())),
    dBaum_dq_dot_mul_(memorymanager::NewVector(robot_ptr->dimv())),
    dBaum_dv_dot_mul_(memorymanager::NewVector(robot_ptr->dimv())),
    dBaum_da_dot_mul_(memorymanager::NewVector(robot_ptr->dimv())) {
}

ZeroHorizonOCP::~ZeroHorizonOCP() {
  memorymanager::DeleteVector(u_);
  memorymanager::DeleteVector(hu_);
  memorymanager::DeleteVector(phiv_);
  memorymanager::DeleteVector(dRNEA_da_dot_hu_);
  memorymanager::DeleteVector(dRNEA_dfext_dot_hu_);
  memorymanager::DeleteVector(dBaum_dq_dot_mul_);
  memorymanager::DeleteVector(dBaum_dv_dot_mul_);
  memorymanager::DeleteVector(dBaum_da_dot_mul_);
}

void ZeroHorizonOCP::computeOptimalityResidual(const double t, const double* q, 
                                               const double* v,  
                                               const double* solution, 
                                               double* optimality_residual) {
  // If there are contacts, the OCP with contacts is called.
  if (dimf_ > 0) {
    computeOptimalityResidualWithContacts(
        t, q, v, solution, &(solution[dimv_]), &(solution[dimv_+dimf_]), 
        &(solution[dimv_+2*dimf_]), &(solution[dimv_+2*dimf_+dim_passive_]), 
        optimality_residual, &(optimality_residual[dimv_]), 
        &(optimality_residual[dimv_+dimf_]), 
        &(optimality_residual[dimv_+2*dimf_]), 
        &(optimality_residual[dimv_+2*dimf_+dim_passive_]));
  } 
  // If there are no contacts, the OCP without contacts is called.
  else {
    computeOptimalityResidualWithoutContacts(
        t, q, v, solution, &(solution[dimv_]), &(solution[dimv_+dim_passive_]), 
        optimality_residual, &(optimality_residual[dimv_]), 
        &(optimality_residual[dimv_+dim_passive_]));
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
    const double* mul_passive, const double* mul_constraints, double* res_ha, 
    double* res_passive, double* res_constraints) {
  // Compute the generalized torque for fully actuated robot u_
  robot_->RNEA(q, v, a, u_);
  // Compute hu and store it in beta_
  cost_function_->lu(robot_, t, q, v, a, u_, hu_);
  robot_->addVecToPassiveIndices(mul_passive, hu_);
  constraints_->addCuDotVec(robot_, t, q, v, a, u_, mul_constraints, hu_);
  // Residuals with respect to the generalized acceleration
  robot_->RNEADerivativesTransDotVec(q, hu_, dRNEA_da_dot_hu_);
  cost_function_->la(robot_, t, q, v, a, u_, res_ha);
  constraints_->addCaDotVec(robot_, t, q, v, a, u_, mul_constraints, res_ha);
  cost_function_->phiv(robot_, t, q, v, phiv_);
  for (int i=0; i<dimv_; ++i) {
    res_ha[i] += phiv_[i] + dRNEA_da_dot_hu_[i];
  }
  // Compute the residual of the passive joints constraints.
  robot_->passiveTorqueViolation(u_, res_passive);
  // Compute the residual of the equality constraints.
  constraints_->residual(robot_, t, q, v, a, u_, res_constraints);
}

inline void ZeroHorizonOCP::computeOptimalityResidualWithContacts(
    const double t, const double* q, const double* v, const double* a, 
    const double* fext, const double* mul_contacts, const double* mul_passive, 
    const double* mul_constraints, double* res_ha, double* res_hfext, 
    double* res_contacts, double* res_passive, double* res_constraints) {
  // Compute the generalized torque u_
  robot_->setFext(fext);
  robot_->RNEA(q, v, a, fext, u_);
  // Compute hu and store in beta_
  cost_function_->lu(robot_, t, q, v, a, u_, fext, hu_);
  constraints_->addCuDotVec(robot_, t, q, v, a, u_, fext, mul_constraints, hu_);
  robot_->addVecToPassiveIndices(mul_passive, hu_);
  // Residuals of acceleration, fext, and contact constraints.
  robot_->updateKinematics(q, v, a);
  robot_->RNEADerivativesTransDotVec(q, hu_, dRNEA_da_dot_hu_, 
                                     dRNEA_dfext_dot_hu_);
  cost_function_->la(robot_, t, q, v, a, u_, fext, res_ha);
  constraints_->addCaDotVec(robot_, t, q, v, a, u_, fext, mul_constraints, 
                            res_ha);
  robot_->baumgarteDerivativesDotVec(mul_contacts, dBaum_dq_dot_mul_, 
                                     dBaum_dv_dot_mul_, dBaum_da_dot_mul_);
  cost_function_->phiv(robot_, t, q, v, phiv_);
  for (int i=0; i<dimv_; ++i) {
    res_ha[i] += phiv_[i] + dRNEA_da_dot_hu_[i] + dBaum_da_dot_mul_[i];
  }
  cost_function_->lf(robot_, t, q, v, a, u_, fext, res_hfext);
  constraints_->addCfDotVec(robot_, t, q, v, a, u_, fext, mul_constraints, 
                            res_hfext);
  for (int i=0; i<dimf_; ++i) {
    res_hfext[i] += dRNEA_dfext_dot_hu_[i];
  }
  // Compute the residual of the contact.
  robot_->baumgarteResidual(res_contacts);
  // Compute the residual of the passive joints constraints.
  robot_->passiveTorqueViolation(u_, res_passive);
  // Compute the residual of the equality constraints.
  constraints_->residual(robot_, t, q, v, a, u_, fext, res_constraints);
}

} // namespace robotcgmres