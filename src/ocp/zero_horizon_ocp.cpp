#include "ocp/zero_horizon_ocp.hpp"


namespace robotcgmres {

ZeroHorizonOCP::ZeroHorizonOCP(Robot* robot_ptr) 
  : robot_(robot_ptr),
    cost_function_(),
    constraints_(),
    dimq_(robot_ptr->dimq()),
    dimv_(robot_ptr->dimv()),
    dimtau_(robot_ptr->dimtau()),
    dimf_(robot_ptr->dimf()),
    dimC_(constraints_.num_constraints()),
    dim_solution_(robot_ptr->dimv()+robot_ptr->dimf() +robot_ptr->dimf()
                  +constraints_.dimC()),
    u_vec_(linearalgebra::NewVector(robot_ptr->dim_velocity())),
    beta_vec_(linearalgebra::NewVector(robot_ptr->dim_velocity())),
    gamma_vec_(linearalgebra::NewVector(robot_ptr->dim_velocity())),
    dRNEA_da_dot_beta_(linearalgebra::NewVector(robot_ptr->dim_velocity())),
    dRNEA_dfext_dot_beta_(linearalgebra::NewVector(robot_ptr->dim_fext())),
    dcontact_dq_dot_alpha_(
        linearalgebra::NewVector(robot_ptr->dim_velocity())),
    dcontact_dv_dot_alpha_(
        linearalgebra::NewVector(robot_ptr->dim_velocity())),
    dcontact_da_dot_alpha_(
        linearalgebra::NewVector(robot_ptr->dim_velocity())) {
}

ZeroHorizonOCP::~ZeroHorizonOCP() {
}

void ZeroHorizonOCP::computeOptimalityResidual(const double t, const double* q, 
                                               const double* v, 
                                               const double* solution,
                                               double* optimality_residual) {
  cost_function_.phiv(t, q, v, gamma_vec_);
  robot_->updateKinematics(q, v, solution);
  robot_->updateContactsState();

  robot_->setFext(&(solution[dimv_]));
  robot_->RNEA(q, v, solution, u_vec_);
  computeLuCu(time, state_vec, &(state_vec[dim_configuration_]), solution_vec, 
              &(solution_vec[dim_velocity_]), u_vec_,
              &(solution_vec[dim_velocity_+2*dim_fext_]), beta_vec_);
  robot_->CRBADotVec(state_vec, beta_vec_, dRNEA_da_dot_beta_);
  robot_->RNEAPartialFextDotVec(state_vec, beta_vec_, dRNEA_dfext_dot_beta_);
  robot_->contactResidualsAndDerivatives(
      state_vec, &(state_vec[dim_configuration_]), solution_vec,
      &(solution_vec[dim_velocity_+dim_fext_]),
      &(optimality_residual[dim_velocity_+dim_fext_]),
      dcontact_dq_dot_alpha_, dcontact_dv_dot_alpha_, 
      dcontact_da_dot_alpha_);


  computeLaCa(time, state_vec, &(state_vec[dim_configuration_]), solution_vec, 
              &(solution_vec[dim_velocity_]), u_vec_,
              &(solution_vec[dim_velocity_+2*dim_fext_]), optimality_residual);
  for (int i=0; i<dim_velocity_; ++i) {
    optimality_residual[i] += gamma_vec_[i] + dRNEA_da_dot_beta_[i] 
                              + dcontact_da_dot_alpha_[i];
  }

  computeLfCf(time, state_vec, &(state_vec[dim_configuration_]), solution_vec, 
              &(solution_vec[dim_velocity_]), u_vec_,
              &(solution_vec[dim_velocity_+2*dim_fext_]), 
              &(optimality_residual[dim_velocity_]));
  for (int i=0; i<dim_fext_; ++i) {
    optimality_residual[dim_velocity_+i] += dRNEA_dfext_dot_beta_[i];
  }

  constraints_.constraintsResidual(
      time, state_vec, &(state_vec[dim_configuration_]), solution_vec,
      &(solution_vec[dim_velocity_]), u_vec_, 
      &(optimality_residual[dim_velocity_+2*dim_fext_]));
}

}

} // namespace robotcgmres