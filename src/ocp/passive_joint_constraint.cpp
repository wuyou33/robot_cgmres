#include "ocp/passive_joint_constraint.hpp"


namespace robotocgmres {

PassiveJointConstraint::PassiveJointConstraint(const Robot* robot_ptr)
  : passive_torques_(robot_ptr.passive_torques()) {
}

void PassiveJointConstraint::setPassiveTorques(double* tau) const {
  for (int i=0; i<passive_torques_.size(); ++i) {
    tau[i] = 0.0;
  }
}

void PassiveJointConstraint::constraintsResidual(const double* tau, 
                                                 double* residual) const {
  for (int i=0; i<passive_torques_.size(); ++i) {
    residual[passive_torques_[i]] = tau[i];
  }
}
           
void PassiveJointConstraint::constraintsDerivativeDotVec(
    const double* vec, double* derivative_dot_vec) const {
  for (int i=0; i<passive_torques_.size(); ++i) {
    derivative_dot_vec[i] = vec[i];
  }
}

int PassiveJointConstraint::dim_passive() const {
  return derivative_dot_vec.size();
}

} // namespace robotocgmres