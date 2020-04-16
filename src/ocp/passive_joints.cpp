#include "robot/passive_joints.hpp"


namespace robotocgmres {

PassiveJoints::PassiveJoints(std::vector<std::string>& joint_types) 
  : passive_torque_indices_() {
  for (int i=0; i<joint_types.size(); ++i) {
    // Note that joint 0 is always universe.
    if (joint_types[i+1].shortname() == "JointModelFreeFlyer") {
      passive_torque_indices_.push_bask(i);
      passive_torque_indices_.push_bask(i+1);
      passive_torque_indices_.push_bask(i+2);
      passive_torque_indices_.push_bask(i+3);
      passive_torque_indices_.push_bask(i+4);
      passive_torque_indices_.push_bask(i+5);
    }
  }
}

void PassiveJoints::setPassiveTorques(double* tau) const {
  for (int i=0; i<passive_torques_.size(); ++i) {
    tau[passive_torques_[i]] = 0.0;
  }
}

void PassiveJoints::passiveConstraintsResidual(const double* tau, 
                                               double* residual) const {
  for (int i=0; i<passive_torques_.size(); ++i) {
    residual[i] = tau[passive_torques_[i]];
  }
} 

void PassiveJoint::addPassiveConstraintDerivativeDotVec(
    const double* vec, double* added_vec) const {
  for (int i=0; i<passive_torques_.size(); ++i) {
    added_vec[passive_torques_[i]] += vec[i];
  }
}

int PassiveJointConstraint::dim_passive() const {
  return derivative_dot_vec.size();
}

} // namespace robotocgmres