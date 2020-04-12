#ifndef ROBOTCGMRES_PASSIVE_JOINT_CONSTRAINTS_HPP_
#define ROBOTCGMRES_PASSIVE_JOINT_CONSTRAINTS_HPP_

#include <string>
#include <vector>

#include "robot/robot.hpp"


namespace robotcgmres {

class PassiveJointConstraint {
public:

  // Constructor. Extract the passive joints from robot_ptr.
  PassiveJointConstraint(const Robot* robot_ptr);

  // Destructor. 
  ~PassiveJointConstraint();

  // Substitutes zero in the generalized torques tau corresponding to the 
  // passive joints.
  void setPassiveTorques(double* tau) const;

  // Calculates the residual of the constrains on the passive joints under given
  // generalized torques tau.
  void constraintsResidual(const double* tau, double* residual) const;

  // Calculates the product of a vector and the derivative of the residual of 
  // the constrains on the passive joints. In usually, the vector must be 
  // the corresponding Lagrange multiplier.
  void constraintsDerivativeDotVec(const double* vec, 
                                   double* derivative_dot_vec) const;

  // Prohibits copy constructor.
  PassiveJointConstraint(const PassiveJointConstraint&) = delete;

  // Prohibits copy operator.
  PassiveJointConstraint& operator=(const PassiveJointConstraint&) = delete;

  // Returns dim_passive_, the dimensiton of the generalized torques 
  // corresponding to the passive joints.
  int dim_passive() const;

private:
  std::vector<int> passive_torques_;

};

} // namespace robotcgmres


#endif // ROBOTCGMRES_PASSIVE_JOINT_CONSTRAINTS_HPP_