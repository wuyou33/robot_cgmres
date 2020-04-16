#ifndef ROBOTCGMRES_PASSIVE_JOINTS_HPP_
#define ROBOTCGMRES_PASSIVE_JOINTS_HPP_

#include <string>
#include <vector>


namespace robotcgmres {

class PassiveJoints {
public:

  // Constructor. Extract the passive joints from robot_ptr.
  // Argments:
  //   robot_ptr: Pointer of the robot model.
  PassiveJoints(std::vector<std::string>& joint_types);

  // Destructor. 
  ~PassiveJoints();

  // Substitutes zero in the generalized torques tau corresponding to the 
  // passive joints.
  // Argments:
  //   tau: The generalized torque for fully actuated system. The size is dimv.
  void setPassiveTorques(double* tau) const;

  // Calculates the residual of the constrains on the passive joints under given
  // generalized torques tau.
  // Argments:
  //   tau: The generalized torque for fully actuated system. The size is dimv.
  //   residual: The residual of the constraints of the zero torques. The size
  //      is dim_passive.
  void passiveConstraintsResidual(const double* tau, double* residual) const;

  // Calculates the product of a vector and the derivative of the residual of 
  // the constrains on the passive joints. In usually, the vector must be 
  // the corresponding Lagrange multiplier.
  // Argments:
  //   vec: A vector whose size is dim_passive.
  //   added_vec: A vector which vec is added. The size must be Robot::dimv().
  void addPassiveConstraintsDerivativeDotVec(const double* vec, 
                                            double* added_vec) const;

  // Returns the dimensiton of the generalized torques corresponding to the 
  // passive joints.
  int dim_passive() const;

  // Prohibits copy constructor.
  PassiveJoints(const PassiveJoints&) = delete;

  // Prohibits copy operator.
  PassiveJoints& operator=(const PassiveJoints&) = delete;

private:
  std::vector<int> passive_torque_indices_;

};

} // namespace robotcgmres


#endif // ROBOTCGMRES_PASSIVE_JOINTS_HPP_