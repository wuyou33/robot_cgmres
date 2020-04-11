#ifndef ROBOTCGMRES_FIXED_BASE_ROBOT_HPP_
#define ROBOTCGMRES_FIXED_BASE_ROBOT_HPP_

#include <string>
#include <vector>
#include <algorithm>

#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"

#include "robot/point_contact.hpp"


namespace robotcgmres {

class Robot {
public:

  // Constructor. Build workspace for pinocchio.
  Robot(const std::string& urdf_file_name);

  // Destructor. 
  ~Robot();

  // Sets the stack of the generalized forces represented in world frame.
  // The array is converted into joint forces.
  void setFext(const double* fext);

  // Computes generalized torques tau corresponding to given q, v, and a.
  // No external forces are assumed.
  void RNEA(const double* q, const double* v, const double* a, double* tau);

  // Computes generalized torques tau corresponding to given q, v, a, and fext.
  // External forces are set by setFext() before calling this function.
  void RNEA(const double* q, const double* v, const double* a,  
            const double* fext, double* tau);

  // Computes the partial dervatives of the function of inverse dynamics with 
  // respect to q, v, and a.
  // No external forces are assumed.
  void RNEADerivativesTransDotVector(const double* q, const double* v, 
                                     const double* a, const double* vec, 
                                     double* dRNEA_dq_dot_vec, 
                                     double* dRNEA_dv_dot_vec, 
                                     double* dRNEA_da_dot_vec);

  // Computes the partial dervatives of the function of inverse dynamics with 
  // respect to q, v, a, and fext.
  // External forces are set by setFext() before calling this function.
  void RNEADerivativesTransDotVector(const double* q, const double* v, 
                                     const double* a, const double* vec, 
                                     double* dRNEA_dq_dot_vec, 
                                     double* dRNEA_dv_dot_vec, 
                                     double* dRNEA_da_dot_vec,
                                     double* dRNEA_dfext_dot_vec);

  // Computes the joint inertia matrix.
  void CRBADotVec(const double* q, const double* vec, double* result);

  // Add the point contact to the robot.
  void addPointContact(const int contact_frame_id, const double baumgarte_alpha, 
                       const double baumgarte_beta);

  // Remove the point contact from the robot.
  void removePointContact(const int contact_frame_id);

  // Returns dimq, the dimensiton of the generalized configuration.
  int dimq() const;

  // Returns dimv, the dimensiton of the generalized velocity.
  int dimv() const;

  // Returns dimv, the dimensiton of the actuated generalized torques.
  int dimtau() const;

  // Returns dimf, the dimensiton of the contact and the contact forces.
  int dimf() const;

  // Prohibits copy constructor.
  Robot(const Robot&) = delete;

  // Prohibits copy operator.
  Robot& operator=(const Robot&) = delete;

private:
  pinocchio::Model model_;
  pinocchio::Data data_;
  std::vector<PointContact> contacts_;
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint_;
  Eigen::MatrixXd dtau_dfext_;
  int dimq_, dimv_, dimtau_, dimf_;
};

} // namespace robotcgmres


#endif // ROBOTCGMRES_FIXED_BASE_ROBOT_HPP_