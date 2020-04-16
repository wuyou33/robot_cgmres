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
#include "pinocchio/algorithm/kinematics-derivatives.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"

#include "robot/point_contact.hpp"


namespace robotcgmres {

class Robot {
public:

  // Constructor. Builds workspace for pinocchio.
  Robot(const std::string& urdf_file_name);

  // Destructor. 
  ~Robot();

  // Updates the kinematics of the robot. The frame placements, frame velocity,
  // frame acceleration, and the relevant Jacobians are calculated. After that, 
  // the each contact residual is updated.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  void updateKinematics(const double* q, const double* v, const double* a);

  // Computes the residual of the contact constriants represented by 
  // Baumgarte's stabilization method. Before calling this function, 
  // updateKinematics() must be called.
  // Argments: 
  //   residual: The array where the result is stored. The must be dimf.
  void baumgarteResidual(double* residual);

  // Computes the product of a vector and the derivatives of the contact 
  // constriants represented by Baumgarte's stabilization method. 
  // Before calling this function, updateKinematics() must be called.
  // Argments: 
  //   vec: A vector multiplied to partial derivatives of the contact 
  //      contstraints represented by Baumgarte's method. 
  //      The size is assumed to be dimf.
  //   dBaum_dq_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  //   dBaum_dv_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  //   dBaum_da_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  void baumgarteDerivativesDotVec(const double* vec, double* dBaum_dq_dot_vec, 
                                  double* dBaum_dv_dot_vec, 
                                  double* dBaum_da_dot_vec);

  // Sets the stack of the generalized forces represented in the world frame.
  // The array is converted into joint forces.
  //   fext: The stack of the contact forces represented in the world frame.
  //      The size is assumed to be dimf.
  void setFext(const double* fext);

  // Computes generalized torques tau corresponding to given q, v, and a.
  // No external forces are assumed.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   tau: Generalized torques for fully actuated system. Size must be dimv.
  void RNEA(const double* q, const double* v, const double* a, double* tau);

  // Computes generalized torques tau corresponding to given q, v, a, and fext.
  // External forces are set by setFext() before calling this function.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   fext: The stack of the contact forces. Size must be dimf.
  //   tau: Generalized torques for fully actuated system. Size must be dimv.
  void RNEA(const double* q, const double* v, const double* a,  
            const double* fext, double* tau);

  // Computes the partial dervatives of the function of inverse dynamics with 
  // respect to q, v, and a, and returns products of them and vec.
  // No external forces are assumed.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   vec: A vector multiplied to partial derivatives of RNEA. 
  //      The size is assumed to be dimv.
  //   dRNEA_dq_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  //   dRNEA_dv_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  //   dRNEA_da_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  void RNEADerivativesTransDotVec(const double* q, const double* v, 
                                  const double* a, const double* vec, 
                                  double* dRNEA_dq_dot_vec, 
                                  double* dRNEA_dv_dot_vec, 
                                  double* dRNEA_da_dot_vec);

  // Computes the partial dervatives of the function of inverse dynamics with 
  // respect to a, and returns products of them and vec.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq.
  //   vec: A vector multiplied to partial derivatives of RNEA. 
  //      The size is assumed to be dimv.
  //   dRNEA_da_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  void RNEADerivativesTransDotVec(const double* q, const double* vec, 
                                  double* dRNEA_da_dot_vec);

  // Computes the partial dervatives of the function of inverse dynamics with 
  // respect to q, v, a, and fext, and returns products of them and vec.
  // External forces are set by setFext() and kinematics are updated by 
  // updateKinematics() before calling this function.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq.
  //   v: Generalized velocity. Size must be dimv.
  //   a: Generalized acceleration. Size must be dimv.
  //   vec: A vector multiplied to partial derivatives of RNEA. 
  //      The size is assumed to be dimv.
  //   dRNEA_dq_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  //   dRNEA_dv_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  //   dRNEA_da_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  //   dRNEA_dfext_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimf.
  void RNEADerivativesTransDotVec(const double* q, const double* v, 
                                  const double* a, const double* vec, 
                                  double* dRNEA_dq_dot_vec, 
                                  double* dRNEA_dv_dot_vec, 
                                  double* dRNEA_da_dot_vec,
                                  double* dRNEA_dfext_dot_vec);

  // Computes the partial dervatives of the function of inverse dynamics with 
  // respect to a and fext, and returns products of them and vec. 
  // updateKinematics() must be called before calling this function.
  // Argments: 
  //   q: Generalized configuration. Size must be dimq.
  //   vec: A vector multiplied to partial derivatives of RNEA. 
  //      The size is assumed to be dimv.
  //   dRNEA_da_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimv.
  //   dRNEA_dfext_dot_vec: The array where the result is stored. The size is
  //      assumed to be dimf.
  void RNEADerivativesTransDotVec(const double* q, const double* vec, 
                                  double* dRNEA_da_dot_vec,
                                  double* dRNEA_dfext_dot_vec);

  // Adds the point contact to the robot.
  // If there is a contact that have contact_frame_id, the contact does not 
  //    increased. Otherwise, a contact that has contact_frame_id is added.
  // Argments:
  //    contact_frame_id: The frame index of the contact. 
  //    baumgarte_alpha: The weight parameter of the Baumgrate's stabilization
  //      method
  //    baumgarte_beta: The weight parameter of the Baumgrate's stabilization
  //      method
  void addPointContact(const int contact_frame_id, const double baumgarte_alpha, 
                       const double baumgarte_beta);

  // Removes the point contact from the robot. If there is no contact that has 
  //    contact_frame_id, this function does not do anything.
  // Argments:
  //    contact_frame_id: The frame index of the contact. 
  void removePointContact(const int contact_frame_id);

  // Returns the dimensiton of the generalized configuration.
  int dimq() const;

  // Returns the dimensiton of the generalized velocity.
  int dimv() const;

  // Returns the dimensiton of the contacts that is equal to the dimension of 
  //    the contact forces.
  int dimf() const;

  // Returns a reference of a vector of indices of generalized torques 
  // corresponding to the passive joints.
  const std::vector<std::string>& joint_types() const;

  // Prohibits copy constructor.
  Robot(const Robot&) = delete;

  // Prohibits copy operator.
  Robot& operator=(const Robot&) = delete;

private:
  pinocchio::Model model_;
  pinocchio::Data data_;
  std::vector<PointContact> contacts_;
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint_;
  std::vector<std::string> joint_types_;
  Eigen::MatrixXd dtau_dfext_, dBaum_dq_, dBaum_dv_, dBaum_da_;
  int dimq_, dimv_, dimf_;
};

} // namespace robotcgmres


#endif // ROBOTCGMRES_FIXED_BASE_ROBOT_HPP_