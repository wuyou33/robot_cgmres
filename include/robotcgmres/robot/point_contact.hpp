#ifndef ROBOTCGMRES_POINT_CONTACT_HPP_
#define ROBOTCGMRES_POINT_CONTACT_HPP_

#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"
#include "pinocchio/container/aligned-vector.hpp"
#include "pinocchio/spatial/force.hpp"


namespace robotcgmres {

class PointContact {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Constructor. Not allocate matrices.
  PointContact(const int contact_frame_id, const double baumgarte_alpha, 
               const double baumgarte_beta);
 
  // Destructor. 
  ~PointContact();

  // Move constructor.
  PointContact(PointContact&& other) noexcept;

  // Move assignment operator. 
  PointContact& operator=(PointContact&& other) noexcept;

  // Loads pinocchio::Model and allocate matrices.
  void loadPinocchioModel(const pinocchio::Model& model);

  // Resets the parameters for the Baumgarte's stabilization method.
  void resetBaugrarteParameters(const double alpha, const double beta);

  // Resets the contact point.
  void resetContactPoint(const Eigen::Vector3d& contact_point);

  // Resets the contact point by current kinematics of the robot. The kinematics
  // is passed through pinocchio::Data. Before calling this function, you have 
  // to update the kinematics (only related to configuraion) in pinocchio::Data.
  void resetContactPointByCurrentState(const pinocchio::Data& data);

  // Converts the double array of the contact force to joint forces.
  void contactForceToJointForce(
      const double* contact_force, 
      pinocchio::container::aligned_vector<pinocchio::Force>& joint_forces);

  // Computes the 3xdimv contact Jacobian in world frame. The size of J_contact 
  // must be 3xdimv.
  void contactJacobian(const pinocchio::Model& model, 
                       const pinocchio::Data& data, Eigen::MatrixXd& J_contact);

  // Computes the 3xdimv contact Jacobian in world frame. The size of J_contact 
  // must be at least 3xdimv. The contact Jacobian is set in the block part of 
  // J_contacts with the indices of row and column begin at row_begin and 
  // column_begin.
  void contactJacobian(const pinocchio::Model& model, 
                       const pinocchio::Data& data, Eigen::MatrixXd& J_contacts,
                       const int row_begin, const int column_begin);

  // Computes the residual of the contact constraints considered by the 
  // Baumgarte's stabilization method.
  void baumgarteResidual(const pinocchio::Model& model, 
                         const pinocchio::Data& data, 
                         double* baumgarte_residual);

  // Computes the the partial derivatives of the contact constraints
  // considered by the Baumgarte's stabilization method. The size of each matrix 
  // must be 3xdimv. 
  void baumgarteDerivatives(const pinocchio::Model& model, 
                            pinocchio::Data& data, 
                            Eigen::MatrixXd& baumgarte_partial_dq, 
                            Eigen::MatrixXd& baumgarte_partial_dv, 
                            Eigen::MatrixXd& baumgarte_partial_da);

  // Computes the the partial derivatives of the contact constraints
  // considered by the Baumgarte's stabilization method. The size of each matrix 
  // must be at least 3xdimv. The derivatives is set in the block part of 
  // baumgarte_partial_dx with the indices of row and column begin at 
  // row_begin and column_begin.
  void baumgarteDerivatives(const pinocchio::Model& model, 
                            pinocchio::Data& data, 
                            Eigen::MatrixXd& baumgarte_partial_dq, 
                            Eigen::MatrixXd& baumgarte_partial_dv, 
                            Eigen::MatrixXd& baumgarte_partial_da,
                            const int row_begin, const int column_begin);

  // Returns contact_frame_id, the index of the contact frame.
  int contact_frame_id() const;

  // Returns parent_joint_id, the index of the parent joint of the contact 
  // frame.
  int parent_joint_id() const;

  // Returns dimv, the dimensiton of the generalized velocity.
  int dimv() const;

  // Returns baumgarte_alpha, the weight parameter on the contact velocity in 
  // the Baumgarte's stabilization method.
  double baumgarte_alpha() const;

  // Returns baumgarte_alpha, the weight parameter on the contact position in 
  // the Baumgarte's stabilization method.
  double baumgarte_beta() const;
  
  // Returns the contact point.
  Eigen::Vector3d contact_point() const;

  // Returns SE(3) translation from the parent joint to the contact frame.
  pinocchio::SE3 jXf() const;

  // Prohibits copy constructor.
  PointContact(const PointContact&) = delete;

  // Prohibits copy operator.
  PointContact& operator=(const PointContact&) = delete;


private:
  int contact_frame_id_, parent_joint_id_, dimv_;
  double baumgarte_alpha_, baumgarte_beta_;
  Eigen::Vector3d contact_point_;

  pinocchio::SE3 jXf_;
  pinocchio::SE3::ActionMatrixType fXj_;
  Eigen::MatrixXd J_frame_, joint_v_partial_dq_, joint_a_partial_dq_, 
                  joint_a_partial_dv_, joint_a_partial_da_;
};

} // namespace robotcgmres


#endif // ROBOTCGMRES_POINT_CONTACT_HPP_