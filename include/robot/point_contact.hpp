#ifndef POINT_CONTACT_HPP_
#define POINT_CONTACT_HPP_

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

  PointContact(const int contact_frame_id, const double baumgarte_alpha, 
               const double baumgarte_beta);

  ~PointContact();

  void loadPinocchioModel(pinocchio::Model& model);

  void resetBaugrarteParameters(const double alpha, const double beta);

  void resetContactPoint(const Eigen::Vector3d& contact_point);

  void resetContactPointByCurrentState(pinocchio::Data& data);

  void contactForceToJointForce(
      const double* contact_force, 
      pinocchio::container::aligned_vector<pinocchio::Force>& joint_forces);

  void contactForceToJointForce(
      const Eigen::Vector3d& contact_force, 
      pinocchio::container::aligned_vector<pinocchio::Force>& joint_forces);

  void contactJacobian(const pinocchio::Model& model, 
                       const pinocchio::Data& data, Eigen::MatrixXd& J_contact);

  void contactJacobian(const pinocchio::Model& model, 
                       const pinocchio::Data& data, Eigen::MatrixXd& J_contacts,
                       const int block_start_row, const int block_start_column);

  void baumgarteResidual(const pinocchio::Model& model, 
                         const pinocchio::Data& data, 
                         double* baumgarte_residual);
  void baumgarteResidual(const pinocchio::Model& model, 
                         const pinocchio::Data& data, 
                         Eigen::Vector3d& baumgarte_residual);

  void baumgarteDerivatives(const pinocchio::Model& model, 
                            pinocchio::Data& data, 
                            Eigen::MatrixXd& baumgarte_partial_dq, 
                            Eigen::MatrixXd& baumgarte_partial_dv, 
                            Eigen::MatrixXd& baumgarte_partial_da);

  void baumgarteDerivatives(const pinocchio::Model& model, 
                            pinocchio::Data& data, 
                            Eigen::MatrixXd& baumgarte_partial_dq, 
                            Eigen::MatrixXd& baumgarte_partial_dv, 
                            Eigen::MatrixXd& baumgarte_partial_da,
                            const int block_start_row, 
                            const int block_start_column);

  int contact_frame_id() const;

  int parent_joint_id() const;

  double baumgarte_alpha() const;

  double baumgarte_beta() const;
  
  Eigen::Vector3d contact_point() const;

  PointContact(const PointContact&) = delete;
  PointContact& operator=(const PointContact&) = delete;

private:
  const int contact_frame_id_;
  int parent_joint_id_, dimv_;
  double baumgarte_alpha_, baumgarte_beta_;
  Eigen::Vector3d contact_point_;

  pinocchio::SE3 jXf_;
  pinocchio::SE3::ActionMatrixType fXj_;
  Eigen::Matrix3d v_skew_, w_skew_;
  Eigen::MatrixXd J_local_, J_frame_, joint_v_partial_dq_, 
                  joint_a_partial_dq_, joint_a_partial_dv_, joint_a_partial_da_;
};

} // namespace robotcgmres


#endif // POINT_CONTACT_HPP_