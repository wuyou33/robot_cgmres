#include "robot/point_contact.hpp"


namespace robotcgmres {

PointContact::PointContact(const int contact_frame_id, 
                           const double baumgarte_alpha, 
                           const double baumgarte_beta) 
  : contact_frame_id_(contact_frame_id),
    parent_joint_id_(0), 
    dimv_(0),
    baumgarte_alpha_(baumgarte_alpha),
    baumgarte_beta_(baumgarte_beta),
    contact_point_(Eigen::Vector3d::Zero()),
    jXf_(pinocchio::SE3::Identity()),
    fXj_(jXf_.inverse().toActionMatrix().transpose()),
    J_frame_(),
    joint_v_partial_dq_(),
    joint_a_partial_dq_(),
    joint_a_partial_dv_(),
    joint_a_partial_da_() {
}

PointContact::~PointContact() {
}

PointContact::PointContact(PointContact&& other) noexcept
  : contact_frame_id_(other.contact_frame_id()),
    parent_joint_id_(other.parent_joint_id()), 
    dimv_(other.dimv()),
    baumgarte_alpha_(other.baumgarte_alpha()),
    baumgarte_beta_(other.baumgarte_beta()),
    contact_point_(other.contact_point()),
    jXf_(other.jXf()),
    fXj_(other.jXf().inverse().toActionMatrix().transpose()),
    J_frame_(),
    joint_v_partial_dq_(),
    joint_a_partial_dq_(),
    joint_a_partial_dv_(),
    joint_a_partial_da_() {
  if (dimv_> 0) {
    J_frame_.resize(6, dimv_); 
    joint_v_partial_dq_.resize(6, dimv_); 
    joint_a_partial_dq_.resize(6, dimv_); 
    joint_a_partial_dv_.resize(6, dimv_); 
    joint_a_partial_da_.resize(6, dimv_); 
    J_frame_.fill(0.0); 
    joint_v_partial_dq_.fill(0.0); 
    joint_a_partial_dq_.fill(0.0); 
    joint_a_partial_dv_.fill(0.0); 
    joint_a_partial_da_.fill(0.0); 
  }
}

PointContact& PointContact::operator=(PointContact&& other) noexcept {
  contact_frame_id_ = other.contact_frame_id();
  parent_joint_id_ = other.parent_joint_id();
  dimv_ = other.dimv();
  baumgarte_alpha_ = other.baumgarte_alpha();
  baumgarte_beta_ = other.baumgarte_beta();
  contact_point_ = other.contact_point();
  jXf_ = other.jXf();
  fXj_ = other.jXf().inverse().toActionMatrix().transpose();
  if (dimv_> 0) {
    J_frame_.resize(6, dimv_); 
    joint_v_partial_dq_.resize(6, dimv_); 
    joint_a_partial_dq_.resize(6, dimv_); 
    joint_a_partial_dv_.resize(6, dimv_); 
    joint_a_partial_da_.resize(6, dimv_); 
    J_frame_.fill(0.0); 
    joint_v_partial_dq_.fill(0.0); 
    joint_a_partial_dq_.fill(0.0); 
    joint_a_partial_dv_.fill(0.0); 
    joint_a_partial_da_.fill(0.0); 
  }
}

void PointContact::loadPinocchioModel(const pinocchio::Model& model) {
  parent_joint_id_ = model.frames[contact_frame_id_].parent;
  dimv_ = model.nv;
  jXf_ = model.frames[contact_frame_id_].placement;
  fXj_ = jXf_.inverse().toActionMatrix().transpose();
  J_frame_.resize(6, dimv_); 
  joint_v_partial_dq_.resize(6, dimv_); 
  joint_a_partial_dq_.resize(6, dimv_); 
  joint_a_partial_dv_.resize(6, dimv_); 
  joint_a_partial_da_.resize(6, dimv_); 
  J_frame_.fill(0.0); 
  joint_v_partial_dq_.fill(0.0); 
  joint_a_partial_dq_.fill(0.0); 
  joint_a_partial_dv_.fill(0.0); 
  joint_a_partial_da_.fill(0.0); 
}

void PointContact::resetBaugrarteParameters(const double alpha, 
                                            const double beta) {
  baumgarte_alpha_ = alpha;
  baumgarte_beta_ = beta;
}

void PointContact::resetContactPoint(const Eigen::Vector3d& contact_point) {
  contact_point_ = contact_point;
}

void PointContact::resetContactPointByCurrentState(
    const pinocchio::Data& data) {
  contact_point_ = data.oMf[contact_frame_id_].translation();
}

void PointContact::contactForceToJointForce(
    const double* contact_force, 
    pinocchio::container::aligned_vector<pinocchio::Force>& joint_forces) {
  joint_forces[parent_joint_id_] 
      = jXf_.act(
          pinocchio::Force(Eigen::Map<const Eigen::Vector3d>(contact_force), 
                           Eigen::Vector3d::Zero()));
}

void PointContact::contactJacobian(const pinocchio::Model& model, 
                                   const pinocchio::Data& data, 
                                   Eigen::MatrixXd& J_contact) {
  pinocchio::getFrameJacobian(model, data, contact_frame_id_, pinocchio::LOCAL, 
                              J_contact);
}

void PointContact::contactJacobian(const pinocchio::Model& model, 
                                   const pinocchio::Data& data, 
                                   Eigen::MatrixXd& J_contacts,
                                   const int row_begin, 
                                   const int column_begin) {
  pinocchio::getFrameJacobian(model, data, contact_frame_id_, pinocchio::LOCAL, 
                              J_contacts.block(row_begin, column_begin, 
                                               6, dimv_));
}

void PointContact::baumgarteResidual(const pinocchio::Model& model, 
                                     const pinocchio::Data& data, 
                                     double* baumgarte_residual) {
  Eigen::Map<Eigen::VectorXd>(baumgarte_residual, 3) 
      = (data.oMi[parent_joint_id_].act(data.a[parent_joint_id_])).linear()
          + baumgarte_alpha_ 
              * (data.oMi[parent_joint_id_].act(data.v[parent_joint_id_])).linear()
          + baumgarte_beta_  
              * (data.oMf[contact_frame_id_].translation()-contact_point_);
}

void PointContact::baumgarteDerivatives(const pinocchio::Model& model, 
                                        pinocchio::Data& data, 
                                        Eigen::MatrixXd& baumgarte_partial_dq, 
                                        Eigen::MatrixXd& baumgarte_partial_dv, 
                                        Eigen::MatrixXd& baumgarte_partial_da) {
  pinocchio::getJointAccelerationDerivatives(model, data, parent_joint_id_, 
                                             pinocchio::WORLD,
                                             joint_v_partial_dq_, 
                                             joint_a_partial_dq_, 
                                             joint_a_partial_dv_, 
                                             joint_a_partial_da_);
  pinocchio::getFrameJacobian(model, data, contact_frame_id_, 
                              pinocchio::LOCAL_WORLD_ALIGNED, J_frame_);
  baumgarte_partial_dq
      = joint_a_partial_dq_.template topRows<3>()
          + baumgarte_alpha_ * joint_v_partial_dq_.template topRows<3>()
          + baumgarte_beta_ * J_frame_.template topRows<3>();
  baumgarte_partial_dv 
      = joint_a_partial_dv_.template topRows<3>()
          + baumgarte_alpha_ * joint_a_partial_da_.template topRows<3>();
  baumgarte_partial_da 
      = joint_a_partial_da_.template topRows<3>();
}

void PointContact::baumgarteDerivatives(const pinocchio::Model& model, 
                                        pinocchio::Data& data, 
                                        Eigen::MatrixXd& baumgarte_partial_dq, 
                                        Eigen::MatrixXd& baumgarte_partial_dv, 
                                        Eigen::MatrixXd& baumgarte_partial_da,
                                        const int row_begin, 
                                        const int column_begin) {
  pinocchio::getJointAccelerationDerivatives(model, data, parent_joint_id_, 
                                             pinocchio::WORLD,
                                             joint_v_partial_dq_, 
                                             joint_a_partial_dq_, 
                                             joint_a_partial_dv_, 
                                             joint_a_partial_da_);
  pinocchio::getFrameJacobian(model, data, contact_frame_id_, 
                              pinocchio::LOCAL_WORLD_ALIGNED, J_frame_);
  baumgarte_partial_dq.block(row_begin, column_begin, 3, dimv_) 
      = joint_a_partial_dq_.template topRows<3>()
          + baumgarte_alpha_ * joint_v_partial_dq_.template topRows<3>()
          + baumgarte_beta_ * J_frame_.template topRows<3>();
  baumgarte_partial_dv.block(row_begin, column_begin, 3, dimv_) 
      = joint_a_partial_dv_.template topRows<3>()
          + baumgarte_alpha_ * joint_a_partial_da_.template topRows<3>();
  baumgarte_partial_da.block(row_begin, column_begin, 3, dimv_) 
      = joint_a_partial_da_.template topRows<3>();
}

int PointContact::contact_frame_id() const {
  return contact_frame_id_;
}

int PointContact::parent_joint_id() const {
  return parent_joint_id_;
}

int PointContact::dimv() const {
  return dimv_;
}

double PointContact::baumgarte_alpha() const {
  return baumgarte_alpha_;
}

double PointContact::baumgarte_beta() const {
  return baumgarte_beta_;
}

Eigen::Vector3d PointContact::contact_point() const {
  return contact_point_;
}

pinocchio::SE3 PointContact::jXf() const {
  return jXf_;
}

} // namespace robotcgmres