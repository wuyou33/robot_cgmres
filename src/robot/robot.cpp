#include "robot/robot.hpp"


namespace robotcgmres {

Robot::Robot(const std::string& urdf_file_name)
  : model_(),
    data_(model_),
    contacts_(),
    fjoint_(),
    passive_torques_(),
    dtau_dfext_(),
    dBaum_dq_(),
    dBaum_dv_(),
    dBaum_da_(),
    dimq_(0),
    dimv_(0),
    dimf_(0) {
  // Build Pinocchio model from URDF.
  pinocchio::urdf::buildModel(urdf_file_name, model_);
  data_ = pinocchio::Data(model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  dimq_ = model_.nq;
  dimv_ = model_.nv;
  // Set passive joints.
  int total_dim_torques = 0;
  for (const auto& joint : model_.joints) {
    if (joint.shortname() == "JointModelFreeFlyer") {
      for (int j=0; j<6; ++j) {
        passive_torques_.push_back(total_dim_torques+j);
      }
      total_dim_torques += 6;
    }
    else {
      total_dim_torques += 1;
    }
  }
}

Robot::~Robot() {
}

void Robot::updateKinematics(const double* q, const double* v, 
                             const double* a) {
  pinocchio::forwardKinematics(model_, data_, 
                               Eigen::Map<const Eigen::VectorXd>(q, dimq_),
                               Eigen::Map<const Eigen::VectorXd>(v, dimv_),
                               Eigen::Map<const Eigen::VectorXd>(a, dimv_));
  pinocchio::computeForwardKinematicsDerivatives(
      model_, data_, Eigen::Map<const Eigen::VectorXd>(q, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(v, dimv_), 
      Eigen::Map<const Eigen::VectorXd>(a, dimv_));
  pinocchio::updateFramePlacements(model_, data_);
  for (int i=0; i<contacts_.size(); ++i) {
    contacts_[i].resetContactPointByCurrentState(data_);
  }
}

void Robot::baumgarteResidual(double* residual) {
  for (int i=0; i<contacts_.size(); ++i) {
    contacts_[i].baumgarteResidual(model_, data_, &(residual[3*i]));
  }
}

void Robot::baumgarteDerivativesDotVec(const double* vec,
                                       double* dBaum_dq_dot_vec, 
                                       double* dBaum_dv_dot_vec, 
                                       double* dBaum_da_dot_vec) {
  for (int i=0; i<contacts_.size(); ++i) {
    contacts_[i].baumgarteDerivatives(model_, data_, dBaum_dq_, dBaum_dv_, 
                                      dBaum_da_, 0, 3*i);
  }
  Eigen::Map<Eigen::VectorXd>(dBaum_dq_dot_vec, dimv_) 
      = dBaum_dq_.transpose() * Eigen::Map<const Eigen::VectorXd>(vec, dimf_);
  Eigen::Map<Eigen::VectorXd>(dBaum_dv_dot_vec, dimv_) 
      = dBaum_dv_.transpose() * Eigen::Map<const Eigen::VectorXd>(vec, dimf_);
  Eigen::Map<Eigen::VectorXd>(dBaum_da_dot_vec, dimv_) 
      = dBaum_da_.transpose() * Eigen::Map<const Eigen::VectorXd>(vec, dimf_);
}

void Robot::setFext(const double* fext) {
  for (int i=0; i<contacts_.size(); ++i) {
    contacts_[i].contactForceToJointForce(&(fext[i*3]), fjoint_);
  }
}

void Robot::RNEA(const double* q, const double* v, const double* a, 
                 double* tau) {
  Eigen::Map<Eigen::VectorXd>(tau, dimv_)
      = pinocchio::rnea(model_, data_, 
                        Eigen::Map<const Eigen::VectorXd>(q, dimq_), 
                        Eigen::Map<const Eigen::VectorXd>(v, dimv_), 
                        Eigen::Map<const Eigen::VectorXd>(a, dimv_));
}

void Robot::RNEA(const double* q, const double* v, const double* a, 
                 const double* fext, double* tau) {
  Eigen::Map<Eigen::VectorXd>(tau, dimv_)
      = pinocchio::rnea(model_, data_, 
                        Eigen::Map<const Eigen::VectorXd>(q, dimq_), 
                        Eigen::Map<const Eigen::VectorXd>(v, dimv_), 
                        Eigen::Map<const Eigen::VectorXd>(a, dimv_),
                        fjoint_);
}

void Robot::RNEADerivativesTransDotVec(const double* q, const double* v, 
                                       const double* a, const double* vec, 
                                       double* dRNEA_dq_dot_vec, 
                                       double* dRNEA_dv_dot_vec, 
                                       double* dRNEA_da_dot_vec) {
  pinocchio::computeRNEADerivatives(model_, data_, 
                                    Eigen::Map<const Eigen::VectorXd>(q, dimq_), 
                                    Eigen::Map<const Eigen::VectorXd>(v, dimv_), 
                                    Eigen::Map<const Eigen::VectorXd>(a, dimv_));
  data_.M.triangularView<Eigen::StrictlyLower>() 
      = data_.M.transpose().triangularView<Eigen::StrictlyLower>();
  Eigen::Map<Eigen::VectorXd>(dRNEA_dq_dot_vec, dimv_) 
      = data_.dtau_dq.transpose() 
          * Eigen::Map<const Eigen::VectorXd>(vec, dimv_);
  Eigen::Map<Eigen::VectorXd>(dRNEA_dv_dot_vec, dimv_) 
      = data_.dtau_dv.transpose() 
          * Eigen::Map<const Eigen::VectorXd>(vec, dimv_);
  Eigen::Map<Eigen::VectorXd>(dRNEA_da_dot_vec, dimv_) 
      = data_.M * Eigen::Map<const Eigen::VectorXd>(vec, dimv_);
}

void Robot::RNEADerivativesTransDotVec(const double* q, const double* vec, 
                                       double* dRNEA_da_dot_vec) {
  pinocchio::crba(model_, data_, Eigen::Map<const Eigen::VectorXd>(q, dimq_));
  data_.M.triangularView<Eigen::StrictlyLower>() 
      = data_.M.transpose().triangularView<Eigen::StrictlyLower>();
  Eigen::Map<Eigen::VectorXd>(dRNEA_da_dot_vec, dimv_) 
      = data_.M * Eigen::Map<const Eigen::VectorXd>(vec, dimv_);
}

void Robot::RNEADerivativesTransDotVec(const double* q, const double* v, 
                                       const double* a, const double* vec, 
                                       double* dRNEA_dq_dot_vec, 
                                       double* dRNEA_dv_dot_vec, 
                                       double* dRNEA_da_dot_vec,
                                       double* dRNEA_dfext_dot_vec) {
  pinocchio::computeRNEADerivatives(model_, data_, 
                                    Eigen::Map<const Eigen::VectorXd>(q, dimq_), 
                                    Eigen::Map<const Eigen::VectorXd>(v, dimv_), 
                                    Eigen::Map<const Eigen::VectorXd>(a, dimv_),
                                    fjoint_);
  data_.M.triangularView<Eigen::StrictlyLower>() 
      = data_.M.transpose().triangularView<Eigen::StrictlyLower>();
  Eigen::Map<Eigen::VectorXd>(dRNEA_dq_dot_vec, dimv_) 
      = data_.dtau_dq.transpose() 
          * Eigen::Map<const Eigen::VectorXd>(vec, dimv_);
  Eigen::Map<Eigen::VectorXd>(dRNEA_dv_dot_vec, dimv_) 
      = data_.dtau_dv.transpose() 
          * Eigen::Map<const Eigen::VectorXd>(vec, dimv_);
  Eigen::Map<Eigen::VectorXd>(dRNEA_da_dot_vec, dimv_) 
      = data_.M * Eigen::Map<const Eigen::VectorXd>(vec, dimv_);
  for (int i=0; i<contacts_.size(); ++i) {
    contacts_[i].contactJacobian(model_, data_, dtau_dfext_, 0, 3*i);
  }
  Eigen::Map<Eigen::VectorXd>(dRNEA_dfext_dot_vec, dimf_) 
      = dtau_dfext_.topRows(dimf_) * Eigen::Map<const Eigen::VectorXd>(vec,   
                                                                       dimv_);
}

void Robot::RNEADerivativesTransDotVec(const double* q, const double* vec, 
                                       double* dRNEA_da_dot_vec,
                                       double* dRNEA_dfext_dot_vec) {
  pinocchio::crba(model_, data_, Eigen::Map<const Eigen::VectorXd>(q, dimq_));
  data_.M.triangularView<Eigen::StrictlyLower>() 
      = data_.M.transpose().triangularView<Eigen::StrictlyLower>();
  Eigen::Map<Eigen::VectorXd>(dRNEA_da_dot_vec, dimv_) 
      = data_.M * Eigen::Map<const Eigen::VectorXd>(vec, dimv_);
  for (int i=0; i<contacts_.size(); ++i) {
    contacts_[i].contactJacobian(model_, data_, dtau_dfext_, 0, 3*i);
  }
  Eigen::Map<Eigen::VectorXd>(dRNEA_dfext_dot_vec, dimf_) 
      = dtau_dfext_.topRows(dimf_) * Eigen::Map<const Eigen::VectorXd>(vec,   
                                                                       dimv_);
}

void Robot::addPointContact(const int contact_frame_id, 
                            const double baumgarte_alpha, 
                            const double baumgarte_beta) {
  auto find
      = std::find_if(contacts_.begin(), contacts_.end(), 
                     [contact_frame_id](PointContact& contact)
                     { return contact.contact_frame_id()==contact_frame_id; });
  if (find == contacts_.end()) {
    contacts_.push_back(PointContact(contact_frame_id, baumgarte_alpha, 
                                     baumgarte_beta));
    contacts_.back().loadPinocchioModel(model_);
    dimf_ += 3;
  }
  dtau_dfext_.resize(dimf_+3, dimv_);
  dBaum_dq_.resize(dimf_, dimv_);
  dBaum_dv_.resize(dimf_, dimv_);
  dBaum_da_.resize(dimf_, dimv_);
}

void Robot::removePointContact(const int contact_frame_id) {
  auto find
      = std::find_if(contacts_.begin(), contacts_.end(), 
                     [contact_frame_id](PointContact& contact)
                     { return contact.contact_frame_id()==contact_frame_id; });
  if (find != contacts_.end()) {
    std::swap<PointContact>(*find, contacts_.back());
    contacts_.pop_back();
    dimf_ -= 3;
  }
  dtau_dfext_.resize(dimf_+3, dimv_);
  dBaum_dq_.resize(dimf_, dimv_);
  dBaum_dv_.resize(dimf_, dimv_);
  dBaum_da_.resize(dimf_, dimv_);
}

int Robot::dimq() const {
  return dimq_;
}

int Robot::dimv() const {
  return dimv_;
}

int Robot::dimf() const {
  return dimf_;
}

const std::vector<int>& Robot::passive_torques() const {
  return passive_torques_;
}

} // namespace robotcgmres