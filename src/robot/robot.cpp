#include "robot/robot.hpp"


namespace robotcgmres {

Robot::Robot(const std::string& urdf_file_name)
  : model_(),
    data_(model_),
    contacts_(),
    fjoint_(),
    dtau_dfext_(),
    dimq_(0),
    dimv_(0),
    dimtau_(0),
    dimf_(0) {
  pinocchio::urdf::buildModel(urdf_file_name, model_);
  data_ = pinocchio::Data(model_);
  fjoint_ = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  dimq_ = model_.nq;
  dimv_ = model_.nv;
  dimtau_ = model_.nv;
}

Robot::~Robot() {
}

void Robot::updateKinematics(const double* q) {
  pinocchio::forwardKinematics(model_, data_, 
                               Eigen::Map<const Eigen::VectorXd>(q, dimq_));
  pinocchio::computeJointJacobians(model_, data_, 
                                   Eigen::Map<const Eigen::VectorXd>(q, dimq_));
  pinocchio::updateFramePlacements(model_, data_);
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

void Robot::RNEADerivativesTransDotVector(const double* q, const double* v, 
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

void Robot::RNEADerivativesTransDotVector(const double* q, const double* v, 
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
      = dtau_dfext_.topRows(dimf_) * Eigen::Map<const Eigen::VectorXd>(vec, dimv_);
}

void Robot::CRBADotVec(const double* q, const double* vec, double* result) {
  pinocchio::crba(model_, data_, Eigen::Map<const Eigen::VectorXd>(q, dimq_));
  data_.M.triangularView<Eigen::StrictlyLower>() 
      = data_.M.transpose().triangularView<Eigen::StrictlyLower>();
  Eigen::Map<Eigen::VectorXd>(result, dimv_) 
      = data_.M * Eigen::Map<const Eigen::VectorXd>(vec, dimv_);
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
  dtau_dfext_.fill(0.0);
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
}

int Robot::dimq() const {
  return dimq_;
}

int Robot::dimv() const {
  return dimv_;
}

int Robot::dimtau() const {
  return dimv_;
}

int Robot::dimf() const {
  return dimf_;
}

} // namespace robotcgmres