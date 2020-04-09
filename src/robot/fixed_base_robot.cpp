#include "robot/fixed_base_robot.hpp"


namespace robotcgmres {

FixedBaseRobot::FixedBaseRobot(const std::string& urdf_file_name)
  : model_(),
    data_(model_),
    dimq_(0) {
  pinocchio::urdf::buildModel(urdf_file_name, model_);
  data_ = pinocchio::Data(model_);
  dimq_ = model_.nq;
}

FixedBaseRobot::~FixedBaseRobot() {
}

void FixedBaseRobot::RNEA(const double* q, const double* v, const double* a, 
                          double* tau) {
  Eigen::Map<Eigen::VectorXd>(tau, dimq_)
      = pinocchio::rnea(model_, data_, 
                        Eigen::Map<const Eigen::VectorXd>(q, dimq_), 
                        Eigen::Map<const Eigen::VectorXd>(v, dimq_), 
                        Eigen::Map<const Eigen::VectorXd>(a, dimq_));
}

void FixedBaseRobot::RNEADerivativesTransDotVector(const double* q, 
                                                   const double* v, 
                                                   const double* a, 
                                                   const double* vec, 
                                                   double* dRNEA_dq_dot_vec, 
                                                   double* dRNEA_dv_dot_vec, 
                                                   double* dRNEA_da_dot_vec) {
  pinocchio::computeRNEADerivatives(model_, data_, 
                                    Eigen::Map<const Eigen::VectorXd>(q, dimq_), 
                                    Eigen::Map<const Eigen::VectorXd>(v, dimq_), 
                                    Eigen::Map<const Eigen::VectorXd>(a, dimq_));
  data_.M.triangularView<Eigen::StrictlyLower>() 
      = data_.M.transpose().triangularView<Eigen::StrictlyLower>();
  Eigen::Map<Eigen::VectorXd>(dRNEA_dq_dot_vec, dimq_) 
      = data_.dtau_dq.transpose() 
          * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  Eigen::Map<Eigen::VectorXd>(dRNEA_dv_dot_vec, dimq_) 
      = data_.dtau_dv.transpose() 
          * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  Eigen::Map<Eigen::VectorXd>(dRNEA_da_dot_vec, dimq_) 
      = data_.M * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
}

void FixedBaseRobot::CRBADotVec(const double* q, const double* vec, 
                                double* result) {
  pinocchio::crba(model_, data_, Eigen::Map<const Eigen::VectorXd>(q, dimq_));
  data_.M.triangularView<Eigen::StrictlyLower>() 
      = data_.M.transpose().triangularView<Eigen::StrictlyLower>();
  Eigen::Map<Eigen::VectorXd>(result, dimq_) 
      = data_.M * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
}

int FixedBaseRobot::dimq() const {
  return dimq_;
}

int FixedBaseRobot::dimv() const {
  return dimq_;
}

int FixedBaseRobot::dimtau() const {
  return dimq_;
}

} // namespace robotcgmres