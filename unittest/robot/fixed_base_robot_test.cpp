#include <string>
#include <random>

#include <gtest/gtest.h>
#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"

#include "common/memory_manager.hpp"
#include "robot/point_contact.hpp"
#include "robot/robot.hpp"


namespace robotcgmres {

class FixedBaseRobotTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_file_name_ = "../../../examples/iiwa14/iiwa14.urdf";
    pinocchio::urdf::buildModel(urdf_file_name_, model_);
    data_ = pinocchio::Data(model_);
    dimq_ = model_.nq;
    contact_frame_id_ = 18;
    q_ = memorymanager::NewVector(dimq_);
    v_ = memorymanager::NewVector(dimq_);
    a_ = memorymanager::NewVector(dimq_);
    Eigen::Map<Eigen::VectorXd>(q_, dimq_) 
        = pinocchio::randomConfiguration(model_, -Eigen::VectorXd::Ones(dimq_), 
                                         Eigen::VectorXd::Ones(dimq_));
    Eigen::Map<Eigen::VectorXd>(v_, dimq_) = Eigen::VectorXd::Random(dimq_);
    Eigen::Map<Eigen::VectorXd>(a_, dimq_) = Eigen::VectorXd::Random(dimq_);
    baumgarte_alpha_ = Eigen::VectorXd::Random(2)[0];
    baumgarte_beta_ = Eigen::VectorXd::Random(2)[0];
  }

  virtual void TearDown() {
    memorymanager::DeleteVector(q_);
    memorymanager::DeleteVector(v_);
    memorymanager::DeleteVector(a_);
  }

  std::string urdf_file_name_;
  pinocchio::Model model_;
  pinocchio::Data data_;
  int dimq_, contact_frame_id_;
  double *q_, *v_, *a_;
  double baumgarte_alpha_, baumgarte_beta_;
};


TEST_F(FixedBaseRobotTest, dim) {
  Robot robot(urdf_file_name_);
  EXPECT_EQ(robot.dimq(), dimq_);
  EXPECT_EQ(robot.dimv(), dimq_);
  EXPECT_EQ(robot.dimtau(), dimq_);
}


TEST_F(FixedBaseRobotTest, RNEAWithoutFext) {
  Robot robot(urdf_file_name_);
  double *tau = memorymanager::NewVector(dimq_);
  robot.RNEA(q_, v_, a_, tau);
  Eigen::VectorXd tau_ref = pinocchio::rnea(
      model_, data_,
      Eigen::Map<const Eigen::VectorXd>(q_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(v_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(a_, dimq_));
  EXPECT_TRUE(tau_ref.isApprox(Eigen::Map<Eigen::VectorXd>(tau, dimq_)));
  memorymanager::DeleteVector(tau);
}


TEST_F(FixedBaseRobotTest, RNEADerivativesTransDotVecWithoutFext) {
  Robot robot(urdf_file_name_);
  double *vec = memorymanager::NewVector(dimq_);
  Eigen::Map<Eigen::VectorXd>(vec, dimq_) = Eigen::VectorXd::Random(dimq_);
  double *dRNEA_dq_dot_vec = memorymanager::NewVector(dimq_);
  double *dRNEA_dv_dot_vec = memorymanager::NewVector(dimq_);
  double *dRNEA_da_dot_vec = memorymanager::NewVector(dimq_);
  robot.RNEADerivativesTransDotVector(q_, v_, a_, vec, dRNEA_dq_dot_vec, 
                                      dRNEA_dv_dot_vec, dRNEA_da_dot_vec);
  pinocchio::computeRNEADerivatives(
      model_, data_,
      Eigen::Map<const Eigen::VectorXd>(q_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(v_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(a_, dimq_));
  data_.M.triangularView<Eigen::StrictlyLower>() 
      = data_.M.transpose().triangularView<Eigen::StrictlyLower>();
  Eigen::VectorXd dRNEA_dq_dot_vec_ref 
      = data_.dtau_dq.transpose() 
          * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  Eigen::VectorXd dRNEA_dv_dot_vec_ref 
      = data_.dtau_dv.transpose() 
          * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  Eigen::VectorXd dRNEA_da_dot_vec_ref 
      = data_.M * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  EXPECT_TRUE(dRNEA_dq_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_dq_dot_vec, dimq_)));
  EXPECT_TRUE(dRNEA_dv_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_dv_dot_vec, dimq_)));
  EXPECT_TRUE(dRNEA_da_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_da_dot_vec, dimq_)));
  memorymanager::DeleteVector(dRNEA_dq_dot_vec);
  memorymanager::DeleteVector(dRNEA_dv_dot_vec);
  memorymanager::DeleteVector(dRNEA_da_dot_vec);
  memorymanager::DeleteVector(vec);
}


TEST_F(FixedBaseRobotTest, CRBADotVec) {
  Robot robot(urdf_file_name_);
  double *vec = memorymanager::NewVector(dimq_);
  Eigen::Map<Eigen::VectorXd>(vec, dimq_) = Eigen::VectorXd::Random(dimq_);
  double *result = memorymanager::NewVector(dimq_);
  robot.CRBADotVec(q_, vec, result);
  pinocchio::crba(model_, data_, Eigen::Map<const Eigen::VectorXd>(q_, dimq_));
  data_.M.triangularView<Eigen::StrictlyLower>() 
      = data_.M.transpose().triangularView<Eigen::StrictlyLower>();
  Eigen::VectorXd result_ref
      = data_.M * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  EXPECT_TRUE(result_ref.isApprox(Eigen::Map<Eigen::VectorXd>(result, dimq_)));
  memorymanager::DeleteVector(vec);
  memorymanager::DeleteVector(result);
}


TEST_F(FixedBaseRobotTest, RNEAWithFext) {
  Robot robot(urdf_file_name_);
  double *tau = memorymanager::NewVector(dimq_);
  double *tau0 = memorymanager::NewVector(dimq_);
  double *fext = memorymanager::NewVector(3*1);
  Eigen::Map<Eigen::VectorXd>(fext, 3*1) = Eigen::VectorXd::Random(3*1);
  robot.addPointContact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  EXPECT_EQ(robot.dimf(), 3*1);
  robot.RNEA(q_, v_, a_, tau0);
  robot.RNEA(q_, v_, a_, fext, tau);
  EXPECT_TRUE((Eigen::Map<Eigen::VectorXd>(tau0, dimq_)).isApprox(Eigen::Map<Eigen::VectorXd>(tau, dimq_)));
  robot.setFext(fext);
  robot.RNEA(q_, v_, a_, fext, tau);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint 
      = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  PointContact contact_ref(contact_frame_id_, baumgarte_alpha_, 
                           baumgarte_beta_);
  contact_ref.loadPinocchioModel(model_);
  contact_ref.contactForceToJointForce(fext, fjoint);
  Eigen::VectorXd tau_ref = pinocchio::rnea(
      model_, data_,
      Eigen::Map<const Eigen::VectorXd>(q_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(v_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(a_, dimq_),
      fjoint);
  EXPECT_TRUE(tau_ref.isApprox(Eigen::Map<Eigen::VectorXd>(tau, dimq_)));
  memorymanager::DeleteVector(tau);
  memorymanager::DeleteVector(tau0);
  memorymanager::DeleteVector(fext);
}


TEST_F(FixedBaseRobotTest, RNEADerivativesTransDotVecWithFext) {
  Robot robot(urdf_file_name_);
  double *vec = memorymanager::NewVector(dimq_);
  Eigen::Map<Eigen::VectorXd>(vec, dimq_) = Eigen::VectorXd::Random(dimq_);
  double *dRNEA_dq_dot_vec = memorymanager::NewVector(dimq_);
  double *dRNEA_dv_dot_vec = memorymanager::NewVector(dimq_);
  double *dRNEA_da_dot_vec = memorymanager::NewVector(dimq_);
  double *dRNEA_dfext_dot_vec = memorymanager::NewVector(3*1);
  double *dRNEA_dq_dot_vec0 = memorymanager::NewVector(dimq_);
  double *dRNEA_dv_dot_vec0 = memorymanager::NewVector(dimq_);
  double *dRNEA_da_dot_vec0 = memorymanager::NewVector(dimq_);
  double *fext = memorymanager::NewVector(3*1);
  Eigen::Map<Eigen::VectorXd>(fext, 3*1) = Eigen::VectorXd::Random(3*1);
  robot.addPointContact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  robot.updateKinematics(q_);
  robot.RNEADerivativesTransDotVector(q_, v_, a_, vec, dRNEA_dq_dot_vec0, 
                                      dRNEA_dv_dot_vec0, dRNEA_da_dot_vec0);
  robot.RNEADerivativesTransDotVector(q_, v_, a_, vec, dRNEA_dq_dot_vec, 
                                      dRNEA_dv_dot_vec, dRNEA_da_dot_vec, 
                                      dRNEA_dfext_dot_vec);
  EXPECT_TRUE(Eigen::Map<Eigen::VectorXd>(dRNEA_dq_dot_vec, dimq_).isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_dq_dot_vec0, dimq_)));
  EXPECT_TRUE(Eigen::Map<Eigen::VectorXd>(dRNEA_dv_dot_vec, dimq_).isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_dv_dot_vec0, dimq_)));
  EXPECT_TRUE(Eigen::Map<Eigen::VectorXd>(dRNEA_da_dot_vec, dimq_).isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_da_dot_vec0, dimq_)));
  robot.setFext(fext);
  robot.RNEADerivativesTransDotVector(q_, v_, a_, vec, dRNEA_dq_dot_vec, 
                                      dRNEA_dv_dot_vec, dRNEA_da_dot_vec, 
                                      dRNEA_dfext_dot_vec);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint 
      = pinocchio::container::aligned_vector<pinocchio::Force>(
                 model_.joints.size(), pinocchio::Force::Zero());
  PointContact contact_ref(contact_frame_id_, baumgarte_alpha_, 
                           baumgarte_beta_);
  contact_ref.loadPinocchioModel(model_);
  robot.updateKinematics(q_);
  contact_ref.contactForceToJointForce(fext, fjoint);
  pinocchio::computeRNEADerivatives(
      model_, data_,
      Eigen::Map<const Eigen::VectorXd>(q_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(v_, dimq_), 
      Eigen::Map<const Eigen::VectorXd>(a_, dimq_), fjoint);
  data_.M.triangularView<Eigen::StrictlyLower>() 
      = data_.M.transpose().triangularView<Eigen::StrictlyLower>();
  Eigen::VectorXd dRNEA_dq_dot_vec_ref 
      = data_.dtau_dq.transpose() 
          * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  Eigen::VectorXd dRNEA_dv_dot_vec_ref 
      = data_.dtau_dv.transpose() 
          * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  Eigen::VectorXd dRNEA_da_dot_vec_ref 
      = data_.M * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  Eigen::MatrixXd dtau_dfext = Eigen::MatrixXd::Zero(6, dimq_);
  pinocchio::forwardKinematics(model_, data_, 
                               Eigen::Map<const Eigen::VectorXd>(q_, dimq_));
  pinocchio::computeJointJacobians(model_, data_, 
                                   Eigen::Map<const Eigen::VectorXd>(q_, dimq_));
  pinocchio::updateFramePlacements(model_, data_);
  contact_ref.contactJacobian(model_, data_, dtau_dfext, 0, 0);
  Eigen::VectorXd dRNEA_dfext_dot_vec_ref 
      = dtau_dfext.topRows(3) * Eigen::Map<const Eigen::VectorXd>(vec, dimq_);
  EXPECT_TRUE(dRNEA_dq_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_dq_dot_vec, dimq_)));
  EXPECT_TRUE(dRNEA_dv_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_dv_dot_vec, dimq_)));
  EXPECT_TRUE(dRNEA_da_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_da_dot_vec, dimq_)));
  EXPECT_TRUE(dRNEA_dfext_dot_vec_ref.isApprox(Eigen::Map<Eigen::VectorXd>(dRNEA_dfext_dot_vec, 3)));

  memorymanager::DeleteVector(fext);
  memorymanager::DeleteVector(dRNEA_dq_dot_vec0);
  memorymanager::DeleteVector(dRNEA_dv_dot_vec0);
  memorymanager::DeleteVector(dRNEA_da_dot_vec0);
  memorymanager::DeleteVector(dRNEA_dq_dot_vec);
  memorymanager::DeleteVector(dRNEA_dv_dot_vec);
  memorymanager::DeleteVector(dRNEA_da_dot_vec);
  memorymanager::DeleteVector(dRNEA_dfext_dot_vec);
  memorymanager::DeleteVector(vec);
}


TEST_F(FixedBaseRobotTest, numContacts) {
  Robot robot(urdf_file_name_);
  int num_contact = 5;
  EXPECT_EQ(robot.dimf(), 0);
  for (int i=0; i<num_contact; ++i) {
    robot.addPointContact(2*i+2, 0, 0);
    EXPECT_EQ(robot.dimf(), (i+1)*3);
  }
  for (int i=0; i<num_contact; ++i) {
    robot.removePointContact(2*i+1);
    EXPECT_EQ(robot.dimf(), num_contact*3);
  }
  for (int i=0; i<num_contact; ++i) {
    robot.removePointContact(2*i+2);
    EXPECT_EQ(robot.dimf(), (num_contact-i-1)*3);
  }
}

} // namespace robot_cgmres


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}