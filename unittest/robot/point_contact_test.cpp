#include <string>
#include <random>
#include <utility>
#include <vector>

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


namespace robotcgmres {

class PointContactTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    std::random_device rnd;
    urdf_file_name_ = "../../../examples/iiwa14/iiwa14.urdf";
    pinocchio::urdf::buildModel(urdf_file_name_, model_);
    data_ = pinocchio::Data(model_);
    dimq_ = model_.nq;
    dimv_ = model_.nv;
    njoints_ = model_.njoints;
    contact_frame_id_ = 2 + (int)((rnd()%(model_.nframes-2)) / 2) * 2;
    q_ = pinocchio::randomConfiguration(model_, -Eigen::VectorXd::Ones(dimq_), 
                                         Eigen::VectorXd::Ones(dimq_));
    v_ = Eigen::VectorXd::Random(dimv_);
    a_ = Eigen::VectorXd::Random(dimv_);
    baumgarte_alpha_ = Eigen::VectorXd::Random(2)[1];
    baumgarte_beta_ = Eigen::VectorXd::Random(2)[1];
  }

  virtual void TearDown() {
  }

  std::string urdf_file_name_;
  pinocchio::Model model_;
  pinocchio::Data data_;
  int dimq_, dimv_, njoints_, contact_frame_id_;
  double baumgarte_alpha_, baumgarte_beta_;
  Eigen::VectorXd q_, v_, a_;
};


TEST_F(PointContactTest, parameters) {
  PointContact contact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  contact.loadPinocchioModel(model_);

  // Check the constructor works well.
  EXPECT_EQ(contact.contact_frame_id(), contact_frame_id_);
  EXPECT_EQ(contact.parent_joint_id(), model_.frames[contact_frame_id_].parent);
  EXPECT_EQ(contact.dimv(), model_.nv);
  EXPECT_DOUBLE_EQ(contact.baumgarte_alpha(), baumgarte_alpha_);
  EXPECT_DOUBLE_EQ(contact.baumgarte_beta(), baumgarte_beta_);
  EXPECT_TRUE(contact.contact_point().isApprox(Eigen::Vector3d::Zero()));
  EXPECT_TRUE(contact.jXf().isApprox(model_.frames[contact_frame_id_].placement));

  // Check reset parameters works well.
  double alpha_tmp = Eigen::VectorXd::Random(2)[1];
  double beta_tmp = Eigen::VectorXd::Random(2)[1];
  Eigen::Vector3d contact_point = Eigen::Vector3d::Random();
  contact.resetBaugrarteParameters(alpha_tmp, beta_tmp);
  contact.resetContactPoint(contact_point);
  EXPECT_EQ(alpha_tmp, contact.baumgarte_alpha());
  EXPECT_EQ(beta_tmp, contact.baumgarte_beta());
  EXPECT_TRUE(contact_point.isApprox(contact.contact_point()));

  // Check the contact point assignment by kinematics.
  pinocchio::forwardKinematics(model_, data_, q_);
  pinocchio::updateFramePlacement(model_, data_, contact_frame_id_);
  contact.resetContactPointByCurrentState(data_);
  EXPECT_TRUE(contact.contact_point().isApprox(data_.oMf[contact_frame_id_].translation()));
}


TEST_F(PointContactTest, moveConstructor) {
  PointContact contact1(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  contact1.loadPinocchioModel(model_);

  // Ckeck the contact1 parameters.
  EXPECT_EQ(contact1.contact_frame_id(), contact_frame_id_);
  EXPECT_EQ(contact1.parent_joint_id(), model_.frames[contact_frame_id_].parent);
  EXPECT_EQ(contact1.dimv(), model_.nv);
  EXPECT_DOUBLE_EQ(contact1.baumgarte_alpha(), baumgarte_alpha_);
  EXPECT_DOUBLE_EQ(contact1.baumgarte_beta(), baumgarte_beta_);
  EXPECT_TRUE(contact1.contact_point().isApprox(Eigen::Vector3d::Zero()));
  EXPECT_TRUE(contact1.jXf().isApprox(model_.frames[contact_frame_id_].placement));

  PointContact contact2(std::move(contact1));
  // Ckeck that the contact2 parameters are these was assinged in the contact1.
  EXPECT_EQ(contact2.contact_frame_id(), contact_frame_id_);
  EXPECT_EQ(contact2.parent_joint_id(), model_.frames[contact_frame_id_].parent);
  EXPECT_EQ(contact2.dimv(), model_.nv);
  EXPECT_DOUBLE_EQ(contact2.baumgarte_alpha(), baumgarte_alpha_);
  EXPECT_DOUBLE_EQ(contact2.baumgarte_beta(), baumgarte_beta_);
  EXPECT_TRUE(contact2.contact_point().isApprox(Eigen::Vector3d::Zero()));
  EXPECT_TRUE(contact2.jXf().isApprox(model_.frames[contact_frame_id_].placement));
}


TEST_F(PointContactTest, vector) {
  std::vector<PointContact> contacts;
  // Check the size at default.
  EXPECT_TRUE(contacts.empty());
  EXPECT_EQ(contacts.size(), 0);
  contacts.push_back(PointContact(contact_frame_id_, baumgarte_alpha_, 
                                  baumgarte_beta_));
  // Check the size when an element is appended.
  EXPECT_FALSE(contacts.empty());
  EXPECT_EQ(contacts.size(), 1);

  contacts.push_back(PointContact((int)(contact_frame_id_/2), 
                                  0.5*baumgarte_alpha_, 
                                  0.5*baumgarte_beta_));
  // Check each contact's parameters.
  EXPECT_EQ(contacts.size(), 2);
  EXPECT_EQ(contacts[0].contact_frame_id(), contact_frame_id_);
  EXPECT_EQ(contacts[0].baumgarte_alpha(), baumgarte_alpha_);
  EXPECT_EQ(contacts[0].baumgarte_beta(), baumgarte_beta_);
  EXPECT_EQ(contacts[1].contact_frame_id(), (int)(contact_frame_id_/2));
  EXPECT_EQ(contacts[1].baumgarte_alpha(), 0.5*baumgarte_alpha_);
  EXPECT_EQ(contacts[1].baumgarte_beta(), 0.5*baumgarte_beta_);
  int remove_idx = 0;
  std::swap<PointContact>(contacts[remove_idx], contacts.back());
  contacts.pop_back();
  // Check whether swap-and-pop_back works correctly.
  EXPECT_EQ(contacts.size(), 1);
  EXPECT_EQ(contacts[0].contact_frame_id(), (int)(contact_frame_id_/2));
  EXPECT_EQ(contacts[0].baumgarte_alpha(), 0.5*baumgarte_alpha_);
  EXPECT_EQ(contacts[0].baumgarte_beta(), 0.5*baumgarte_beta_);

  contacts.push_back(PointContact(contact_frame_id_, baumgarte_alpha_, 
                                  baumgarte_beta_));
  EXPECT_EQ(contacts.size(), 2);
  // Check each contact's parameters.
  EXPECT_EQ(contacts[0].contact_frame_id(), (int)(contact_frame_id_/2));
  EXPECT_EQ(contacts[0].baumgarte_alpha(), 0.5*baumgarte_alpha_);
  EXPECT_EQ(contacts[0].baumgarte_beta(), 0.5*baumgarte_beta_);
  EXPECT_EQ(contacts[1].contact_frame_id(), contact_frame_id_);
  EXPECT_EQ(contacts[1].baumgarte_alpha(), baumgarte_alpha_);
  EXPECT_EQ(contacts[1].baumgarte_beta(), baumgarte_beta_);
  remove_idx = 1;
  contacts.erase(contacts.begin()+remove_idx);
  // Check whether erase works correctly.
  EXPECT_EQ(contacts.size(), 1);
  EXPECT_EQ(contacts[0].contact_frame_id(), (int)(contact_frame_id_/2));
  EXPECT_EQ(contacts[0].baumgarte_alpha(), 0.5*baumgarte_alpha_);
  EXPECT_EQ(contacts[0].baumgarte_beta(), 0.5*baumgarte_beta_);

  contacts.push_back(PointContact(contact_frame_id_, baumgarte_alpha_, 
                                  baumgarte_beta_));
  contacts.push_back(PointContact(contact_frame_id_, baumgarte_alpha_, 
                                  baumgarte_beta_));
  // Check whether erase works correctly.
  EXPECT_EQ(contacts.size(), 3);
  contacts.erase(contacts.begin()+remove_idx, contacts.end());
  EXPECT_EQ(contacts.size(), 1);
  contacts.pop_back();
  EXPECT_TRUE(contacts.empty());
}


TEST_F(PointContactTest, contactForceToJointForce) {
  PointContact contact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  contact.loadPinocchioModel(model_);
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint
      = pinocchio::container::aligned_vector<pinocchio::Force>(
            model_.joints.size(), pinocchio::Force::Zero());
  pinocchio::container::aligned_vector<pinocchio::Force> fjoint_ref = fjoint;
  double *fext = memorymanager::NewVector(3);
  Eigen::Map<Eigen::VectorXd>(fext, 3) = Eigen::Vector3d::Random();
  contact.contactForceToJointForce(fext, fjoint);
  const int parent_joint_id = contact.parent_joint_id();
  fjoint_ref[parent_joint_id] 
      = model_.frames[contact_frame_id_].placement.act(
          pinocchio::Force(Eigen::Map<const Eigen::Vector3d>(fext), 
                           Eigen::Vector3d::Zero()));
  EXPECT_TRUE(fjoint[parent_joint_id].isApprox(fjoint_ref[parent_joint_id]));
  memorymanager::DeleteVector(fext);
}


TEST_F(PointContactTest, contactJacobianAndItsBlock) {
  PointContact contact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  contact.loadPinocchioModel(model_);
  pinocchio::forwardKinematics(model_, data_, q_);
  pinocchio::computeJointJacobians(model_, data_, q_);

  Eigen::MatrixXd J_contact = Eigen::MatrixXd::Zero(6, dimv_);
  Eigen::MatrixXd J_contact_ref = Eigen::MatrixXd::Zero(6, dimv_);
  contact.contactJacobian(model_, data_, J_contact);
  pinocchio::getFrameJacobian(model_, data_, contact_frame_id_, 
                              pinocchio::LOCAL, J_contact_ref);
  EXPECT_TRUE(J_contact.isApprox(J_contact_ref));

  Eigen::MatrixXd J_contact_block = Eigen::MatrixXd::Zero(2*6, 2*dimv_);
  Eigen::MatrixXd J_contact_block_ref = Eigen::MatrixXd::Zero(2*6, 2*dimv_);
  const int row_begin = rand() % dimv_;
  const int column_begin = rand() % 6;
  contact.contactJacobian(model_, data_, J_contact_block, row_begin, 
                          column_begin);
  pinocchio::getFrameJacobian(model_, data_, contact_frame_id_, 
                              pinocchio::LOCAL, 
                              J_contact_block_ref.block(row_begin, column_begin, 
                                                        6, dimv_));
  EXPECT_TRUE(J_contact_block.isApprox(J_contact_block_ref));
  EXPECT_TRUE(J_contact.isApprox(J_contact_block.block(row_begin, column_begin, 
                                                       6, dimv_)));
}


TEST_F(PointContactTest, baumgarteResidual) {
  PointContact contact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  contact.loadPinocchioModel(model_);
  pinocchio::forwardKinematics(model_, data_, q_, v_, a_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q_, v_, a_);
  pinocchio::updateFramePlacements(model_, data_);
  const int parent_joint_id = contact.parent_joint_id();
  double residual[3], residual_ref[3];
  contact.baumgarteResidual(model_, data_, residual);
  Eigen::Map<Eigen::VectorXd>(residual_ref, 3) 
      = (data_.oMi[parent_joint_id].act(data_.a[parent_joint_id])).linear()
          + baumgarte_alpha_ 
              * (data_.oMi[parent_joint_id].act(data_.v[parent_joint_id])).linear()
          + baumgarte_beta_  
              * (data_.oMf[contact_frame_id_].translation()-contact.contact_point());
  for (int i=0; i<3; ++i) {
    EXPECT_DOUBLE_EQ(residual[i], residual_ref[i]);
  } 
}


TEST_F(PointContactTest, baumgarteDerivatives) {
  PointContact contact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  contact.loadPinocchioModel(model_);
  pinocchio::forwardKinematics(model_, data_, q_, v_, a_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q_, v_, a_);
  pinocchio::updateFramePlacement(model_, data_, contact_frame_id_);
  const int parent_joint_id = contact.parent_joint_id();
  Eigen::MatrixXd baum_partial_dq = Eigen::MatrixXd::Zero(3, dimv_);
  Eigen::MatrixXd baum_partial_dv = Eigen::MatrixXd::Zero(3, dimv_);
  Eigen::MatrixXd baum_partial_da = Eigen::MatrixXd::Zero(3, dimv_);
  contact.baumgarteDerivatives(model_, data_, baum_partial_dq, baum_partial_dv, 
                               baum_partial_da);
  Eigen::MatrixXd baum_partial_dq_ref = Eigen::MatrixXd::Zero(3, dimv_);
  Eigen::MatrixXd baum_partial_dv_ref = Eigen::MatrixXd::Zero(3, dimv_);
  Eigen::MatrixXd baum_partial_da_ref = Eigen::MatrixXd::Zero(3, dimv_);
  Eigen::MatrixXd J_frame = Eigen::MatrixXd::Zero(6, dimv_);
  Eigen::MatrixXd joint_v_partial_dq = Eigen::MatrixXd::Zero(6, dimv_);
  Eigen::MatrixXd joint_a_partial_dq = Eigen::MatrixXd::Zero(6, dimv_);
  Eigen::MatrixXd joint_a_partial_dv = Eigen::MatrixXd::Zero(6, dimv_);
  Eigen::MatrixXd joint_a_partial_da = Eigen::MatrixXd::Zero(6, dimv_);
  pinocchio::getJointAccelerationDerivatives(model_, data_, parent_joint_id, 
                                             pinocchio::WORLD,
                                             joint_v_partial_dq, 
                                             joint_a_partial_dq, 
                                             joint_a_partial_dv, 
                                             joint_a_partial_da);
  pinocchio::getFrameJacobian(model_, data_, contact_frame_id_, 
                              pinocchio::LOCAL_WORLD_ALIGNED, J_frame);
  baum_partial_dq_ref
      = joint_a_partial_dq.template topRows<3>()
          + baumgarte_alpha_ * joint_v_partial_dq.template topRows<3>()
          + baumgarte_beta_ * J_frame.template topRows<3>();
  baum_partial_dv_ref
      = joint_a_partial_dv.template topRows<3>()
          + baumgarte_alpha_ * joint_a_partial_da.template topRows<3>();
  baum_partial_da_ref
      = joint_a_partial_da.template topRows<3>();
  EXPECT_TRUE(baum_partial_dq_ref.isApprox(baum_partial_dq));
  EXPECT_TRUE(baum_partial_dv_ref.isApprox(baum_partial_dv));
  EXPECT_TRUE(baum_partial_da_ref.isApprox(baum_partial_da));
}


TEST_F(PointContactTest, baumgarteDerivativesBlock) {
  PointContact contact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  contact.loadPinocchioModel(model_);
  pinocchio::forwardKinematics(model_, data_, q_, v_, a_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q_, v_, a_);
  pinocchio::updateFramePlacement(model_, data_, contact_frame_id_);
  const int parent_joint_id = contact.parent_joint_id();
  Eigen::MatrixXd baum_partial_dq = Eigen::MatrixXd::Zero(2*3, 2*dimv_);
  Eigen::MatrixXd baum_partial_dv = Eigen::MatrixXd::Zero(2*3, 2*dimv_);
  Eigen::MatrixXd baum_partial_da = Eigen::MatrixXd::Zero(2*3, 2*dimv_);
  const int row_begin = rand() % dimv_;
  const int column_begin = rand() % 3;
  contact.baumgarteDerivatives(model_, data_, baum_partial_dq, baum_partial_dv, 
                               baum_partial_da, row_begin, 
                               column_begin);
  Eigen::MatrixXd baum_partial_dq_ref = Eigen::MatrixXd::Zero(2*3, 2*dimv_);
  Eigen::MatrixXd baum_partial_dv_ref = Eigen::MatrixXd::Zero(2*3, 2*dimv_);
  Eigen::MatrixXd baum_partial_da_ref = Eigen::MatrixXd::Zero(2*3, 2*dimv_);
  Eigen::MatrixXd J_frame = Eigen::MatrixXd::Zero(6, dimv_);
  Eigen::MatrixXd joint_v_partial_dq = Eigen::MatrixXd::Zero(6, dimv_);
  Eigen::MatrixXd joint_a_partial_dq = Eigen::MatrixXd::Zero(6, dimv_);
  Eigen::MatrixXd joint_a_partial_dv = Eigen::MatrixXd::Zero(6, dimv_);
  Eigen::MatrixXd joint_a_partial_da = Eigen::MatrixXd::Zero(6, dimv_);
  pinocchio::getJointAccelerationDerivatives(model_, data_, parent_joint_id, 
                                             pinocchio::WORLD,
                                             joint_v_partial_dq, 
                                             joint_a_partial_dq, 
                                             joint_a_partial_dv, 
                                             joint_a_partial_da);
  pinocchio::getFrameJacobian(model_, data_, contact_frame_id_, 
                              pinocchio::LOCAL_WORLD_ALIGNED, J_frame);
  baum_partial_dq_ref.block(row_begin, column_begin, 3, dimv_)
      = joint_a_partial_dq.template topRows<3>()
          + baumgarte_alpha_ * joint_v_partial_dq.template topRows<3>()
          + baumgarte_beta_ * J_frame.template topRows<3>();
  baum_partial_dv_ref.block(row_begin, column_begin, 3, dimv_)
      = joint_a_partial_dv.template topRows<3>()
          + baumgarte_alpha_ * joint_a_partial_da.template topRows<3>();
  baum_partial_da_ref.block(row_begin, column_begin, 3, dimv_)
      = joint_a_partial_da.template topRows<3>();
  EXPECT_TRUE(baum_partial_dq_ref.isApprox(baum_partial_dq));
  EXPECT_TRUE(baum_partial_dv_ref.isApprox(baum_partial_dv));
  EXPECT_TRUE(baum_partial_da_ref.isApprox(baum_partial_da));
}


TEST_F(PointContactTest, baumgarteDerivativesAndItsBlock) {
  PointContact contact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  contact.loadPinocchioModel(model_);
  pinocchio::forwardKinematics(model_, data_, q_, v_, a_);
  pinocchio::computeForwardKinematicsDerivatives(model_, data_, q_, v_, a_);
  pinocchio::updateFramePlacement(model_, data_, contact_frame_id_);
  const int parent_joint_id = contact.parent_joint_id();
  Eigen::MatrixXd baum_partial_dq = Eigen::MatrixXd::Zero(3, dimv_);
  Eigen::MatrixXd baum_partial_dv = Eigen::MatrixXd::Zero(3, dimv_);
  Eigen::MatrixXd baum_partial_da = Eigen::MatrixXd::Zero(3, dimv_);
  contact.baumgarteDerivatives(model_, data_, baum_partial_dq, baum_partial_dv, 
                               baum_partial_da);
  Eigen::MatrixXd baum_partial_dq_block = Eigen::MatrixXd::Zero(2*3, 2*dimv_);
  Eigen::MatrixXd baum_partial_dv_block = Eigen::MatrixXd::Zero(2*3, 2*dimv_);
  Eigen::MatrixXd baum_partial_da_block = Eigen::MatrixXd::Zero(2*3, 2*dimv_);
  const int row_begin = rand() % 3;
  const int column_begin = rand() % dimv_;
  contact.baumgarteDerivatives(model_, data_, baum_partial_dq_block, 
                               baum_partial_dv_block, baum_partial_da_block, 
                               row_begin, column_begin);
  EXPECT_TRUE(baum_partial_dq.
              isApprox(baum_partial_dq_block.block(row_begin, 
                                                   column_begin, 
                                                   3, dimv_)));
  EXPECT_TRUE(baum_partial_dv.
              isApprox(baum_partial_dv_block.block(row_begin, 
                                                   column_begin, 
                                                   3, dimv_)));
  EXPECT_TRUE(baum_partial_da.
              isApprox(baum_partial_da_block.block(row_begin, 
                                                   column_begin, 
                                                   3, dimv_)));
}

} // namespace robot_cgmres


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}