#include <gtest/gtest.h>
#include <time.h>
#include <string>
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
    urdf_file_name_ = "../../../examples/iiwa14/iiwa14.urdf";
    pinocchio::urdf::buildModel(urdf_file_name_, model_);
    contact_frame_id_ = 6;
    data_ = pinocchio::Data(model_);
    dimq_ = model_.nq;
    dimv_ = model_.nv;
    q0_ = pinocchio::randomConfiguration(model_, -Eigen::VectorXd::Ones(dimq_), 
                                         Eigen::VectorXd::Ones(dimq_));
    v0_ = Eigen::VectorXd::Random(dimv_);
    a0_ = Eigen::VectorXd::Random(dimv_);
    q_ = memorymanager::NewVector(dimq_);
    v_ = memorymanager::NewVector(dimv_);
    a_ = memorymanager::NewVector(dimv_);
    Eigen::Map<Eigen::VectorXd>(q_, dimq_)
        = pinocchio::randomConfiguration(model_, -Eigen::VectorXd::Ones(dimq_), 
                                         Eigen::VectorXd::Ones(dimq_));
    Eigen::Map<Eigen::VectorXd>(v_, dimv_) = Eigen::VectorXd::Random(dimv_);
    Eigen::Map<Eigen::VectorXd>(a_, dimv_) = Eigen::VectorXd::Random(dimv_);
  }

  virtual void TearDown() {
    memorymanager::DeleteVector(q_);
    memorymanager::DeleteVector(v_);
    memorymanager::DeleteVector(a_);
  }

  std::string urdf_file_name_;
  pinocchio::Model model_;
  pinocchio::Data data_;
  int dimq_, dimv_, contact_frame_id_;
  double baumgarte_alpha_, baumgarte_beta_;
  double *q_, *v_, *a_;
  Eigen::VectorXd q0_, v0_, a0_;
};


TEST_F(PointContactTest, parameters) {
  PointContact contact(contact_frame_id_, baumgarte_alpha_, baumgarte_beta_);
  contact.loadPinocchioModel(model_);

  EXPECT_EQ(contact_frame_id_, contact.contact_frame_id());
  EXPECT_EQ(model_.frames[contact_frame_id_].parent, contact.parent_joint_id());
  EXPECT_EQ(baumgarte_alpha_, contact.baumgarte_alpha());
  EXPECT_EQ(baumgarte_beta_, contact.baumgarte_beta());

  double alpha_tmp = Eigen::VectorXd::Random(2)[1];
  double beta_tmp = Eigen::VectorXd::Random(2)[1];
  Eigen::Vector3d contact_point = Eigen::Vector3d::Random();
  contact.resetBaugrarteParameters(alpha_tmp, beta_tmp);
  contact.resetContactPoint(contact_point);
  EXPECT_EQ(alpha_tmp, contact.baumgarte_alpha());
  EXPECT_EQ(beta_tmp, contact.baumgarte_beta());
  EXPECT_TRUE(contact_point.isApprox(contact.contact_point()));

  pinocchio::forwardKinematics(model_, data_, q0_);
  contact.resetContactPointByCurrentState(data_);
  EXPECT_TRUE(contact.contact_point().isApprox(data_.oMf[contact_frame_id_].translation()));
}

} // namespace robot_cgmres


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}