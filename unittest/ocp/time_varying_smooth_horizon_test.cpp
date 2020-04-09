#include <gtest/gtest.h>
#include <time.h>
#include <Eigen/Core>

#include "ocp/time_varying_smooth_horizon.hpp"


namespace robotcgmres {

class TimeVaryingSmoothHorizonTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand((unsigned int) time(0));
    T_f = std::abs(Eigen::VectorXd::Random(2)[1]);
    alpha = std::abs(Eigen::VectorXd::Random(2)[1]);
    initial_time = std::abs(Eigen::VectorXd::Random(2)[1]); 
  }

  double T_f, alpha, initial_time;
};


TEST_F(TimeVaryingSmoothHorizonTest, getLength) {
  const double current_time = std::abs(Eigen::VectorXd::Random(2)[1]) + initial_time;
  const double expect_length = T_f * (1.0 - std::exp(-alpha*(current_time - initial_time)));
  TimeVaryingSmoothHorizon horizon(T_f, alpha, initial_time);
  EXPECT_DOUBLE_EQ(horizon.getLength(current_time), expect_length);
  TimeVaryingSmoothHorizon horizon0(T_f, alpha);
  const double expect_length0 = T_f * (1.0 - std::exp(-alpha*current_time));
  EXPECT_DOUBLE_EQ(horizon0.getLength(current_time), expect_length0);
}

TEST_F(TimeVaryingSmoothHorizonTest, resetLength) {
  double current_time = std::abs(Eigen::VectorXd::Random(2)[1]) + initial_time;
  TimeVaryingSmoothHorizon horizon(T_f, alpha, initial_time);

  double reset_time = std::abs(Eigen::VectorXd::Random(2)[1]);
  horizon.resetLength(reset_time);
  double expect_length = T_f * (1.0 - std::exp(-alpha*(current_time - reset_time)));
  EXPECT_DOUBLE_EQ(horizon.getLength(current_time), expect_length);

  double reset_T = std::abs(Eigen::VectorXd::Random(2)[1]);
  double reset_alpha = std::abs(Eigen::VectorXd::Random(2)[1]);
  reset_time = std::abs(Eigen::VectorXd::Random(2)[1]);
  horizon.resetLength(reset_T, reset_alpha, reset_time);
  current_time = std::abs(Eigen::VectorXd::Random(2)[1]) + reset_time;
  expect_length = reset_T * (1.0 - std::exp(-reset_alpha*(current_time - reset_time)));
  EXPECT_DOUBLE_EQ(horizon.getLength(current_time), expect_length);
}

} // namespace robot_cgmres


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}