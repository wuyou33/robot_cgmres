#include <gtest/gtest.h>
#include <Eigen/Core>

#include "common/memory_manager.hpp"
#include "common/linear_algebra.hpp"


namespace robotcgmres {

class LinearAlgebraTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    dim = 50;
    vec1 = memorymanager::NewVector(dim);
    vec2 = memorymanager::NewVector(dim);
    Eigen::Map<Eigen::VectorXd>(vec1, dim) = Eigen::VectorXd::Random(dim);
    Eigen::Map<Eigen::VectorXd>(vec2, dim) = Eigen::VectorXd::Random(dim);
  }

  virtual void TearDown() {
    memorymanager::DeleteVector(vec1);
    memorymanager::DeleteVector(vec2);
  }

  int dim;
  double *vec1, *vec2;
};


TEST_F(LinearAlgebraTest, InnerProduct) {
  double ref = 0;
  for (int i=0; i<dim; ++i) {
    ref += vec1[i] * vec2[i];
  }
  double res = linearalgebra::InnerProduct(dim, vec1, vec2);
  EXPECT_DOUBLE_EQ(res, ref);
}

TEST_F(LinearAlgebraTest, SquaredNorm) {
  double ref = 0;
  for (int i=0; i<dim; ++i) {
    ref += vec1[i] * vec1[i];
  }
  double res = linearalgebra::SquaredNorm(dim, vec1);
  EXPECT_DOUBLE_EQ(res, ref);
}

} // namespace robotcgmres


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}