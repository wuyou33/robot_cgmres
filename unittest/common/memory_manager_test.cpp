#include <random>

#include <gtest/gtest.h>

#include "common/memory_manager.hpp"


namespace robotcgmres {
namespace memorymanager {

class MemoryManagerTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    int max_dim = 100;
    std::random_device rnd;
    dim = rnd() % max_dim + 1;
  }

  int dim;
};


TEST_F(MemoryManagerTest, Vector) {
  double *vec = NewVector(dim);
  for (int i=0; i<dim; ++i) {
    EXPECT_DOUBLE_EQ(vec[i], 0);
  }
  DeleteVector(vec);
}


TEST_F(MemoryManagerTest, Matrix) {
  double **mat = NewMatrix(dim, dim);
  for (int i=0; i<dim; ++i) {
    for (int j=0; j<dim; ++j) {
      EXPECT_DOUBLE_EQ(mat[i][j], 0);
    }
  }
  DeleteMatrix(mat);
}

} // namespace memorymanager
} // namespace robotcgmres


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}