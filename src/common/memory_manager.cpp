#include "common/memory_manager.hpp"


namespace robotcgmres {
namespace memorymanager {

double* NewVector(const int dim) {
  double* vec = new double[dim];
  for (int i=0; i<dim; ++i) {
    vec[i] = 0;
  }
  return vec;
}

void DeleteVector(double* vec) {
  delete[] vec;
}

double** NewMatrix(const int dim_row, const int dim_column) {
  double** mat = new double*[dim_row];
  mat[0] = new double[dim_row*dim_column];
  for (int i=1; i<dim_row; ++i) {
    mat[i] = mat[i-1] + dim_column;
  }
  for (int i=0; i<dim_row*dim_column; ++i) {
    mat[0][i] = 0;
  }
  return mat;
}

void DeleteMatrix(double** mat) {
  delete[] mat[0];
  delete[] mat;
}

} // namespace memorymanager
} // namespace robotcgmres