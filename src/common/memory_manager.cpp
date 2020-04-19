#include "common/memory_manager.hpp"


namespace robotcgmres {
namespace memorymanager {

double* NewVector(const int dim) {
  if (dim > 0) {
    double* vec = new double[dim];
    for (int i=0; i<dim; ++i) {
      vec[i] = 0;
    }
    return vec;
  } 
  else {
    return nullptr;
  }
}

void DeleteVector(double* vec) {
  delete[] vec;
}

double** NewMatrix(const int dim_row, const int dim_column) {
  if (dim_row > 0 && dim_column > 0) {
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
  else {
    return nullptr;
  }
}

void DeleteMatrix(double** mat) {
  delete[] mat[0];
  delete[] mat;
}

} // namespace memorymanager
} // namespace robotcgmres