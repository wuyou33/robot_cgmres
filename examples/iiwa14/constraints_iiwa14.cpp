#include "constraints_iiwa14.hpp"


namespace robotcgmres {

Constraints_iiwa14::Constraints_iiwa14(const Robot* robot_ptr) 
  : ConstraintsInterface() {
}

Constraints_iiwa14::~Constraints_iiwa14() {
}

void Constraints_iiwa14::residual(const Robot* robot_ptr, const double t, 
                                  const double* q, const double* v, 
                                  const double* a, const double* u, 
                                  const double* f, double* residual) {
}

void Constraints_iiwa14::addCqDotVec(const Robot* robot_ptr, const double t, 
                                     const double* q, const double* v, 
                                     const double* a, const double* u, 
                                     const double* f, const double* vec, 
                                     double* hq) {

}

void Constraints_iiwa14::addCvDotVec(const Robot* robot_ptr, const double t, 
                                     const double* q, const double* v, 
                                     const double* a, const double* u, 
                                     const double* f, const double* vec, 
                                     double* hv) {

}

void Constraints_iiwa14::addCaDotVec(const Robot* robot_ptr, const double t, 
                                     const double* q, const double* v, 
                                     const double* a, const double* u, 
                                     const double* f, const double* vec, 
                                     double* ha) {

}

void Constraints_iiwa14::addCuDotVec(const Robot* robot_ptr, const double t, 
                                     const double* q, const double* v, 
                                     const double* a, const double* u, 
                                     const double* f, const double* vec, 
                                     double* hu) {

}

void Constraints_iiwa14::addCfDotVec(const Robot* robot_ptr, const double t, 
                                     const double* q, const double* v, 
                                     const double* a, const double* u, 
                                     const double* f, const double* vec, 
                                     double* hf) {

}

void Constraints_iiwa14::residual(const Robot* robot_ptr, const double t, 
                                  const double* q, const double* v, 
                                  const double* a, const double* u, 
                                  double* residual) {

}

void Constraints_iiwa14::addCqDotVec(const Robot* robot_ptr, const double t, 
                                     const double* q, const double* v, 
                                     const double* a, const double* u, 
                                     const double* vec, double* hq) {
}

void Constraints_iiwa14::addCvDotVec(const Robot* robot_ptr, const double t, 
                                     const double* q, const double* v, 
                                     const double* a, const double* u, 
                                     const double* vec, double* hv) {
}

void Constraints_iiwa14::addCaDotVec(const Robot* robot_ptr, const double t, 
                                     const double* q, const double* v, 
                                     const double* a, const double* u, 
                                     const double* vec, double* ha) {
}

void Constraints_iiwa14::addCuDotVec(const Robot* robot_ptr, const double t, 
                                     const double* q, const double* v, 
                                     const double* a, const double* u, 
                                     const double* vec, double* hu) {
}

void Constraints_iiwa14::addCfDotVec(const Robot* robot_ptr, const double t, 
                                     const double* q, const double* v, 
                                     const double* a, const double* u, 
                                     const double* vec, double* hf) {
}

int Constraints_iiwa14::dim_constraints() const {
  return 0;
}

} // namespace robotcgmres