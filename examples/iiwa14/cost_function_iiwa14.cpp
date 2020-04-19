#include "cost_function_iiwa14.hpp"


namespace robotcgmres {

CostFunction_iiwa14::CostFunction_iiwa14(const Robot* robot_ptr) 
  : CostFunctionInterface(),
    joint_space_cost_(robot_ptr) {
  std::vector<int> all_joints = {0, 1, 2, 3, 4, 5, 6};
  joint_space_cost_.set_q_weight(all_joints, 1.0);
  joint_space_cost_.set_v_weight(all_joints, 0.1);
  joint_space_cost_.set_a_weight(all_joints, 0.01);
  joint_space_cost_.set_u_weight(all_joints, 0.001);
}

CostFunction_iiwa14::~CostFunction_iiwa14() {
}

void CostFunction_iiwa14::lq(const Robot* robot_ptr, const double t, 
                             const double* q, const double* v, const double* a, 
                             const double* u, const double* f, double* lq) {
  joint_space_cost_.lq(q, lq);
}

void CostFunction_iiwa14::lv(const Robot* robot_ptr, const double t, 
                             const double* q, const double* v, const double* a, 
                             const double* u, const double* f, double* lv) {
  joint_space_cost_.lv(v, lv);
}

void CostFunction_iiwa14::la(const Robot* robot_ptr, const double t, 
                             const double* q, const double* v, const double* a, 
                             const double* u, const double* f, double* la) {
  joint_space_cost_.la(a, la);
}

void CostFunction_iiwa14::lu(const Robot* robot_ptr, const double t, 
                             const double* q, const double* v, const double* a, 
                             const double* u, const double* f, double* lu) {
  joint_space_cost_.lu(u, lu);
}

void CostFunction_iiwa14::lf(const Robot* robot_ptr, const double t, 
                             const double* q, const double* v, const double* a, 
                             const double* u, const double* f, double* lf) {
}

void CostFunction_iiwa14::lq(const Robot* robot_ptr, const double t, 
                             const double* q, const double* v, const double* a, 
                             const double* u, double* lq) {
  joint_space_cost_.lq(q, lq);
}

void CostFunction_iiwa14::lv(const Robot* robot_ptr, const double t, 
                             const double* q, const double* v, const double* a, 
                             const double* u, double* lv) {
  joint_space_cost_.lv(v, lv);
}

void CostFunction_iiwa14::la(const Robot* robot_ptr, const double t, 
                             const double* q, const double* v, const double* a, 
                             const double* u, double* la) {
  joint_space_cost_.la(a, la);
}

void CostFunction_iiwa14::lu(const Robot* robot_ptr, const double t, 
                             const double* q, const double* v, const double* a, 
                             const double* u, double* lu) {
  joint_space_cost_.lu(u, lu);
}

void CostFunction_iiwa14::phiq(const Robot* robot_ptr, const double t, 
                               const double* q, const double* v, double* phiq) {
  joint_space_cost_.lq(q, phiq);
}

void CostFunction_iiwa14::phiv(const Robot* robot_ptr, const double t, 
                               const double* q, const double* v, double* phiv) {
  joint_space_cost_.lv(v, phiv);
}

} // namespace robotcgmres