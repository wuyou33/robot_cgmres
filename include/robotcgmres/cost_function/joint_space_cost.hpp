#ifndef ROBOTCGMRES_JOINT_SPACE_COST_HPP_
#define ROBOTCGMRES_JOINT_SPACE_COST_HPP_

#include <vector>

#include "robot/robot.hpp"
#include "common/memory_manager.hpp"


namespace robotcgmres {

class JointSpaceCost {
public:

  JointSpaceCost(const Robot* robot_ptr);

  ~JointSpaceCost();

  void lq(const double* q, double* lq) const;

  void lv(const double* v, double* lv) const;

  void la(const double* a, double* la) const;

  void lu(const double* u, double* lu) const;

  void set_q_ref(const std::vector<int>& joint_indices, const double q_ref);

  void set_v_ref(const std::vector<int>& joint_indices, const double q_ref);

  void set_a_ref(const std::vector<int>& joint_indices, const double q_ref);

  void set_u_ref(const std::vector<int>& joint_indices, const double q_ref);

  void set_q_weight(const std::vector<int>& joint_indices, 
                    const double q_weight);

  void set_v_weight(const std::vector<int>& joint_indices, 
                    const double v_weight);

  void set_a_weight(const std::vector<int>& joint_indices, 
                    const double a_weight);

  void set_u_weight(const std::vector<int>& joint_indices, 
                    const double u_weight);

  // Prohibits copy constructor.
  JointSpaceCost(const JointSpaceCost&) = delete;

  // Prohibits copy operator.
  JointSpaceCost& operator=(const JointSpaceCost&) = delete;

private:
  int dimq_, dimv_, dimf_;
  double *q_ref_, *v_ref_, *a_ref_, *u_ref_, *f_ref_;
  double *wq_, *wv_, *wa_, *wu_, *wf_;
};
  
} // namespace robotcgmres


#endif // ROBOTCGMRES_JOINT_SPACE_COST_HPP_