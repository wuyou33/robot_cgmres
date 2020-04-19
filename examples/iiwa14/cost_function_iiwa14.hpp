#ifndef ROBOTCGMRES_COST_FUNCTION_IIWA14_HPP_
#define ROBOTCGMRES_COST_FUNCTION_IIWA14_HPP_

#include <vector>

#include "cost_function/cost_function_interface.hpp"
#include "cost_function/joint_space_cost.hpp"
#include "robot/robot.hpp"


namespace robotcgmres {

class CostFunction_iiwa14 : public CostFunctionInterface {
public:

  CostFunction_iiwa14(const Robot* robot_ptr);

  ~CostFunction_iiwa14();

  void lq(const Robot* robot_ptr, const double t, const double* q, 
          const double* v, const double* a, const double* u, const double* f, 
          double* lq) override;

  void lv(const Robot* robot_ptr, const double t, const double* q, 
          const double* v, const double* a, const double* u, const double* f, 
          double* lv) override;

  void la(const Robot* robot_ptr, const double t, const double* q, 
          const double* v, const double* a, const double* u, const double* f, 
          double* la) override;

  void lu(const Robot* robot_ptr, const double t, const double* q, 
          const double* v, const double* a, const double* u, const double* f, 
          double* lu) override;

  void lf(const Robot* robot_ptr, const double t, const double* q, 
          const double* v, const double* a, const double* u, const double* f, 
          double* lf) override;

  void lq(const Robot* robot_ptr, const double t, const double* q, 
          const double* v, const double* a, const double* u, 
          double* lq) override;

  void lv(const Robot* robot_ptr, const double t, const double* q, 
          const double* v, const double* a, const double* u, 
          double* lv) override;

  void la(const Robot* robot_ptr, const double t, const double* q, 
          const double* v, const double* a, const double* u, 
          double* la) override;

  void lu(const Robot* robot_ptr, const double t, const double* q, 
          const double* v, const double* a, const double* u, 
          double* lu) override;

  void phiq(const Robot* robot_ptr, const double t, const double* q, 
            const double* v, double* phiq) override;

  void phiv(const Robot* robot_ptr, const double t, const double* q, 
            const double* v, double* phiv) override;

private:
  JointSpaceCost joint_space_cost_;

};

} // namespace robotcgmres

#endif // ROBOTCGMRES_COST_FUNCTION_IIWA14_HPP_