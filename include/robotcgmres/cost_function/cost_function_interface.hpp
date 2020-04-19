#ifndef ROBOTCGMRES_COST_FUNCTION_INTERFACE_HPP_
#define ROBOTCGMRES_COST_FUNCTION_INTERFACE_HPP_

#include "robot/robot.hpp"


namespace robotcgmres {

class CostFunctionInterface {
public:
  CostFunctionInterface() {}

  virtual ~CostFunctionInterface() {}

  virtual void lq(const Robot* robot_ptr, const double t, const double* q, 
                  const double* v, const double* a, const double* u, 
                  const double* f, double* lq) = 0;

  virtual void lv(const Robot* robot_ptr, const double t, const double* q, 
                  const double* v, const double* a, const double* u, 
                  const double* f, double* lv) = 0;

  virtual void la(const Robot* robot_ptr, const double t, const double* q, 
                  const double* v, const double* a, const double* u, 
                  const double* f, double* la) = 0;

  virtual void lu(const Robot* robot_ptr, const double t, const double* q, 
                  const double* v, const double* a, const double* u, 
                  const double* f, double* lu) = 0;

  virtual void lf(const Robot* robot_ptr, const double t, const double* q, 
                  const double* v, const double* a, const double* u, 
                  const double* f, double* lf) = 0;

  virtual void lq(const Robot* robot_ptr, const double t, const double* q, 
                  const double* v, const double* a, const double* u, 
                  double* lq) = 0;

  virtual void lv(const Robot* robot_ptr, const double t, const double* q, 
                  const double* v, const double* a, const double* u, 
                  double* lv) = 0;

  virtual void la(const Robot* robot_ptr, const double t, const double* q, 
                  const double* v, const double* a, const double* u, 
                  double* la) = 0;

  virtual void lu(const Robot* robot_ptr, const double t, const double* q, 
                  const double* v, const double* a, const double* u, 
                  double* lu) = 0;

  virtual void phiq(const Robot* robot_ptr, const double t, const double* q, 
                    const double* v, double* phiq) = 0;

  virtual void phiv(const Robot* robot_ptr, const double t, const double* q, 
                    const double* v, double* phiv) = 0;

  // Prohibits copy constructor.
  CostFunctionInterface(const CostFunctionInterface&) = delete;

  // Prohibits copy operator.
  CostFunctionInterface& operator=(const CostFunctionInterface&) = delete;

};

} // namespace robotcgmres

#endif // ROBOTCGMRES_COST_FUNCTION_INTERFACE_HPP_