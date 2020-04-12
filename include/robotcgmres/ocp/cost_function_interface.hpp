#ifndef ROBOTCGMRES_COST_FUNCTION_INTERFACE_HPP_
#define ROBOTCGMRES_COST_FUNCTION_INTERFACE_HPP_

#include "robot/robot.hpp"


class CostFunctionInterface {
public:
  CostFunctionInterface(const Robot* robot_ptr);
  ~CostFunctionInterface();

  virtual void lq(const double t, const double* q, const double* v, 
                  const double* a, const double* u, const double* f, 
                  double* lq);

  virtual void lv(const double t, const double* q, const double* v, 
                  const double* a, const double* u, const double* f, 
                  double* lv);

  virtual void la(const double t, const double* q, const double* v, 
                  const double* a, const double* u, const double* f, 
                  double* la);

  virtual void lu(const double t, const double* q, const double* v, 
                  const double* a, const double* u, const double* f, 
                  double* lu);

  virtual void lf(const double t, const double* q, const double* v, 
                  const double* a, const double* u, const double* f, 
                  double* lf);


  virtual void phiq(const double t, const double* q, const double* v, 
                    double* phiq);

  virtual void phiv(const double t, const double* q, const double* v, 
                    double* phiv);

private:

};


#endif // ROBOTCGMRES_COST_FUNCTION_INTERFACE_HPP_