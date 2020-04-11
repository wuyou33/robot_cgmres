#ifndef ROBOTCGMRES_ZERO_HORIZON_OCP_HPP_
#define ROBOTCGMRES_ZERO_HORIZON_OCP_HPP_

#include "common/memory_manager.hpp"
#include "robot/robot.hpp"
#include "robot/cost_function_interface.hpp"
#include "robot/constraints_interface.hpp"


namespace robotcgmres {

class ZeroHorizonOCP {
public:
  ZeroHorizonOCP(Robot* robot_ptr);
  ~ZeroHorizonOCP();

  void computeOptimalityResidual(const double t, const double* q, 
                                 const double* v, const double* solution,
                                 double* optimality_residual);

  int dimq() const;

  int dimv() const;

  int dimtau() const;

  int dimtau() const;

private:
  Robot *robot_;
  CostFunctionInterface cost_function_;
  ConstraintsInterface constraints_;
  int dimq_, dimv_, dimtau_, dimf_, dimC_, dim_solution_;
  double *u_vec_, *beta_vec_, *gamma_vec_, *dRNEA_da_dot_beta_;
};

} // namespace robotcgmres


#endif // ROBOTCGMRES_ZERO_HORIZON_OCP_HPP_