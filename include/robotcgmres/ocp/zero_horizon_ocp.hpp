#ifndef ROBOTCGMRES_ZERO_HORIZON_OCP_HPP_
#define ROBOTCGMRES_ZERO_HORIZON_OCP_HPP_

#include "common/memory_manager.hpp"
#include "robot/robot.hpp"
#include "ocp/passive_joint_constraint.hpp"
#include "ocp/cost_function_interface.hpp"
#include "ocp/constraints_interface.hpp"


namespace robotcgmres {

class ZeroHorizonOCP {
public:
  ZeroHorizonOCP(const Robot* robot_ptr);
  ~ZeroHorizonOCP();

  void computeOptimalityResidual(const double t, const double* q, 
                                 const double* v, const double* solution,
                                 double* optimality_residual);

  int dimq() const;

  int dimv() const;

  int dimtau() const;

  int dim_constraints() const;

  int dim_solution() const;

  ZeroHorizonOCP(const ZeroHorizonOCP&) = delete;

  ZeroHorizonOCP& operator=(const ZeroHorizonOCP&) = delete;

private:
  Robot *robot_;
  PassiveJointConstraint passive_joint_constraint_;
  CostFunctionInterface cost_function_;
  ConstraintsInterface constraints_;
  int dimq_, dimv_, dimf_, dim_passive_, dim_constraints_, dim_solution_;
  double *u_, *beta_, *gamma_, *dRNEA_da_dot_beta_, *dRNEA_dfext_dot_beta_,
          *dBaum_dq_dot_alpha_, *dBaum_dv_dot_alpha_, *dBaum_da_dot_alpha_;
};

} // namespace robotcgmres


#endif // ROBOTCGMRES_ZERO_HORIZON_OCP_HPP_