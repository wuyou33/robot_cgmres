#ifndef ROBOTCGMRES_ZERO_HORIZON_OCP_HPP_
#define ROBOTCGMRES_ZERO_HORIZON_OCP_HPP_

#include "robot/robot.hpp"
#include "robot/passive_joints.hpp"
#include "cost_function/cost_function_interface.hpp"
#include "constraints/constraints_interface.hpp"
#include "common/memory_manager.hpp"


namespace robotcgmres {

class ZeroHorizonOCP {
public:

  // Allocates arrays used in the OCP.
  // Argments:
  //   robot_ptr: Pointer of the robot model.
  ZeroHorizonOCP(const Robot* robot_ptr);

  // Frees arrays used in the OCP.
  ~ZeroHorizonOCP();

  // Computes optimality residual. If there are contacts, the solution of 
  // the OCP is given by [a, multipliers], whose dimensions are
  // dimv and dim_passive+dim_constraints, respectively.
  // The optimality_residual must have the size of 
  // dimv+dim_passive+dim_constraints.
  // Argments: 
  //   t: Time.
  //   q: Generalized configuration. The size must be dimq.
  //   v: Generalized velocity. The size must be dimv.
  //   a: Generalized acceleration. The size must be dimv.
  //   multipliers: The Lagrange multipliers for equality constraints. The size
  //      must be dim_passive+dim_constraints.
  //   optimality_residual: The residual of the otimality condition. The size
  //      must be dimv+dim_passive+dim_constraints.
  void computeOptimalityResidual(const double t, const double* q, 
                                 const double* v, const double* a, 
                                 const double* mutipliers, 
                                 double* optimality_residual);

  // Computes optimality residual. If there are contacts, the solution of 
  // the OCP is given by [a, fext, multipliers], whose dimensions are
  // dimv_ , dimf_, dimf_+dim_passive_+dim_constraints_, respectively.
  // The optimality_residual must have the size of 
  // dimv+2*dimf_+dim_passive+dim_constraints.
  // Argments: 
  //   t: Time.
  //   q: Generalized configuration. The size must be dimq.
  //   v: Generalized velocity. The size must be dimv.
  //   a: Generalized acceleration. The size must be dimv.
  //   fext: Generalized acceleration. The size must be dimf.
  //   multipliers: The Lagrange multipliers for equality constraints. The size
  //      must be dimf+dim_passive+dim_constraints.
  //   optimality_residual: The residual of the otimality condition. The size
  //      must be dimv+2*dimf+dim_passive+dim_constraints.
  void computeOptimalityResidual(const double t, const double* q, 
                                 const double* v, const double* a, 
                                 const double* fext, const double* mutipliers, 
                                 double* optimality_residual);

  // Integrates the time variation of the solution of the OCP.
  // Argments: 
  //   solution_update: The time variation of the solution.  
  //   integration_length: The length of the integration. 
  //   solution: The solution of the OCP.
  void integrateSolution(const double* solution_update, 
                         const double integration_length, double* solution);

  // Returns the dimension of the generalized configuration.
  int dimq() const;

  // Returns the dimension of the generalized velocity.
  int dimv() const;

  // Returns the total dimension of the contacts.
  int dimf() const;

  // Returns the dimension of the passive torques.
  int dim_passive() const;

  // Returns the dimension of the constrainpassiveConstraintsResidualts.
  int dim_constraints() const;

  // Returns the dimension of the solution. If there are contacts, returns 
  // dimv_+2*dimf_+dim_passive_+dim_constraints_. If there are no  contacts, 
  // returns dimv_+dim_passive_+dim_constraints_. 
  int dim_solution() const;

  // Prohibits copy constructor.
  ZeroHorizonOCP(const ZeroHorizonOCP&) = delete;

  // Prohibits copy operator.
  ZeroHorizonOCP& operator=(const ZeroHorizonOCP&) = delete;

private:
  Robot *robot_;
  PassiveJoints passive_joints_;
  CostFunctionInterface cost_function_;
  ConstraintsInterface constraints_;
  int dimq_, dimv_, dimf_, dim_passive_, dim_constraints_, dim_solution_;
  double *u_, *beta_, *gamma_, *dRNEA_da_dot_beta_, *dRNEA_dfext_dot_beta_,
         *dBaum_dq_dot_alpha_, *dBaum_dv_dot_alpha_, *dBaum_da_dot_alpha_;
};

} // namespace robotcgmres


#endif // ROBOTCGMRES_ZERO_HORIZON_OCP_HPP_