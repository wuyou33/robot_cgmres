#ifndef ROBOTCGMRES_ZERO_HORIZON_OCP_HPP_
#define ROBOTCGMRES_ZERO_HORIZON_OCP_HPP_

#include "robot/robot.hpp"
#include "cost_function/cost_function_interface.hpp"
#include "constraints/constraints_interface.hpp"
#include "common/memory_manager.hpp"


namespace robotcgmres {

class ZeroHorizonOCP {
public:

  // Allocates arrays used in the OCP.
  // Argments:
  //   robot_ptr: Pointer of the robot model.
  ZeroHorizonOCP(const Robot* robot_ptr,  
                 const CostFunctionInterface* cost_function,
                 const ConstraintsInterface* constraints);

  // Frees arrays used in the OCP.
  ~ZeroHorizonOCP();

  // Computes optimality residual. If there are contacts, the solution of 
  // the OCP is given by [a, multipliers], whose dimensions are
  // dimv and dim_passive+dim_constraints, respectively.
  // The solutin is restricted to the 
  // The optimality_residual must have the size of 
  // dimv+dim_passive+dim_constraints.
  // Argments: 
  //   t: Time.
  //   q: Generalized configuration. The size must be dimq.
  //   v: Generalized velocity. The size must be dimv.
  //   solution: composed of the acceleration, the contact forces, and 
  //      the Lagrange multipliers for equality constraints. The size
  //      must be dimv+2*dimf+dim_passive+dim_constraints.
  //   optimality_residual: The residual of the otimality condition. The size
  //      must be dimv+2*dimf+dim_passive+dim_constraints.
  void computeOptimalityResidual(const double t, const double* q, 
                                 const double* v,  const double* solution, 
                                 double* optimality_residual);

  // Integrates the time variation of the solution of the OCP.
  // Argments: 
  //   solution_update: The time variation of the solution.  
  //   integration_length: The length of the integration. 
  //   solution: The solution of the OCP.
  void integrateSolution(const double* solution_update, 
                         const double integration_length, double* solution);

  // Extracts the control input torques of the system from the solution of the 
  // OCP.
  // Argments: 
  //   solution: The given solution. The size must be dim_solution.
  //   control_input_torques: The memory where the resultant control input is 
  //      stored.
  void getControlInputFromSolution(const double* q, const double* v, 
                                   const double* solution, 
                                   double* control_input_torques);

  // Returns the dimension of the generalized configuration.
  int dimq() const;

  // Returns the dimension of the generalized velocity.
  int dimv() const;

  // Returns the total dimension of the contacts.
  int dimf() const;

  // Returns the dimension of the constraints on the passive jonits.
  int dim_passive() const;

  // Returns the dimension of the equality constraints.
  int dim_constraints() const;

  // Returns the dimension of the solution.  
  // dimv_+2*dimf_+dim_passive_+dim_constraints_. 
  int dim_solution() const;

  // Prohibits copy constructor.
  ZeroHorizonOCP(const ZeroHorizonOCP&) = delete;

  // Prohibits copy operator.
  ZeroHorizonOCP& operator=(const ZeroHorizonOCP&) = delete;

private:
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
  //   nu: The Lagrange multipliers with respect to the passive joint 
  //      constraints. The size must be dim_passive.
  //   mu: The Lagrange multipliers with respect to the user-defined 
  //      equality constraints. The size must be dim_constraints.
  //   optimality_residual: The residual of the otimality condition. The size
  //      must be dimv+dim_passive+dim_constraints.
  inline void computeOptimalityResidualWithoutContacts(
      const double t, const double* q, const double* v, const double* a, 
      const double* nu, const double* mu, double* optimality_residual);

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
  //   alpha: The Lagrange multipliers with respect to the contact constraints. 
  //      The size must be dim_contacts.
  //   nu: The Lagrange multipliers with respect to the passive joint 
  //      constraints. The size must be dim_passive.
  //   mu: The Lagrange multipliers with respect to the user-defined 
  //      equality constraints. The size must be dim_constraints.
  //   optimality_residual: The residual of the otimality condition. The size
  //      must be dimv+2*dimf+dim_passive+dim_constraints.
  inline void computeOptimalityResidualWithContacts(
      const double t, const double* q, const double* v, const double* a, 
      const double* fext, const double* alpha, const double* nu, 
      const double* mu, double* optimality_residual);

  Robot *robot_;
  CostFunctionInterface *cost_function_;
  ConstraintsInterface *constraints_;
  int dimq_, dimv_, dimf_, dim_passive_, dim_constraints_, dim_solution_;
  double *u_, *beta_, *gamma_, *dRNEA_da_dot_beta_, *dRNEA_dfext_dot_beta_,
         *dBaum_dq_dot_alpha_, *dBaum_dv_dot_alpha_, *dBaum_da_dot_alpha_;
};

} // namespace robotcgmres


#endif // ROBOTCGMRES_ZERO_HORIZON_OCP_HPP_