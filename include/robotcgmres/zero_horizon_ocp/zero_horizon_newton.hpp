#ifndef ROBOTCGMRES_ZERO_HORIZON_NEWTON_HPP_
#define ROBOTCGMRES_ZERO_HORIZON_NEWTON_HPP_

#include "robot/robot.hpp"
#include "zero_horizon_ocp/zero_horizon_ocp.hpp"
#include "solver/newton_system_for_ocp.hpp"
#include "solver/fdgmres.hpp"
#include "cost_function/cost_function_interface.hpp"
#include "common/memory_manager.hpp"


namespace robotcgmres {

class ZeroHorizonNewton {
public:
  // Allocates arrays used in Newton iteration.
  // Argments:
  //   robot_ptr: Pointer of the robot model.
  //   finite_difference_increment: The increment of the finite difference 
  //      approximation.
  //   kmax: The maximum number of GMRES iteration.
  ZeroHorizonNewton(const Robot* robot_ptr, 
                    const CostFunctionInterface* cost_function,
                    const double finite_difference_increment, 
                    const unsigned int kmax);

  // Frees arrays used in Newton iteration.
  ~ZeroHorizonNewton();

  // Sets parameters for Newton iteration.
  // Argments:
  //   robot_ptr: Pointer of the robot model.
  //   finite_difference_increment: The increment of the finite difference.
  void setCriteriaForTermination(const unsigned int max_newton_iteration, 
                                 const double criteria_newton_termination); 

  // Sets initial guess of the solution of the Newton iteration.
  // Argments:
  //   initial_guess_solution: The initial guess of the solutino. The size must 
  //      be dim_solution().
  void setInitialGuessSolution(const double* initial_guess_solution); 

  // Solves the zero-horizon OCP and the result is stored in solution.
  // Argments: 
  //   t: Time.
  //   q: Generalized configuration. The size must be dimq().
  //   v: Generalized velocity. The size must be dimv().
  //   solution: The solution of the zero horizon OCP. The size must be 
  //      dim_solution().
  void solveZeroHorizonOCP(const double t, const double* q, const double* v, 
                           double* solution);

  // Returns the dimension of the generalized configuration.
  int dimq() const;

  // Returns the dimension of the generalized velocity.
  int dimv() const;

  // Returns the dimenstion of the solution.
  int dim_solution() const;

  // Prohibits copy constructor.
  ZeroHorizonNewton(const ZeroHorizonNewton&) = delete;

  // Prohibits copy operator.
  ZeroHorizonNewton& operator=(const ZeroHorizonNewton&) = delete;

private:
  NewtonSystemForOCP<ZeroHorizonOCP> newton_;
  FDGMRES<NewtonSystemForOCP<ZeroHorizonOCP>> fdgmres_;
  double *solution_update_, *initial_guess_solution_;
  double criteria_newton_termination_;
  int max_newton_iteration_;
};

} // namespace robotcgmres


#endif // ROBOTCGMRES_ZERO_HORIZON_NEWTON_HPP_