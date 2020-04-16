#ifndef ROBOTCGMRES_ZERO_HORIZON_NEWTON_HPP_
#define ROBOTCGMRES_ZERO_HORIZON_NEWTON_HPP_

#include "robot/robot.hpp"
#include "ocp/zero_horizon_ocp.hpp"
#include "ocp/newton_system_for_ocp.hpp"
#include "ocp/fdgmres.hpp"
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
                    const double finite_difference_increment, const int kmax);

  // Frees arrays used in Newton iteration.
  ~ZeroHorizonNewton();

  // Sets arrays used in Newton iteration.
  // Argments:
  //   robot_ptr: Pointer of the robot model.
  //   finite_difference_increment: The increment of the finite difference 
  void setCriteriaForTermination(const int max_newton_iteration, 
                                 const double criteria_newton_termination); 

  void setInitialGuessSolution(const double* initial_guess_solution); 


  // Solves the zero-horizon OCP and the result is stored in solution.
  void solveZeroHorizonOCP(const double t, const double* q, const double* v, 
                           const double* solution);

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