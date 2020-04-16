#ifndef ROBOTCGMRES_NEWTON_SYSTEM_FOR_OCP_HPP_
#define ROBOTCGMRES_NEWTON_SYSTEM_FOR_OCP_HPP_

#include <cmath>

#include "robot/robot.hpp"
#include "common/memory_manager.hpp"
#include "common/linear_algebra.hpp"


namespace robotcgmres {

template <class OCPType>
class NewtonSystemForOCP {
public:

  // Allocates arrays used in Newton iteration.
  // Argments:
  //   robot_ptr: Pointer of the robot model.
  //   finite_difference_increment: The increment of the finite difference 
  //      approximation.
  NewtonSystemForOCP(const Robot* robot_ptr, 
                     const double finite_difference_increment)
    : ocp_(robot_ptr),
      finite_difference_increment_(finite_difference_increment),
      incremented_solution_(memorymanager::NewVector(ocp_.dim_solution())), 
      optimality_residual_(memorymanager::NewVector(ocp_.dim_solution())), 
      incremented_optimality_residual_(memorymanager::NewVector(ocp_.dim_solution())) {
  }

  // Frees arrays used in Newton iteration.
  NewtonSystemForOCP::~NewtonSystemForOCP() {
    memorymanager::DeleteVector(incremented_solution_);
    memorymanager::DeleteVector(optimality_residual_);
    memorymanager::DeleteVector(incremented_optimality_residual_);
  }

  // Compute b - Ax used in the FDGMRES.
  // Argments: 
  //   t: Time.
  //   q: Generalized configuration. The size must be dimq.
  //   v: Generalized velocity. The size must be dimv.
  //   solution: Current iterate solution. The size must be dim_solution.
  //   initial_direction: Initial value of the update of the iterate solution. 
  //      The size must be dim_solution.
  //   initial_residual: The residual corresponding to b-Ax of the linear problem. 
  //      The size must be dim_solution.
  void computeInitialResidual(const double t, const double* q, const double* v, 
                              const double* solution, 
                              const double* initial_direction, 
                              double* initial_residual) {
    for (int i=0; i<ocp_.dim_solution(); ++i) {
      incremented_solution_[i] = solution[i] 
          + finite_difference_increment_ * initial_direction[i];
    }
    // OCP with contacts
    if (ocp_.dimf() > 0) {
      ocp_.computeOptimalityResidual(t, q, v, solution, 
                                    &(solution[ocp_.dimv()]), 
                                    &(solution[ocp_.dimv()+ocp_.dimf()]),
                                    optimality_residual_);
      ocp_.computeOptimalityResidual(t, q, v, incremented_solution_, 
                                    &(incremented_solution_[ocp_.dimv()]), 
                                    &(incremented_solution_[ocp_.dimv()+ocp_.dimf()]),
                                    incremented_optimality_residual_);
    }
    // OCP without contacts
    else {
      ocp_.computeOptimalityResidual(t, q, v, solution, 
                                    &(solution[ocp_.dimv()]), 
                                    optimality_residual_);
      ocp_.computeOptimalityResidual(t, q, v, 
                                    incremented_solution_, 
                                    &(incremented_solution_[ocp_.dimv()]), 
                                    incremented_optimality_residual_);
    }
    for (int i=0; i<ocp_.dim_solution(); ++i) {
      initial_residual[i] = - optimality_residual_[i] 
          - (incremented_optimality_residual_[i]-optimality_residual_[i]) 
              / finite_difference_increment_;
    }
  }

  // Compute Ax used in the FDGMRES.
  // Argments: 
  //   t: Time.
  //   q: Generalized configuration. The size must be dimq.
  //   v: Generalized velocity. The size must be dimv.
  //   solution: Current iterate solution. The size must be dim_solution.
  //   direction: The update of the iterate solution. The size must be 
  //      dim_solution.
  //   Ax: Ax of the linear problem. The size must be dim_solution.
  void computeAx(const double t, const double* q, const double* v, 
                  const double* solution, const double* direction, double* Ax) {
    for (int i=0; i<ocp_.dim_solution(); ++i) { 
      incremented_solution_[i] = solution[i] 
                                  + finite_difference_increment_ * direction[i];
    }
    // OCP with contacts
    if (ocp_.dimf() > 0) {
      ocp_.computeOptimalityResidual(
          t, q, v, incremented_solution_, &(incremented_solution_[ocp_.dimv()]), 
          &(incremented_solution_[ocp_.dimv()+ocp_.dimf()]),
          incremented_optimality_residual_);
    }
    // OCP without contacts
    else {
      ocp_.computeOptimalityResidual(
          t, q, v, incremented_solution_, &(incremented_solution_[ocp_.dimv()]), 
          incremented_optimality_residual_);
    }
    for (int i=0; i<ocp_.dim_solution(); ++i) {
      Ax[i] = (incremented_optimality_residual_[i]-optimality_residual_[i]) 
                  / finite_difference_increment_;
    }
  }

  // Computes and returns the residual of the OCP.
  // Argments: 
  //   t: Time.
  //   q: Generalized configuration. The size must be dimq.
  //   v: Generalized velocity. The size must be dimv.
  //   solution: The solution of the Newton system. The size must be 
  //      dim_solution.
  double residualNorm(const double t, const double* q, const double* v, 
                      const double* solution) {
    // OCP with contacts
    if (ocp_.dimf() > 0) {
      ocp_.computeOptimalityResidual(t, q, v, solution, 
                                     &(solution[ocp_.dimv()]), 
                                     &(solution[ocp_.dimv()+ocp_.dimf()]), 
                                     optimality_residual_);
    }
    // OCP without contacts
    else {
      ocp_.computeOptimalityResidual(t, q, v, solution, 
                                     &(solution[ocp_.dimv()]), 
                                     optimality_residual_);
    }
    return std::sqrt(linearalgebra::SquaredNorm(ocp_.dim_solution(), 
                                                optimality_residual_));
  }

  // Integrates the time variation of the solution of the OCP.
  // Argments: 
  //   solution_update: The time variation of the solution.  
  //   integration_length: The length of the integration. 
  //   solution: The solution of the OCP.
  void integrateSolution(const double* solution_update, 
                         const double integration_length, double* solution) {
    ocp_.integrateSolution(solution_update, integration_length, solution);
  }

  // Returns the dimenstion of the solution.
  int dim_solution() const {
    return ocp_.dim_solution();
  }

  // Prohibits copy constructor.
  NewtonSystemForOCP(const NewtonSystemForOCP&) = delete;

  // Prohibits copy operator.
  NewtonSystemForOCP& operator=(const NewtonSystemForOCP&) = delete;

private:
  OCPType ocp_;
  int dim_solution_;
  double finite_difference_increment_;
  double *incremented_solution_, *optimality_residual_, 
         *incremented_optimality_residual_;

};
  
} // namespace robotcgmres


#endif // ROBOTCGMRES_NEWTON_SYSTEM_FOR_OCP_HPP_