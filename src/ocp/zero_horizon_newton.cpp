#include "ocp/zero_horizon_newton.hpp"


namespace robotcgmres {

ZeroHorizonNewton::ZeroHorizonNewton(const Robot* robot_ptr, 
                                     const double finite_difference_increment, 
                                     const int kmax)
  : newton_(robot_ptr, finite_difference_increment),
    fdgmres_(newton_.dim_solution(), kmax),
    solution_update_(memorymanager::NewVector(newton_.dim_solution())),
    initial_guess_solution_(memorymanager::NewVector(newton_.dim_solution())),
    criteria_newton_termination_(1e-08),
    max_newton_iteration_(50) {
}

ZeroHorizonNewton::~ZeroHorizonNewton() {
  memorymanager::DeleteVector(solution_update_);
  memorymanager::DeleteVector(initial_guess_solution_);
}

void ZeroHorizonNewton::setCriteriaForTermination(
    const int max_newton_iteration, const double criteria_newton_termination) {
  max_newton_iteration_ = max_newton_iteration; 
  criteria_newton_termination_ = criteria_newton_termination; 
}

void ZeroHorizonNewton::setInitialGuessSolution(
    const double* initial_guess_solution) {
  for (int i=0; i<newton_.dim_solution(); ++i) {
    solution_[i] = initial_guess_solution[i];
  }
}

void ZeroHorizonNewton::solveZeroHorizonOCP(const double t, const double* q, 
                                            const double* v, double* solution) {
  for (int i=0; i<newton_.dim_solution(); ++i) {
    solution[i] = initial_guess_solution_[i];
  }
  int num_newton = 0;
  double error = newton_.residualNorm(t, q, v, solution);
  while (num_newton < max_newton_iteration_   
             && error > criteria_newton_termination_) {
    fdgmres_.computeSolutionUpdate(newton_, t, q, v, solution, solution_update_);
    newton_.integrateSolution(solution_update_, 1.0, solution);
  }
}

int ZeroHorizonNewton::dim_solution() const {
  return newton_.dim_solution();
}

} // namespace robotcgmres