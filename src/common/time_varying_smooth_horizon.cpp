#include "common/time_varying_smooth_horizon.hpp"


namespace robotcgmres {

TimeVaryingSmoothHorizon::TimeVaryingSmoothHorizon(const double T_f, 
                                                   const double alpha, 
                                                   const double initial_time) 
  : T_f_(std::abs(T_f)),
    alpha_(std::abs(alpha)),
    initial_time_(std::abs(initial_time)) {
}

TimeVaryingSmoothHorizon::TimeVaryingSmoothHorizon(const double T_f, 
                                                   const double alpha)
  : T_f_(std::abs(T_f)),
    alpha_(std::abs(alpha)),
    initial_time_(double(0.0)) {
}

TimeVaryingSmoothHorizon::~TimeVaryingSmoothHorizon() {
}

double TimeVaryingSmoothHorizon::getLength(const double time) const {
  return T_f_ * (1.0-std::exp(-alpha_*(time-initial_time_)));
}

// Resets the parameters of the horizon.
void TimeVaryingSmoothHorizon::resetLength(const double initial_time) {
  initial_time_ = std::abs(initial_time);
}

// Resets the parameters of the horizon.
void TimeVaryingSmoothHorizon::resetLength(const double T_f, const double alpha, 
                                           const double initial_time) {
  T_f_ = std::abs(T_f);
  alpha_ = std::abs(alpha);
  initial_time_ = std::abs(initial_time);
}

} // namespace robot_cgmres