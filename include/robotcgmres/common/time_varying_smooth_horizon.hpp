#ifndef ROBOTCGMRES_TIME_VARYING_SMOOTH_HORIZON_HPP_
#define ROBOTCGMRES_TIME_VARYING_SMOOTH_HORIZON_HPP_

#include <cmath>


namespace robotcgmres {

// Provides time varying smooth horizon. The length of the horizon is given by
// T_f_ * (1.0-std::exp(-alpha_*(time-initial_time_))).
class TimeVaryingSmoothHorizon {
public:

  // Sets the parameters of the horizon.
  // The length of the horizon is set by 
  // T_f_ * (1.0-std::exp(-alpha_*(time-initial_time_))).
  // Argments:
  //    T_f: Maxmum length of the horizon. The horizon length converges to this 
  //        value.
  //    alpha: Determines the speed of the convergence of the length of the 
  //        horizon.
  //    initial_time: Criteria of the time-varying parameter of the horizon. 
  TimeVaryingSmoothHorizon(const double T_f, const double alpha, 
                           const double initial_time);

  // Sets the parameters of the horizon. The initial_time is set by zeto.
  // The length of the horizon is set by T_f_ * (1.0-std::exp(-alpha_*time)).
  // Argments:
  //    T_f: Maxmum length of the horizon. The horizon length converges to this 
  //        value.
  //    alpha: Determines the speed of the convergence of the length of the 
  //        horizon.
  TimeVaryingSmoothHorizon(const double T_f, const double alpha);

  ~TimeVaryingSmoothHorizon();

  // Returns the length of the horizon given by 
  // T_f_ * (1.0-std::exp(-alpha_*(time-initial_time_))).
  // Argments:
  //    time: The time parameter of the horizon.
  double getLength(const double time) const;

  // Resets the parameters of the horizon.
  // Argments:
  //    initial_time: Criteria of the time-varying parameter of the horizon. 
  void resetLength(const double initial_time);

  // Resets the parameters of the horizon.
  // Argments:
  //    T_f: Maxmum length of the horizon. The horizon length converges to this 
  //        value.
  //    alpha: Determines the speed of the convergence of the length of the 
  //        horizon.
  //    initial_time: Criteria of the time-varying parameter of the horizon. 
  void resetLength(const double T_f, const double alpha, 
                   const double initial_time);

  // Prohibits copy.
  TimeVaryingSmoothHorizon(const TimeVaryingSmoothHorizon&) = delete;

  // Prohibits copy operator.
  TimeVaryingSmoothHorizon& operator=(const TimeVaryingSmoothHorizon&) = delete;

private:
  double T_f_, alpha_, initial_time_;
};

} // namespace robot_cgmres


#endif // ROBOTCGMRES_TIME_VARYING_SMOOTH_HORIZON_HPP_