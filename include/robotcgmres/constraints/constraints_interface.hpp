#ifndef ROBOTCGMRES_CONSTRAINTS_INTERFACE_HPP_
#define ROBOTCGMRES_CONSTRAINTS_INTERFACE_HPP_

#include "robot/robot.hpp"


namespace robotcgmres {

class ConstraintsInterface {
public:
  ConstraintsInterface() {}

  virtual ~ConstraintsInterface() {}

  virtual void residual(const Robot* robot_ptr, const double t, const double* q, 
                        const double* v, const double* a, const double* u, 
                        const double* f, double* residual) = 0;

  virtual void addCqDotVec(const Robot* robot_ptr, const double t, 
                           const double* q, const double* v, const double* a, 
                           const double* u, const double* f, const double* vec, 
                           double* hq) = 0;

  virtual void addCvDotVec(const Robot* robot_ptr, const double t, 
                           const double* q, const double* v, const double* a, 
                           const double* u, const double* f, const double* vec, 
                           double* hv) = 0;

  virtual void addCaDotVec(const Robot* robot_ptr, const double t, 
                           const double* q, const double* v, const double* a, 
                           const double* u, const double* f, const double* vec, 
                           double* ha) = 0;

  virtual void addCuDotVec(const Robot* robot_ptr, const double t, 
                           const double* q, const double* v, const double* a, 
                           const double* u, const double* f, const double* vec, 
                           double* hu) = 0;

  virtual void addCfDotVec(const Robot* robot_ptr, const double t, 
                           const double* q, const double* v, const double* a, 
                           const double* u, const double* f, const double* vec, 
                           double* hf) = 0;

  virtual void residual(const Robot* robot_ptr, const double t, const double* q, 
                        const double* v, const double* a, const double* u, 
                        double* residual) = 0;

  virtual void addCqDotVec(const Robot* robot_ptr, const double t, 
                           const double* q, const double* v, const double* a, 
                           const double* u, const double* vec, double* hq) = 0;

  virtual void addCvDotVec(const Robot* robot_ptr, const double t, 
                           const double* q, const double* v, const double* a, 
                           const double* u, const double* vec, double* hv) = 0;

  virtual void addCaDotVec(const Robot* robot_ptr, const double t, 
                           const double* q, const double* v, const double* a, 
                           const double* u, const double* vec, double* ha) = 0;

  virtual void addCuDotVec(const Robot* robot_ptr, const double t, 
                           const double* q, const double* v, const double* a, 
                           const double* u, const double* vec, double* hu) = 0;

  virtual void addCfDotVec(const Robot* robot_ptr, const double t, 
                           const double* q, const double* v, const double* a, 
                           const double* u, const double* vec, double* hf) = 0;

  virtual int dim_constraints() const = 0;

  // Prohibits copy constructor.
  ConstraintsInterface(const ConstraintsInterface&) = delete;

  // Prohibits copy operator.
  ConstraintsInterface& operator=(const ConstraintsInterface&) = delete;

};

} // namespace robotcgmres

#endif // ROBOTCGMRES_CONSTRAINTS_INTERFACE_HPP_