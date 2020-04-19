#ifndef ROBOTCGMRES_CONSTRAINTS_IIWA14_HPP_
#define ROBOTCGMRES_CONSTRAINTS_IIWA14_HPP_

#include "constraints/constraints_interface.hpp"
#include "robot/robot.hpp"


namespace robotcgmres {

class Constraints_iiwa14 : public ConstraintsInterface {
public:

  Constraints_iiwa14(const Robot* robot_ptr);

  ~Constraints_iiwa14();

  void residual(const Robot* robot_ptr, const double t, const double* q, 
                const double* v, const double* a, const double* u, 
                const double* f, double* residual) override;

  void addCqDotVec(const Robot* robot_ptr, const double t, 
                   const double* q, const double* v, const double* a, 
                   const double* u, const double* f, const double* vec, 
                   double* hq) override;

  void addCvDotVec(const Robot* robot_ptr, const double t, 
                   const double* q, const double* v, const double* a, 
                   const double* u, const double* f, const double* vec, 
                   double* hv) override;

  void addCaDotVec(const Robot* robot_ptr, const double t, 
                   const double* q, const double* v, const double* a, 
                   const double* u, const double* f, const double* vec, 
                   double* ha) override;

  void addCuDotVec(const Robot* robot_ptr, const double t, 
                   const double* q, const double* v, const double* a, 
                   const double* u, const double* f, const double* vec, 
                   double* hu) override;

  void addCfDotVec(const Robot* robot_ptr, const double t, 
                   const double* q, const double* v, const double* a, 
                   const double* u, const double* f, const double* vec, 
                   double* hf) override;

  void residual(const Robot* robot_ptr, const double t, const double* q, 
                const double* v, const double* a, const double* u, 
                double* residual) override;

  void addCqDotVec(const Robot* robot_ptr, const double t, 
                   const double* q, const double* v, const double* a, 
                   const double* u, const double* vec, double* hq) override;

  void addCvDotVec(const Robot* robot_ptr, const double t, 
                   const double* q, const double* v, const double* a, 
                   const double* u, const double* vec, double* hv) override;

  void addCaDotVec(const Robot* robot_ptr, const double t, 
                   const double* q, const double* v, const double* a, 
                   const double* u, const double* vec, double* ha) override;

  void addCuDotVec(const Robot* robot_ptr, const double t, 
                   const double* q, const double* v, const double* a, 
                   const double* u, const double* vec, double* hu) override;

  void addCfDotVec(const Robot* robot_ptr, const double t, 
                   const double* q, const double* v, const double* a, 
                   const double* u, const double* vec, double* hf) override;

  int dim_constraints() const override;

private:

};

} // namespace robotcgmres

#endif // ROBOTCGMRES_CONSTRAINTS_IIWA14_HPP_