#ifndef ROBOTCGMRES_CONSTRAINTS_INTERFACE_HPP_
#define ROBOTCGMRES_CONSTRAINTS_INTERFACE_HPP_

#include "robot/robot.hpp"


class ConstraintsInterface {
public:
  ConstraintsInterface(const Robot* robot_ptr);
  virtual ~ConstraintsInterface();


  virtual void C(const double t, const double* q, const double* v, 
                 const double* a, const double* u, double* C);

  virtual void addCqDotVec(const double t, const double* q, const double* v, 
                           const double* a, const double* u,  
                           const double* vec, double* hq);

  virtual void addCvDotVec(const double t, const double* q, const double* v, 
                           const double* a, const double* u,  
                           const double* vec, double* hv);

  virtual void addCaDotVec(const double t, const double* q, const double* v, 
                           const double* a, const double* u,  
                           const double* vec, double* ha);

  virtual void addCuDotVec(const double t, const double* q, const double* v, 
                           const double* a, const double* u,  
                           const double* vec, double* hu);


  virtual void C(const double t, const double* q, const double* v, 
                 const double* a, const double* u, const double* f, double* C);

  virtual void addCqDotVec(const double t, const double* q, const double* v, 
                           const double* a, const double* u, const double* f, 
                           const double* vec, double* hq);

  virtual void addCvDotVec(const double t, const double* q, const double* v, 
                           const double* a, const double* u, const double* f, 
                           const double* vec, double* hv);

  virtual void addCaDotVec(const double t, const double* q, const double* v, 
                           const double* a, const double* u, const double* f, 
                           const double* vec, double* ha);

  virtual void addCuDotVec(const double t, const double* q, const double* v, 
                           const double* a, const double* u, const double* f, 
                           const double* vec, double* hu);

  virtual void addCfDotVec(const double t, const double* q, const double* v, 
                           const double* a, const double* u, const double* f, 
                           const double* vec, double* hf);


  virtual int dim_constraints() const;

private:

};


#endif // ROBOTCGMRES_CONSTRAINTS_INTERFACE_HPP_