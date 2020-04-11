#ifndef ROBOTCGMRES_CONSTRAINTS_INTERFACE_HPP_
#define ROBOTCGMRES_CONSTRAINTS_INTERFACE_HPP_

#include "robot/robot.hpp"


class ConstraintsInterface {
public:
  ConstraintsInterface();
  ~ConstraintsInterface();

  virtual void C(const double t, const double* q, const double* v, 
                 const double* a, const double* u, double* C);

  virtual void addCq(const double t, const double* q, const double* v, 
                     const double* a, const double* u, double* lq);

  virtual void addCv(const double t, const double* q, const double* v, 
                     const double* a, const double* u, double* lv);

  virtual void addCa(const double t, const double* q, const double* v, 
                     const double* a, const double* u, double* la);

  virtual void addCu(const double t, const double* q, const double* v, 
                     const double* a, const double* u, double* lu);


  virtual void C(const double t, const double* q, const double* v, 
                 const double* a, const double* u, const double* f, double* C);

  virtual void addCq(const double t, const double* q, const double* v, 
                     const double* a, const double* u, const double* f, 
                     double* lq);

  virtual void addCv(const double t, const double* q, const double* v, 
                     const double* a, const double* u, const double* f, 
                     double* lv);

  virtual void addCa(const double t, const double* q, const double* v, 
                     const double* a, const double* u, const double* f, 
                     double* la);

  virtual void addCu(const double t, const double* q, const double* v, 
                     const double* a, const double* u, const double* f, 
                     double* lu);

  virtual void addCf(const double t, const double* q, const double* v, 
                     const double* a, const double* u, const double* f, 
                     double* lu);


  virtual int dimC() const;

private:

};


#endif // ROBOTCGMRES_CONSTRAINTS_INTERFACE_HPP_