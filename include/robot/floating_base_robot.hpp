#ifndef FLOATING_BASE_ROBOT_HPP_
#define FLOATING_BASE_ROBOT_HPP_

#include <string>
#include <cassert>

#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/rnea-derivatives.hpp"


namespace robotcgmres {

class FloatingBaseRobot {
public:
  FloatingBaseRobot(const std::string& urdf_file_name);

  ~FloatingBaseRobot();

  void RNEA(const double* q, const double* v, const double* a, double* tau);

  void RNEADerivativesTransDotVector(const double* q, const double* v, 
                                     const double* a, const double* f,
                                     const double* vec, 
                                     double* dRNEA_dq_dot_vec, 
                                     double* dRNEA_dv_dot_vec, 
                                     double* dRNEA_da_dot_vec);

  void CRBADotVec(const double* q, const double* vec, double* result);

  int dimq() const;

  int dimv() const;

  int dimtau() const;

  FloatingBaseRobot(const FloatingBaseRobot&) = delete;
  FloatingBaseRobot& operator=(const FloatingBaseRobot&) = delete;


private:
  pinocchio::Model model_;
  pinocchio::Data data_;
  int dimq_, dimv_, dimtau_;
};

} // namespace robotcgmres


#endif // FLOATING_BASE_ROBOT_HPP_