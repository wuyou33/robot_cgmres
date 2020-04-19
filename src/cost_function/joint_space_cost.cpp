#include "cost_function/joint_space_cost.hpp"


namespace robotcgmres {

JointSpaceCost::JointSpaceCost(const Robot* robot_ptr)
  : dimq_(robot_ptr->dimq()),
    dimv_(robot_ptr->dimv()),
    q_ref_(memorymanager::NewVector(dimq_)),
    v_ref_(memorymanager::NewVector(dimv_)),
    a_ref_(memorymanager::NewVector(dimv_)),
    u_ref_(memorymanager::NewVector(dimv_)),
    wq_(memorymanager::NewVector(dimq_)),
    wv_(memorymanager::NewVector(dimv_)),
    wa_(memorymanager::NewVector(dimv_)),
    wu_(memorymanager::NewVector(dimv_)) {
}

JointSpaceCost::~JointSpaceCost() {
  memorymanager::DeleteVector(q_ref_);
  memorymanager::DeleteVector(v_ref_);
  memorymanager::DeleteVector(a_ref_);
  memorymanager::DeleteVector(u_ref_);
  memorymanager::DeleteVector(wq_);
  memorymanager::DeleteVector(wv_);
  memorymanager::DeleteVector(wa_);
  memorymanager::DeleteVector(wu_);
}

void JointSpaceCost::lq(const double* q, double* lq) const {
  for (int i=0; i<dimq_; ++i) {
    lq[i] = wq_[i] * (q[i] - q_ref_[i]);
  }
}

void JointSpaceCost::lv(const double* v, double* lv) const {
  for (int i=0; i<dimv_; ++i) {
    lv[i] = wv_[i] * (v[i] - v_ref_[i]);
  }
}

void JointSpaceCost::la(const double* a, double* la) const {
  for (int i=0; i<dimv_; ++i) {
    la[i] = wa_[i] * (a[i] - a_ref_[i]);
  }
}

void JointSpaceCost::lu(const double* u, double* lu) const {
  for (int i=0; i<dimv_; ++i) {
    lu[i] = wu_[i] * (u[i] - u_ref_[i]);
  }
}

void JointSpaceCost::set_q_ref(const std::vector<int>& joint_indices, 
                               const double q_ref) {
  for (int i : joint_indices) {
    if (i >= 0 && i < dimq_) {
      q_ref_[i] = q_ref;
    }
  }
}

void JointSpaceCost::set_v_ref(const std::vector<int>& joint_indices, 
                               const double v_ref) {
  for (int i : joint_indices) {
    if (i >= 0 && i < dimq_) {
      v_ref_[i] = v_ref;
    }
  }
}

void JointSpaceCost::set_a_ref(const std::vector<int>& joint_indices, 
                               const double a_ref) {
  for (int i : joint_indices) {
    if (i >= 0 && i < dimq_) {
      a_ref_[i] = a_ref;
    }
  }
}

void JointSpaceCost::set_u_ref(const std::vector<int>& joint_indices, 
                               const double u_ref) {
  for (int i : joint_indices) {
    if (i >= 0 && i < dimq_) {
      u_ref_[i] = u_ref;
    }
  }
}

void JointSpaceCost::set_q_weight(const std::vector<int>& joint_indices, 
                                  const double q_weight) {
  for (int i : joint_indices) {
    if (i >= 0 && i < dimq_) {
      wq_[i] = q_weight;
    }
  }
}

void JointSpaceCost::set_v_weight(const std::vector<int>& joint_indices, 
                                  const double v_weight) {
  for (int i : joint_indices) {
    if (i >= 0 && i < dimv_) {
      wv_[i] = v_weight;
    }
  }
}

void JointSpaceCost::set_a_weight(const std::vector<int>& joint_indices, 
                                  const double a_weight) {
  for (int i : joint_indices) {
    if (i >= 0 && i < dimv_) {
      wa_[i] = a_weight;
    }
  }
}

void JointSpaceCost::set_u_weight(const std::vector<int>& joint_indices, 
                                  const double u_weight) {
  for (int i : joint_indices) {
    if (i >= 0 && i < dimv_) {
      wu_[i] = u_weight;
    }
  }
}

} // namespace robotcgmres