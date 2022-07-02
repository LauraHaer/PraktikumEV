#ifndef LANCZOS_SAAD_HH_
#define LANCZOS_SAAD_HH_

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <cassert>
#include <cmath>
#include <iostream>

/// Computes `m` Ritz values for a given **Hermitian** matrix.
/// @param A: Target Matrix
/// @param m: desired number of Eigenvalues
/// @param aV1: start vector av1
/// @returns m eigenvalues of aMat
template <class DeMat>
struct result_lanczos_saad {
  DeMat a;
  DeMat b;
  DeMat v;
  Eigen::VectorXcd ev;
};

template <class aMat, class aVec>
result_lanczos_saad<Eigen::Matrix<typename aVec::value_type, -1, -1>> lanczos_saad(
    const aMat& A, const int& m, const aVec& aV1) {
  // We define a Dense Matrix of Dynamik Size used in the calculations based on
  // the Template Datatype
  typedef Eigen::Matrix<typename aVec::value_type, -1, -1> tmpMat;
  assert(A.cols() == A.rows());
  assert(A.cols() == aV1.rows());
  //  static_assert(typeid(typename aMat::value_type) ==
  //                typeid(typename aVec::value_type));
  assert(m <= aV1.rows());
  assert(std::abs(sqrt(aV1.dot(aV1))) >= 1e-16);

  tmpMat TMat = tmpMat::Zero(m, m);
  tmpMat v = tmpMat::Zero(aV1.rows(), m + 2);
  v.col(1) = aV1 / sqrt(aV1.dot(aV1));
  typename aVec::value_type b_i = 0;

  // the algorithm
  for (int i = 1; i <= m; i++) {
    aVec w_i = A * v.col(i) - b_i * v.col(i - 1);  // line2
    // aVec w_i = v.col(i) - b_i * v.col(i - 1);              // line2
    typename aVec::value_type a_i = w_i.dot(v.col(i));     // line3
    w_i = w_i - a_i * v.col(i);                            // line4
    typename aVec::value_type b_ip1 = sqrt(w_i.dot(w_i));  // line5

    TMat(i - 1, i - 1) = a_i;
    if (i < m) {
      TMat(i, i - 1) = b_ip1;
      TMat(i - 1, i) = b_ip1;
    }
    if (std::abs(b_ip1 - 0) <= 1e-16) break;  // line 6
    v.col(i + 1) = w_i / b_ip1;
  }
  std::cout << "TMat = " << TMat << std::endl;
  result_lanczos_saad<tmpMat> res;
  res.v = v;
  Eigen::EigenSolver<tmpMat> es(TMat);
   std::cout << "the evs are " << es.eigenvalues() << std::endl;
  res.ev = es.eigenvalues();
  return res;
}
#endif
