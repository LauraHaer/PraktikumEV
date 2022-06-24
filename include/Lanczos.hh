#ifndef LANCZOS_HH_
#define LANCZOS_HH_

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <cassert>
#include <cmath>
#include <iostream>

#include "Lanczos_saad.hh"

/// Computes `m` Ritz values for a given **Hermitian** matrix.
/// @param A: Target Matrix
/// @param m: desired number of Eigenvalues
/// @param aV1: start vector av1
/// @returns m eigenvalues of aMat
template <class DeMat>
struct result_lanczos {
  DeMat a;
  DeMat b;
  DeMat v;
  Eigen::VectorXcd ev;
};

template <class aMat, class aVec>
result_lanczos<Eigen::Matrix<typename aVec::value_type, -1, -1>> lanczos_factorization(
    const aMat& A, const int& m, const aVec& aR, const double& aRho) {
  // We define a Dense Matrix of Dynamik Size used in the calculations based on
  // the Template Datatype
  typedef typename aVec::value_type aData;
  typedef Eigen::Matrix<aData, -1, -1> tmpMat;
//  assert("Matrix must be quadratic", A.cols() == A.rows());
//  assert("aR and A must have the same number of rows", A.cols() == aR.rows());
  //  static_assert(typeid(typename aMat::value_type) ==
  //                typeid(typename aVec::value_type));
//  assert("Number of calculated Eigenvalues must be <= dim(A)", m <= aR.rows());

  aVec r = aR;
  //Eigen::Matrix<typename aVec::value_type, -1, 1> beta =
  Eigen::Matrix<aData, -1, 1> beta =
      Eigen::Matrix<aData, -1, 1>::Zero(m + 1);
  Eigen::Matrix<aData, -1, 1> alpha =
      Eigen::Matrix<aData, -1, 1>::Zero(m + 1);
  beta(0) = r.norm();
//  assert("startvector aR should have norm >> 0", beta(0) >= 1e-16);
  tmpMat v = tmpMat::Zero(aR.rows(), m + 1);

  // the algorithm
  for (int i = 1; i <= m; ++i) {
    v.col(i) = r / r.norm();                        // Step (3)
    r = A * v.col(i) - v.col(i - 1) * beta(i - 1);  // Step (4)
    alpha(i) = v.col(i).dot(r);                     // Step (5)
    r = r - v.col(i) * alpha(i);                    // Step (6)
    if (r.norm() <
        aRho * sqrt(std::pow(alpha(i), 2) + std::pow(beta(i - 1), 2))) {
      Eigen::Matrix<aData, -1, 1> s =
          v.transpose() * r;  // Step (8)
      r = r - v * s;          // Step (9)
      alpha(i) = alpha(i) + s(i);
      beta(i) = beta(i) + s(i - 1);
    }
  }
  // std::cout << "TMat = " << TMat << std::endl;
  result_lanczos<tmpMat> res;
  //res.v = v;
  //Eigen::EigenSolver<tmpMat> es(TMat);
  //// std::cout << "the evs are " << es.eigenvalues() << std::endl;
  //res.ev = es.eigenvalues();
  return res;
}

template <class aMat, class aVec>
result_lanczos<Eigen::Matrix<typename aVec::value_type, -1, -1>> lanczos(
    const aMat& A, const int& m, const aVec& aR, const double& aRho=1.0) {
    result_lanczos<Eigen::Matrix<typename aVec::value_type, -1, -1>> res;
    auto res1 = lanczos_saad(A, m , aR);
    res.v = res1.v;
    res.ev = res1.ev;
  return res;
}
#endif
