#ifndef LANCZOS_HH_
#define LANCZOS_HH_

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <cassert>
#include <cmath>
#include <iostream>
#include <algorithm>

#include "Lanczos_saad.hh"
#include "tridiag_ev_solver.hh"
#include "standard_include.hh"
#include "helpfunctions.hh"

/// Computes `m` Ritz values for a given **Hermitian** matrix.
/// @param A: Target Matrix
/// @param m: desired number of Eigenvalues
/// @param aV1: start vector av1
/// @returns m eigenvalues of aMat
template <class DeMat>
struct result_lanczos {
  DeMat v;
  Eigen::VectorXcd ev;
};

/// Computes the Lanczos factorisation of A, starting at the k+1 step
/// @param A: Target Matrix
/// @param m: desired number of Eigenvalues
/// @param aR: start vector aR
/// @param aRho: factor for checking if reorthogonalization is necessary
/// @param aV: start Matrix V
/// @param k: number of already computed steps
/// @returns alpha, beta, V, r
template <class aMat, class aVec>
auto lanczos_factorization(const aMat& A, const int& m, const aVec& aR,
                           const double& aRho, Eigen::Matrix<typename aVec::value_type, -1, -1> aV = Eigen::Matrix<typename aVec::value_type, -1, -1>::Zero(0,0), const int& k=0) {
  //Type Definitions and Asserts
  typedef typename aVec::value_type aData;
  typedef Eigen::Matrix<aData, -1, -1> tmpMat; // Dense Matrix with dynamic size
  if (aV.rows() == 0) aV = tmpMat::Zero(aR.rows(), m + 1);
  assert("Matrix must be quadratic" &&  A.cols() == A.rows());
  assert("aR and A must have the same number of rows" &&  A.cols() == aR.rows());
  //static_assert(std::is_same(typeid(typename aMat::value_type),
  //               typeid(typename aVec::value_type)));
  assert("Number of calc Eigenvalues must be <= dim(A)" && m <= aR.rows());

  aVec r = aR;
  Eigen::Matrix<aData, -1, 1> beta = Eigen::Matrix<aData, -1, 1>::Zero(m + 1);
  Eigen::Matrix<aData, -1, 1> alpha = Eigen::Matrix<aData, -1, 1>::Zero(m + 1);
  beta(k) = r.norm();
  tmpMat v = aV;

  // the algorithm
  for (int i = k+1; i <= m; ++i) {
//    mos << "Start of loop" << std::endl;
//    mos << PRINT_REFLECTION(r) << std::endl;
    v.col(i) = r / r.norm();                        // Step (3)
    r = A * v.col(i) - v.col(i - 1) * beta(i - 1);  // Step (4)
//    mos << PRINT_REFLECTION(r) << std::endl;
    alpha(i) = v.col(i).dot(r);                     // Step (5)
    r = r - v.col(i) * alpha(i);                    // Step (6)
//    mos << PRINT_REFLECTION(r) << std::endl;
    beta(i) = r.norm();
    for(int ii = 0; ii < 6; ++ii) {
      if (ii == 5) std::cout << "Reorthogonalisation failed" << std::endl;
      if (r.norm() <
          aRho * sqrt(std::pow(alpha(i), 2) + std::pow(beta(i - 1), 2))) {
        //mos << "Reorthogonalisation required:" << PRINT_REFLECTION(ii) << std::endl;
        Eigen::Matrix<aData, -1, 1> s = v.transpose() * r;  // Step (8)
        r = r - v * s;                                      // Step (9)
        alpha(i) = alpha(i) + s(i);
        beta(i) = beta(i) + s(i - 1);
      } else {
        break;
      }
    }
  }
  return std::make_tuple(alpha, beta, v, r);
}

template <class aMat, class aVec>
auto simple_lanczos(const aMat& A, const int& m, const aVec& aR,
                    const double& aRho = 1.0) {
  auto [alpha, beta, v, r] = lanczos_factorization(A, m, aR, aRho);  // Step (2)
  result_lanczos<Eigen::Matrix<typename aVec::value_type, -1, -1>> res;
  res.v = v;
  res.ev = tridiag_ev_solver(alpha, beta);
  return res;
}


template <class aMat, class aVec>
// result_lanczos<Eigen::Matrix<typename aVec::value_type, -1, -1>> lanczos(
auto lanczos_ir(const aMat& A, const int& m, const aVec& aR, const int& k,
             const double& aRho = 1.0, const double& aEps = 1e-4) {
  typedef typename aVec::value_type aData;
  typedef Eigen::Matrix<aData, -1, -1> tmpMat;

  assert("Kept number of eigenvalues must be <= calculated eigenvalues" && k < m);
  auto [alpha, beta, vm, rm] = lanczos_factorization(A, m, aR, aRho);  // Step (2)
  //vm = vm(Eigen::all, Eigen::lastN(m));
  auto TMat = createTMatrix(alpha, beta);
  //mos << PRINT_REFLECTION(TMat) << std::endl;
  while (TMat(Eigen::seqN(0,k), Eigen::seqN(0,k)).diagonal(1).maxCoeff() >= aEps) {  // Step (3)
    auto eigenvalues = tridiag_ev_solver(alpha, beta); // Select last p eigenvalues
    std::sort(eigenvalues.begin(), eigenvalues.end(), greaterEigenvalue);
    Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(m,m);
    for (int j=k; j < m; ++j) {
      Eigen::MatrixXd tmp_input = TMat - eigenvalues[j] * Eigen::MatrixXd::Identity(m,m);
      //auto Qj = qr.householderQ();
      Eigen::MatrixXd Qj = givens_q(tmp_input);
      TMat = Qj.transpose() * TMat * Qj;
      Q = Q * Qj;
    }
    aVec rk = vm.col(k+1) * TMat(k, k-1) + rm*eigenvalues[k-1] * Q(m-1,k-1); // Step (10)
    //mos << PRINT_REFLECTION(rk) << std::endl;
    tmpMat vk = vm(Eigen::all, Eigen::lastN(m)) * Q(Eigen::all, Eigen::seqN(0, k));
    tmpMat tk = TMat(Eigen::seqN(0, k), Eigen::seqN(0, k));
    vm = tmpMat::Zero(aR.rows(), m + 1);
    vm(Eigen::all, Eigen::seqN(1, k)) = vk;
    std::tie(alpha, beta, vm, rm) = lanczos_factorization(A, m, rk, aRho, vm, k);  // Step (2)
    TMat = createTMatrix(alpha, beta);
    TMat(Eigen::seqN(0, k), Eigen::seqN(0, k)) = tk;
    //mos << PRINT_REFLECTION(TMat) << std::endl;
    //std::getchar();
  }

  result_lanczos<Eigen::Matrix<typename aVec::value_type, -1, -1>> res;
  res.v = vm;
  alpha(Eigen::lastN(m)) = TMat.diagonal();
  beta(Eigen::lastN(m-1)) = TMat.diagonal(1);
  res.ev = tridiag_ev_solver(alpha, beta);
  return res;
  }
#endif
