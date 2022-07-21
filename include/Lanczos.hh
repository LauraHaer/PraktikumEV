#ifndef LANCZOS_HH_
#define LANCZOS_HH_

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <type_traits>  // std::is_same_v

#include "Lanczos_saad.hh"
#include "helpfunctions.hh"
#include "standard_include.hh"
#include "tridiag_ev_solver.hh"

// Trace how often the reorthogonalization fails
#define trace_reorthogonalization

/// Computes `m` Ritz values for a given **Hermitian** matrix.
/// @param A: Target Matrix
/// @param m: desired number of Eigenvalues
/// @param aV1: start vector av1
/// @returns m eigenvalues of aMat
template <class DeMat>
struct result_lanczos {
  Eigen::VectorXd ev;
  Eigen::MatrixXd vec;
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
auto lanczos_factorization(
    const aMat& A, const int& m, const aVec& aR, const double& aRho,
    Eigen::Matrix<typename aVec::value_type, -1, -1> aV =
        Eigen::Matrix<typename aVec::value_type, -1, -1>::Zero(0, 0),
    const int& k = 0) {
  // Type Definitions and Asserts
  typedef typename aVec::value_type aData;
  typedef Eigen::Matrix<aData, -1, -1>
      tmpMat;  // Dense Matrix with dynamic size
  if (aV.rows() == 0) aV = tmpMat::Zero(aR.rows(), m + 1);
  assert("Matrix must be quadratic" && A.cols() == A.rows());
  assert("aR and A must have the same number of rows" && A.cols() == aR.rows());
  static_assert(
      "A and aR must have the same data type" &&
      std::is_same_v<typename aMat::value_type, typename aVec::value_type>);
  assert("Number of calc Eigenvalues must be <= dim(A)" && m <= aR.rows());
  assert("aRho must be positive" && aRho >= 0);

  aVec r = aR;
  aVec beta = aVec::Zero(m + 1);
  aVec alpha = aVec::Zero(m + 1);
  beta(k) = r.norm();
  tmpMat v = aV;

#ifdef trace_reorthogonalization
  int number_of_reorthorthogonalization_fails = 0;
#endif

  // the algorithm
  for (int i = k + 1; i <= m; ++i) {
    v.col(i) = r / r.norm();                        // Step (3)
    r = A * v.col(i) - v.col(i - 1) * beta(i - 1);  // Step (4)
    alpha(i) = v.col(i).dot(r);                     // Step (5)
    r = r - v.col(i) * alpha(i);                    // Step (6)
    beta(i) = r.norm();
    if (r.norm() <
        aRho * sqrt(std::pow(alpha(i), 2) + std::pow(beta(i - 1), 2))) {
      int ii = 0;
      aVec s;
      do {
        if (ii == 5) {
#ifdef trace_reorthogonalization
          ++number_of_reorthorthogonalization_fails;
#endif
          beta(i) = 0;
          r = aVec::Random(aR.rows());
          //++i;
          if (i >= m) break;
          ii = 0;
        }
        s = v.transpose() * r;  // Step (8)
        r = r - v * s;          // Step (9)
        alpha(i) = alpha(i) + s(i);
        beta(i) = beta(i) + s(i - 1);
        ++ii;
      } while (r.norm() <= aRho * s.norm());
      //               aRho * sqrt(std::pow(alpha(i), 2) + std::pow(beta(i - 1),
      //               2) +
      //                           std::pow(s.norm(), 2)));
    }
  }
#ifdef trace_reorthogonalization
  std::cout << "Reorthogonalisation failed: "
            << number_of_reorthorthogonalization_fails << " times."
            << std::endl;
#endif
  return std::make_tuple(alpha, beta, v, r);
}

template <class aMat, class aVec>
auto simple_lanczos(const aMat& A, const int& m, const aVec& aR,
                    const double& aRho = 1.0) {
  auto [alpha, beta, v, r] = lanczos_factorization(A, m, aR, aRho);  // Step (2)
  result_lanczos<Eigen::Matrix<typename aVec::value_type, -1, -1>> res;
  mos << PRINT_REFLECTION(alpha) << std::endl;
  mos << PRINT_REFLECTION(beta) << std::endl;
  res.ev = tridiag_ev_solver(alpha, beta);
  // res.vec = eigenvectorsA(createTMatrix(alpha, beta), res.ev);
  return res;
}

template <class aMat, class aVec>
// result_lanczos<Eigen::Matrix<typename aVec::value_type, -1, -1>> lanczos(
auto lanczos_ir(const aMat& A, const int& m, const aVec& aR, const int& k,
                const double& aRho = 1.0, const double& aEps = 1e-4) {
  typedef typename aVec::value_type aData;
  typedef Eigen::Matrix<aData, -1, -1> tmpMat;

  assert("Kept number of eigenvalues must be <= calculated eigenvalues" &&
         k < m);
  assert("aEps must be positive" && aEps >= 0);
  auto [alpha, beta, vm, rm] =
      lanczos_factorization(A, m, aR, aRho);  // Step (2)
  // vm = vm(Eigen::all, Eigen::lastN(m));
  auto TMat = createTMatrix(alpha, beta);

  while (TMat(Eigen::seqN(0, k), Eigen::seqN(0, k)).diagonal(1).maxCoeff() >=
         aEps) {  // Step (3)
    auto eigenvalues =
        tridiag_ev_solver(alpha, beta);  // Select last p eigenvalues
    std::sort(eigenvalues.begin(), eigenvalues.end(), greaterEigenvalue);
    tmpMat Q = tmpMat::Identity(m, m);
    for (int j = k; j < m; ++j) {
      tmpMat tmp_input = TMat - eigenvalues[j] * tmpMat::Identity(m, m);
      tmpMat Qj = givens_q(tmp_input);
      TMat = Qj.transpose() * TMat * Qj;
      Q = Q * Qj;
    }
    aVec rk = vm.col(k + 1) * TMat(k, k - 1) +
              rm * eigenvalues[k - 1] * Q(m - 1, k - 1);  // Step (10)
    tmpMat vk =
        vm(Eigen::all, Eigen::lastN(m)) * Q(Eigen::all, Eigen::seqN(0, k));
    tmpMat tk = TMat(Eigen::seqN(0, k), Eigen::seqN(0, k));
    vm = tmpMat::Zero(aR.rows(), m + 1);
    vm(Eigen::all, Eigen::seqN(1, k)) = vk;
    std::tie(alpha, beta, vm, rm) =
        lanczos_factorization(A, m, rk, aRho, vm, k);  // Step (2)
    TMat = createTMatrix(alpha, beta);
    TMat(Eigen::seqN(0, k), Eigen::seqN(0, k)) = tk;

    // std::getchar();
  }

  result_lanczos<Eigen::Matrix<typename aVec::value_type, -1, -1>> res;
  alpha(Eigen::lastN(m)) = TMat.diagonal();
  beta(Eigen::lastN(m - 1)) = TMat.diagonal(1);
  res.ev = tridiag_ev_solver(alpha, beta);
  // res.vec = eigenvectorsA(createTMatrix(alpha, beta), res.ev);
  return res;
}
#endif
