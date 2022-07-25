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
//#define trace_reorthogonalization

struct result_lanczos {
  Eigen::VectorXd ev;
  Eigen::MatrixXd vec;
};

template <class aVec>
struct factorization_result {
  aVec alpha;
  aVec beta;
  Eigen::Matrix<typename aVec::value_type, -1, -1> V;
  aVec r;
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
factorization_result<aVec> lanczos_factorization(
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
    if (i == A.rows()) break;  // TODO Check if this is correkt
    if (r.norm() <
        aRho * sqrt(std::pow(alpha(i), 2) + std::pow(beta(i - 1), 2))) {
      int ii = 0;
      aVec s;
      do {
        if (ii == 5) {  // Reorthogonalization failed
#ifdef trace_reorthogonalization
          ++number_of_reorthorthogonalization_fails;
#endif
          beta(i) = 0;
          do {
            createRandomVector(r);
            aVec old_r = r;
            for (int iii = 1; iii <= i; ++iii) {
              r = r - v.col(iii) * (v.col(iii).dot(old_r) / v.col(iii).norm());
            }
          } while (r.norm() <= 1);
          break;
        }
        s = v.transpose() * r;  // Step (8)
        r = r - v * s;          // Step (9)
        alpha(i) = alpha(i) + s(i);
        beta(i) = beta(i) + s(i - 1);
        ++ii;
      } while (r.norm() <= aRho * s.norm());  // TODO Check alternatives
      //} while (r.norm() <
      //         aRho * sqrt(std::pow(alpha(i), 2) + std::pow(beta(i - 1), 2)
      //         +
      //                     std::pow(s.norm(), 2)));
    }
  }
#ifdef trace_reorthogonalization
  std::cout << "Reorthogonalisation failed: "
            << number_of_reorthorthogonalization_fails << " times."
            << std::endl;
#endif
  // mos << PRINT_REFLECTION(v) << std::endl;
  //
  factorization_result<aVec> res = {(aVec)alpha(Eigen::lastN(m)),
                                    (aVec)beta(Eigen::seqN(1, m - 1)), v, r};
  return res;
  //  return std::make_tuple((aVec)alpha(Eigen::lastN(m)),
  //                         (aVec)beta(Eigen::seqN(1, m - 1)), v, r);
}

template <class aMat, class aVec>
auto simple_lanczos(const aMat& A, const int& m, const aVec& aR,
                    const double& aRho = 1.0) {
  factorization_result<aVec> step =
      lanczos_factorization(A, m, aR, aRho);  // Step (2)
  result_lanczos res;
  res.ev = tridiag_ev_solver(step.alpha, step.beta);
  // res.vec = eigenvectorsA(createTMatrix(alpha, beta), res.ev);
  return res;
}

template <class aMat, class aVec>
result_lanczos lanczos_ir(const aMat& A, const int& m, const aVec& aR,
                          const int& k, const double& aRho = 1.0,
                          const double& aEps = 1e-4) {
  typedef typename aVec::value_type aData;
  typedef Eigen::Matrix<aData, -1, -1> tmpMat;

  assert("Kept number of eigenvalues must be <= calculated eigenvalues" &&
         k < m);
  assert("aEps must be positive" && aEps >= 0);

  factorization_result<aVec> step =
      lanczos_factorization(A, m, aR, aRho);  // Step (2)
  // vm = vm(Eigen::all, Eigen::lastN(m));
  Eigen::VectorXd diag_k;
  Eigen::VectorXd sdiag_k;

  while (step.beta.maxCoeff() >= aEps) {  // Step (3)
    Eigen::MatrixXd TMat = new_createTMatrix(step.alpha, step.beta);
    Eigen::EigenSolver<Eigen::MatrixXd> es(TMat);

    // Select last p eigenvalues
    Eigen::VectorXd eigenvalues = tridiag_ev_solver(step.alpha, step.beta);
    std::sort(eigenvalues.begin(), eigenvalues.end(), greaterEigenvalue);
    // mos << PRINT_REFLECTION(eigenvalues) << std::endl;

    tmpMat Q = tmpMat::Identity(m, m);
    for (int j = k; j < m; ++j) {
      tmpMat tmp_input = TMat - eigenvalues[j] * tmpMat::Identity(m, m);
      tmpMat Qj = givens_q(tmp_input);
      TMat = Qj.transpose() * TMat * Qj;
      Q = Q * Qj;
    }

    tmpMat V = step.V(Eigen::all, Eigen::lastN(m)) * Q;
    step.V(Eigen::all, Eigen::lastN(m)) = V;
    aVec rk = step.V.col(k + 1) * TMat(k, k - 1) +
              step.r * Q(m - 1, k - 1);  // Step (10)
    step.V(Eigen::all, Eigen::lastN(m-k)) = tmpMat::Zero(step.V.rows(), m-k);

    diag_k = TMat.diagonal()(Eigen::seqN(0, k));
    sdiag_k = TMat.diagonal(1)(Eigen::seqN(0, k - 1));

    if (rk.maxCoeff() < 1e-8) break;
    step = lanczos_factorization(A, m, rk, aRho, step.V, k);  // Step (2)
    step.alpha(Eigen::seqN(0, k)) = diag_k;
    step.beta(Eigen::seqN(0, k - 1)) = sdiag_k;
  }

  result_lanczos res;

  res.ev = tridiag_ev_solver(diag_k, sdiag_k);
  //  res.vec = eigenvectorsA(createTMatrix(alpha, beta), res.ev);
  return res;
}
#endif
