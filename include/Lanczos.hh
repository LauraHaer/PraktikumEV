#ifndef LANCZOS_HH_
#define LANCZOS_HH_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>    // std::sort
#include <cassert>      // assert
#include <type_traits>  // std::is_same_v

#include "helpfunctions.hh"

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

template <class aMat, class aVec>
factorization_result<aVec> lanczos_factorization(
    const aMat& A, const int& m, const aVec& aR, const double& aRho,
    Eigen::Matrix<typename aVec::value_type, -1, -1> aV =
        Eigen::Matrix<typename aVec::value_type, -1, -1>(0, 0),
    const int& k = 0) {
  //
  // Type Definitions and Asserts
  typedef typename aVec::value_type aData;
  typedef Eigen::Matrix<aData, -1, -1> tmpMat;

  assert("Matrix must be quadratic" && A.cols() == A.rows());
  assert("aR and A must have the same number of rows" && A.cols() == aR.rows());
  static_assert(
      "A and aR must have the same data type" &&
      std::is_same_v<typename aMat::value_type, typename aVec::value_type>);
  assert("Number of calc Eigenvalues must be <= dim(A)" && m <= aR.rows());
  assert("aRho must be positive" && aRho >= 0);
  assert("aV must have m columns" && (aV.cols() == m || aV.cols() == 0));
  assert("aV must have same number of rows as A" &&
         (aV.rows() == A.rows() || aV.rows() == 0));

  aVec r = aR;
  aVec beta = aVec::Zero(m + 1);
  aVec alpha = aVec::Zero(m + 1);
  beta(k) = r.norm();
  tmpMat v = tmpMat::Zero(aR.rows(), m + 1);
  if (aV.cols() != 0) {
    v(Eigen::all, Eigen::lastN(m)) = aV;
  }

  // the algorithm
  for (int i = k + 1; i <= m; ++i) {
    v.col(i) = r / r.norm();                        // Step (3)
    r = A * v.col(i) - v.col(i - 1) * beta(i - 1);  // Step (4)
    alpha(i) = v.col(i).dot(r);                     // Step (5)
    r = r - v.col(i) * alpha(i);                    // Step (6)
    beta(i) = r.norm();
    if (i == A.rows()) break;  // v already spans R^(A.rows())
    if (r.norm() <
        aRho * sqrt(std::pow(alpha(i), 2) + std::pow(beta(i - 1), 2))) {
      int ii = 0;
      aVec s;
      do {
        if (ii == 5) {  // Reorthogonalization failed
          beta(i) = 0;
          do {
            createRandomVector(r);
            // Orthogonalize r using modified gram-schmidt
            for (int iii = 1; iii <= i; ++iii) {
              r = r - v.col(iii) * (v.col(iii).dot(r) / r.dot(r));
            }
          } while (r.norm() <= 1);
          break;
        }
        s = v.transpose() * r;  // Step (8)
        r = r - v * s;          // Step (9)
        alpha(i) = alpha(i) + s(i);
        beta(i) = beta(i) + s(i - 1);
        ++ii;
      } while (r.norm() <= aRho * s.norm());  // could include alpha and beta
    }
  }
  return factorization_result<aVec>{(aVec)alpha(Eigen::lastN(m)),
                                    (aVec)beta(Eigen::seqN(1, m - 1)),
                                    (tmpMat)v(Eigen::all, Eigen::lastN(m)), r};
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
  Eigen::VectorXd alpha_k;
  Eigen::VectorXd beta_k;

  while (step.beta(Eigen::seqN(0, k - 1)).maxCoeff() >= aEps) {  // Step (3)
    Eigen::MatrixXd TMat = createTMatrix(step.alpha, step.beta);

    // Reduce to solution of k largest eigenvalues
    Eigen::VectorXd eigenvalues = tridiag_ev_solver(step.alpha, step.beta);
    std::sort(eigenvalues.begin(), eigenvalues.end(), greaterEigenvalue);
    tmpMat Q = tmpMat::Identity(m, m);
    for (int j = k; j < m; ++j) {
      tmpMat tmp_input = TMat - eigenvalues[j] * tmpMat::Identity(m, m);
      tmpMat Qj = givens_q(tmp_input);
      TMat = Qj.transpose() * TMat * Qj;
      Q = Q * Qj;
    }
    step.V = step.V * Q;
    aVec rk =
        step.V.col(k) * TMat(k, k - 1) + step.r * Q(m - 1, k - 1);  // Step (10)
    step.V(Eigen::all, Eigen::lastN(m - k)).setZero();

    alpha_k = TMat.diagonal()(Eigen::seqN(0, k));
    beta_k = TMat.diagonal(1)(Eigen::seqN(0, k - 1));

    if (rk.maxCoeff() < 1e-8) break;
    step = lanczos_factorization(A, m, rk, aRho, step.V, k);  // Step (2)
    step.alpha(Eigen::seqN(0, k)) = alpha_k;
    step.beta(Eigen::seqN(0, k - 1)) = beta_k;
  }

  result_lanczos res;

  res.ev = tridiag_ev_solver(alpha_k, beta_k);
  // res.vec = eigenvectorsA(createTMatrix(alpha, beta), res.ev);
  return res;
}
#endif
