#ifndef HELPFUNCTIONS_HH
#define HELPFUNCTIONS_HH

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>  // sqrt, pow etc.
#include <complex>
#include <ctime>  // For initialization for random

bool greaterEigenvalue(std::complex<double> a, std::complex<double> b);

template <class aVec>
void createRandomVector(aVec& aV, int aSeed = std::time(nullptr)) {
  std::srand(aSeed);
  for (auto iter = aV.begin(); iter != aV.end(); iter++) {
    *iter = std::rand() % 20;
  }
  aV / aV.norm();
  return;
}

template <class aVec>
Eigen::MatrixXd createTMatrix(aVec diag, aVec sdiag) {
  int n = diag.rows();

  // build T with given diagonal and subdiagonal
  Eigen::MatrixXd TMat = Eigen::MatrixXd::Zero(n, n);
  TMat.diagonal(1) = sdiag;
  TMat.diagonal(-1) = sdiag;
  TMat.diagonal() = diag;

  return TMat;
}

template <class aMat>
Eigen::MatrixXd givens_q(aMat A) {
  int n = A.cols();
  Eigen::MatrixXd q = Eigen::MatrixXd::Identity(n, n);
  for (int i = 0; i < n - 1; ++i) {
    // get givens matrix entries c and s
    double c = 1;
    double s = 0;
    double x = A(i, i);
    double y = A(i + 1, i);
    double r = std::sqrt((std::pow(x, 2) + std::pow(y, 2)));
    if (r > 1e-8) {  // 1e-8 is chosen arbitarily
      c = x / r;
      s = y / r;
    }

    // q_i contains givens matrix inside identity matrix
    aMat q_i = aMat::Identity(n, n);
    q_i(i, i) = c;
    q_i(i + 1, i + 1) = c;
    q_i(i, i + 1) = s;
    q_i(i + 1, i) = -s;

    A = q_i * A * q_i.transpose();
    q = q_i * q;
  }
  return q.transpose();
}

template <class aVec>
aVec tridiag_ev_solver(aVec diag, aVec sdiag) {
  Eigen::MatrixXd TMat = createTMatrix(diag, sdiag);
  Eigen::EigenSolver<Eigen::MatrixXd> es(TMat);
  aVec res = es.eigenvalues().real();
  return res;
}

#endif
