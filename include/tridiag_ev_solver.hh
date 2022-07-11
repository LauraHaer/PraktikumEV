#ifndef TRIDIAG_EV_SOLVER_HH
#define TRIDIAG_EV_SOLVER_HH

//#include "Lanczos.hh"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>

#include "standard_include.hh"

// typedef Eigen::VectorXd aVec;  // change into template
// typedef Eigen::MatrixXd aMat;

typedef Eigen::SparseMatrix<double>
    SpMat;  // declares a column-major sparse matrix type of double

template <class aVec>
Eigen::MatrixXd createTMatrix(aVec diag, aVec sdiag) {
  int n = diag.rows() - 1;

  // build T with given diagonal and subdiagonal
  Eigen::MatrixXd TMat = Eigen::MatrixXd::Zero(n, n);
  TMat.diagonal(1) = sdiag(Eigen::seqN(1, n - 1));
  TMat.diagonal(-1) = sdiag(Eigen::seqN(1, n - 1));
  TMat.diagonal() = diag(Eigen::lastN(n));

  return TMat;
}

template<class aMat>
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
    if (r > 1e-8) {
      c = x / r;
      s = y / r;
    }

    // q_i has givens matrix inside identity matrix
    SpMat q_i(n, n);
    for(int ii=0; ii < n; ++ii) q_i.coeffRef(ii,ii) = 1;
    q_i.coeffRef(i, i) = c;
    q_i.coeffRef(i + 1, i + 1) = c;
    q_i.coeffRef(i, i + 1) = s;
    q_i.coeffRef(i + 1, i) = -s;

    A = q_i * A * q_i.transpose();

    q = q_i * q;
  }
  return q.transpose();
}

template <class aVec>
aVec tridiag_ev_solver(aVec diag, aVec sdiag) {
  // build T with given diagonal and subdiagonal
  Eigen::MatrixXd TMat = createTMatrix(diag, sdiag);
  // document error as maximum of subdiagonal,
  double e = 1;

  // qr iteration until T converged to diagonal matrix
  while (e > 1e-8) {
    Eigen::MatrixXd q = givens_q(TMat);
    TMat = q.transpose() * TMat * q;
    // update error
    e = TMat.diagonal<-1>().maxCoeff();
  }
  return TMat.diagonal();
}


template <class aMat>
Eigen::VectorXd inverseIteration(const aMat T, const double mu, const double eps = 1e-8) {
  Eigen::VectorXd x_prev = Eigen::VectorXd::Zero(T.cols());
  double e;
  Eigen::VectorXd x;
  do {
    x = (T - mu * Eigen::MatrixXd::Identity(T.cols(), T.cols())).inverse() * x;
    e = (x-x_prev).norm();
    x_prev = x;
  } while( e > eps );
  return x;
}

//template <class aMat, class aVec>
//Eigen::MatrixXd eigenvectorsA(aMat T, aVec ev, aMat v) {
//  v = v(Eigen::all, Eigen::lastN(ev.cols()));
//  Eigen::MatrixXd vec(ev.cols(), ev.cols());
//
//  for (int i = 0; i < res.ev.cols(); i++) {
////    res.vec.col(i) = inverseIteration(TMat, res.ev(i).real());
//    //compute eigenvectors of T
//    vec.col(i) = inverseIteration(t, ev(i));
//    //compute eigenvectors of A from eigenvectors of T
//    vec.col(i) = v*vec.col(i);
//  }
//
//  for
//}

#endif
