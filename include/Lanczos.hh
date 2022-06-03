#ifndef LANCZOS_HH_
#define LANCZOS_HH_

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <cassert>
#include <cmath>
#include <iostream>

typedef Eigen::SparseMatrix<double>
    SpMat;  // declares a column-major sparse matrix type of double
typedef Eigen::SparseVector<double>
    SpVec;  // declares a column-major sparse matrix type of double

typedef Eigen::MatrixXd
    DeMat;  // declares a column-major sparse matrix type of double
typedef Eigen::VectorXd
    DeVec;  // declares a column-major sparse matrix type of double

struct result {
  DeMat a;
  DeMat b;
  DeMat v;
  DeVec ev;
};

//  give matrix A and number of eigenvalues you want m
// template <typename Data, template <typename, int, int> class Matrix>
// void lanczos(Matrix<Data, int, int> A, int m) {
//  std::cout << A << std::endl;
result lanczos(DeMat A, int m, DeVec aV1) {
  int n = A.cols();  // dimension of A
  assert(A.cols() == A.rows());
  assert(m <= n);
  assert(sqrt(aV1.dot(aV1)));

  DeMat TMat = DeMat::Zero(m, m);
  DeMat v = DeMat::Zero(n, m + 2);
  v.col(1) = aV1;
  double b_i = 0;

  // the algorithm
  for (int i = 1; i <= m; i++) {
    DeVec w_i = A * v.col(i) - b_i * v.col(i - 1);  // line2
    double a_i = w_i.dot(v.col(i));                 // line3
    w_i = w_i - a_i * v.col(i);                     // line4
    double b_ip1 = sqrt(w_i.dot(w_i));              // line5

    TMat(i - 1, i - 1) = a_i;
    if (i < m) {
      TMat(i, i - 1) = b_ip1;
      TMat(i - 1, i) = b_ip1;
    }
    if (b_ip1 == 0) break;  // line 6
    v.col(i + 1) = w_i / b_ip1;
  }
  std::cout << "TMat = " << TMat << std::endl;
  result res;
  res.v = v;
  Eigen::EigenSolver<DeMat> es(TMat);
  std::cout << "the evs are " << es.eigenvalues() << std::endl;
  // res.ev = es.eigenvalues();
  // res.ev = TMat.eigenvalues();
  return res;
}
#endif
