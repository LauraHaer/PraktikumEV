#include <gtest/gtest.h>

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <ctime>
#include <cstdlib>

#include "Lanczos.hh"
#include "helpfunctions.hh"
#include "standard_include.hh"
#include "tridiag_ev_solver.hh"

TEST(INVERSE_ITERATION, CalculateFromRandomFullDense) {
  std::srand(std::time(nullptr));
  int m = std::rand() % 5 + 5;
  int n = m * 2;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
  A = A + A.transpose().eval();
  Eigen::EigenSolver<Eigen::MatrixXd> es(A);
  auto eigenvalues = es.eigenvalues();
  auto eigenvectors = es.eigenvectors();

  Eigen::MatrixXd vecs(A.cols(),A.cols());
  for (int i = 0; i < A.cols(); i++) {
    vecs.col(i) = inverseIteration(A, eigenvalues(i).real());
  }
  
  EXPECT_LE((vecs.col(0) - eigenvectors.col(0).real()).norm(), 1e-8 );

}

TEST(INVERSE_ITERATION, CalculateFromDenseDiagonal) {
  std::srand(std::time(nullptr));
  int m = std::rand() % 5 + 5;
  int n = m * 2;
  Eigen::MatrixXd A = Eigen::MatrixXd::Identity(n,n);
  A.diagonal().setRandom();
  Eigen::EigenSolver<Eigen::MatrixXd> es(A);
  auto eigenvalues = es.eigenvalues();;
  auto eigenvectors = es.eigenvectors();

  Eigen::MatrixXd vecs(A.cols(),A.cols());
  for (int i = 0; i < A.cols(); i++) {
    vecs.col(i) = inverseIteration(A, eigenvalues(i).real());
  }
  EXPECT_LE((vecs.col(0) - eigenvectors.col(0).real()).norm(), 1e-8 );

}
