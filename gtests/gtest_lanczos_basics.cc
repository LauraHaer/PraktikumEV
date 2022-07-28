#include <gtest/gtest.h>
#include <stdlib.h>

#include <Eigen/Sparse>
#include <iostream>
#include <typeinfo>
#include <vector>

#include "Lanczos.hh"
#include "gtest/gtest.h"
#include "gtest_helpfunctions.hh"
#include "helpfunctions.hh"

TEST(LANCZOS_BASICS, SparseAndDenseGenerateSameValues) {
  int n = 5;
  int k = 15;
  int m = 10;
  Eigen::VectorXd v1 = Eigen::VectorXd::Random(n * n);
  Eigen::SparseMatrix<double> ASparse =
      CreateLaplaceMatrix<Eigen::SparseMatrix<double>>(n);
  Eigen::MatrixXd ADense = ASparse * Eigen::MatrixXd::Identity(n * n, n * n);
  std::vector<double> errorSparse =
      runLanczos(ASparse, v1, k, m, 1, 1e-10, false, true);
  std::vector<double> errorDense =
      runLanczos(ADense, v1, k, m, 1, 1e-10, false, true);
  for (int i = 0; i < 3; ++i) {
    EXPECT_DOUBLE_EQ(errorSparse.at(i), errorDense.at(i));
  }
}

TEST(LANCZOS_BASICS, TestSizeMismatchVector) {
  int n = 5;
  Eigen::MatrixXd A = CreateRandomDense(n);
  Eigen::VectorXd v1 = Eigen::VectorXd::Zero(n + 1);
  v1(0) = 1;
  v1 << -1.0, 0.0, 1.0, -1.0, 1.0, 100.0;
  ASSERT_DEBUG_DEATH(lanczos_ir(A, n, v1, 1),
                     "aR and A must have the same number of rows");
  ASSERT_DEBUG_DEATH(simple_lanczos(A, n, v1),
                     "aR and A must have the same number of rows");
  v1 = v1(Eigen::seq(2, 4));
  ASSERT_DEBUG_DEATH(lanczos_ir(A, n, v1, 1),
                     "aR and A must have the same number of rows");
  ASSERT_DEBUG_DEATH(simple_lanczos(A, n, v1),
                     "aR and A must have the same number of rows");
}

TEST(LANCZOS_BASICS, NonSquareMatrix) {
  int n = 5;
  Eigen::MatrixXd A = CreateRandomDense(n);
  Eigen::VectorXd v1 = Eigen::VectorXd::Zero(n);
  v1(0) = 1;
  // Lanczos_ir
  ASSERT_DEBUG_DEATH(
      lanczos_ir(A(Eigen::all, Eigen::lastN(n - 1)).eval(), n, v1, 1),
      "Matrix must be quadratic");
  ASSERT_DEBUG_DEATH(
      lanczos_ir(A(Eigen::lastN(n - 1), Eigen::all).eval(), n, v1, 1),
      "Matrix must be quadratic");
  // Simple Lanczos
  ASSERT_DEBUG_DEATH(
      simple_lanczos(A(Eigen::all, Eigen::lastN(n - 1)).eval(), n, v1, 1),
      "Matrix must be quadratic");
  ASSERT_DEBUG_DEATH(
      simple_lanczos(A(Eigen::lastN(n - 1), Eigen::all).eval(), n, v1, 1),
      "Matrix must be quadratic");
}

TEST(LANCZOS_BASICS, NonConformingParameters) {
  int n = 5;
  int m = 4;
  Eigen::MatrixXd A = CreateRandomDense(n);
  Eigen::VectorXd v1 = Eigen::VectorXd::Random(n);
  Eigen::MatrixXd V = Eigen::MatrixXd::Zero(A.rows(), m + 1);
  // Lanczos_ir
  ASSERT_DEBUG_DEATH(lanczos_ir(A, n + 1, v1, 1),
                     "Number of calc Eigenvalues must be <= dim[(A)]");
  ASSERT_DEBUG_DEATH(
      lanczos_ir(A, 2, v1, 3),
      "Kept number of eigenvalues must be <= calculated eigenvalues");
  ASSERT_DEBUG_DEATH(lanczos_ir(A, 5, v1, 3, -1.0), "aRho must be positive");
  ASSERT_DEBUG_DEATH(lanczos_ir(A, 5, v1, 3, 1.0, -2), "aEps must be positive");
  // simple_lanczos
  ASSERT_DEBUG_DEATH(simple_lanczos(A, n + 1, v1),
                     "Number of calc Eigenvalues must be <= dim[(A)]");
  ASSERT_DEBUG_DEATH(simple_lanczos(A, 5, v1, -1.0), "aRho must be positive");
  ASSERT_DEBUG_DEATH(lanczos_factorization(A, m, v1, 1.0, V),
                     "aV must have m columns");
  ASSERT_DEBUG_DEATH(
      lanczos_factorization(A, m, v1, 1.0,
                            V(Eigen::lastN(A.rows() - 1), Eigen::lastN(m))),
      "aV must have same number of rows of A");
}
