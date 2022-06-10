#include <gtest/gtest.h>

#include <iostream>
#include <Eigen/Sparse>

#include "create_laplace_matrix.hh"
#include "Lanczos.hh"

TEST(LANCZOS, InputIsValid) {
  Eigen::VectorXd v1 = Eigen::VectorXd::Zero(16);
  v1(0) = 1;
  Eigen::SparseMatrix<double> LMat;
  LMat = CreateLaplaceMatrix<Eigen::SparseMatrix<double>>(4);
  lanczos(LMat, 10, v1);
}


TEST(LANCZOS, TetsStartWithEigenvector) {
  int n = 5;
  Eigen::MatrixXd A{{4, -1, 0, 0, -1},
                    {-1, 4, -1, 0, 0},
                    {0, -1, 4, -1, 0},
                    {0, 0, -1, 4, 0},
                    {-1, 0, 0, 0, 4}};
  DeVec v1 = DeVec::Zero(n);
  v1(0) = 1;
  v1 << -1, 0, 1, -1, 1;
  v1 = v1 / sqrt(v1.dot(v1));
  result res = lanczos(A, n, v1);
}
