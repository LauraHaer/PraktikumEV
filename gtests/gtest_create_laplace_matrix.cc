#include <gtest/gtest.h>

#include <Eigen/Sparse>
#include <iostream>

#include "gtest_helpfunctions.hh"

TEST(LaplaceMatrix, DISABLED_HasCorrectSize) {
  Eigen::SparseMatrix<double> LMat;
  LMat = CreateLaplaceMatrix<Eigen::SparseMatrix<double>>(4);
  EXPECT_EQ(LMat.rows(), 16);
  EXPECT_EQ(LMat.cols(), 16);
  LMat = CreateLaplaceMatrix<Eigen::SparseMatrix<double>>(6);
  EXPECT_EQ(LMat.rows(), 36);
  EXPECT_EQ(LMat.cols(), 36);
}

TEST(LaplaceMatrix, DISABLED_CheckZeroDiagonals) {
  int n = 6;
  Eigen::SparseMatrix<double> LMat;
  LMat = CreateLaplaceMatrix<Eigen::SparseMatrix<double>>(n);
  for (int k = 0; k < LMat.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(LMat, k); it; ++it) {
      EXPECT_TRUE(std::abs(it.row() - it.col()) <= 1 ||
                  std::abs(it.row() - it.col()) == n);
    }
  }
}

TEST(LaplaceMatrix, DISABLED_CorrectNumberNonZeros) {
  int n = 12;
  Eigen::SparseMatrix<double> LMat;
  LMat = CreateLaplaceMatrix<Eigen::SparseMatrix<double>>(n);
  int nonzeros = 5 * n * n - 4 * n;
  EXPECT_EQ(LMat.nonZeros(), nonzeros);
}
