#include <gtest/gtest.h>

#include <Eigen/Sparse>
#include <iostream>

#include "create_laplace_matrix.hh"

TEST(EIGEN, CreateLaplaceMatrix) {
  Eigen::SparseMatrix<double> LMat;
  LMat = CreateLaplaceMatrix<Eigen::SparseMatrix<double>>(4);
}
