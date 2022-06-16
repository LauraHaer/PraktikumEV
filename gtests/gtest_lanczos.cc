#include <gtest/gtest.h>

#include <iostream>
#include <Eigen/Sparse>
#include <vector>

#include "create_laplace_matrix.hh"
#include "Lanczos.hh"
#include "gtest/gtest.h"
#include <typeinfo>

TEST(LANCZOS, InputIsValid) {
  int n = 4;
  Eigen::VectorXd v1 = Eigen::VectorXd::Zero(16);
  v1(0) = 1;
  Eigen::SparseMatrix<double> LMat;
  LMat = CreateLaplaceMatrix<Eigen::SparseMatrix<double>>(n);
  result_lanczos<Eigen::MatrixXd> res = lanczos(LMat, 10, v1);
}


TEST(LANCZOS, TetsStartWithEigenvector) {
  int n = 5;
  Eigen::MatrixXd A{{4.0, -1.0, 0.0, 0.0, -1.0},
                    {-1.0, 4.0, -1.0, 0.0, 0.0},
                    {0.0, -1.0, 4.0, -1.0, 0.0},
                    {0.0, 0.0, -1.0, 4.0, 0.0},
                    {-1.0, 0.0, 0.0, 0.0, 4.0}};
  Eigen::VectorXd v1 = Eigen::VectorXd::Zero(n);
  v1(0) = 1;
  v1 << -1.0, 0.0, 1.0, -1.0, 1.0;
  v1 = v1 / sqrt(v1.dot(v1));
  result_lanczos<Eigen::MatrixXd> res = lanczos(A, n, v1);
  // std::cout << typeid(Eigen::VectorXcd).name() << std::endl;
  // std::cout << typeid(res.ev[0]).name() << std::endl;
  std::vector<double> ev_values = {5.0,0.0,0.0,0.0,0.0};
  // EXPECT_DOUBLE_EQ(res.ev[0], val2);
  for (int i = 0; i<ev_values.size(); i++){
    EXPECT_DOUBLE_EQ(res.ev[i].real(), ev_values[i]);
  }

  EXPECT_EQ(res.ev.size(), v1.size());

}

TEST(LANCZOS, TestSizeMismatchVector){
  int n = 5;
  Eigen::MatrixXd A{{4.0, -1.0, 0.0, 0.0, -1.0},
                    {-1.0, 4.0, -1.0, 0.0, 0.0},
                    {0.0, -1.0, 4.0, -1.0, 0.0},
                    {0.0, 0.0, -1.0, 4.0, 0.0},
                    {-1.0, 0.0, 0.0, 0.0, 4.0}};
  Eigen::VectorXd v1 = Eigen::VectorXd::Zero(n+1);
  v1(0) = 1;
  v1 << -1.0, 0.0, 1.0, -1.0, 1.0, 100.0;
  v1 = v1 / sqrt(v1.dot(v1));
  // [a-zA-Z0-9//<>)()/:, \][=;._-]*.
  ASSERT_DEBUG_DEATH(lanczos(A, n, v1), "[*]*" );
  // TODO: put better regex :3
}

TEST(LANCZOS, TestSizeMismatchMatrix){
  int n = 5;
  Eigen::MatrixXd A{{4.0, -1.0, 0.0, 0.0, -1.0},
                    {-1.0, 4.0, -1.0, 0.0, 0.0},
                    {0.0, -1.0, 4.0, -1.0, 0.0},
                    {-1.0, 0.0, 0.0, 0.0, 4.0}};
  Eigen::VectorXd v1 = Eigen::VectorXd::Zero(n);
  v1(0) = 1;
  v1 << -1.0, 0.0, 1.0, -1.0, 1.0;
  v1 = v1 / sqrt(v1.dot(v1));
  // [a-zA-Z0-9//<>)()/:, \][=;._-]*.
  ASSERT_DEBUG_DEATH(lanczos(A, n, v1), "[*]*" );
  // TODO: put better regex :3
}

// TODO dense matrix?, non-hermitian matrix?
//
