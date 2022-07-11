/*
#include <gtest/gtest.h>

#include <iostream>
#include <Eigen/Sparse>
#include <vector>

#include <stdlib.h>

#include "create_laplace_matrix.hh"
#include "Lanczos.hh"
#include "Lanczos_saad.hh"
#include "helpfunctions.hh"
#include "gtest/gtest.h"
#include <typeinfo>


TEST(LANCZOS, DISABLED_InputIsValid) {
  int n = 4;
  Eigen::VectorXd v1 = Eigen::VectorXd::Zero(16);
  v1(0) = 1;
  Eigen::SparseMatrix<double> LMat;
  result_lanczos<Eigen::MatrixXd> res;
  LMat = CreateLaplaceMatrix<Eigen::SparseMatrix<double>>(n);
  res = lanczos_ir(LMat, 10, v1, 5);
  Eigen::SparseMatrix<int> L2Mat;
//  L2Mat = CreateLaplaceMatrix<Eigen::SparseMatrix<int>>(10*n);
//  Eigen::VectorXi v11;
//  v11 = Eigen::VectorXi::Zero(L2Mat.cols());
//  v11(0) = 1;
//  res = lanczos_ir(L2Mat, 10, v11, 5);
}

//TEST(LANCZOS, IncompatibleDataType) {
//  int n = 4;
//  result_lanczos<Eigen::MatrixXd> res;
//  Eigen::VectorXd v1 = Eigen::VectorXd::Zero(n*n);
//  v1(0) = 1;
//  Eigen::SparseMatrix<float> L1Mat;
//  L1Mat = CreateLaplaceMatrix<Eigen::SparseMatrix<float>>(n);
//  L1Mat * (double)3.0;
//  res = lanczos_ir(L1Mat, 10, v1, 5);
//}


TEST(LANCZOS, DISABLED_TestStartWithEigenvector) {
  int n = 5;
  Eigen::MatrixXd A{{4.0, -1.0, 0.0, 0.0, -1.0},
                    {-1.0, 4.0, -1.0, 0.0, 0.0},
                    {0.0, -1.0, 4.0, -1.0, 0.0},
                    {0.0, 0.0, -1.0, 4.0, 0.0},
                    {-1.0, 0.0, 0.0, 0.0, 4.0}};
  Eigen::VectorXd v1 = Eigen::VectorXd::Zero(n);
  v1(0) = 1;
  v1 << -1.0, 0.0, 1.0, -1.0, 1.0;
  //result_lanczos<Eigen::MatrixXd> res = lanczos_ir(A, n, v1, 1);
  result_lanczos<Eigen::MatrixXd> res = simple_lanczos(A, n, v1,0);
  auto res2 = lanczos_saad(A, n, v1);
  EXPECT_DOUBLE_EQ(res.ev[0].real(), 5.0);

  EXPECT_EQ(res.ev.size(), n);

}


TEST(LANCZOS, DISABLED_TestStartWithoutEigenvector) {
  int n = 5;
  Eigen::MatrixXd A{{4.0, -1.0, 0.0, 0.0, -1.0},
                    {-1.0, 4.0, -1.0, 0.0, 0.0},
                    {0.0, -1.0, 4.0, -1.0, 0.0},
                    {0.0, 0.0, -1.0, 4.0, 0.0},
                    {-1.0, 0.0, 0.0, 0.0, 4.0}};
  Eigen::VectorXd v1 = Eigen::VectorXd::Zero(n);
  v1(0) = 1;
  v1 << 1.0, 1.0, 1.0, 1.0, 1.0;
  result_lanczos<Eigen::MatrixXd> res = lanczos_ir(A, n, v1, 4, 1.0, 1e-4);
  //result_lanczos<Eigen::MatrixXd> res = simple_lanczos(A, n, v1, 1.0);
  //auto res2 = lanczos_saad(A, n, v1);
  //EXPECT_DOUBLE_EQ(res.ev[0].real(), 5.73, 1e-);
  EXPECT_LE(std::abs(res.ev[0].real() -  5.73), 1e-5);

  EXPECT_EQ(res.ev.size(), n);

}



TEST(LANCZOS, TestRandomMatrixWithoutEigenvector) {
  for(int i = 0; i < 5; ++i) {
//    int n = rand() % 15 + 5;
//    int m = rand() % 16 + 10;
    int n = rand() % 5 + 5;
    int m = rand() % 6 + 10;
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(n*n, n*n);
    A = A + A.transpose().eval();
    Eigen::EigenSolver<Eigen::MatrixXd> es(A);
    //Eigen::VectorXd v1 = es.eigenvectors().col(0).real();
    Eigen::VectorXd v1 = Eigen::VectorXd::Zero(n*n);
    v1(0) = 1;
    //result_lanczos<Eigen::MatrixXd> res = lanczos_ir(A, m, v1, std::ceil(m/2), 1.0, 1e-2);
    auto eigenvalues = es.eigenvalues();
    std::sort(eigenvalues.begin(), eigenvalues.end(), greaterEigenvalue);
    result_lanczos<Eigen::MatrixXd> res = simple_lanczos(A, m, v1, 1.0);
    mos << "Eigen calculated eigenvalues" << std::endl;
    for( int i =0; i < m; ++i) std::cout << eigenvalues[i] << ", ";
    std::cout << std::endl;
    mos << "Lanczos calculated eigenvalues" << std::endl;
    for( int i =0; i < m; ++i) std::cout << res.ev[i] << ", ";
    std::cout << std::endl;
    EXPECT_LE(std::abs(res.ev[0].real() - es.eigenvalues()[0].real()), 1e-8 );
  }
}

TEST(LANCZOS, DISABLED_TestRandomMatrixWithEigenvector) {
  for(int i = 0; i < 10; ++i) {
//    int n = rand() % 15 + 5;
//    int m = rand() % 16 + 10;
    int n = rand() % 5 + 5;
    int m = rand() % 6 + 10;
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(n*n, n*n);
    A = A + A.transpose().eval();
    Eigen::EigenSolver<Eigen::MatrixXd> es(A);
    Eigen::VectorXd v1 = es.eigenvectors().col(0).real();
    result_lanczos<Eigen::MatrixXd> res = lanczos_ir(A, m, v1, std::ceil(m/2), 1.0, 1e-17);
    EXPECT_LE(std::abs(res.ev[0].real() - es.eigenvalues()[0].real()), 1e-8 );
  }
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
  ASSERT_DEBUG_DEATH(lanczos_ir(A, n, v1, 1), "[*]*" );
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
  ASSERT_DEBUG_DEATH(lanczos_ir(A, n, v1, 1), "[*]*" );
  // TODO: put better regex :3
}

// NOTTODO
// TODO dense matrix?, non-hermitian matrix?
//
*/
