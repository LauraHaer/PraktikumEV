#include <gtest/gtest.h>

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <ctime>
#include <cstdlib>

#include "Lanczos.hh"
#include "helpfunctions.hh"
#include "standard_include.hh"


TEST(LANCZOS_IR, CalculateFromFixedDenseDiagonal) {
  int n = 5;
  int m = 4;
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n,n);
  Eigen::VectorXd diag{{1,2,3,4,5}};
  A.diagonal() = diag;
  Eigen::VectorXd v1 = Eigen::VectorXd::Ones(n);
  Eigen::EigenSolver<Eigen::MatrixXd> es(A);
  auto eigenvalues = es.eigenvalues();
  std::sort(eigenvalues.begin(), eigenvalues.end(), greaterEigenvalue);
  result_lanczos<Eigen::MatrixXd> res = lanczos_ir(A, m, v1, std::ceil(m/2), 1.0, 1e-7);
  mos << "Eigen calculated eigenvalues" << std::endl;
  for( int i =0; i < m; ++i) std::cout << eigenvalues[i] << ", ";
  std::cout << std::endl;
  mos << "Lanczos calculated eigenvalues" << std::endl;
  for( int i =0; i < m; ++i) std::cout << res.ev[i] << ", ";
  std::cout << std::endl;
  EXPECT_LE(std::abs(res.ev[0] - eigenvalues[0].real()), 1e-8 );
}


TEST(LANCZOS_IR, CalculateFromDenseDiagonal) {
  std::srand(std::time(nullptr));
  int m = std::rand() % 5 + 5;
  int n = m * 2;
  Eigen::MatrixXd A = Eigen::MatrixXd::Identity(n,n);
  A.diagonal().setRandom();
  Eigen::EigenSolver<Eigen::MatrixXd> es(A);
  auto eigenvalues = es.eigenvalues();
  std::sort(eigenvalues.begin(), eigenvalues.end(), greaterEigenvalue);
  Eigen::VectorXd v1 = Eigen::VectorXd::Random(n);
  result_lanczos<Eigen::MatrixXd> res = lanczos_ir(A, m, v1, std::ceil(m/2), 1.0, 1e-7);
  mos << "Eigen calculated eigenvalues" << std::endl;
  for( int i =0; i < m; ++i) std::cout << eigenvalues[i] << ", ";
  std::cout << std::endl;
  mos << "Lanczos calculated eigenvalues" << std::endl;
  for( int i =0; i < m; ++i) std::cout << res.ev[i] << ", ";
  std::cout << std::endl;
  EXPECT_LE(std::abs(res.ev[0] - eigenvalues[0].real()), 1e-8 );

}

TEST(LANCZOS_IR, CalculateFromRandomFullDense) {
  std::srand(std::time(nullptr));
  int m = std::rand() % 5 + 5;
  int n = m * 2;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
  A = A + A.transpose().eval();
  Eigen::EigenSolver<Eigen::MatrixXd> es(A);
  auto eigenvalues = es.eigenvalues();
  std::sort(eigenvalues.begin(), eigenvalues.end(), greaterEigenvalue);
  Eigen::VectorXd v1 = Eigen::VectorXd::Random(n);
  result_lanczos<Eigen::MatrixXd> res = lanczos_ir(A, m, v1, std::ceil(m/2), 1.0, 1e-7);
  mos << "Eigen calculated eigenvalues" << std::endl;
  for( int i =0; i < m; ++i) std::cout << eigenvalues[i] << ", ";
  std::cout << std::endl;
  mos << "Lanczos calculated eigenvalues" << std::endl;
  for( int i =0; i < m; ++i) std::cout << res.ev[i] << ", ";
  std::cout << std::endl;
  EXPECT_LE(std::abs(res.ev[0] - eigenvalues[0].real()), 1e-8 );

}

TEST(LANCZOS_IR, CalculateFromRandomSparse) {
  std::srand(std::time(nullptr));
  int m = std::rand() % 5 + 5;
  int n = m * 2;
  Eigen::SparseMatrix<double> A(n,n); //= Eigen::SparseMatrix(n,n);
  for (int i=0; i<2*n; i++) {
    int j = std::rand() % n;
    int k = std::rand() % n;
    double x = std::rand();
    A.coeffRef(j,k) = x;
    A.coeffRef(k,j) = x;
  }
  Eigen::MatrixXd ADense = Eigen::MatrixXd::Identity(n,n) * A;
  Eigen::EigenSolver<Eigen::MatrixXd> es(ADense);
  auto eigenvalues = es.eigenvalues();
  std::sort(eigenvalues.begin(), eigenvalues.end(), greaterEigenvalue);
  Eigen::VectorXd v1 = Eigen::VectorXd::Random(n);
  result_lanczos<Eigen::MatrixXd> res = lanczos_ir(A, m, v1, std::ceil(m/2), 1.0, 1e-7);
  mos << "Eigen calculated eigenvalues" << std::endl;
  for( int i =0; i < m; ++i) std::cout << eigenvalues[i] << ", ";
  std::cout << std::endl;
  mos << "Lanczos calculated eigenvalues" << std::endl;
  for( int i =0; i < m; ++i) std::cout << res.ev[i] << ", ";
  std::cout << std::endl;
  EXPECT_LE(std::abs(res.ev[0] - eigenvalues[0].real()), 1e-8 );

}

TEST(SIMPLE_LANCZOS, CalculateFromDenseDiagonal) {
  std::srand(std::time(nullptr));
  int m = std::rand() % 5 + 5;
  int n = m * 2;
  Eigen::MatrixXd A = Eigen::MatrixXd::Identity(n,n);
  A.diagonal().setRandom();
  Eigen::EigenSolver<Eigen::MatrixXd> es(A);
  auto eigenvalues = es.eigenvalues();
  std::sort(eigenvalues.begin(), eigenvalues.end(), greaterEigenvalue);
  Eigen::VectorXd v1 = Eigen::VectorXd::Random(n);
  result_lanczos<Eigen::MatrixXd> res = lanczos_ir(A, m, v1, std::ceil(m/2), 1.0, 1e-7);
  mos << "Eigen calculated eigenvalues" << std::endl;
  for( int i =0; i < m; ++i) std::cout << eigenvalues[i] << ", ";
  std::cout << std::endl;
  mos << "Lanczos calculated eigenvalues" << std::endl;
  for( int i =0; i < m; ++i) std::cout << res.ev[i] << ", ";
  std::cout << std::endl;
  EXPECT_LE(std::abs(res.ev[0] - eigenvalues[0].real()), 1e-8 );

}

TEST(SIMPLE_LANCZOS, CalculateFromDenseDiagonalWithoutReorthogonalization) {
  std::srand(std::time(nullptr));
  int m = std::rand() % 5 + 5;
  int n = m * 2;
  Eigen::MatrixXd A = Eigen::MatrixXd::Identity(n,n);
  A.diagonal().setRandom();
  Eigen::EigenSolver<Eigen::MatrixXd> es(A);
  auto eigenvalues = es.eigenvalues();
  std::sort(eigenvalues.begin(), eigenvalues.end(), greaterEigenvalue);
  Eigen::VectorXd v1 = Eigen::VectorXd::Random(n);
  result_lanczos<Eigen::MatrixXd> res = lanczos_ir(A, m, v1, std::ceil(m/2), 1.0, 1e-7);
  mos << "Eigen calculated eigenvalues" << std::endl;
  for( int i =0; i < m; ++i) std::cout << eigenvalues[i] << ", ";
  std::cout << std::endl;
  mos << "Lanczos calculated eigenvalues" << std::endl;
  for( int i =0; i < m; ++i) std::cout << res.ev[i] << ", ";
  std::cout << std::endl;
  EXPECT_LE(std::abs(res.ev[0] - eigenvalues[0].real()), 1e-8 );

}

TEST(SIMPLE_LANCZOS, CalculateFromRandomFullDense) {
  std::srand(std::time(nullptr));
  int m = std::rand() % 5 + 5;
  int n = m * 2;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(n,n);
  A = A + A.transpose().eval();
  Eigen::EigenSolver<Eigen::MatrixXd> es(A);
  auto eigenvalues = es.eigenvalues();
  std::sort(eigenvalues.begin(), eigenvalues.end(), greaterEigenvalue);
  Eigen::VectorXd v1 = Eigen::VectorXd::Random(n);
  result_lanczos<Eigen::MatrixXd> res = lanczos_ir(A, m, v1, std::ceil(m/2), 1.0, 1e-7);
  mos << "Eigen calculated eigenvalues" << std::endl;
  for( int i =0; i < m; ++i) std::cout << eigenvalues[i] << ", ";
  std::cout << std::endl;
  mos << "Lanczos calculated eigenvalues" << std::endl;
  for( int i =0; i < m; ++i) std::cout << res.ev[i] << ", ";
  std::cout << std::endl;
  EXPECT_LE(std::abs(res.ev[0] - eigenvalues[0].real()), 1e-8 );

}

TEST(SIMPLE_LANCZOS, CalculateFromRandomSparse) {
  std::srand(std::time(nullptr));
  int m = std::rand() % 5 + 5;
  int n = m * 2;
  Eigen::SparseMatrix<double> A(n,n); //= Eigen::SparseMatrix(n,n);
  for (int i=0; i<2*n; i++) {
    int j = std::rand() % n;
    int k = std::rand() % n;
    double x = std::rand();
    A.coeffRef(j,k) = x;
    A.coeffRef(k,j) = x;
  }
  Eigen::MatrixXd ADense = Eigen::MatrixXd::Identity(n,n) * A;
  Eigen::EigenSolver<Eigen::MatrixXd> es(ADense);
  auto eigenvalues = es.eigenvalues();
  std::sort(eigenvalues.begin(), eigenvalues.end(), greaterEigenvalue);
  Eigen::VectorXd v1 = Eigen::VectorXd::Random(n);
  result_lanczos<Eigen::MatrixXd> res = lanczos_ir(A, m, v1, std::ceil(m/2), 1.0, 1e-7);
  mos << "Eigen calculated eigenvalues" << std::endl;
  for( int i =0; i < m; ++i) std::cout << eigenvalues[i] << ", ";
  std::cout << std::endl;
  mos << "Lanczos calculated eigenvalues" << std::endl;
  for( int i =0; i < m; ++i) std::cout << res.ev[i] << ", ";
  std::cout << std::endl;
  EXPECT_LE(std::abs(res.ev[0] - eigenvalues[0].real()), 1e-8 );

}

