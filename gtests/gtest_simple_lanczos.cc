#include <gtest/gtest.h>

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <ctime>
#include <cstdlib>

#include "Lanczos.hh"
#include "Lanczos_saad.hh"
#include "helpfunctions.hh"
#include "standard_include.hh"

TEST(SIMPLE_LANCZOS, CalculateSmallEVWithGivenEigenvector) {
  int n = 5;
  Eigen::MatrixXd A{{4.0, -1.0, 0.0, 0.0, -1.0},
                    {-1.0, 4.0, -1.0, 0.0, 0.0},
                    {0.0, -1.0, 4.0, -1.0, 0.0},
                    {0.0, 0.0, -1.0, 4.0, 0.0},
                    {-1.0, 0.0, 0.0, 0.0, 4.0}};
  Eigen::VectorXd v1{{-1.0, 0.0, 1.0, -1.0, 1.0}};
  //result_lanczos<Eigen::MatrixXd> res = lanczos_ir(A, n, v1, 1);
  std::cout << v1 << std::endl;
  result_lanczos<Eigen::MatrixXd> res = simple_lanczos(A, n, v1, 0.9);

  std::cout << res.ev << std::endl;
  EXPECT_DOUBLE_EQ(res.ev[0], 5.0);
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

TEST(SIMPLE_LANCZOS, DISABLED_CalculateFromDenseDiagonalWithoutReorthogonalization) {
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

TEST(SIMPLE_LANCZOS, DISABLED_CalculateFromRandomSparse) {
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

