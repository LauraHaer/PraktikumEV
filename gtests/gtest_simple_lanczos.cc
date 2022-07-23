#include <gtest/gtest.h>

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <ctime>
#include <cstdlib>

#include "Lanczos.hh"
#include "gtest_helpfunctions.hh"
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
  std::vector<double> error = runLanczos(A,v1, 5, 0, 1, 0, true, true);
  EXPECT_LE(error.at(0), 1e-15);

  //Two EVs
  v1 = Eigen::VectorXd{{-1.0, -1.0, 1.0, 0, 2.0}};
  mos << PRINT_REFLECTION(v1) << std::endl;
  error = runLanczos(A,v1, 5, 0, 1, 0, true, true);
  EXPECT_LE(error.at(0), 1e-15);
  EXPECT_LE(error.at(1), 1e-15);
}

TEST(SIMPLE_LANCZOS, CalculateFromDenseDiagonal) {
  std::srand(std::time(nullptr));
  int m = std::rand() % 5 + 5;
  int n = m * 2;
  Eigen::MatrixXd A = CreateRandomDiagonal(n);
  Eigen::VectorXd v1 = Eigen::VectorXd::Ones(n);

  std::vector<double> error = runLanczos(A,v1, m, 0, 1, 0, true, true);
  EXPECT_LE(error.at(0), 1e-8 );

}

TEST(SIMPLE_LANCZOS, CalculateFromRandomFullDense) {
  std::srand(std::time(nullptr));
  int m = std::rand() % 5 + 5;
  int n = m * 2;
  Eigen::MatrixXd A = CreateStdRandom(n);
  Eigen::VectorXd v1 = Eigen::VectorXd::Random(n);

  std::vector<double> error = runLanczos(A,v1, m, 0, 1, 0, true, true);
  EXPECT_LE(error.at(0), 1e-8 );
}


TEST(SIMPLE_LANCZOS, CalculateFromRandomFullDenseGoodStart) {
  std::srand(std::time(nullptr));
  int m = std::rand() % 5 + 5;
  int n = m * 2;
  Eigen::MatrixXd A = CreateStdRandom(n);
  Eigen::VectorXd v1 = CreateGoodStartVector(A, 4);

  std::vector<double> error = runLanczos(A,v1, m, 0, 1, 0, true, true);
  EXPECT_LE(error.at(0), 1e-8 );
  EXPECT_LE(error.at(1), 1e-8 );
  EXPECT_LE(error.at(2), 1e-8 );
  EXPECT_LE(error.at(3), 1e-8 );
}


TEST(SIMPLE_LANCZOS, CalculateFromLargeRandomFullDense) {
  std::srand(std::time(nullptr));
  Eigen::MatrixXd A = CreateStdRandom(1000);
  Eigen::VectorXd v1 = Eigen::VectorXd::Random(1000);

  std::vector<double> error = runLanczos(A,v1, 30, 0, 1, 0, true, true);
  EXPECT_LE(error.at(0), 1e-8 );
}


TEST(SIMPLE_LANCZOS, CalculateFromLargeRandomFullDenseGoodStart) {
  std::srand(std::time(nullptr));
  Eigen::MatrixXd A = CreateStdRandom(1000);
  Eigen::VectorXd v1 = CreateGoodStartVector(A, 10);

  std::vector<double> error = runLanczos(A,v1, 30, 0, 1, 0, true, true);
  EXPECT_LE(error.at(0), 1e-8 );
}

TEST(SIMPLE_LANCZOS, DISABLED_CalculateFromRandomSparse) {
  std::srand(std::time(nullptr));
  int m = std::rand() % 5 + 5;
  int n = m * 2;
  Eigen::SparseMatrix<double> A = CreateRandomSparse(n, n/3*n);
  Eigen::VectorXd v1 = Eigen::VectorXd::Random(n);

  std::vector<double> error = runLanczos(A,v1, m, 0, 1, 0, true, true);
  EXPECT_LE(error.at(0), 1e-8 );
}

