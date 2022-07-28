#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

#include "Lanczos.hh"
#include "gtest_helpfunctions.hh"
#include "helpfunctions.hh"

TEST(SIMPLE_LANCZOS, CalculateSmallEVWithGivenEigenvector) {
  // int n = 5;
  Eigen::MatrixXd A{{4.0, -1.0, 0.0, 0.0, -1.0},
                    {-1.0, 4.0, -1.0, 0.0, 0.0},
                    {0.0, -1.0, 4.0, -1.0, 0.0},
                    {0.0, 0.0, -1.0, 4.0, 0.0},
                    {-1.0, 0.0, 0.0, 0.0, 4.0}};
  Eigen::VectorXd v1{{-1.0, 0.0, 1.0, -1.0, 1.0}};
  std::vector<double> error = runLanczos(A, v1, 5, 0, 1, 0, true, true);
  EXPECT_LE(error.at(0), 1e-13);

  // Two EVs
  v1 = Eigen::VectorXd{{-1.0, -1.0, 1.0, 0, 2.0}};
  error = runLanczos(A, v1, 5, 0, 1, 0, true, true);
  EXPECT_LE(error.at(0), 1e-13);
  EXPECT_LE(error.at(1), 1e-13);
}

TEST(SIMPLE_LANCZOS, CalculateFromDenseDiagonal) {
  std::srand(std::time(nullptr));
  int m = std::rand() % 5 + 5;
  int n = m * 2;
  Eigen::MatrixXd A = CreateRandomDiagonal(n);
  Eigen::VectorXd v1 = Eigen::VectorXd::Ones(n);

  std::vector<double> error = runLanczosN(A, 60, m, 0, 1, 0, true, true);
  EXPECT_LE(error.at(0), 1e-1);
}

TEST(SIMPLE_LANCZOS, CalculateFromRandomFullDense) {
  std::srand(std::time(nullptr));
  int m = std::rand() % 5 + 5;
  int n = m * 2;
  Eigen::MatrixXd A = CreateStdRandom(n);
  Eigen::VectorXd v1 = Eigen::VectorXd::Random(n);

  std::vector<double> error = runLanczosN(A, 6, m, 0, 1, 0, true, true);
  EXPECT_LE(error.at(0), 1e-6);
}

TEST(SIMPLE_LANCZOS, CalculateFromLargeRandomFullDense) {
  std::srand(std::time(nullptr));
  Eigen::MatrixXd A = CreateStdRandom(1000);
  Eigen::VectorXd v1 = Eigen::VectorXd::Random(1000);

  std::vector<double> error = runLanczosN(A, 4, 30, 0, 1, 0, true, true);
  EXPECT_LE(error.at(0), 1e-8);
}

TEST(SIMPLE_LANCZOS, DISABLED_CalculateFromRandomSparse) {
  std::srand(std::time(nullptr));
  int m = std::rand() % 5 + 5;
  int n = m * 2;
  Eigen::SparseMatrix<double> A = CreateRandomSparse(n, n / 3 * n);

  std::vector<double> error = runLanczosN(A, 4, m, 0, 1, 0, true, true);
  EXPECT_LE(error.at(0), 1e-8);
}
