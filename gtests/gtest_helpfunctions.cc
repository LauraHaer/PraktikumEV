#include "gtest_helpfunctions.hh"
#include <iostream>

Eigen::SparseMatrix<double> CreateRandomSparse(const int aN, const int aEntries,
                                               const int aSeed) {
  std::srand(aSeed);
  Eigen::SparseMatrix<double> res(aN, aN);
  for (int i = 1; i < aEntries; i++) {
    int j = std::rand() % aN;
    int k = std::rand() % aN;
    double x = std::rand();
    res.coeffRef(j, k) = x;
    res.coeffRef(k, j) = x;
  }
  return res;
}

Eigen::MatrixXd CreateRandomDense(const int aN, const int aSeed) {
  std::srand(aSeed);
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(aN, aN);
  A = A + A.transpose().eval();
  return A;
}

Eigen::MatrixXd CreateRandomDiagonal(const int aN, const int aSeed) {
  std::srand(aSeed);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(aN, aN);
  Eigen::VectorXd diag = Eigen::VectorXd::Random(aN);
  A.diagonal() = diag;
  return A;
}

Eigen::MatrixXd CreateStdRandom(const int aN, const int aSeed) {
  std::srand(aSeed);
  Eigen::MatrixXd A(aN, aN);
  for (int i = 0; i<aN; i++) {
    for (int j=0; j<= i; j++) {
      int x = std::rand();
      A(i,j) = x;
      A(j,i) = x;
    }
  }
  return A;
}