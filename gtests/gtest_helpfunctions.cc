#include "gtest_helpfunctions.hh"

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
  auto A = Eigen::MatrixXd::Random(aN, aN);
  return A + A.transpose().eval();
}

Eigen::MatrixXd CreateRandomDiagonal(const int aN, const int aSeed) {
  std::srand(aSeed);
  auto A = Eigen::MatrixXd::Zero(aN, aN);
  auto diag = Eigen::VectorXd::Random(aN);
  A.diagonal() = diag;
  return A;
}
