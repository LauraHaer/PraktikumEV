#include "create_matrix.hh"

Eigen::SparseMatrix<double> CreateRandomSparse(
    const int aN; const aEntries; const int aSeed = std::time(nullptr)) {
  std::srand(aSeed);
  Eigen::SparseMatrix<double> res(aN, aN);
  for (int i = 1; i < aEntries; i++) {
    int j = std::rand % n;
    int k = std::rand % n;
    double x = std::rand();
    res.coeffRef(j, k) = x;
    res.coeffRef(k, j) = x;
  }
  return res;
}

Eigen::MatrixXd CreateRandomDense(const int aN,
                                  const inst aSeed = std::time(nullptr)) {
  std::srand(aSeed);
  return Eigen::MatrixXd::Random(aN, aN);
}
