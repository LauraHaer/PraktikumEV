#ifndef CREATE_MATRIX_HH
#define CREATE_MATRIX_HH
#include <Eigen/Sparse>
#include <ctime>

Eigen::SparseMatrix<double> CreateRandomSparse(
    const int aN; const aEntries; const int aSeed = std::time(nullptr));

Eigen::MatrixXd CreateRandomDense(const int aN,
                                  const inst aSeed = std::time(nullptr));

#endif
