#ifndef CREATE_LAPLACE_MATRIX_HH_
#define CREATE_LAPLACE_MATRIX_HH_
#include <Eigen/Sparse>
#include <vector>

typedef Eigen::Triplet<int> T;

template <class Mat>
Mat CreateLaplaceMatrix(int n) {
  std::vector<T> entries;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      entries.push_back(T(i * n + j, i * n + j, 4));
      if (i > 0) entries.push_back(T(i * n + j, (i - 1) * n + j, -1));
      if (i < n - 1) entries.push_back(T(i * n + j, (i + 1) * n + j, -1));
      if (j > 0) entries.push_back(T(i * n + j, i * n + j - 1, -1));
      if (j < n - 1) entries.push_back(T(i * n + j, i * n + j + 1, -1));
    }
  }
  Mat Matrix(n * n, n * n);
  Matrix.setFromTriplets(entries.begin(), entries.end());

  return Matrix;
}
#endif
