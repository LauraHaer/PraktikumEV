#ifndef GTEST_HELPFUNCTIONS_HH
#define GTEST_HELPFUNCTIONS_HH
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <vector>

#include "Lanczos.hh"

Eigen::SparseMatrix<double> CreateRandomSparse(
    const int aN, const int aEntries, const int aSeed = std::time(nullptr));

Eigen::MatrixXd CreateRandomDense(const int aN,
                                  const int aSeed = std::time(nullptr));

Eigen::MatrixXd CreateStdRandom(const int aN,
                                  const int aSeed = std::time(nullptr));

Eigen::VectorXd CreateStdRandomVector(const int aN,
                                  const int aSeed = std::time(nullptr));

Eigen::MatrixXd CreateRandomDiagonal(const int aN,
                                     const int aSeed = std::time(nullptr));

Eigen::VectorXd CreateGoodStartVector(const Eigen::MatrixXd A, int n = 0);

template <class aMat, class aVec>
std::vector<double> runLanczos(const aMat& A, const aVec& aV1, const int m,
                               const int k = 0, const double aRho = 1,
                               const double aEps = 1e-8,
                               const bool aSimple = false,
                               const bool aPrint_result = false) {
  typedef Eigen::Matrix<typename aMat::value_type, -1, -1> tmpMat;
  result_lanczos<Eigen::MatrixXd> res;
  auto start_time = std::chrono::high_resolution_clock::now();
  if (aSimple) {
    res = simple_lanczos(A, m, aV1, aRho);
  } else {
    res = lanczos_ir(A, m, aV1, k, aRho, aEps);
  }
  auto runtime = std::chrono::high_resolution_clock::now() - start_time;
  tmpMat DenseA = A * tmpMat::Identity(A.cols(), A.cols());
  Eigen::EigenSolver<tmpMat> es(DenseA);
  auto comp_eigenvalues = es.eigenvalues();
  std::vector<double> eigenvalue_error;
  std::vector<double> min_eigenvalue_error;
  std::sort(comp_eigenvalues.begin(), comp_eigenvalues.end(),
            greaterEigenvalue);
  auto eigenvalues = comp_eigenvalues.real();

  for (int i = 0; i < res.ev.rows(); ++i) {
    eigenvalue_error.push_back(std::abs(eigenvalues(i) - res.ev(i)));
    double min = std::abs(eigenvalues(0) - res.ev(i));
    for (auto x : eigenvalues) {
      if (std::abs(x - res.ev(i)) < min) {
        min = std::abs(x - res.ev(i));
      }
    }
    min_eigenvalue_error.push_back(min);
  }
  if (aPrint_result) {
    std::cout
        << "Test successfully executed!"
        << std::endl
        << "Runtime: "
        << std::chrono::duration_cast<std::chrono::seconds>(runtime).count()
        << " sec" << std::endl
        << "Eigenvalues: ";
    for (int i = 0; i < res.ev.rows(); ++i) {
      std::cout << res.ev(i) << ", ";
    }
    std::cout << std::endl << "Eigenvalue Error: ";
    for (auto x : eigenvalue_error) {
      std::cout << x << ", ";
    }
    std::cout << std::endl << "Min Eigenvalue Error: ";
    for (auto x : min_eigenvalue_error) {
      std::cout << x << ", ";
    }
    std::cout << std::endl;
  }
  if (aSimple) {
    return min_eigenvalue_error;
  }

  return eigenvalue_error;
}

template <class Mat>
Mat CreateLaplaceMatrix(int n) {
  typedef Eigen::Triplet<int> T;
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
