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

Eigen::MatrixXd CreateTridiagMatrix(const int n);

Eigen::VectorXd CreateGoodStartVector(const Eigen::MatrixXd A, int n = 0);

bool greaterError(std::vector<double> first, std::vector<double> last);

template <class aMat, class aVec>
std::vector<double> runLanczos(const aMat& A, const aVec& aV1, const int m,
                               const int k = 0, const double aRho = 1,
                               const double aEps = 1e-8,
                               const bool aSimple = false,
                               const bool aPrint_result = false) {
  typedef Eigen::Matrix<typename aMat::value_type, -1, -1> tmpMat;
  result_lanczos res;
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

template <class aMat>
std::vector<double> runLanczosN(const aMat& A, const int n, const int m,
                               const int k = 0, const double aRho = 1,
                               const double aEps = 1e-8,
                               const bool aSimple = false,
                               const bool aPrint_result = false) {
  typedef Eigen::Matrix<typename aMat::value_type, -1, -1> tmpMat;

  // Generate Eigenvalues with Lanczos
  std::vector<result_lanczos> res_vec(n);
  std::vector<std::chrono::duration<double>> runtime(n);
  for (int i = 0; i < n; ++i) {
    Eigen::VectorXd V1 = CreateStdRandomVector(A.rows(), std::time(nullptr));
    auto start_time = std::chrono::high_resolution_clock::now();
    if (aSimple) {
      res_vec.at(i) = simple_lanczos(A, m, V1, aRho);
    } else {
      res_vec.at(i) = lanczos_ir(A, m, V1, k, aRho, aEps);
    }
    runtime.at(i) = std::chrono::high_resolution_clock::now() - start_time;
  }

  //Find Eigenvalues using Eigen Package
  tmpMat DenseA = A * tmpMat::Identity(A.cols(), A.cols());
  auto start_time = std::chrono::high_resolution_clock::now();
  Eigen::EigenSolver<tmpMat> es(DenseA);
  auto comp_eigenvalues = es.eigenvalues();
  std::chrono::duration<double> runtime_eigen = std::chrono::high_resolution_clock::now() - start_time;
  std::sort(comp_eigenvalues.begin(), comp_eigenvalues.end(),
            greaterEigenvalue);
  auto eigenvalues = comp_eigenvalues.real();

  // Calculate the error
  std::vector<std::vector<double>> eigenvalue_error(n);
  std::vector<std::vector<double>> min_eigenvalue_error(n);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < res_vec.at(j).ev.rows(); ++i) {
      eigenvalue_error.at(j).push_back(std::abs(eigenvalues(i) - res_vec.at(j).ev(i)));
      double min = std::abs(eigenvalues(0) - res_vec.at(j).ev(i));
      for (auto x : eigenvalues) {
        if (std::abs(x - res_vec.at(j).ev(i)) < min) {
          min = std::abs(x - res_vec.at(j).ev(i));
        }
      }
      min_eigenvalue_error.at(j).push_back(min);
    }
  }

  std::cout << "NEW MATRIX" << std::endl << "Eigen Eigenvalues: ";
  for(auto x : eigenvalues) {
    std::cout << x << ", ";
  }
  std::vector<double> best_min_error = *std::max_element(min_eigenvalue_error.begin(), min_eigenvalue_error.end(), greaterError);
  std::cout << std::endl;
  std::cout << "Best approximation of any Eigenvalues: ";
  for (auto x : best_min_error) {
    std::cout << x << ", ";
  }
  std::cout << std::endl;

  for (int j = 0; j < n; ++j) {
    if (aPrint_result) {
      std::cout
          << "Test successfully executed!"
          << std::endl
          << "Runtime: "
          << std::chrono::duration_cast<std::chrono::milliseconds>(runtime.at(j)).count()
          << " ms" << std::endl
          << "Runtime Eigen: "
          << std::chrono::duration_cast<std::chrono::milliseconds>(runtime_eigen).count()
          << " ms" << std::endl
          << "Eigenvalues: ";
      for (int i = 0; i < res_vec.at(j).ev.rows(); ++i) {
        std::cout << res_vec.at(j).ev(i) << ", ";
      }
      std::cout << std::endl << "Eigenvalue Error: ";
      for (auto x : eigenvalue_error.at(j)) {
        std::cout << x << ", ";
      }
      std::cout << std::endl << "Min Eigenvalue Error: ";
      for (auto x : min_eigenvalue_error.at(j)) {
        std::cout << x << ", ";
      }
      std::cout << std::endl;
    }
  std::cout << std::endl;
  }

    //if (aSimple) {
    //  return min_eigenvalue_error;
    //}

  return best_min_error;
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
