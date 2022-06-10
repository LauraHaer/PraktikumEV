#include <gtest/gtest.h>

#include <iostream>


#include <Eigen/Sparse>
#include <Eigen/SparseQR>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

double
ConvergenceCheck(SpMat A) {
  double max_val = 0;
  int size_of_matrix = A.rows();
  for ( int i = 0; i < size_of_matrix; i++ ) {
    for ( int j = 0; j < size_of_matrix; j++ ) {
      if ( i == j ) continue;
      max_val = std::max(abs(A.coeff(i, j)), max_val);
    }
  }
  return max_val;
}



SpMat
CalculateEigenvalues(SpMat A) {
  SpMat B = A;
  Eigen::SparseQR<SpMat,  Eigen::NaturalOrdering<int>> QRSolver;
  int i = 0;
  while ( ConvergenceCheck(A) > 1e-10 ) {
    i++;
    QRSolver.compute(A);
    EXPECT_EQ(QRSolver.info(), Eigen::Success);
    Eigen::SparseMatrix<double> Q;
    Q = QRSolver.matrixQ();
    // Q = Eigen::SparseQR<Eigen::SparseMatrix<double> >(A).matrixQ(); // False line the documentation
    A = Q.transpose() * A * Q;
  }
  return A;
}

TEST(EIGEN, TestEigenvaluesNoErrors) {
  int size_of_matrix = 5;
  SpMat  A(size_of_matrix, size_of_matrix);
  std::vector<T> entries;
  entries.push_back(T(0, 0, 2));
  for (int i=1; i < size_of_matrix; i++) {
    entries.push_back(T(i, i, 2));
    entries.push_back(T(i-1 , i, -1));
    entries.push_back(T(i, i-1, -1));
  }
  A.setFromTriplets(entries.begin(), entries.end());

  SpMat D = CalculateEigenvalues(A);
  EXPECT_EQ(size_of_matrix, 5);
}
