#include <gtest/gtest.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Sparse>

#include "gtest_helpfunctions.hh"


TEST(LANCZOS_IR_EV, CalculateFromLaplaceMatrix) {
  //Define the parameters
  std::vector<int> MatrixSize{4,6,13};
  std::vector<int> NumberOfEigenvalues{2,4,7};
  std::vector<int> NumberOfTmpEigenvalues{10,15,20};
  bool print_result = true;

  for(int i= 0; i < (int)MatrixSize.size(); ++i ) {
    Eigen::SparseMatrix<double> A = CreateLaplaceMatrix<Eigen::SparseMatrix<double>>(MatrixSize.at(i));
    std::vector<double> error = runLanczosN(A, 10, NumberOfTmpEigenvalues.at(i), NumberOfEigenvalues.at(i), 1, 1e-8, false, print_result);
    EXPECT_LE(*std::max_element(error.begin(), error.end()), 1e-10);
  }
}


TEST(LANCZOS_IR_EV, CalculateFromRandomSparse) {
  std::vector<int> MatrixSize{15,40,200};
  std::vector<int> NumberOfEntries{15,60,800};
  std::vector<int> NumberOfEigenvalues{2,6,10};
  std::vector<int> NumberOfTmpEigenvalues{10,20,50};
  bool print_result = true;

  for(int i= 0; i < (int)MatrixSize.size(); ++i ) {
    Eigen::SparseMatrix<double> A = CreateRandomSparse(MatrixSize.at(i), NumberOfEntries.at(i), 42*i);
    std::vector<double> error = runLanczosN(A, 10, NumberOfTmpEigenvalues.at(i), NumberOfEigenvalues.at(i), 1, 1e-8, false, print_result);
    EXPECT_LE(*std::max_element(error.begin(), error.end()), 1e-10);
  }
}


TEST(LANCZOS_IR_EV, CalculateFromRandomense) {
  std::vector<int> MatrixSize{15,40,200};
  std::vector<int> NumberOfEigenvalues{2,6,6};
  std::vector<int> NumberOfTmpEigenvalues{10,20,20};
  bool print_result = true;

  for(int i= 0; i < (int)MatrixSize.size(); ++i ) {
    auto A = CreateRandomDense(MatrixSize.at(i), 42*i);
    std::vector<double> error = runLanczosN(A, 10, NumberOfTmpEigenvalues.at(i), NumberOfEigenvalues.at(i), 1, 1e-8, false, print_result);
    EXPECT_LE(*std::max_element(error.begin(), error.end()), 1e-10);
  }
}


TEST(LANCZOS_IR_EV, CalculateFromRandomDiagonal) {
  std::vector<int> MatrixSize{15,40,60};
  std::vector<int> NumberOfEigenvalues{2,6,10};
  std::vector<int> NumberOfTmpEigenvalues{10,20,50};
  bool print_result = true;

  for(int i= 0; i < (int)MatrixSize.size(); ++i ) {
    auto A = CreateRandomDiagonal(MatrixSize.at(i), 42*i);
    std::vector<double> error = runLanczosN(A, 10, NumberOfTmpEigenvalues.at(i), NumberOfEigenvalues.at(i), 1, 1e-8, false, print_result);
    EXPECT_LE(*std::max_element(error.begin(), error.end()), 1e-10);
  }
}

