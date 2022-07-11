#include <gtest/gtest.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Sparse>

#include "gtest_helpfunctions.hh"


TEST(LANCZOS_IR_EV, CalculateFromLaplaceMatrix) {
  //Define the parameters
  std::vector<int> MatrixSize{4,6,8};
  std::vector<int> NumberOfEigenvalues{2,6,10};
  std::vector<int> NumberOfTmpEigenvalues{10,20,20};
  bool print_result = true;

  for(int i= 0; i < (int)MatrixSize.size(); ++i ) {
    Eigen::SparseMatrix<double> A = CreateLaplaceMatrix<Eigen::SparseMatrix<double>>(MatrixSize.at(i));
    Eigen::VectorXd v1 = Eigen::VectorXd::Zero(MatrixSize.at(i) * MatrixSize.at(i));
    v1(0) = 1;
    std::vector<double> error = runLanczos(A, v1, NumberOfTmpEigenvalues.at(i), NumberOfEigenvalues.at(i), 1, 1e-8, false, print_result);
    EXPECT_LE(error.at(0), 1e-10);
    //EXPECT_LE(*std::max_element(error.begin(), error.end()), 1e-5);
  }
}


TEST(LANCZOS_IR_EV, CalculateFromRandomSparse) {
  std::vector<int> MatrixSize{15,40,60};
  std::vector<int> NumberOfEntries{15,60,240};
  std::vector<int> NumberOfEigenvalues{2,6,10};
  std::vector<int> NumberOfTmpEigenvalues{10,20,50};
  bool print_result = true;

  for(int i= 0; i < (int)MatrixSize.size(); ++i ) {
    Eigen::SparseMatrix<double> A = CreateRandomSparse(MatrixSize.at(i), NumberOfEntries.at(i), 42*i);
    Eigen::VectorXd v1 = Eigen::VectorXd::Zero(MatrixSize.at(i));
    v1(0) = 1;
    std::vector<double> error = runLanczos(A, v1, NumberOfTmpEigenvalues.at(i), NumberOfEigenvalues.at(i), 1, 1e-8, false, print_result);
    EXPECT_LE(error.at(0), 1e-10);
    //EXPECT_LE(*std::max_element(error.begin(), error.end()), 1e-5);
  }
}


TEST(LANCZOS_IR_EV, CalculateFromRandomDense) {
  std::vector<int> MatrixSize{15,40,60};
  std::vector<int> NumberOfEigenvalues{2,6,10};
  std::vector<int> NumberOfTmpEigenvalues{10,20,50};
  bool print_result = true;

  for(int i= 0; i < (int)MatrixSize.size(); ++i ) {
    auto A = CreateRandomDense(MatrixSize.at(i), 42*i);
    Eigen::VectorXd v1 = Eigen::VectorXd::Zero(MatrixSize.at(i));
    v1(0) = 1;
    std::vector<double> error = runLanczos(A, v1, NumberOfTmpEigenvalues.at(i), NumberOfEigenvalues.at(i), 1, 1e-8, false, print_result);
    EXPECT_LE(error.at(0), 1e-10);
    //EXPECT_LE(*std::max_element(error.begin(), error.end()), 1e-5);
  }
}


TEST(LANCZOS_IR_EV, CalculateFromRandomDiagonal) {
  std::vector<int> MatrixSize{15,40,60};
  std::vector<int> NumberOfEigenvalues{2,6,10};
  std::vector<int> NumberOfTmpEigenvalues{10,20,50};
  bool print_result = true;

  for(int i= 0; i < (int)MatrixSize.size(); ++i ) {
    auto A = CreateRandomDiagonal(MatrixSize.at(i), 42*i);
    Eigen::VectorXd v1 = Eigen::VectorXd::Zero(MatrixSize.at(i));
    v1(0) = 1;
    std::vector<double> error = runLanczos(A, v1, NumberOfTmpEigenvalues.at(i), NumberOfEigenvalues.at(i), 1, 1e-8, false, print_result);
    //EXPECT_LE(error.at(0), 1e-10);
    //EXPECT_LE(*std::max_element(error.begin(), error.end()), 1e-5);
  }
}

