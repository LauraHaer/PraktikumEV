#include <gtest/gtest.h>

#include <iostream>

#include "dummy.h"
#include <Eigen/Dense>

using Eigen::MatrixXd;

/// Simple test function for confirming a successfull installation
TEST(Eigen, TestFunctioning)
{
  MatrixXd m(2, 2);
  m(0, 0) = 3;
  m(1, 0) = 2.5;
  m(0, 1) = -1;
  m(1, 1) = m(1, 0) + m(0, 1);
  //   EXPECT_EQ((unsigned int)2, i);
}



TEST(Eigen, LibraryFunctioning)
{
  test_eigen();
}
