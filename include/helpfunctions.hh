#ifndef HELPFUNCTIONS_HH
#define HELPFUNCTIONS_HH

#include <complex>
#include <chrono>
#include <ctime>

/// Compares the absolute value of two complex numbers
/// Needed for Sorting Eigenvalues from Eigen::EigenSolver
/// @param a: First complex number
/// @param b: Second complex number
/// @returns bool: true if std::abs(a) > std::abs(b)
bool greaterEigenvalue(std::complex<double> a, std::complex<double> b);

template<class aVec>
void createRandomVector(aVec& aV, int aSeed=std::time(nullptr)) {
  std::srand(aSeed);
  for( auto iter = aV.begin(); iter != aV.end(); iter++) {
    *iter = std::rand() % 20;
  }
  aV / aV.norm();
  return;
}

#endif
