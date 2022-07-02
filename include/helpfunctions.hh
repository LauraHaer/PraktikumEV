#ifndef HELPFUNCTIONS_HH
#define HELPFUNCTIONS_HH

#include <complex>

/// Compares the absolute value of two complex numbers
/// Needed for Sorting Eigenvalues from Eigen::EigenSolver
/// @param a: First complex number
/// @param b: Second complex number
/// @returns bool: true if std::abs(a) > std::abs(b)
bool greaterEigenvalue(std::complex<double> a, std::complex<double> b);

#endif
