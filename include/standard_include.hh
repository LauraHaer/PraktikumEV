#pragma once

#include <cassert>
#include <cinttypes>
#include <cmath>
#include <iostream>
#include <string>

// implemented operator<< for some container classes in std:
#include <vector>
#include <list>
#include <unordered_set>
#include <map>
#include <unordered_map>

// command line parser
// https://github.com/CLIUtils/CLI11
#if __has_include("infra/CLI11.hpp")
#include "infra/CLI11.hpp"
#endif
/*
 * Usage:
 * struct Cb {
 *  int i;
 *  std::string s;
 * };
 *
 * int main (int argc, char** argv) {
 *   CLI::App argparser;
 *
 *   std::cout << argv[0] << std::endl;
 *   Cb cb;
 *
 *   // flag. optional.
 *   bool myflag = false;
 *   argparser.add_flag("-f,--flag", myflag, "it is a flag");
 *
 *   // "Option" parameter. required, expects three arguments
 *   std::vector<std::string> mystrings = {};
 *   argparser.add_option("-s,--stringvar", mystrings, "it is a string vector of length 3")->required()->expected(3);
 *
 *   // positional parameters
 *   argparser.add_option("one", cb.i, "it is a positional int")->required();
 *   int j = 0;
 *   argparser.add_option("two", j, "it is a positional int");
 *
 *   // parse a (lambda) function
 *   int a, b;
 *   argparser.add_option_function<int>("--funcparse", [&a, &b] (const int& arg) {a = b = arg;}, "parse via lambda");
 *
 *   bool funcbool = true;
 *   argparser.add_flag_function("--functionboolparse",
 *                               [&funcbool](int count) {assert(count > 0); funcbool = false;},
 *                               "flag via lambda");
 *
 *   CLI11_PARSE(argparser, argc, argv);
 */


#define mos std::cout << __FILE__ << ":" << __LINE__ << ": "

#define PRINT_REFLECTION(X) #X << ": " << X


using byte_t = uint8_t;


// === Functions to print containers from std library ===

template <typename TFirst, typename TSecond>
std::ostream& operator<< (std::ostream& os, const std::pair<TFirst, TSecond>& P) {
  os << P.first << ":" << P.second;
  return os;
}


template <typename TContainerIter>
std::ostream& __CONTAINER_ITER_PRINT(std::ostream& os, TContainerIter beginiter, TContainerIter enditer) {
  os << '[';
  if (beginiter != enditer) {
    os << *beginiter;
  }
  ++beginiter;
  for (; beginiter != enditer; ++beginiter) {
    os << ", " << *beginiter;
  }
  os << ']';
  return os;
}


template <typename TContainer>
std::ostream& __CONTAINER_PRINT(std::ostream& os, const TContainer& C) {
  return __CONTAINER_ITER_PRINT(os, C.begin(), C.end());
}


#ifndef INFRA_VECTOR_OPERATORS_HH
template <typename T>
std::ostream& operator<< (std::ostream& os, const std::vector<T>& V) {
  return __CONTAINER_PRINT(os, V);
}
#endif


template <typename T>
std::ostream& operator<< (std::ostream& os, const std::list<T>& V) {
  return __CONTAINER_PRINT(os, V);
}


template <typename T>
std::ostream& operator<< (std::ostream& os, const std::unordered_set<T>& S) {
  return __CONTAINER_PRINT(os, S);
}


template <typename T1, typename T2>
std::ostream& operator<< (std::ostream& os, const std::map<T1, T2>& M) {
  return __CONTAINER_PRINT(os, M);
}


template <typename T1, typename T2>
std::ostream& operator<< (std::ostream& os, const std::unordered_map<T1, T2>& M) {
  return __CONTAINER_PRINT(os, M);
}

