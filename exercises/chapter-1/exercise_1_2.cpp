// g++ --std=c++23 -g -Wall -Wextra -Wpedantic exercise_1_2.cpp -lginac -lcln &&
// ./a.out

#include <cassert>
#include <ginac/ginac.h>
#include <ginac/matrix.h>
#include <iostream>
#include <string>

static void basic_example() {
  using namespace GiNaC;

  GiNaC::symbol x("x");
  GiNaC::symbol y("y");

  GiNaC::ex poly;

  for (auto i{0}; i < 3; ++i) {
    poly += GiNaC::factorial(i) * GiNaC::pow(x, i) * GiNaC::pow(y, 2 - i);
  }

  std::cout << poly << '\n';
}
//
// Compute the 5 points centered finite difference coefficients to approximate
// the order 2 deriviative
static void linear_system_example() {
  using namespace GiNaC;

  matrix A{{symbol("a"), 1}, {-1, 1}};

  matrix x(2, 1);
  x(0, 0) = symbol("x");
  x(1, 0) = symbol("y");

  matrix b(2, 1);
  b(0, 0) = 1;
  b(1, 0) = 0;

  std::cout << A.solve(x, b) << '\n';
  std::cout << A.mul(x) << '\n';
  std::cout << b << '\n';
}

// Compute the 5 points centered finite difference coefficients to
// approximate the order 2 deriviative
static void second_derivative_centered_5() {
  using namespace GiNaC;

  constexpr auto order{2};
  constexpr auto size{5};

  assert(size > order);

  symbol h{"h"};
  matrix vandermonde{size, size};

  for (auto i{0}; i < size; ++i) {
    for (auto j{0}; j < size; ++j) {
      vandermonde(i, j) =
          i > 0 ? pow((-size / 2 + j) * h, i) / factorial(i) : 1;
    }
  }

  matrix coeffs(size, 1);
  // TODO: should be doable with symbolic_matrix but it does not work
  for (auto i{0}; i < size; ++i)
    coeffs(i, 0) = symbol("x" + std::to_string(i));

  matrix rhs(size, 1);
  rhs(order, 0) = 1;

  auto res{vandermonde.solve(coeffs, rhs)};

  std::cout << (res * pow(h, order)).evalm() << '\n';
}

int main() {
  // basic_example();
  // linear_system_example();
  second_derivative_centered_5();
}
