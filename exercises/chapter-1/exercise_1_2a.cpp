#include <cassert>
#include <ginac/ginac.h>
#include <ginac/matrix.h>
#include <iostream>
#include <string>

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
  // TODO: should be doable with symbolic_matrix? but it does not work
  for (auto i{0}; i < size; ++i)
    coeffs(i, 0) = symbol("x" + std::to_string(i));

  matrix rhs(size, 1);
  rhs(order, 0) = 1;

  auto res{vandermonde.solve(coeffs, rhs)};

  std::cout << (res * pow(h, order)).evalm() << '\n';
}

int main() {
  second_derivative_centered_5();
}
