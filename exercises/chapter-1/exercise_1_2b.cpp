#include "fornberg.hpp"
#include <cstdlib>
#include <eigen3/Eigen/Core>
#include <iostream>

int main() {
  using Matrix =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using Vector = Eigen::VectorXd;

  constexpr auto size{5};

  Vector grid{Vector::LinSpaced(size, -2.0, 2.0)};
  Matrix coeffs{Matrix::Zero(size, size)};

  if (!fornberg(grid, 0.0, coeffs))
    std::exit(EXIT_FAILURE);

  std::cout << coeffs << '\n';
}
