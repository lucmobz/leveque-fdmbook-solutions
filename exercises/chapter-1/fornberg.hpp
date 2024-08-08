#include <eigen3/Eigen/Dense>

// Fills coeffs with coefficients to approximate all derivatives from 0 to
// points.size() - 1, with best possible order of accuracy using the finite
// difference method. The implementation uses the explicit Fornberg recursive
// algorithm on arbitrary grids
template <typename Scalar, int Storage>
auto fornberg(
    const Eigen::Vector<Scalar, Eigen::Dynamic> &grid, const Scalar center,
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Storage> &coeffs)
    -> bool {
  const auto size{grid.size()};

  if (size <= 0)
    return false;

  coeffs.resize(size - 1, size - 1);
  coeffs.setZero();
  coeffs(0, 0) = Scalar{1};

  auto prev_denom{1.0};

  // The algorithm works by updating a table, each value influences up to 4
  // cells of the next iteration of the table
  for (auto j{1}; j < size; ++j) {
    const auto xbar_m_xj{center - grid(j)};

    // Loop over the inner table by column to fill in place the outer table
    for (auto i{0}; i < j - 1; ++i) {
      auto prev_coeff{coeffs(0, i)};
      coeffs(0, i) = Scalar{0};
      const auto xi_m_xj{grid(i) - grid(j)};

      for (auto k{0}; k < j; ++k) {
        const auto coeff{prev_coeff};

        coeffs(k, i) += k * coeff / xi_m_xj;
        coeffs(k + 1, i) = xbar_m_xj * coeff / xi_m_xj;

        prev_coeff = coeffs(k + 1, i);
      }
    }

    // Last column influences also the adjacent column
    auto prev_coeff{coeffs(0, j - 1)};
    coeffs(0, j - 1) = Scalar{0};
    const auto xjm1_m_xj{grid(j - 1) - grid(j)};
    const auto xbar_m_xjm1{center - grid(j - 1)};

    for (auto k{0}; k < j; ++k) {
      const auto coeff{prev_coeff};
      const auto denom{1.0}; // TODO: compute the denominator
      const auto ratio{prev_denom / denom};

      coeffs(k, j - 1) += k * coeff / xjm1_m_xj;
      coeffs(k + 1, j - 1) = xbar_m_xj * coeff / xjm1_m_xj;
      // These two will already be 0
      coeffs(k, j) = ratio * k * coeff;
      coeffs(k + 1, j) = ratio * xbar_m_xjm1 * coeff;

      prev_coeff = coeffs(k + 1, j - 1);
      prev_denom = denom;
    }
  }
}
