# Exercise 1.1 - Solution

An alternative way to derive finite difference formulas is to interprete the function whose derivative is discretized at the stencil points as a grid function to interpolate. The derivative of the resutling polynomial then takes the place of the derivative of the original function.

To approximate $u'(\bar{x})$ with this method one first needs to find the interpolating polynomial coefficients, this polynomial is then derived once and its derivative at the point $\bar{x}$ defines the finite difference scheme.

The interpolating polynomial can be expressed with respect to a basis such as the power basis, the Lagrange basis or the Newton basis. The latter is the most convenient because both finding the polynomial coefficients and evaluating the derivatives is easier, in fact one can use the divided difference algorithm (see, https://en.wikipedia.org/wiki/Newton_polynomial).

Let $p(x)$ be the (quadratic) polynomial that interpolates $u$ at $\bar{x}$, $\bar{x} - h$, $\bar{x} - 2h$, and let $u_0$, $u_{-1}$, $u_{-2}$ be the values of the function at those points. By using divided difference one gets the polynomial coefficients $a_0$, $a_1$, $a_2$ with respect to the Newton basis:

$$
  \begin{equation*}
  \begin{split}
  a_0 &= [u_{-2}] = u(\bar{x} - 2h) \\
  a_1 &= [u_{-2}, u_{-1}] = \frac{u_{-1} - u_{-2}}{h} \\
  a_2 &= [u_{-2}, u_{-1}, u_0] = \frac{\frac{u_0 - u_{-1}}{h} - a_1}{2h} 
      = \frac{u_{0} - 2u_{-1} + u_{-2}}{2h^2}
  \end{split}
  \end{equation*}
$$

Therefore the polynomial and its derivative are:

$$
  \begin{equation*}
  \begin{split}
  p(x) &= u_{-2} + \frac{u_{-1} - u_{-2}}{h}(x - \bar{x} + 2h)
       + \frac{u_{0} - 2u_{-1} + u_{-2}}{2h^2}(x - \bar{x} + h)(x - \bar{x} + 2h) \\
  p'(x) &= \frac{u_{-1} - u_{-2}}{h} 
        + \frac{u_{0} - 2u_{-1} + u_{-2}}{2h^2}(2x - 2\bar{x} + 3h) \\
  p'(\bar{x}) &= \frac{u_{-1} - u_{-2}}{h} 
        + \frac{3}{2h}(u_{0} - 2u_{-1} + u_{-2}) 
        = \frac{3u(\bar{x}) - 4u(\bar{x} - h) + u(\bar{x} - 2h)}{2h}
  \end{split}
  \end{equation*}
$$

The last formula is known as BDF2. This is the expected result. However when using this method one loses the information on the order of accuracy that comes out from using Taylor approximations
