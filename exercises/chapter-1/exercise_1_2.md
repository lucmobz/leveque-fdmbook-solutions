# Exercise 1.2 - Solution


The method of undeterminate coefficients is based on Taylor approximations of $u(x_j)$ centered at the evaluation point $\bar{x}$, where $x_j$ for $j = 0 \dots n-1$, are the $n$ grid points (distinct and for simplicity in increasing order). Consider the difference (also the error) between the $k$-th derivative and the linear combination that approximates it:

$$
  u^{(k)}(\bar{x}) - \sum_{j=0}^{n - 1} c_j u(x_j)
$$

using Taylor expansions up to order $n-1$ (to have convergence the above difference must depend on some function of $h_j = x_j - \bar{x}$, so at least $n=k+1$ points are needed to have a chance to cancel the $k$-th derivative at $\bar{x}$)

$$
  u(x_j) = \sum_{i=0}^{n-1} u^{(i)}(\bar{x}) \frac{(h_j)^i}{i!} + O((h_j)^n)
$$

Substituing and rearranging one gets the following error:

$$
  e = u^{(k)}(\bar{x}) - \sum_{i=0}^{n-1} u^{(i)}(\bar{x}) \sum_{j=0}^{n-1} c_j \frac{(h_j)^i}{i!} + \sum_{j=0}^{n-1} c_j O((h_j)^n)
$$

Normally it turns out (except on some symmetric cases, where extra cancellation may happen and one needs to expand the Taylor series further) that $c_j = O(h^{-k})$ (assuming for simplicity the equispaced case where $x_j=\bar{x}$, for some $j$ and $h$ is the interval length) so the error will be of order $e=O(h^p)$, with $p=n-k$ if these conditions are met:

$$
  \sum_{j=0}^{n-1} c_j \frac{(h_j)^i}{i!} =
  \begin{cases}
    0, \quad \forall i \neq k \\
    1, \quad i = k
  \end{cases}
$$

So, the vector $c_j$ needs to solve an $n \times n$ linear system where the matrix is a Vandermonde matrix $a_{i,j}=\frac{(h_j)^i}{i!}$, and the right hand is all zeros except $b_k=1$.

## Part A

The *GiNaC* C++ library (https://www.ginac.de/) implements a computer algebra system (CAS) that can solve the linear system symbolically and give exact results in terms of rational numbers. On *WSL*, using *Ubuntu 22.04LTS*, it is available with `sudo apt install libginac-dev` (which also installs the dependency `libcln-dev` which *GiNaC* uses).

The program `exercise_1_2a.cpp`, compiled with `g++ --std=c++23 -g -Wall -Wextra -Wpedantic exercise_1_2a.cpp -lginac -lcln && ./a.out`, finds the following coefficients (multiplied by $h^2$): 
```shell
[[-1/12],[4/3],[-5/2],[4/3],[-1/12]]
```

## Fornberg Algorithm

The Vandermonde system is not well-conditioned when the number of points $n$ is large. A better algorithm is the Fornberg algorithm appearing in the 1998 paper (https://www.jstor.org/stable/2653239) which is a recursive algorithm.

There are two main intuitions behind this algorithm:

First, express the interpolating polynomial $p$, over the first $j+1$ nodes using the Lagrange basis $L_{i,j}$ to obtain the following formula for the $k$-th derivative at $\bar{x}$:

$$ p^{(k)}(\bar{x}) = \sum_{i=0}^n L_{i,j}^{(k)}(\bar{x}) u_i = \sum_{i=0}^n c_{i,j}^k u_i $$

$u_i$ is the function value at the $x_i$ node, and $c_{i,j}^k$ is the $i$-th coefficient when approximating the $k$-th derivative on $j+1$ nodes. Since $L_{i,j}$ are polynomials the Taylor expansion holds exactly (at all points) and so: 

$$ L_{i,j}(x) = \sum_{k=0}^j c_{i,j}^k \frac{(x - \bar{x})^k}{k!} $$

Second, the Lagrange basis over a set of $j+1$ nodes is related to the one over the same nodes minus the last one ($j$ nodes in total). Defining the denominator of $L_{j,j}$ as $D_j = \prod_{\nu=0}^{j-1} (x_j - x_{\nu})$ one gets: 

$$
 \begin{split}
    L_{i,j}(x) &= \frac{x - x_j}{x_i - x_j} L_{i, j - 1}(x), \quad \forall i < j \\
    L_{j,j}(x) &= \frac{D_{j-1}}{D_j}(x - x_{j-1}) L_{j-1, j-1}(x), \quad i = j 
  \end{split} 
$$ 

Substituting the Taylor expansions in the recursive relation, equating the polynomial members term by term and assuming:

* $c_{0, 0}^0 = 1$ (the coefficient to interpolate the function with the Lagrange basis at a single point)
* $c_{i,j}^k$ is defined if $0 \leq i \leq j$ and $0 \leq k \leq j$, else it is zero
* $D_0 = 1$ (denominator of the Lagrange polynomial referred to a single node)

one gets (remember to sum and subtract $\bar{x}$):

$$
  \begin{split}
    c_{i,j}^k &= \frac{(\bar{x} - x_j) c_{i,j-1}^k + k c_{i,j-1}^{k-1}}{x_i-x_j}, \quad \forall i < j \\
    c_{j,j}^k &= \frac{D_{j-1}}{D_j} \left((\bar{x} - x_{j-1}) c_{j-1,j-1}^k + k c_{j-1,j-1}^{k-1} \right), \quad i = j
  \end{split}
$$

## Part B
  
The program `fornberg.hpp` implements the Fornberg algorithm in the version where, given a grid of size $n+1$ and a point $\bar{x}$, it returns the coefficients to approximate all the derivatives of order $k$ for $k \leq n$ using the grid as stencil.

Compiling and running the file `exercise_1_2b.cpp` with `g++ --std=c++23 -g -Wall -Wextra -Wpedantic exercise_1_2b.cpp && ./a.out` (it requires the *Eigen* library which can be installed with `sudo apt install libeigen3-dev`) gives the following:

```shell
         0          0          1          0          0
 0.0833333  -0.666667          0   0.666667 -0.0833333
-0.0833333    1.33333       -2.5    1.33333 -0.0833333
      -0.5          1          0         -1        0.5
         1         -4          6         -4          1
```

The third row gives the same coefficients as found with the method of undetermined coefficient using the CAS.
