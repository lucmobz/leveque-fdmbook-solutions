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

## A

The *GiNaC* C++ library (https://www.ginac.de/) implements a computer algebra system that can solve the linear system symbolically and give exact results in terms of rational numbers. On *WSL*, using *Ubuntu 22.04LTS*, it is available with `sudo apt install libginac-dev` (which also installs the dependency `libcln-dev` which *GiNaC* uses).

The program `exercise_1_2a.cpp`, compiled with `g++ --std=c++23 -g -Wall -Wextra -Wpedantic exercise_1_2a.cpp -lginac -lcln && ./a.out`, finds the following coefficients: 
```shell
[[-1/12],[4/3],[-5/2],[4/3],[-1/12]]
```

## B
