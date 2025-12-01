
This part of the report summarizes the implementation and testing of an automatic differentiation (AutoDiff) class applied to two problems: solving the pendulum ODE and evaluating Legendre polynomials with their derivatives. Automatic differentiation provides exact derivatives without numerical approximation errors.

The AutoDiff class template stores both values and derivatives, with comprehensive operator overloading for basic arithmetic:

-Addition, subtraction, multiplication, and division operators
-Support for AutoDiff-AutoDiff and constant-AutoDiff combinations
-Proper implementation of derivative rules (product rule, quotient rule)

Moreover we extended for the eelementary functions: sin, cos, exp, log.

In the other hand The pendulum equation $\alpha'' = -\frac{g}{l}\sin(\alpha)$ was reformulated as a first-order system:
$y_0' &= y_1$
$y_1' &= -\frac{g}{l}\sin(y_0)$

Testing with initial condition $[\pi/4, 0.1]$ gave:

Function: $f(x) = [0.1, -0.567669]$
Jacobian: $\begin{bmatrix} 0 & 1 \\ -0.800633 & 0 \end{bmatrix}$

The Jacobian matches the analytical solution, confirming correct derivative computation.

Legendre polynomials were evaluated using their recursive definition as it was defined.

The AutoDiff implementation automatically computes both polynomial values and their derivatives. Polynomials up to order 5 were evaluated over $[-1, 1]$, with derivatives computed accurately through the automatic differentiation mechanism.

```{image} ../demos/legendre_plot.png