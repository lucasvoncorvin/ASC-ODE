# General Runge Kutta

There are two timestepper available for general Runge-Kutta methods. Since they are not trivial to use we provide a short explanation. Besides the two general Runge Kutta solvers, there is also a explicit RK2 solver, which can be used the same as most other provided timesteppers.

The special timesteppers in question are ImplicitRungeKutta and ExplicitRungeKutta. Both take the same input parameters.

```cpp
ImplicitRungeKutta( // ExplicitRungeKutta takes the same parameters
    std::shared_ptr<NonlinearFunction> rhs, 
    const Matrix<> &a, 
    const Vector<> &b, 
    const Vector<> &c)
```

The parameter a, b, and c represent the Butcher tableau.

$$
\begin{array}{c|ccc}
c_0 & a_{0,0} & \cdots & a_{0,s-1} \\\\
\vdots & \vdots & \ddots & \vdots \\\\
c_{s-1} & a_{s-1,0} & \cdots & a_{s-1,s-1} \\\\
\hline
& b_0 & \cdots & b_{s-1}
\end{array}
$$

It completely defines any Runge Kutta method via
$$
    k_j = f(y_i + \tau \sum_{l=0}^{s-1}{a_{jl}k_l})
$$ 
and
$$
    y_{i+1} = y_i + \tau \sum_{l=0}^{s-1}{b_l k_l}
$$


## Explicit Runge Kutta

The general runge kutta method requires the solution of a system of equation to get the factors $k_l$. This is what the class ImplicitRungeKutta does.
However, the class ExplicitRungeKutta, does not solve a system of equations.
Instead it explicitly computes the factors $k_l$ by supposing that $a_{jl} = 0$ for $l \ge j$

Therefore when using this class the **matrix must be strictly lower triangular**. Otherwise it will not yield the expected results
