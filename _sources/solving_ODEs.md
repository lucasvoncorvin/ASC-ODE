# Solving ODEs

We solve ODEs of the form
$$
    y' = f(x, y)
$$
Higher order ODEs must be converted into a system of first order ODEs. For example an undamped mass spring system is described by the second order ODE
$$
\begin{aligned}
    m y''(t) + k y(t) &= 0 \\
    y''(t) &= - \frac{k}{m} y(t) 
\end{aligned}
$$

which can be converted to a system of first order ODEs by setting
$$
\begin{aligned}
    y &= y_0 \\
    y' = y_0' &= y_1 \\
    y'' = y_1' &= - \frac{k}{m} y_0
\end{aligned}
$$

which can be written in matrix form as

$$
\begin{pmatrix}
    y_0' \\
    y_1'
\end{pmatrix}
=
\begin{pmatrix}
    0 & 1 \\
    -\frac{k}{m} & 0
\end{pmatrix}
\begin{pmatrix}
    y_0 \\
    y_1
\end{pmatrix}
$$

## Right-Hand Side Function

The right hand side $f(t, x)$ is represented as a derived class from the abstract NonlinearFunction base class.
We have to implement the following member functions
```cpp
    class NonlinearFunction
    {
    public:
        virtual ~NonlinearFunction() = default; // optional destructor
        virtual size_t dimX() const = 0; // returns the dimension of the input
        virtual size_t dimF() const = 0; // returns the dimension of the output
        // evaluates the function at input x and writes the result to output f
        virtual void evaluate (VectorView<double> x, VectorView<double> f) const = 0; 
        // evaluate the derivative at input x and writes the result to output f
        virtual void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const = 0; 
    };
```

For the mass spring example we implement 

```cpp
    class MassSpring : public NonlinearFunction
    {
    private:
    // the function has two parameter, which are assigned in the constructor
    double mass;
    double stiffness;

    public:
    MassSpring(double m, double k) : mass(m), stiffness(k) {}

    size_t dimX() const override { return 2; } // inputs are 2d vectors [y_0, y_1] = [y, y']
    size_t dimF() const override { return 2; } // outputs are 2d vector [f_0, f_1]
    
    void evaluate (VectorView<double> x, VectorView<double> f) const override
    {   
        // Here the function is defined as seen above
        f(0) = x(1);
        f(1) = -stiffness/mass*x(0);
    }
    
    void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
        df = 0.0;
        df(0,1) = 1;
        df(1,0) = -stiffness/mass;
    }
    };
```

## The TimerStepper class

With the right hand side defined as a NonlinearFunction, we can now solve the ODE system using any of the provided implementations of the Timestepper class
The provided implementations are:
- Explicit Euler
- Implicit Euler
- Explicit Runge-Kutta
  - Runge-Kutta 2
- Implicit Runge-Kutta
- Crank-Nicolson

Depending on the problem a different method may be suitable. The user should choose the method that best fits their requirements.

The timestepper is initialized with (at least) a shared pointer to the right hand side function. It implements the doStep method, which does one solution step. It should hence be used in a loop until reaching the desired time.

```cpp
class MyStepper : public TimeStepper {

    ...
public:
    // This function does one time step of size tau and writes it to the y parameter
    // @param
    //  - tau   size of the timestep
    //  - y     solution at current time
    void doStep(double tau, VectorView<double> y) override
    ...
};
```

## Example

We provide an example using the ExplicitEuler implementation

```cpp
    // define the timestep
    double t = 10;          // we seek the solution a time 10
    int steps = 100;        // number of steps the solver should take
    double tau = t/steps;   // size of the timestep

    // define initial condition
    Vector<> y = { 1, 0 };

    // define the right hand side function (Note that the MassSpring function implementation is shown above)
    shared_ptr<NonlinearFunction> rhs = std::make_shared<MassSpring>(1, 1);

    // define the solution method (here ExplicitEuler)
    ExplicitEuler stepper(rhs);

    // solve
    for (int i = 0; i < steps; i++)
    {
        stepper.doStep(tau, y);
    }

    // the solution is now stored in the vector y
    std::cout << "Solution: y[0] = " << y(0) << std::endl;
```
