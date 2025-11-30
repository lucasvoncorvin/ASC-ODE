#include <iostream>
#include <autodiff.hpp>
#include <nonlinfunc.hpp>
#include <cmath>

using namespace ASC_ode;


template <typename T>
T func1 (T x, T y)
{
  return x * sin(y);
  // return 1e6 + y;
}



int main()
{
  double x = 1, y = 2;
  AutoDiff<2> adx = Variable<0>(x);
  AutoDiff<2> ady = Variable<1>(y);

  std::cout << "adx = " << adx << std::endl;
  std::cout << "ady = " << ady << std::endl;

  AutoDiff<2> prod = adx * ady;
  std::cout << "prod = " << prod << std::endl;

  std::cout << "func1(adx, ady) = " << func1(adx, ady) << std::endl;

  double eps = 1e-8;
  std::cout << "numdiff df/dx = " << (func1(x + eps, y) - func1(x-eps, y)) / (2*eps) << std::endl;
  std::cout << "numdiff df/dy = " << (func1(x, y + eps) - func1(x, y-eps)) / (2*eps) << std::endl;


  {
    // we can do second derivatives:
    AutoDiff<1, AutoDiff<1>> addx{Variable<0>(2)};
    std::cout << "addx = " << addx << std::endl;
    // func = x*x
    // func' = 2*x
    // func'' = 2
    std::cout << "addx*addx = " << addx * addx << std::endl;

    // std::cout << "sin(addx) = " << sin(addx) << std::endl;
  }

  //Test for the pendulum
  PendulumAD pendulum(1.0); //Our pendulum via our class PendulumAD

  std::vector<double> x={M_PI/4, 0.1}; //Our test (x,x') (ie. (pi/4) rad with angular speed 0.1)
    std::vector<double> f(2); //Where we are going to store our evaluation
    std::vector<std::vector<double>> df(2, std::vector<double>(2)); //for storing the Jacobian Matrix
  
    pendulum.evaluate(x, f); //testing the evaluation of our defined f for the ODE system
    std::cout << "x = [" << x[0] << ", " << x[1] << "]\n";
    std::cout << "f(x) = [" << f[0] << ", " << f[1] << "]\n";
    
    pendulum.evaluateDeriv(x, df); //Jacobian tsting
    std::cout << "\nJacobian f'(x):\n";
    std::cout << "[" << df[0][0] << ", " << df[0][1] << "]\n";
    std::cout << "[" << df[1][0] << ", " << df[1][1] << "]\n";
  
  
    return 0;
}