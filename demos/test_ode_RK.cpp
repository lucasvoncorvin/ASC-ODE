#include <iostream>
#include <fstream> 
#include <numbers>
#include <nonlinfunc.hpp>
#include <timestepper.hpp>
#include <implicitRK.hpp>
#include <ExplicitRK.hpp>

using namespace ASC_ode;


class MassSpring : public NonlinearFunction
{
private:
  double mass;
  double stiffness;

public:
  MassSpring(double m, double k) : mass(m), stiffness(k) {}

  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }
  
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
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


int main()
{
  double tend = 4 * std::numbers::pi;
  int steps = 10;
  double tau = tend/steps;

  auto rhs = std::make_shared<MassSpring>(1.0, 1.0);

  std::ofstream file_RK4("output_RK4.txt");
  std::ofstream file_Gauss("output_Gauss4.txt");
  std::ofstream file_Radau("output_Radau4.txt");


// 1) EXPLICIT RK4
    {
        Matrix<> RK4a {
            {0,   0,   0,   0},
            {1.0/2, 0,   0,   0},
            {0,   1.0/2, 0,   0},
            {0,   0,   1.0, 0}
        };

        Vector<> RK4b { 1.0/6, 1.0/3, 1.0/3, 1.0/6 };
        Vector<> RK4c { 0.0, 0.5, 0.5, 1.0 };

        ExplicitRungeKutta stepper(rhs, RK4a, RK4b, RK4c);

        Vector<> y = {1, 0};   // initial condition
        file_RK4 << 0 << "  " << y(0) << " " << y(1) << "\n";

        for (int i = 0; i < steps; i++)
        {
            stepper.DoStep(tau, y);
            file_RK4 << (i+1)*tau << "  " << y(0) << " " << y(1) << "\n";
        }
    }
  


  // 2) GAUSS–LEGENDRE 4-stage implicit RK (order 8)

    {
        int s = 4;
        Vector<> c(s), b(s);
        GaussLegendre(c, b);

        auto [a, b2] = ComputeABfromC(c);
        ImplicitRungeKutta stepper(rhs, a, b2, c);

        Vector<> y = {1, 0};   // reset initial condition
        file_Gauss << 0 << "  " << y(0) << " " << y(1) << "\n";

        for (int i = 0; i < steps; i++)
        {
            stepper.DoStep(tau, y);
            file_Gauss << (i+1)*tau << "  " << y(0) << " " << y(1) << "\n";
        }
    }
  
    // 3) RADAU IIA 4-stage implicit RK (order 7)
   
    {
        int s = 4;
        Vector<> c(s), b(s);
        GaussRadau(c, b);

        auto [a, b2] = ComputeABfromC(c);
        ImplicitRungeKutta stepper(rhs, a, b2, c);

        Vector<> y = {1, 0};   // reset initial condition
        file_Radau << 0 << "  " << y(0) << " " << y(1) << "\n";

        for (int i = 0; i < steps; i++)
        {
            stepper.DoStep(tau, y);
            file_Radau << (i+1)*tau << "  " << y(0) << " " << y(1) << "\n";
        }
    }


 std::cout << "All methods completed. Files generated:\n";
    std::cout << "  - output_RK4.txt\n";
    std::cout << "  - output_Gauss4.txt\n";
    std::cout << "  - output_Radau4.txt\n";

    return 0;
}
