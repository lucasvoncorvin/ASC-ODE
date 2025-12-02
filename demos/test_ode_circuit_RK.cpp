#include <iostream>
#include <fstream>
#include <cmath>
#include <numbers>
#include <memory>
#include <nonlinfunc.hpp>
#include <timestepper.hpp>
#include <implicitRK.hpp>
#include <ExplicitRK.hpp>

using namespace ASC_ode;

// 
// RC CIRCUIT MODEL
// 

class RCCircuit : public NonlinearFunction
{
private:
    double R, C;

public:
    RCCircuit(double R_, double C_) : R(R_), C(C_) {}

    size_t dimX() const override { return 2; }
    size_t dimF() const override { return 2; }

    void evaluate(VectorView<double> x, VectorView<double> f) const override
    {
        double UC = x(0);
        double t  = x(1);

        double U0 = std::cos(100.0 * std::numbers::pi * t);

        f(0) = (U0 - UC) / (R*C);    // dUC/dt
        f(1) = 1.0;                  // dt/dt
    }

    void evaluateDeriv(VectorView<double> x, MatrixView<double> df) const override
    {
        df = 0.0;

        double t = x(1);

        df(0,0) = -1.0 / (R*C);
        df(0,1) = (-100.0 * std::numbers::pi *
                   std::sin(100.0 * std::numbers::pi * t)) / (R*C);
    }
};


// 
// MAIN PROGRAM
// 

int main()
{
    double R = 100.0;
    double C = 1.0e-6;

    double t_end = 1.0;
    int steps = 1000;
    double tau = t_end / steps;

    auto rhs = std::make_shared<RCCircuit>(R, C);

    // Output files
    std::ofstream file_RK4("output_RK4.txt");
    std::ofstream file_Gauss("output_Gauss4.txt");
    std::ofstream file_Radau("output_Radau4.txt");

    // 
    // 1) EXPLICIT RK4
    // 

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

        Vector<> y = {0.0, 0.0};

        file_RK4 << 0 << "  " << y(0) << " " << y(1) << "\n";

        for (int i = 0; i < steps; i++)
        {
            stepper.DoStep(tau, y);
            file_RK4 << (i+1)*tau << "  " << y(0) << " " << y(1) << "\n";
        }
    }

    // 
    // 2) GAUSS–LEGENDRE (4 stages, order 8)
    // 

    {
        int s = 4;
        Vector<> c(s), b(s);
        GaussLegendre(c, b);

        auto [a, b2] = ComputeABfromC(c);
        ImplicitRungeKutta stepper(rhs, a, b2, c);

        Vector<> y = {0.0, 0.0};

        file_Gauss << 0 << "  " << y(0) << " " << y(1) << "\n";

        for (int i = 0; i < steps; i++)
        {
            stepper.DoStep(tau, y);
            file_Gauss << (i+1)*tau << "  " << y(0) << " " << y(1) << "\n";
        }
    }

    // 
    // 3) RADAU IIA (4 stages, order 7)
    // 

    {
        int s = 4;
        Vector<> c(s), b(s);
        GaussRadau(c, b);

        auto [a, b2] = ComputeABfromC(c);
        ImplicitRungeKutta stepper(rhs, a, b2, c);

        Vector<> y = {0.0, 0.0};

        file_Radau << 0 << "  " << y(0) << " " << y(1) << "\n";

        for (int i = 0; i < steps; i++)
        {
            stepper.DoStep(tau, y);
            file_Radau << (i+1)*tau << "  " << y(0) << " " << y(1) << "\n";
        }
    }

    std::cout << "All methods completed. Files written.\n";

    return 0;
}
