#include <iostream>
#include <fstream>
#include <cmath>
#include <numbers>

#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;

// autonomous ODE for RC Circuit
// y(0) = UC 
// y(1) = t  
// f0 = (cos(100*pi*t) - UC) / (R*C)
// f1 = 1

class RCCircuit : public NonlinearFunction
{
private:
    double R, C;

public:
    RCCircuit(double R_, double C_) : R(R_), C(C_) {}

    size_t dimX() const override { return 2; }
    size_t dimF() const override { return 2; }

    // f(y)
    void evaluate(VectorView<double> x, VectorView<double> f) const override
    {
        double UC = x(0);
        double t  = x(1);

        double U0 = std::cos(100.0 * std::numbers::pi * t);

        f(0) = (U0 - UC) / (R * C);  // dUC/dt
        f(1) = 1.0;                  // dt/dt = 1
    }

    // Jacobian df/dy
    void evaluateDeriv(VectorView<double> x, MatrixView<double> df) const override
    {
        df = 0.0;

        double t = x(1);

        df(0,0) = -1.0 / (R*C);
        df(0,1) = (-100.0 * std::numbers::pi * std::sin(100.0 * std::numbers::pi * t)) / (R*C);

        df(1,1) = 0.0;
    }
};


// MAIN

int main()
{
    double R = 100.0;
    double C = 1e-6;

    double t_end = 1;
    int steps = 100;
    double tau = t_end / steps;

    // RHS
    auto rhs = std::make_shared<RCCircuit>(R, C);

    // Initial condition
    Vector<> y_exp = { 0.0, 0.0 };
    Vector<> y_imp = y_exp;
    Vector<> y_cn  = y_exp;

    // Time steppers
    ExplicitEuler step_exp(rhs);
    ImplicitEuler step_imp(rhs);
    CrankNicolson step_cn(rhs);

    std::ofstream fE("rc_explicit.txt");
    std::ofstream fI("rc_implicit.txt");
    std::ofstream fC("rc_crank.txt");

    fE << y_exp(1) << " " << y_exp(0) << "\n";
    fI << y_imp(1) << " " << y_imp(0) << "\n";
    fC << y_cn(1)  << " " << y_cn(0)  << "\n";

    // Time loop
    for (int i = 0; i < steps; i++)
    {
        step_exp.DoStep(tau, y_exp);
        step_imp.DoStep(tau, y_imp);
        step_cn.DoStep(tau, y_cn);

        fE << y_exp(1) << " " << y_exp(0) << "\n";
        fI << y_imp(1) << " " << y_imp(0) << "\n";
        fC << y_cn(1)  << " " << y_cn(0)  << "\n";
    }

    std::cout << "Simulation finished. Files written:\n"
              << "  rc_explicit.txt\n"
              << "  rc_implicit.txt\n"
              << "  rc_crank.txt\n";

    return 0;
}
