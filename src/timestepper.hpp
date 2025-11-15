#ifndef TIMERSTEPPER_HPP
#define TIMERSTEPPER_HPP

#include <functional>
#include <exception>
#include <memory>
#include "Newton.hpp"


namespace ASC_ode
{
  
  class TimeStepper
  { 
  protected:
    std::shared_ptr<NonlinearFunction> m_rhs;
  public:
    TimeStepper(std::shared_ptr<NonlinearFunction> rhs) : m_rhs(rhs) {}
    virtual ~TimeStepper() = default;
    virtual void DoStep(double tau, VectorView<double> y) = 0;
  };

  class ExplicitEuler : public TimeStepper
  {
    Vector<> m_vecf;
  public:
    ExplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_vecf(rhs->dimF()) {}
    void DoStep(double tau, VectorView<double> y) override
    {
      this->m_rhs->evaluate(y, m_vecf);
      y += tau * m_vecf;
    }
  };

class RK2 : public TimeStepper
{
private:
    Vector<> k1, k2;

public:
    RK2(std::shared_ptr<NonlinearFunction> rhs)
        : TimeStepper(rhs),
          k1(rhs->dimF()),
          k2(rhs->dimF())
    {}

    void DoStep(double tau, VectorView<> y) override
    {
        this->m_rhs->evaluate(y, k1);
        Vector<> ytilde(y.size());
        ytilde = y + 0.5 * tau * k1;
        this->m_rhs->evaluate(ytilde, k2);
        y += tau * k2;
    }
};

 
class ImplicitEuler : public TimeStepper
{
private:
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<ConstantFunction> m_yold;
    std::shared_ptr<Parameter> m_tau;

public:
    ImplicitEuler(std::shared_ptr<NonlinearFunction> rhs)
        : TimeStepper(rhs)
    {
        m_yold = std::make_shared<ConstantFunction>(
                    Vector<>(rhs->dimX()) );
        m_tau = std::make_shared<Parameter>(0.0);
        auto Id = std::make_shared<IdentityFunction>(rhs->dimX());
        auto term1 = std::make_shared<SumFunction>(
            Id, m_yold, 1.0, -1.0
        );

        auto scaled_rhs = (-1.0) * (m_tau * rhs);
        m_equ = std::make_shared<SumFunction>(
            term1, scaled_rhs, 1.0, 1.0
        );
    }

    void DoStep(double tau, VectorView<> y) override
    {
        m_yold->set(y);
        m_tau->set(tau);
        NewtonSolver(m_equ, y, 1e-10, 30, nullptr);
    }
};

class CrankNicolson : public TimeStepper
{
private:
    std::shared_ptr<NonlinearFunction> G;
    std::shared_ptr<ConstantFunction> y_old;
    std::shared_ptr<ConstantFunction> f_old;
    std::shared_ptr<Parameter> tau_param;

public:
    CrankNicolson(std::shared_ptr<NonlinearFunction> rhs)
        : TimeStepper(rhs)
    {
        size_t n = rhs->dimX();

        y_old  = std::make_shared<ConstantFunction>(Vector<>(n));
        f_old  = std::make_shared<ConstantFunction>(Vector<>(n));
        tau_param = std::make_shared<Parameter>(0.0);

        auto Id = std::make_shared<IdentityFunction>(n);

        // (y - y_old)
        auto term1 = std::make_shared<SumFunction>(Id, y_old,
                                                   1.0, -1.0);

        // -(tau/2)*f(y)
        auto half_tau = std::make_shared<Parameter>(0.5);
        auto scaled_rhs = std::make_shared<ScaleFunction>(rhs, tau_param);
        auto minus_half_scaled_rhs =
            std::make_shared<ScaleFunction>(scaled_rhs, std::make_shared<Parameter>(-0.5));

        // -(tau/2)*f_old
        auto scaled_f_old = std::make_shared<ScaleFunction>(f_old, tau_param);
        auto minus_half_f_old =
            std::make_shared<ScaleFunction>(scaled_f_old, std::make_shared<Parameter>(-0.5));

        // G = (y - y_old) - (tau/2)(f(y) + f_old)
        G = std::make_shared<SumFunction>(
            term1,
            std::make_shared<SumFunction>(minus_half_scaled_rhs, minus_half_f_old,
                                          1.0, 1.0),
            1.0, 1.0
        );
    }

    void DoStep(double tau, VectorView<> y) override
    {
     
        y_old->set(y);
        Vector<> ft(y.size());
        m_rhs->evaluate(y, ft);
        f_old->set(ft);

        tau_param->set(tau);

        double t_old = y(1);

        NewtonSolver(G, y, 1e-10, 1000, nullptr);

    
        y(1) = t_old + tau;
    }
};



}

#endif
