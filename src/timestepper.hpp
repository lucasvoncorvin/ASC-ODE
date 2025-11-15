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
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<ConstantFunction> m_yold;
    std::shared_ptr<ConstantFunction> m_fold;
    std::shared_ptr<Parameter> m_tau_half;

    Vector<> m_tmp_f;

public:
    CrankNicolson(std::shared_ptr<NonlinearFunction> rhs)
        : TimeStepper(rhs),
          m_tmp_f(rhs->dimF())
    {
        m_yold = std::make_shared<ConstantFunction>(Vector<>(rhs->dimX()));
        m_fold = std::make_shared<ConstantFunction>(Vector<>(rhs->dimF()));
        m_tau_half = std::make_shared<Parameter>(0.0);
        auto Id = std::make_shared<IdentityFunction>(rhs->dimX());
        auto term1 = std::make_shared<SumFunction>(
            Id, m_yold, 1.0, -1.0
        );
        auto term2 = (-1.0) * (m_tau_half * this->m_rhs);
        m_equ = std::make_shared<SumFunction>(
            term1, term2, 1.0, 1.0
        );
        m_equ = std::make_shared<SumFunction>(
            m_equ, m_fold, 1.0, 1.0
        );
    }

    void DoStep(double tau, VectorView<> y) override
    {
        double half = 0.5 * tau;
        m_tau_half->set(half);
        m_yold->set(y);
        this->m_rhs->evaluate(y, m_tmp_f);
        Vector<> tmp = -half * m_tmp_f;
        m_fold->set(tmp);
        NewtonSolver(m_equ, y, 1e-10, 30, nullptr);
    }
};


}

#endif
