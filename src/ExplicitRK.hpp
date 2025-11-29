#ifndef EXPLICITRK_HPP
#define EXPLICITRK_HPP

#include <vector.hpp>
#include <matrix.hpp>
#include <inverse.hpp>
#include <nonlinfunc.hpp>
#include <timestepper.hpp>
#include <memory>

namespace ASC_ode {
    using namespace nanoblas;

    class ExplicitRungeKutta : public TimeStepper
    {
        Matrix<> m_a;
        Vector<> m_b, m_c;
        int m_stages;
    public:
        ExplicitRungeKutta(std::shared_ptr<NonlinearFunction> rhs, 
            const Matrix<> &a, const Vector<> &b, const Vector<> &c) : 
            TimeStepper(rhs), m_a(a), m_b(b), m_c(c),
            m_stages(c.size())
        {
        }

        void DoStep(double tau, VectorView<double> y) override
        {
            // Implementation of the explicit Runge-Kutta step

            // compute k explicitly
            Matrix<> ks(m_stages, this->m_rhs->dimF());

            // store k values, initialize to 0
            Vector<> ksum(this->m_rhs->dimF());
            ksum = 0.0;

            // compute all k values
            for (size_t l = 0; l < m_stages; ++l) {
                Vector<> y_temp = y;

                // compute y_i + tau * sum_j a(l,j)*k_j (input for f to get k_l)
                for (size_t j = 0; j < l; ++j) { // only previous ks are used since all a(l,j) with j>=l are zero (explicit RK)  
                    y_temp = y_temp + (tau * m_a(l, j) * ks.row(j));
                }

                // evaluate the right-hand side function to get k_l
                this->m_rhs->evaluate(y_temp, ks.row(l)); // k_l = f(y_i + tau * sum_j a(l,j)*k_j)

                // accumulate ksum for the final update
                ksum = ksum + m_b(l) * ks.row(l); // accumulate weighted sum of ks
            }

            // update solution y_{i+1} = y_i + tau * sum_l b(l)*k_l
            y = y + tau * ksum; 
          
        }
    };

} // namespace ASC_ode 

#endif // EXPLICITRK_HPP