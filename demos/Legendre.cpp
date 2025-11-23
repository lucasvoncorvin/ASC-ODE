#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <autodiff.hpp>


using namespace ASC_ode;

template <typename T>
void LegendrePolynomials(int n, T x, std::vector<T>& P) {
    if (n < 0) {
        P.clear();
        return;
    }
    P.resize(n + 1);
    P[0] = T(1);
    if (n == 0) return;
    P[1] = x;
    for (int k = 2; k <= n; ++k) {
        P[k] = ((T(2 * k - 1) * x * P[k - 1]) - T(k - 1) * P[k - 2]) / T(k);
    }
}

int main(){
    int max_order=5; //MaxOrd Legendre Polynomials
    int steps=100; //Personal choice

    //Files for values and derivatives, respectively so we can plot easier (I think)
    std::vector<std::ofstream> value_files(max_order + 1); //+1 because we also consider the 0 order Legendre polynomial
    std::vector<std::ofstream> deriv_files(max_order + 1);

    for (int order = 0; order <= max_order; ++order) { //Just opening files
        std::string value_filename = "legendre_P" + std::to_string(order) + "_values.txt";
        std::string deriv_filename = "legendre_P" + std::to_string(order) + "_derivatives.txt";
        
        value_files[order].open(value_filename);
        deriv_files[order].open(deriv_filename);
    }

    for (int i = 0; i <= steps; i++) {
        double x_val= -1.0 + (2.0 * i) / steps;  //x in [-1, 1]
        AutoDiff<1> x= Variable<0>(x_val);
        std::vector<AutoDiff<1>> P; //Initiating Empty Vector of polynomials
        LegendrePolynomials(max_order, x, P); //Calling it and changing it inside by reference
        //Writing the evaluations in the respective files
        for (int order = 0; order <= max_order; ++order) {
            value_files[order] << " " << x_val << " " << P[order].value() << "\n";
            deriv_files[order] << " " << x_val << " " << P[order].deriv()[0] << "\n";
        }
    }

    //Closing fles
    for (int order = 0; order <= max_order; order++) {
        value_files[order].close();
        deriv_files[order].close();
    }

    return 0;
}