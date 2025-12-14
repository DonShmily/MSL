#include <iostream>
#include <vector>

#include "polynomial.hpp"
#include "polynomial/polynomial.hpp"

using namespace msl;

int test_polynomial()
{
    std::vector<double> x{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    std::vector<double> y{1.0, 2.0, 0.0, 5.0, 4.0, 3.0};

    polynomial::Polynomial poly1(x, y, 3);
    polynomial::Polynomial poly2(x, y, 4);

    std::vector<double> test_x{0.5, 1.5, 2.5, 3.5};
    try
    {
        std::vector<double> results1 = poly1(test_x);
        std::cout << "Polynomial 1 results:" << std::endl;
        std::cout << "Coefficients: ";
        for (const auto &c : poly1.coefficients())
        {
            std::cout << c << " ";
        }
        std::cout << std::endl;
        std::cout << "Test X: ";
        for (const auto &val : test_x)
        {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        std::cout << "Results: ";
        for (const auto &val : results1)
        {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::vector<double> results2 = poly2(test_x);
        std::cout << "Polynomial 2 results:" << std::endl;
        std::cout << "Coefficients: ";
        for (const auto &c : poly2.coefficients())
        {
            std::cout << c << " ";
        }
        std::cout << std::endl;
        std::cout << "Test X: ";
        for (const auto &val : test_x)
        {
            std::cout << val << " ";
        }
        std::cout << std::endl;
        std::cout << "Results: ";
        for (const auto &val : results2)
        {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Polynomial test failed: " << e.what() << std::endl;
        return -1;
    }
    return 0;
}

int main() { return test_polynomial(); }
