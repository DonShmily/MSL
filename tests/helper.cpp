#include "helper.hpp"

void print_coefficients(const FilterCoefficients &coeffs,
                        const std::string &name)
{
    std::cout << name << " Coefficients:\n";
    std::cout << "Order: " << coeffs.order() << "\n\n";

    std::cout << "Numerator (b):\n  [";
    for (size_t i = 0; i < coeffs.b.size(); ++i)
    {
        std::cout << std::setw(12) << std::setprecision(6) << coeffs.b[i];
        if (i < coeffs.b.size() - 1)
            std::cout << ",";
    }
    std::cout << "]\n\n";

    std::cout << "Denominator (a):\n  [";
    for (size_t i = 0; i < coeffs.a.size(); ++i)
    {
        std::cout << std::setw(12) << std::setprecision(6) << coeffs.a[i];
        if (i < coeffs.a.size() - 1)
            std::cout << ",";
    }
    std::cout << "]\n\n";
}

void print_vector(const std::vector<double> &v,
                  const std::string &name,
                  size_t max_print)
{
    std::cout << name << " (size=" << v.size() << "): [";
    size_t n = std::min(v.size(), max_print);
    for (size_t i = 0; i < n; ++i)
    {
        std::cout << std::setw(8) << std::setprecision(4) << v[i];
        if (i < n - 1)
            std::cout << ", ";
    }
    if (v.size() > max_print)
        std::cout << ", ...";
    std::cout << "]\n";
}

void print_complex_vector(const std::vector<std::complex<double>> &v,
                          const std::string &name,
                          size_t max_print)
{
    std::cout << name << " (size=" << v.size() << "):\n";
    size_t n = std::min(v.size(), max_print);
    for (size_t i = 0; i < n; ++i)
    {
        std::cout << "  [" << i << "] = " << std::setw(10)
                  << std::setprecision(4) << v[i].real() << " + "
                  << std::setw(10) << std::setprecision(4) << v[i].imag()
                  << "i\n";
    }
    if (v.size() > max_print)
        std::cout << "  ...\n";
}

// Read CSV file (simple, single column)
std::vector<double> read_csv(const std::string &filename)
{
    std::vector<double> data;
    std::ifstream file(filename);

    if (!file.is_open())
    {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::string line;
    while (std::getline(file, line))
    {
        if (!line.empty())
        {
            data.push_back(std::stod(line));
        }
    }

    return data;
}

// Compute maximum absolute error
double max_abs_error(const std::vector<double> &a, const std::vector<double> &b)
{
    if (a.size() != b.size())
    {
        throw std::runtime_error("Size mismatch in error calculation");
    }

    double max_err = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
    {
        max_err = std::max(max_err, std::abs(a[i] - b[i]));
    }
    return max_err;
}

// Compute maximum absolute error for complex vectors
double max_abs_error_complex(const std::vector<std::complex<double>> &a,
                             const std::vector<std::complex<double>> &b)
{
    if (a.size() != b.size())
    {
        throw std::runtime_error("Size mismatch in error calculation");
    }

    double max_err = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
    {
        max_err = std::max(max_err, std::abs(a[i] - b[i]));
    }
    return max_err;
}

// Compute RMS error
double rms_error(const std::vector<double> &a, const std::vector<double> &b)
{
    if (a.size() != b.size())
    {
        throw std::runtime_error("Size mismatch in RMS calculation");
    }

    double sum_sq = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
    {
        double diff = a[i] - b[i];
        sum_sq += diff * diff;
    }
    return std::sqrt(sum_sq / a.size());
}

// Print test result
void print_result(const std::string &test_name, double error, double threshold)
{
    std::cout << std::setw(50) << std::left << test_name;
    std::cout << " Error: " << std::scientific << std::setprecision(2) << error;

    if (error < threshold)
    {
        std::cout << " ✓ PASS\n";
    }
    else
    {
        std::cout << " ✗ FAIL (threshold: " << threshold << ")\n";
    }
}