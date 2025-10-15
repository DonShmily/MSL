// tests/compare_with_matlab.cpp
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "signal/filter_design.hpp"

namespace fs = std::filesystem;
using namespace msl::signal;

static void write_vector_csv(const fs::path &p, const std::vector<double> &v)
{
    std::ofstream ofs(p);
    ofs << std::setprecision(15);
    for (size_t i = 0; i < v.size(); ++i)
    {
        if (i)
            ofs << ",";
        ofs << v[i];
    }
    ofs << "\n";
}

int test_filter_design()
{
    fs::path outdir = "msl_out";
    if (!fs::exists(outdir))
        fs::create_directory(outdir);

    struct Case
    {
        std::string tag;
        std::string type;
        int N;
        std::vector<double> Wn;
    };
    std::vector<Case> cases = {{"case01", "low", 2, {0.10}},
                               {"case02", "low", 4, {0.20}},
                               {"case03", "low", 6, {0.40}},
                               {"case04", "low", 4, {0.01}},
                               {"case05", "high", 2, {0.10}},
                               {"case06", "high", 4, {0.40}},
                               {"case07", "band", 2, {0.10, 0.15}},
                               {"case08", "band", 3, {0.05, 0.40}},
                               {"case09", "stop", 4, {0.20, 0.25}},
                               {"case10", "low", 8, {0.48}}};

    std::cout << "MSL coefficient generator — output directory: " << outdir
              << "\n\n";

    for (const auto &c : cases)
    {
        std::cout << "Generating " << c.tag << " (" << c.type << ", N=" << c.N
                  << ") ...\n";

        FilterCoefficients coeffs;
        try
        {
            if (c.type == "low")
            {
                ButterworthFilter f(c.N, c.Wn[0], FilterType::lowpass);
                coeffs = f.coefficients();
            }
            else if (c.type == "high")
            {
                ButterworthFilter f(c.N, c.Wn[0], FilterType::highpass);
                coeffs = f.coefficients();
            }
            else if (c.type == "band")
            {
                ButterworthFilter f(
                    c.N, c.Wn[0], c.Wn[1], FilterType::bandpass);
                coeffs = f.coefficients();
            }
            else if (c.type == "stop")
            {
                ButterworthFilter f(
                    c.N, c.Wn[0], c.Wn[1], FilterType::bandstop);
                coeffs = f.coefficients();
            }
            else
            {
                std::cerr << "Unknown case type: " << c.type << "\n";
                continue;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception designing filter: " << e.what() << "\n";
            continue;
        }

        // Normalize if a[0] != 1
        if (!coeffs.a.empty())
        {
            double a0 = coeffs.a[0];
            if (std::abs(a0 - 1.0) > 1e-15 && std::abs(a0) > 0.0)
            {
                for (auto &x : coeffs.b)
                    x /= a0;
                for (auto &x : coeffs.a)
                    x /= a0;
            }
        }

        // write files
        auto bfile = outdir / (c.tag + "_b.csv");
        auto afile = outdir / (c.tag + "_a.csv");
        write_vector_csv(bfile, coeffs.b);
        write_vector_csv(afile, coeffs.a);

        // print to stdout for quick look
        std::cout << "  b (" << coeffs.b.size() << "): ";
        for (auto v : coeffs.b)
            std::cout << v << " ";
        std::cout << "\n  a (" << coeffs.a.size() << "): ";
        for (auto v : coeffs.a)
            std::cout << v << " ";
        std::cout << "\n  files: " << bfile << " , " << afile << "\n\n";
    }

    std::cout << "All cases written into " << outdir << "\n";
    return 0;
}
