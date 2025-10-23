#include <vector>

#include "integral.hpp"
#include "integral/cumtrapz.hpp"
#include "matrix.hpp"
#include "utils/data_io.hpp"

using namespace msl;
int test_integral()
{
    auto ori_data = utils::ReadData("KunmingSSJY.txt", 6, 3e4);
    auto data_1d =
        std::vector<double>(ori_data.begin(), ori_data.begin() + 3e4);
    matrixd ori_matrix(3e4, 6, std::span<const double>(ori_data));
    try
    {
        // test 1d integral
        auto integral_1d = integral::cumtrapz(data_1d, 1.0);
        utils::WriteData(
            "test_result/integral/integral_1d.txt", integral_1d, 1, 3e4);

        // test 1d integral with non-uniform x
        std::vector<double> x_vals(3e4);
        for (size_t i = 0; i < x_vals.size(); ++i)
        {
            x_vals[i] = static_cast<double>(i) * 0.1;
        }
        auto integral_1d_non_uniform = integral::cumtrapz(x_vals, data_1d);
        utils::WriteData("test_result/integral/integral_1d_non_uniform.txt",
                         integral_1d_non_uniform,
                         1,
                         3e4);

        // test matrix integral
        auto integral_matrix = integral::cumtrapz(ori_matrix, 1.0);
        utils::WriteData("test_result/integral/integral_matrix.txt",
                         std::vector<double>(integral_matrix.data(),
                                             integral_matrix.data()
                                                 + integral_matrix.size()),
                         integral_matrix.cols(),
                         integral_matrix.rows());

        // test matrix integral with non-uniform x
        auto integral_matrix_non_uniform =
            integral::cumtrapz(x_vals, ori_matrix);
        utils::WriteData(
            "test_result/integral/integral_matrix_non_uniform.txt",
            std::vector<double>(integral_matrix_non_uniform.data(),
                                integral_matrix_non_uniform.data()
                                    + integral_matrix_non_uniform.size()),
            integral_matrix_non_uniform.cols(),
            integral_matrix_non_uniform.rows());
    }
    catch (const std::exception &e)
    {
        std::cerr << "Integral test failed: " << e.what() << std::endl;
        return -1;
    }

    std::cout << "Integral test completed successfully." << std::endl;
    return 0;
}


int main() { return test_integral(); }