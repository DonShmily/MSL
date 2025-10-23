#include <vector>

#include "interp.hpp"
#include "interp/interp_1d_akima.hpp"
#include "interp/interp_1d_linear.hpp"
#include "interp/interp_1d_spline.hpp"
#include "matrix.hpp"
#include "utils/data_io.hpp"

using namespace msl;

int test_interp()
{
    auto ori_data =
        utils::ReadData("test_result/interp/test_interp_data.txt", 4, 11);
    auto data_x = std::vector<double>(ori_data.begin(), ori_data.begin() + 11);
    auto data_y1 =
        std::vector<double>(ori_data.begin() + 11, ori_data.begin() + 22);
    auto data_y2 =
        std::vector<double>(ori_data.begin() + 22, ori_data.begin() + 33);
    auto data_y3 =
        std::vector<double>(ori_data.begin() + 33, ori_data.begin() + 44);
    std::vector<double> new_x = {
        -1, 0, 0.5, 2, 2.2, 5, 5.6, 8, 8.7, 9.1, 10, 11.5};

    try
    {
        // test linear interp
        auto lin_interp_y1 = interp::interp1_linear(data_x, data_y1, new_x);
        auto lin_interp_y2 = interp::interp1_linear(data_x, data_y2, new_x);
        auto lin_interp_y3 = interp::interp1_linear(data_x, data_y3, new_x);
        std::vector<double> lin_interp_all(4 * new_x.size());
        std::copy(new_x.begin(), new_x.end(), lin_interp_all.begin());
        std::copy(lin_interp_y1.begin(),
                  lin_interp_y1.end(),
                  lin_interp_all.begin() + new_x.size());
        std::copy(lin_interp_y2.begin(),
                  lin_interp_y2.end(),
                  lin_interp_all.begin() + 2 * new_x.size());
        std::copy(lin_interp_y3.begin(),
                  lin_interp_y3.end(),
                  lin_interp_all.begin() + 3 * new_x.size());
        utils::WriteData("test_result/interp/linear_interp_all.txt",
                         lin_interp_all,
                         4,
                         new_x.size());

        // test spline interp
        auto spline_interp_y1 = interp::interp1_cubic(data_x, data_y1, new_x);
        auto spline_interp_y2 = interp::interp1_cubic(data_x, data_y2, new_x);
        auto spline_interp_y3 = interp::interp1_cubic(data_x, data_y3, new_x);
        std::vector<double> spline_interp_all(4 * new_x.size());
        std::copy(new_x.begin(), new_x.end(), spline_interp_all.begin());
        std::copy(spline_interp_y1.begin(),
                  spline_interp_y1.end(),
                  spline_interp_all.begin() + new_x.size());
        std::copy(spline_interp_y2.begin(),
                  spline_interp_y2.end(),
                  spline_interp_all.begin() + 2 * new_x.size());
        std::copy(spline_interp_y3.begin(),
                  spline_interp_y3.end(),
                  spline_interp_all.begin() + 3 * new_x.size());
        utils::WriteData("test_result/interp/spline_interp_all.txt",
                         spline_interp_all,
                         4,
                         new_x.size());

        // test akima interp
        auto akima_interp_y1 = interp::interp1_akima(data_x, data_y1, new_x);
        auto akima_interp_y2 = interp::interp1_akima(data_x, data_y2, new_x);
        auto akima_interp_y3 = interp::interp1_akima(data_x, data_y3, new_x);
        std::vector<double> akima_interp_all(4 * new_x.size());
        std::copy(new_x.begin(), new_x.end(), akima_interp_all.begin());
        std::copy(akima_interp_y1.begin(),
                  akima_interp_y1.end(),
                  akima_interp_all.begin() + new_x.size());
        std::copy(akima_interp_y2.begin(),
                  akima_interp_y2.end(),
                  akima_interp_all.begin() + 2 * new_x.size());
        std::copy(akima_interp_y3.begin(),
                  akima_interp_y3.end(),
                  akima_interp_all.begin() + 3 * new_x.size());
        utils::WriteData("test_result/interp/akima_interp_all.txt",
                         akima_interp_all,
                         4,
                         new_x.size());
    }
    catch (const std::exception &e)
    {
        std::cerr << "Interp test failed: " << e.what() << std::endl;
        return -1;
    }

    std::cout << "Interp test completed successfully." << std::endl;
    return 0;
}

int main() { return test_interp(); }
