#include <vector>

#include "difference.hpp"
#include "matrix.hpp"
#include "utils/data_io.hpp"


using namespace msl;
int test_difference()
{
    auto ori_data = utils::ReadData("KunmingSSJY.txt", 6, 3e4);
    auto data_1d =
        std::vector<double>(ori_data.begin(), ori_data.begin() + 3e4);
    matrix::matrixd ori_matrix(3e4, 6, std::span<const double>(ori_data));
    try
    {
        // 1d differentiation
        auto diff_1d = difference::diff(data_1d);
        utils::WriteData(
            "test_result/difference/diff_1d.txt", diff_1d, 1, diff_1d.size());

        // 2d differentiation along rows
        auto diff_2d_row = difference::diff(ori_matrix, 0);
        auto diff_2d_row_vec = std::vector<double>(
            diff_2d_row.data(), diff_2d_row.data() + diff_2d_row.size());
        utils::WriteData("test_result/difference/diff_2d_row.txt",
                         diff_2d_row_vec,
                         diff_2d_row.cols(),
                         diff_2d_row.rows());

        // 2d differentiation along columns
        auto diff_2d_col = difference::diff(ori_matrix, 1);
        auto diff_2d_col_vec = std::vector<double>(
            diff_2d_col.data(), diff_2d_col.data() + diff_2d_col.size());
        utils::WriteData("test_result/difference/diff_2d_col.txt",
                         diff_2d_col_vec,
                         diff_2d_col.cols(),
                         diff_2d_col.rows());

        // 1d central gradient
        auto cent_grad_1d = difference::central_gradient(data_1d, 1.0);
        utils::WriteData("test_result/difference/cent_grad_1d.txt",
                         cent_grad_1d,
                         1,
                         cent_grad_1d.size());

        // 2d central gradient along rows
        auto cent_grad_2d_row =
            difference::central_gradient(ori_matrix, 1.0, 0);
        auto cent_grad_2d_row_vec = std::vector<double>(
            cent_grad_2d_row.data(),
            cent_grad_2d_row.data() + cent_grad_2d_row.size());
        utils::WriteData("test_result/difference/cent_grad_2d_row.txt",
                         cent_grad_2d_row_vec,
                         cent_grad_2d_row.cols(),
                         cent_grad_2d_row.rows());

        // 2d central gradient
        auto cent_grad_2d_col = difference::central_gradient2d(ori_matrix, 1.0);
        auto cent_grad_2d_col_vec1 = std::vector<double>(
            cent_grad_2d_col.first.data(),
            cent_grad_2d_col.first.data() + cent_grad_2d_col.first.size());
        utils::WriteData("test_result/difference/cent_grad_2d_col_x.txt",
                         cent_grad_2d_col_vec1,
                         cent_grad_2d_col.first.cols(),
                         cent_grad_2d_col.first.rows());
        auto cent_grad_2d_col_vec2 = std::vector<double>(
            cent_grad_2d_col.second.data(),
            cent_grad_2d_col.second.data() + cent_grad_2d_col.second.size());
        utils::WriteData("test_result/difference/cent_grad_2d_col_y.txt",
                         cent_grad_2d_col_vec2,
                         cent_grad_2d_col.second.cols(),
                         cent_grad_2d_col.second.rows());
    }
    catch (const std::exception &e)
    {
        std::cerr << "Differentiation test failed: " << e.what() << std::endl;
        return -1;
    }

    std::cout << "Differentiation test completed successfully." << std::endl;
    return 0;
}


int main() { return test_difference(); }