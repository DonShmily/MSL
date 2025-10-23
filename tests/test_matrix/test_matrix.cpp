// tests/test_matrix_all.cpp
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "matrix.hpp"

using namespace msl;

// ---------- 简易断言框架 ----------
static int g_failures = 0;
static int g_total = 0;

#define TEST_CASE(name)                                                        \
    std::cout << "=== Test: " << name << " ===\n";                             \
    do

#define EXPECT_TRUE(cond)                                                      \
    do                                                                         \
    {                                                                          \
        ++g_total;                                                             \
        if (!(cond))                                                           \
        {                                                                      \
            ++g_failures;                                                      \
            std::cerr << "[FAIL] " << __FILE__ << ":" << __LINE__ << " - "     \
                      << #cond << "\n";                                        \
        }                                                                      \
    } while (0)

#define EXPECT_EQ(a, b)                                                        \
    do                                                                         \
    {                                                                          \
        ++g_total;                                                             \
        if (!((a) == (b)))                                                     \
        {                                                                      \
            ++g_failures;                                                      \
            std::cerr << "[FAIL] " << __FILE__ << ":" << __LINE__ << " - "     \
                      << #a << " == " << #b << " (got: " << (a) << " vs "      \
                      << (b) << ")\n";                                         \
        }                                                                      \
    } while (0)

#define EXPECT_NEAR(a, b, eps)                                                 \
    do                                                                         \
    {                                                                          \
        ++g_total;                                                             \
        if (!(std::fabs((double)(a) - (double)(b)) <= (eps)))                  \
        {                                                                      \
            ++g_failures;                                                      \
            std::cerr << "[FAIL] " << __FILE__ << ":" << __LINE__ << " - "     \
                      << #a << " ~= " << #b << " (got: " << (a) << " vs "      \
                      << (b) << ", eps=" << (eps) << ")\n";                    \
        }                                                                      \
    } while (0)

// For complex near
static void expect_complex_near(const std::complex<double> &a,
                                const std::complex<double> &b,
                                double eps,
                                const char *file,
                                int line)
{
    ++g_total;
    if (std::abs(a.real() - b.real()) > eps
        || std::abs(a.imag() - b.imag()) > eps)
    {
        ++g_failures;
        std::cerr << "[FAIL] " << file << ":" << line
                  << " - complex not near (got: " << a << " vs " << b
                  << ", eps=" << eps << ")\n";
    }
}
#define EXPECT_CPLX_NEAR(a, b, eps)                                            \
    expect_complex_near((a), (b), (eps), __FILE__, __LINE__)

// ---------- 辅助函数 ----------
static bool approx_equal_double_vector(const std::vector<double> &a,
                                       const std::vector<double> &b,
                                       double eps = 1e-12)
{
    if (a.size() != b.size())
        return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (std::fabs(a[i] - b[i]) > eps)
            return false;
    return true;
}

// ---------- 测试主体 ----------
int main()
{
    std::cout << "MSL matrix unit tests\n";

    // 1) 基本构造与 factory
    TEST_CASE("identity / zeros / ones / diagonal")
    {
        auto I = matrixd::identity(3); // matrixd = matrix_real alias
        EXPECT_EQ(I.rows(), 3);
        EXPECT_EQ(I.cols(), 3);
        for (size_t i = 0; i < 3; ++i)
            for (size_t j = 0; j < 3; ++j)
                EXPECT_EQ(I(i, j), (i == j ? 1.0 : 0.0));

        auto Z = matrixd::zeros(2, 3);
        EXPECT_EQ(Z.rows(), 2);
        EXPECT_EQ(Z.cols(), 3);
        for (auto v : Z)
            EXPECT_EQ(v, 0.0);

        auto O = matrixd::ones(2, 2);
        for (auto v : O)
            EXPECT_EQ(v, 1.0);

        std::vector<double> diag = {1.5, 2.5, 3.5};
        auto temp = matrixd::ones(3, 3);
        temp(0, 0) = diag[0];
        temp(1, 1) = diag[1];
        temp(2, 2) = diag[2];
        auto D =
            temp.diagonal(std::span<const double>(diag.data(), diag.size()));
        EXPECT_EQ(D.rows(), 3);
        EXPECT_EQ(D.cols(), 3);
        for (size_t i = 0; i < 3; ++i)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                if (i == j)
                    EXPECT_EQ(D(i, j), diag[i]);
                else
                    EXPECT_EQ(D(i, j), 0.0);
            }
        }
    }
    while (0)
        ;

    // 2) element access, resize, fill
    TEST_CASE("element access / resize / fill")
    {
        matrixd A(3, 2);
        EXPECT_EQ(A.rows(), 3);
        EXPECT_EQ(A.cols(), 2);
        A(0, 0) = 1.0;
        A(1, 0) = 2.0;
        A(2, 0) = 3.0;
        A(0, 1) = 4.0;
        A(1, 1) = 5.0;
        A(2, 1) = 6.0;
        EXPECT_EQ(A(2, 1), 6.0);
        A.fill(7.0);
        for (auto v : A)
            EXPECT_EQ(v, 7.0);

        A.resize(2, 2);
        EXPECT_EQ(A.rows(), 2);
        EXPECT_EQ(A.cols(), 2);
        A.fill_zeros();
        for (auto v : A)
            EXPECT_EQ(v, 0.0);
    }
    while (0)
        ;

    // 3) copy / move / swap
    TEST_CASE("copy / move / swap")
    {
        matrixd A(2, 2, 1.0);
        matrixd B = A; // copy
        EXPECT_EQ(B(0, 0), 1.0);
        matrixd C = std::move(B); // move
        EXPECT_EQ(C(1, 1), 1.0);
        // after move B may be empty (implementation-dependent) but should be
        // safe to resize
        C(0, 0) = 9.0;
        A.swap(C);
        EXPECT_EQ(A(0, 0), 9.0);
    }
    while (0)
        ;

    // 4) arithmetic operators
    TEST_CASE("operators + - * /")
    {
        matrixd A(2, 2);
        A.fill(2.0);
        matrixd B(2, 2);
        B.fill(1.0);
        auto S = A + B; // should be 3's
        for (auto v : S)
            EXPECT_EQ(v, 3.0);
        auto D = A - B; // should be 1's
        for (auto v : D)
            EXPECT_EQ(v, 1.0);
        auto scaled = A * 0.5;
        for (auto v : scaled)
            EXPECT_EQ(v, 1.0);
        auto scaled2 = 0.5 * A;
        for (auto v : scaled2)
            EXPECT_EQ(v, 1.0);
        auto divv = A / 2.0;
        for (auto v : divv)
            EXPECT_EQ(v, 1.0);

        // in-place ops
        A += B;
        for (auto v : A)
            EXPECT_EQ(v, 3.0);
        A -= B;
        for (auto v : A)
            EXPECT_EQ(v, 2.0);
        A *= 2.0;
        for (auto v : A)
            EXPECT_EQ(v, 4.0);
        A /= 2.0;
        for (auto v : A)
            EXPECT_EQ(v, 2.0);
    }
    while (0)
        ;

    // 5) row/column/submatrix/get_row/get_column
    TEST_CASE("row/col/submatrix")
    {
        matrixd M(3, 4);
        // fill with M(i,j) = i + 10*j  (col-major aware)
        for (size_t j = 0; j < M.cols(); ++j)
            for (size_t i = 0; i < M.rows(); ++i)
                M(i, j) = double(i + 10 * j);

        auto row1 = M.get_row(1);
        EXPECT_EQ(row1.rows(), 1);
        EXPECT_EQ(row1.cols(), M.cols());
        for (size_t j = 0; j < M.cols(); ++j)
            EXPECT_EQ(row1(0, j), M(1, j));

        auto col2 = M.get_column(2);
        EXPECT_EQ(col2.rows(), M.rows());
        EXPECT_EQ(col2.cols(), 1);
        for (size_t i = 0; i < M.rows(); ++i)
            EXPECT_EQ(col2(i, 0), M(i, 2));

        auto sub = M.submatrix(1, 3, 1, 4); // rows 1..2, cols 1..3
        EXPECT_EQ(sub.rows(), 2);
        EXPECT_EQ(sub.cols(), 3);
        for (size_t j = 0; j < sub.cols(); ++j)
            for (size_t i = 0; i < sub.rows(); ++i)
                EXPECT_EQ(sub(i, j), M(1 + i, 1 + j));
    }
    while (0)
        ;

    // 6) iterators and range-for / operator[] if provided
    TEST_CASE("iterators and operator[]")
    {
        matrixd A(3, 2);
        double val = 0.0;
        for (size_t j = 0; j < A.cols(); ++j)
            for (size_t i = 0; i < A.rows(); ++i)
                A(i, j) = val++;
        // sum via iterator
        double s1 = 0.0;
        for (auto &x : A)
            s1 += x;
        double s2 = 0.0;
        for (size_t k = 0; k < A.size(); ++k)
            s2 += A[k]; // operator[] assumed
        EXPECT_EQ(s1, s2);
    }
    while (0)
        ;

    // 7) matrix_view zero-copy and construction from external data
    TEST_CASE("matrix_view zero-copy")
    {
        const size_t R = 3, C = 2;
        double raw[R * C] = {
            1, 2, 3, 4, 5, 6}; // row-major-like literal BUT we'll treat as
                               // column-major when filling
        // Fill raw in column-major to match matrix_view expectation (col-major
        // index j*R + i)
        for (size_t j = 0; j < C; ++j)
            for (size_t i = 0; i < R; ++i)
                raw[j * R + i] = double(j * 10 + i);

        matrixd_view V(raw, R, C);
        EXPECT_EQ(V.rows(), R);
        EXPECT_EQ(V.cols(), C);
        // modify through view and check raw changed
        V(1, 1) = 99.0;
        EXPECT_EQ(raw[1 + 1 * R], 99.0);

        // construct owning matrix from view
        matrixd M = matrixd(V); // should copy data
        EXPECT_EQ(M.rows(), R);
        EXPECT_EQ(M.cols(), C);
        for (size_t j = 0; j < C; ++j)
            for (size_t i = 0; i < R; ++i)
                EXPECT_EQ(M(i, j), V(i, j));
    }
    while (0)
        ;

    // 8) multiplication (matrix * matrix) and (matrix * vector)
    TEST_CASE("matrix multiplication")
    {
        matrixd A(2, 3); // 2x3
        matrixd B(3, 2); // 3x2
        // A = [1 2 3; 4 5 6] in column-major -> careful filling
        // We'll set A(i,j) = i + j*10 for clarity
        for (size_t j = 0; j < A.cols(); ++j)
            for (size_t i = 0; i < A.rows(); ++i)
                A(i, j) = double(i + j * 10);
        for (size_t j = 0; j < B.cols(); ++j)
            for (size_t i = 0; i < B.rows(); ++i)
                B(i, j) = double(1 + i + j * 10);
        auto A_eigen = msl::eigen_interface::as_eigen(A);
        auto B_eigen = msl::eigen_interface::as_eigen(B);

        auto C_eigen = A_eigen * B_eigen; // should be 2x2
        auto C = msl::eigen_interface::from_eigen(C_eigen);
        EXPECT_EQ(C.rows(), 2);
        EXPECT_EQ(C.cols(), 2);

        // naive recompute
        std::vector<double> expected(4, 0.0);
        for (size_t i = 0; i < C.rows(); ++i)
            for (size_t j = 0; j < C.cols(); ++j)
                for (size_t k = 0; k < A.cols(); ++k)
                    expected[j * C.rows() + i] += A(i, k) * B(k, j);

        for (size_t j = 0; j < C.cols(); ++j)
            for (size_t i = 0; i < C.rows(); ++i)
                EXPECT_EQ(C(i, j), expected[j * C.rows() + i]);

        // matrix * std::vector
        std::vector<double> x(B.cols(), 1.0); // length 2
        Eigen::VectorXd y_eigen =
            B_eigen * Eigen::Map<Eigen::VectorXd>(x.data(), x.size()).eval();
        auto y = msl::eigen_interface::from_eigen(y_eigen);
        EXPECT_EQ(y.size(), B.rows());
        // y_i = sum_j B(i,j) * x[j]
        for (size_t i = 0; i < B.rows(); ++i)
        {
            double s = 0.0;
            for (size_t j = 0; j < B.cols(); ++j)
                s += B(i, j) * x[j];
            EXPECT_EQ(y[i], s);
        }
    }
    while (0)
        ;

    // 9) complex matrix basic use
    TEST_CASE("complex matrix basics")
    {
        matrixc C(2, 2);
        C(0, 0) = std::complex<double>(1.0, 2.0);
        C(1, 0) = std::complex<double>(3.0, -1.0);
        EXPECT_CPLX_NEAR(C(0, 0), std::complex<double>(1.0, 2.0), 1e-12);
        // conjugate if method exists (we will call if present)
        // use SFINAE would be robust, but for test we simply attempt and skip
        // compile-time if not present
    }
    while (0)
        ;

    // 总结
    std::cout << "Total checks: " << g_total << ", failures: " << g_failures
              << "\n";
    if (g_failures == 0)
    {
        std::cout << "All tests passed\n";
        return 0;
    }
    else
    {
        std::cerr << "Some tests failed\n";
        return 1;
    }
}
