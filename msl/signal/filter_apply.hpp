/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2025, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: filter_apply.hpp
** -----
** File Created: Tuesday, 14th October 2025 11:12:27
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Sunday, 14th December 2025 17:05:23
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_FILTER_APPLY_HPP
#define MSL_FILTER_APPLY_HPP

#include <span>
#include <vector>

#include "filter.hpp"
#include "filter_design.hpp"
#include "filtfilt.hpp"
#include "matrix/real_matrix_owned.hpp"

namespace msl::signal
{
// ============================================================================
// Apply filter to signal (span overload)
// ============================================================================
/**
 * @brief Apply filter to signal (vector overload)
 */
inline std::vector<double> filter(const std::vector<double> &signal,
                                  const FilterCoefficients &coeffs)
{
    return filter(std::span<const double>(signal), coeffs);
}

/**
 * @brief Zero-phase filtering (vector overload)
 */
inline std::vector<double> filtfilt(const std::vector<double> &signal,
                                    const FilterCoefficients &coeffs)
{
    return filtfilt(std::span<const double>(signal), coeffs);
}

// ============================================================================
// Matrix filtering (column-wise)
// ============================================================================

/**
 * @brief Apply filter to each column of a matrix
 *
 * @param signals Matrix where each column is a signal
 * @param coeffs Filter coefficients
 * @return Filtered matrix
 */
inline matrix::matrixd filter_columns(const matrix::matrixd &signals,
                                      const FilterCoefficients &coeffs)
{
    matrix::matrixd output(signals.rows(), signals.cols());

    for (size_t j = 0; j < signals.cols(); ++j)
    {
        auto col_span = signals.column(j);
        auto filtered = filter(col_span, coeffs);

        for (size_t i = 0; i < signals.rows(); ++i)
        {
            output(i, j) = filtered[i];
        }
    }

    return output;
}

/**
 * @brief Apply zero-phase filter to each column of a matrix
 */
inline matrix::matrixd filtfilt_columns(const matrix::matrixd &signals,
                                        const FilterCoefficients &coeffs)
{
    matrix::matrixd output(signals.rows(), signals.cols());

    for (size_t j = 0; j < signals.cols(); ++j)
    {
        auto col_span = signals.column(j);
        auto filtered = filtfilt(col_span, coeffs);

        for (size_t i = 0; i < signals.rows(); ++i)
        {
            output(i, j) = filtered[i];
        }
    }

    return output;
}

} // namespace msl::signal

#endif // MSL_FILTER_APPLY_HPP