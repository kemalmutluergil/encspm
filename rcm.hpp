#pragma once 

#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <random>
#include <set>
#include <limits>
#include <cmath>
#include <string>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <sstream>
#include <tuple>
#include <stdexcept>
#include <sstream>

using namespace std;

#include "csr.hpp"

// A struct to hold calculated statistics for a set of values (e.g., row degrees).
struct DistributionStats {
    double min = 0.0;
    double max = 0.0;
    double mean = 0.0;
    double median = 0.0;
    double stdev = 0.0;
    size_t non_zero_count = 0;
};

// Contains all functions related to the Reverse Cuthill-McKee algorithm.
namespace rcm {

// --- Core RCM and Helper Functions ---
void find_rcm_ordering(const CSRMatrix& A,
                       std::vector<size_t>& row_perm,
                       std::vector<size_t>& col_perm) 
{
    size_t n = A.rows;
    row_perm.assign(n, std::numeric_limits<size_t>::max());
    col_perm.assign(A.cols, std::numeric_limits<size_t>::max());

    std::vector<bool> visited(n, false);
    std::vector<size_t> degree(n, 0);

    // compute degree of each row (number of nonzeros)
    for (size_t i = 0; i < n; i++) {
        degree[i] = A.row_ptr[i+1] - A.row_ptr[i];
    }

    std::vector<size_t> ordering;
    ordering.reserve(n);

    // BFS from the lowest degree node each time (RCM style)
    for (size_t start = 0; start < n; start++) {
        if (visited[start]) continue;

        // BFS queue
        std::queue<size_t> q;
        q.push(start);
        visited[start] = true;

        std::vector<size_t> component;

        while (!q.empty()) {
            size_t u = q.front();
            q.pop();
            component.push_back(u);

            // Collect neighbors
            std::vector<size_t> neighbors;
            for (size_t k = A.row_ptr[u]; k < A.row_ptr[u+1]; k++) {
                size_t v = A.col_idx[k];
                if (v < n && !visited[v]) { // assuming square or symmetric structure
                    neighbors.push_back(v);
                }
            }

            // Sort neighbors by degree before pushing
            std::sort(neighbors.begin(), neighbors.end(),
                      [&](size_t a, size_t b) { return degree[a] < degree[b]; });

            for (size_t v : neighbors) {
                if (!visited[v]) {
                    visited[v] = true;
                    q.push(v);
                }
            }
        }

        // Reverse component order (RCM trick)
        std::reverse(component.begin(), component.end());
        ordering.insert(ordering.end(), component.begin(), component.end());
    }

    // Build forward permutation (old -> new)
    for (size_t new_idx = 0; new_idx < ordering.size(); new_idx++) {
        size_t old_idx = ordering[new_idx];
        row_perm[old_idx] = new_idx;
        col_perm[old_idx] = new_idx; // same permutation for symmetric matrix
    }
}

/**
 * @brief Applies RCM to a rectangular matrix.
 * @param[out] row_P_inv The calculated row permutation (new->old).
 * @param[out] col_P_inv The calculated column permutation (new->old).
 */
// void rcm_rectangular(const CSRMatrix &A,
//                      std::vector<size_t> &row_P_inv,
//                      std::vector<size_t> &col_P_inv) {
//     CSRMatrix B = build_bipartite_sym(A);
//     auto P_inv = rcm_ordering(B);
    
//     row_P_inv.clear();
//     col_P_inv.clear();
    
//     size_t m = A.rows;
//     for (auto p : P_inv) {
//         if (p < m) {
//             row_P_inv.push_back(p);
//         } else {
//             col_P_inv.push_back(p - m);
//         }
//     }
// }

// --- Performance Metrics ---

double bandwidth(const CSRMatrix &A) {
    if (A.col_idx.empty()) return 0.0;
    double max_bw = 0.0;
    for (size_t i = 0; i < A.rows; ++i) {
        for (size_t k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
            max_bw = std::max(max_bw, std::abs(static_cast<double>(i) - static_cast<double>(A.col_idx[k])));
        }
    }
    return max_bw;
}

double average_bandwidth(const CSRMatrix &A) {
    if (A.col_idx.empty()) return 0.0;
    double total_bw = 0.0;
    for (size_t i = 0; i < A.rows; ++i) {
        for (size_t k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
            total_bw += std::abs(static_cast<double>(i) - static_cast<double>(A.col_idx[k]));
        }
    }
    return total_bw / A.col_idx.size();
}

size_t profile_size(const CSRMatrix &A) {
    size_t profile = 0;
    for (size_t i = 0; i < A.rows; ++i) {
        if (A.row_ptr[i] < A.row_ptr[i + 1]) {
            size_t min_j = A.col_idx[A.row_ptr[i]];
            size_t max_j = A.col_idx[A.row_ptr[i + 1] - 1];
            profile += (max_j - min_j + 1);
        }
    }
    return profile;
}

size_t max_row_profile(const CSRMatrix &A) {
    size_t max_prof = 0;
    for (size_t i = 0; i < A.rows; ++i) {
        if (A.row_ptr[i] < A.row_ptr[i + 1]) {
            size_t min_j = A.col_idx[A.row_ptr[i]];
            size_t max_j = A.col_idx[A.row_ptr[i + 1] - 1];
            max_prof = std::max(max_prof, max_j - min_j);
        }
    }
    return max_prof;
}

} // namespace rcm

// --- Utility and Test Functions ---

/**
 * @brief Calculates distribution statistics for a given vector of data.
 * The stats (min, max, mean, median, stdev) are computed over the non-zero elements.
 * @param data The input vector (e.g., of row or column degrees).
 * @return A DistributionStats struct containing the results.
 */
DistributionStats calculate_stats(const std::vector<double>& data) {
    DistributionStats stats;
    std::vector<double> non_zero_data;
    non_zero_data.reserve(data.size());
    for(double val : data) {
        if (val > 1e-6) { // Consider non-zero
            non_zero_data.push_back(val);
        }
    }

    stats.non_zero_count = non_zero_data.size();
    if(stats.non_zero_count == 0) return stats;

    std::sort(non_zero_data.begin(), non_zero_data.end());

    stats.max = non_zero_data.back();
    stats.min = non_zero_data.front();
    
    if (stats.non_zero_count % 2 == 0) {
        stats.median = (non_zero_data[stats.non_zero_count / 2 - 1] + non_zero_data[stats.non_zero_count / 2]) / 2.0;
    } else {
        stats.median = non_zero_data[stats.non_zero_count / 2];
    }

    double sum = std::accumulate(non_zero_data.begin(), non_zero_data.end(), 0.0);
    stats.mean = sum / stats.non_zero_count;

    double sq_sum_diff = 0.0;
    for (double val : non_zero_data) {
        sq_sum_diff += (val - stats.mean) * (val - stats.mean);
    }
    stats.stdev = std::sqrt(sq_sum_diff / stats.non_zero_count);

    return stats;
}


CSRMatrix make_random_sparse(size_t m, size_t n, size_t avg_per_row) {
    CSRMatrix A;
    A.rows = m; A.cols = n;
    A.row_ptr.resize(m + 1, 0);
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> dist(0, n - 1);

    for (size_t i = 0; i < m; ++i) {
        std::set<size_t> unique_cols;
        size_t num_elems = std::max(size_t(1), avg_per_row);
        while (unique_cols.size() < num_elems && unique_cols.size() < n) {
            unique_cols.insert(dist(rng));
        }
        A.row_ptr[i + 1] = A.row_ptr[i] + unique_cols.size();
        A.col_idx.insert(A.col_idx.end(), unique_cols.begin(), unique_cols.end());
    }
    A.values.resize(A.col_idx.size(), 1.0);
    return A;
}

CSRMatrix symmetrize_pattern(const CSRMatrix &A) {
    if (A.rows != A.cols) {
        std::cerr << "Warning: symmetrize_pattern called on non-square matrix." << std::endl;
        return A;
    }
    size_t n = A.rows;
    std::vector<std::set<size_t>> rows(n);
    for (size_t i = 0; i < n; ++i) {
        rows[i].insert(i);
        for (size_t k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
            size_t j = A.col_idx[k];
            rows[i].insert(j);
            rows[j].insert(i);
        }
    }
    CSRMatrix S;
    S.rows = S.cols = n;
    S.row_ptr.resize(n + 1, 0);
    for (size_t i = 0; i < n; ++i) {
        S.row_ptr[i + 1] = S.row_ptr[i] + rows[i].size();
        S.col_idx.insert(S.col_idx.end(), rows[i].begin(), rows[i].end());
    }
    S.values.resize(S.col_idx.size(), 1.0);
    return S;
}

void print_dense(const CSRMatrix &A) {
    std::vector<std::vector<char>> M(A.rows, std::vector<char>(A.cols, '.'));
    for (size_t i = 0; i < A.rows; ++i) {
        for (size_t k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
            if (A.col_idx[k] < A.cols) {
                M[i][A.col_idx[k]] = '#';
            }
        }
    }
    for (size_t i = 0; i < A.rows; ++i) {
        for (size_t j = 0; j < A.cols; ++j) {
            std::cout << M[i][j];
        }
        std::cout << "\n";
    }
}

/**
 * @brief Reports all performance metrics for a before and after matrix.
 */
void report_stats(const CSRMatrix &orig, const CSRMatrix &perm) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "-------------------|----------|----------\n";
    std::cout << "Metric             |   Before |    After\n";
    std::cout << "-------------------|----------|----------\n";
    std::cout << "Max Bandwidth      | " << std::setw(8) << rcm::bandwidth(orig) << " | " << std::setw(8) << rcm::bandwidth(perm) << "\n";
    std::cout << "Average Bandwidth  | " << std::setw(8) << rcm::average_bandwidth(orig) << " | " << std::setw(8) << rcm::average_bandwidth(perm) << "\n";
    std::cout << "Profile Size       | " << std::setw(8) << rcm::profile_size(orig) << " | " << std::setw(8) << rcm::profile_size(perm) << "\n";
    std::cout << "Max Row Profile    | " << std::setw(8) << rcm::max_row_profile(orig) << " | " << std::setw(8) << rcm::max_row_profile(perm) << "\n";
    
    // --- Detailed Degree Distribution ---
    auto get_degrees = [](const CSRMatrix& A, bool by_row) -> std::vector<double> {
        size_t dim = by_row ? A.rows : A.cols;
        std::vector<double> degrees(dim, 0.0);
        if (by_row) {
            for(size_t i = 0; i < A.rows; ++i) degrees[i] = static_cast<double>(A.row_ptr[i+1] - A.row_ptr[i]);
        } else {
            for(size_t idx : A.col_idx) if (idx < degrees.size()) degrees[idx]++;
        }
        return degrees;
    };
    
    auto row_deg_orig = get_degrees(orig, true);
    auto col_deg_orig = get_degrees(orig, false);
    auto row_deg_perm = get_degrees(perm, true);
    auto col_deg_perm = get_degrees(perm, false);
    
    DistributionStats row_stats_orig = calculate_stats(row_deg_orig);
    DistributionStats col_stats_orig = calculate_stats(col_deg_orig);
    DistributionStats row_stats_perm = calculate_stats(row_deg_perm);
    DistributionStats col_stats_perm = calculate_stats(col_deg_perm);
    std::cout << "-------------------|----------|----------\n";
    std::cout << "Statistic (Rows)   |   Before |    After\n";
    std::cout << "-------------------|----------|----------\n";
    std::cout << "Non-Zero Rows      | " << std::setw(8) << row_stats_orig.non_zero_count << " | " << std::setw(8) << row_stats_perm.non_zero_count << "\n";
    std::cout << "Min Degree         | " << std::setw(8) << row_stats_orig.min << " | " << std::setw(8) << row_stats_perm.min << "\n";
    std::cout << "Max Degree         | " << std::setw(8) << row_stats_orig.max << " | " << std::setw(8) << row_stats_perm.max << "\n";
    std::cout << "Mean Degree        | " << std::setw(8) << row_stats_orig.mean << " | " << std::setw(8) << row_stats_perm.mean << "\n";
    std::cout << "Median Degree      | " << std::setw(8) << row_stats_orig.median << " | " << std::setw(8) << row_stats_perm.median << "\n";
    std::cout << "Stdev Degree       | " << std::setw(8) << row_stats_orig.stdev << " | " << std::setw(8) << row_stats_perm.stdev << "\n";

    std::cout << "\nStatistic (Cols)   |   Before |    After\n";
    std::cout << "-------------------|----------|----------\n";
    std::cout << "Non-Zero Cols      | " << std::setw(8) << col_stats_orig.non_zero_count << " | " << std::setw(8) << col_stats_perm.non_zero_count << "\n";
    std::cout << "Min Degree         | " << std::setw(8) << col_stats_orig.min << " | " << std::setw(8) << col_stats_perm.min << "\n";
    std::cout << "Max Degree         | " << std::setw(8) << col_stats_orig.max << " | " << std::setw(8) << col_stats_perm.max << "\n";
    std::cout << "Mean Degree        | " << std::setw(8) << col_stats_orig.mean << " | " << std::setw(8) << col_stats_perm.mean << "\n";
    std::cout << "Median Degree      | " << std::setw(8) << col_stats_orig.median << " | " << std::setw(8) << col_stats_perm.median << "\n";
    std::cout << "Stdev Degree       | " << std::setw(8) << col_stats_orig.stdev << " | " << std::setw(8) << col_stats_perm.stdev << "\n";
}

/**
 * @brief Checks if a CSR matrix is structurally symmetric.
 * It assumes the column indices within each row are sorted.
 * @param A The input matrix.
 * @return True if the matrix pattern is symmetric, false otherwise.
 */
bool is_symmetric(const CSRMatrix& A) {
    // A symmetric matrix must be square.
    if (A.rows != A.cols) {
        return false;
    }

    size_t n = A.rows;
    // For every non-zero entry (i, j), we must find a non-zero entry (j, i).
    for (size_t i = 0; i < n; ++i) {
        for (size_t k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
            size_t j = A.col_idx[k];

            // Now, we need to efficiently check if A(j, i) exists.
            // We can do this with a binary search on the j-th row's column indices.
            auto row_j_start = A.col_idx.begin() + A.row_ptr[j];
            auto row_j_end = A.col_idx.begin() + A.row_ptr[j + 1];
            
            bool found_transpose = std::binary_search(row_j_start, row_j_end, i);

            if (!found_transpose) {
                // If we found (i, j) but not (j, i), it's not symmetric.
                return false;
            }
        }
    }

    // If we checked all non-zero elements and found their transpose counterparts, it's symmetric.
    return true;
}

