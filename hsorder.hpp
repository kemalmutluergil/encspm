#pragma once
#include <queue>
#include <unordered_map>
#include <utility>

using namespace std;

#include "csr.hpp"
#include "utils.hpp"

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ (h2 << 1); // or use boost::hash_combine
    }
};

const int max_iter = 20;
void hsorder(const CSRMatrix& A, std::vector<size_t>& row_perm, std::vector<size_t>& col_perm) {
    size_t n = A.rows;
    row_perm.resize(n);
    col_perm.resize(n);

    // Initial identity permutation
    for (size_t i = 0; i < n; i++) {
        row_perm[i] = i;
        col_perm[i] = i;
    }

    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(row_perm.begin(), row_perm.end(), g);
    std::shuffle(col_perm.begin(), col_perm.end(), g);

    size_t best_score = count_nonempty_hs_diagonals(permute(A, row_perm, col_perm));
    bool improved = true;

    while (improved) {
        improved = false;

        // Try row swaps
        for (size_t i = 0; i < n; i++) {
            for (size_t j = i+1; j < n; j++) {
                std::swap(row_perm[i], row_perm[j]);
                size_t score = count_nonempty_hs_diagonals(permute(A, row_perm, col_perm));
                if (score < best_score) {
                    best_score = score;
                    improved = true;
                } else {
                    std::swap(row_perm[i], row_perm[j]); // undo
                }
            }
        }

        // Try col swaps
        for (size_t i = 0; i < n; i++) {
            for (size_t j = i+1; j < n; j++) {
                std::swap(col_perm[i], col_perm[j]);
                size_t score = count_nonempty_hs_diagonals(permute(A, row_perm, col_perm));
                if (score < best_score) {
                    best_score = score;
                    improved = true;
                } else {
                    std::swap(col_perm[i], col_perm[j]); // undo
                }
            }
        }
    }

    std::cout << "HSOrder complete. Active diagonals: " << best_score << "/" << n << std::endl;
}
