#pragma once
#include <queue>
#include <unordered_map>
#include <utility>

using namespace std;

#include "csr.hpp"
#include "utils.hpp"
#include <climits>

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

void hsorder_kemal(const CSRMatrix& A, std::vector<size_t>& row_perm, std::vector<size_t>& col_perm) {
    // number of nonzeros in each HS diagonal
    unsigned int num_nnz[A.rows] = {0};
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, A.rows - 1);

    // initialize nonzero counts
    for (size_t row = 0; row < A.rows; row++) {
        for (size_t k = A.row_ptr[row]; k < A.row_ptr[row + 1]; k++) {
            num_nnz[(A.col_idx[k] + A.rows - row) % A.rows]++;
        }
    }

    std::vector<size_t> current_row_perm(A.rows);
    std::vector<size_t> current_col_perm(A.cols);
    for (size_t i = 0; i < A.rows; i++) current_row_perm[i] = i;
    for (size_t i = 0; i < A.cols; i++) current_col_perm[i] = i;

    const unsigned int max_steps_without_gain = 20;
    const unsigned int max_steps_without_improvement = 500;
    
    // TODO: use this for some sort of backtracking
    unsigned int steps_without_gain = 0;
    unsigned int steps_without_improvement = 0;

    std::vector<size_t> best_row_perm = current_row_perm;
    std::vector<size_t> best_col_perm = current_col_perm;
    unsigned int best_diagonal_count = count_nonempty_hs_diagonals(A);
    unsigned int current_diagonal_count = best_diagonal_count;

    while (steps_without_improvement < max_steps_without_improvement) {
        int best_move_gain = INT_MIN;
        int best_move_row0, best_move_row1;

        //TODO: This should be adaptive
        const int max_r_iter = 10;

        for (size_t r_iter = 0; r_iter < max_r_iter; r_iter++) {
            // select 2 rows at random
            int row0 = dist(gen);
            int row1 = dist(gen);
            while (row1 == row0) row1 = dist(gen);

            // LeaveGain_i := 1 if we make an HS diagonal fully vacant when row_i is removed from its place
            // ArrivalLoss_i := 1 if we put a nonzero in a previously vacant HS diagonal by inserting row_i into j's place
            int leave_gains = 0;
            int arrival_losses = 0;

            size_t idx_0, idx_1;
            for (size_t i = 0; i < A.rows; i++) {
                if (current_row_perm[i] == row0) idx_0 = i;
                if (current_row_perm[i] == row1) idx_1 = i;
            }

            int num_nnz_copy[A.rows];
            std::copy(num_nnz, num_nnz + A.rows, num_nnz_copy);

            for (size_t k = A.row_ptr[idx_0]; k < A.row_ptr[idx_0 + 1]; k++) {
                if (num_nnz_copy[(A.col_idx[k] + A.rows - row0) % A.rows] == 1) leave_gains++;
                num_nnz_copy[(A.col_idx[k] + A.rows - row0) % A.rows]--;
                if (num_nnz_copy[(A.col_idx[k] + A.rows - row1) % A.rows] == 0) arrival_losses++;
                num_nnz_copy[(A.col_idx[k] + A.rows - row1) % A.rows]++;
            }

            for (size_t k = A.row_ptr[idx_1]; k < A.row_ptr[idx_1+1]; k++) {
                if (num_nnz_copy[(A.col_idx[k] + A.rows - row1) % A.rows] == 1) leave_gains++;
                num_nnz_copy[(A.col_idx[k] + A.rows - row1) % A.rows]--;
                if (num_nnz_copy[(A.col_idx[k] + A.rows - row0) % A.rows] == 0) arrival_losses++;
                num_nnz_copy[(A.col_idx[k] + A.rows - row0) % A.rows]++;
            }

            int move_gain = leave_gains - arrival_losses;

            // only allow negative gain moves if we feel stuck
            if (r_iter < max_r_iter / 2 && move_gain < 0) continue;

            if (move_gain > best_move_gain) {
                best_move_gain = move_gain;
                best_move_row0 = row0;
                best_move_row1 = row1;
            }
        }
        if (best_move_gain == INT_MIN) {
            std::cout << "No beneficial moves found." << std::endl;
            break;
        } else {
            // std::cout << "Best move gain: " << best_move_gain << " (Row " << best_move_row0 << " <-> Row " << best_move_row1 << ")" << std::endl;
            
            if (best_move_gain > 0) steps_without_gain = 0;
            else steps_without_gain++;

            current_diagonal_count -= best_move_gain;

            size_t idx_0, idx_1;
            for (size_t i = 0; i < A.rows; i++) {
                if (current_row_perm[i] == best_move_row0) idx_0 = i;
                if (current_row_perm[i] == best_move_row1) idx_1 = i;
            }
            std::swap(current_row_perm[idx_0], current_row_perm[idx_1]);

            if (current_diagonal_count < best_diagonal_count) {
                best_diagonal_count = current_diagonal_count;
                std::copy(current_row_perm.begin(), current_row_perm.end(), best_row_perm.begin());
                std::copy(current_col_perm.begin(), current_col_perm.end(), best_col_perm.begin());
                steps_without_improvement = 0;
            } else steps_without_improvement++;

            for (size_t k = A.row_ptr[idx_0]; k < A.row_ptr[idx_0 + 1]; k++) {
                num_nnz[(A.col_idx[k] + A.rows - best_move_row0) % A.rows]--;
                num_nnz[(A.col_idx[k] + A.rows - best_move_row1) % A.rows]++;
            }
            for (size_t k = A.row_ptr[idx_1]; k < A.row_ptr[idx_1 + 1]; k++) {
                num_nnz[(A.col_idx[k] + A.rows - best_move_row1) % A.rows]--;
                num_nnz[(A.col_idx[k] + A.rows - best_move_row0) % A.rows]++;
            }

            // // print debug info
            // std::cout << "Current diagonal count: " << current_diagonal_count << std::endl;
            // std::cout << "Best diagonal count: " << best_diagonal_count << std::endl;
            // std::cout << "Current row perm: ";
            // for (size_t i = 0; i < A.rows; i++) {
            //     std::cout << current_row_perm[i] << " ";
            // }
            // std::cout << std::endl;
            // std::cout << "nnz array: ";
            // for (size_t i = 0; i < A.rows; i++) {
            //     std::cout << num_nnz[i] << " ";
            // }
            // std::cout << std::endl;
        }
    }

    // copy best row and col perm to the output
    for (size_t i = 0; i < A.rows; i++) {
        row_perm.push_back(best_row_perm[i]);
    }
    for (size_t i = 0; i < A.cols; i++) {
        col_perm.push_back(best_col_perm[i]);
    }
}

CSRMatrix permute_kemal(const CSRMatrix& A, const std::vector<size_t>& row_perm) {
    std::vector<size_t> inv_row_perm(A.rows);
    for (size_t i = 0; i < A.rows; i++) {
        inv_row_perm[row_perm[i]] = i;
    }
    CSRMatrix A_perm;
    A_perm.rows = A.rows;
    A_perm.cols = A.cols;
    A_perm.row_ptr.resize(A.rows + 1);
    A_perm.col_idx.resize(A.col_idx.size());
    A_perm.values.resize(A.col_idx.size());

    size_t num_nnz = 0;
    for (size_t row = 0; row < A.rows; row++) {
        A_perm.row_ptr[row] = num_nnz;
        for (size_t k = A.row_ptr[inv_row_perm[row]]; k < A.row_ptr[inv_row_perm[row] + 1]; k++) {
            A_perm.col_idx[num_nnz] = A.col_idx[k];
            A_perm.values[num_nnz] = A.values[k];
            num_nnz++;
        }
    }
    A_perm.row_ptr[A.rows] = num_nnz;

    return A_perm;
}