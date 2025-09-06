#pragma once
#include <queue>
#include <unordered_map>
#include <utility>
#include <atomic>

using namespace std;

#include "csr.hpp"
#include "utils.hpp"
#include <climits>
#include <assert.h>


/**
 * @brief Updates a given permutation to take into effect moving a row/col to the
 * given index and displacing all rows/cols in between
 * 
 * @param perm  A vector of permutation (old->new)
 * @param row_to_move   index of the row/col to move
 * @param new_pos   index of the destination
 */
void update_perm_by_move(std::vector<size_t>& perm, size_t row_to_move, size_t new_pos) {
    if (row_to_move < new_pos) {
        for (size_t i = 0; i < perm.size(); i++) {
            if (perm[i] == row_to_move) perm[i] = new_pos;
            else if (perm[i] > row_to_move && perm[i] <= new_pos) perm[i]--;
        }
    } else {
        for (size_t i = 0; i < perm.size(); i++) {
            if (perm[i] == row_to_move) perm[i] = new_pos;
            else if (perm[i] < row_to_move && perm[i] >= new_pos) perm[i]++;
        }

    }
}

/**
 * @brief Helper function for search algorithms to update variables accordingly.
 * Moving a row to desired place
 * 
 * @param A Our CSRMatrix object
 * @param current_row_perm  vector of the row permutation (old->new) 
 * @param current_col_perm vector of the col permutation (old->new)
 * @param best_move_row the index of the row to move
 * @param best_move_new_pos the index of the destination
 * @param num_nnz array holding the number of nonzeros in each HS diagonal
 */
void execute_move_row(const CSRMatrix& A, std::vector<size_t>& current_row_perm, std::vector<size_t>& current_col_perm, 
                      size_t best_move_row, size_t best_move_new_pos, unsigned int* num_nnz) {
    std::vector<size_t> inv_row_perm(A.rows);
    for (size_t i = 0; i < A.rows; i++) inv_row_perm[current_row_perm[i]] = i;

    for (size_t k = A.row_ptr[inv_row_perm[best_move_row]]; k < A.row_ptr[inv_row_perm[best_move_row] + 1]; k++) {
        size_t new_col_idx = current_col_perm[A.col_idx[k]];

        size_t new_diag_idx = (new_col_idx + A.rows - best_move_new_pos) % A.rows;
        size_t old_diag_idx = (new_col_idx + A.rows - best_move_row) % A.rows;

        num_nnz[new_diag_idx]++;
        num_nnz[old_diag_idx]--;
    }

    if (best_move_row > best_move_new_pos) {
        for (size_t current_row = best_move_new_pos; current_row < best_move_row; current_row++) {
            size_t old_row_idx = inv_row_perm[current_row];

            for (size_t k = A.row_ptr[old_row_idx]; k < A.row_ptr[old_row_idx+1]; k++) {
                size_t new_col_idx = current_col_perm[A.col_idx[k]];

                size_t old_diag_idx = (new_col_idx + A.rows - current_row) % A.rows;
                size_t new_diag_idx = (old_diag_idx + A.rows - 1) % A.rows;

                num_nnz[new_diag_idx]++;
                num_nnz[old_diag_idx]--;
            }
        }
    } else {
        for (size_t current_row = best_move_new_pos; current_row > best_move_row; current_row--) {
            size_t old_row_idx = inv_row_perm[current_row];

            for (size_t k = A.row_ptr[old_row_idx]; k < A.row_ptr[old_row_idx+1]; k++) {
                size_t new_col_idx = current_col_perm[A.col_idx[k]];

                size_t old_diag_idx = (new_col_idx + A.rows - current_row) % A.rows;
                size_t new_diag_idx = (old_diag_idx + 1) % A.rows;

                num_nnz[new_diag_idx]++;
                num_nnz[old_diag_idx]--;
            }
        }
    }
    
    update_perm_by_move(current_row_perm, best_move_row, best_move_new_pos);
}

/**
 * @brief Helper function for search algorithms to update variables accordingly.
 * Moving a col to desired place.
 * 
 * @param A Our CSRMatrix object
 * @param current_row_perm  vector of the row permutation (old->new) 
 * @param current_col_perm vector of the col permutation (old->new)
 * @param best_move_col the index of the row to move
 * @param best_move_new_pos the index of the destination
 * @param num_nnz array holding the number of nonzeros in each HS diagonal
 */
void execute_move_col(const CSRMatrix& A, std::vector<size_t>& current_row_perm, std::vector<size_t>& current_col_perm, 
                      size_t best_move_col, size_t best_move_new_pos, unsigned int* num_nnz) {
    if (best_move_col > best_move_new_pos) {
        for (size_t old_row = 0; old_row < A.rows; old_row++) {
            size_t new_row = current_row_perm[old_row];

            for (size_t k = A.row_ptr[old_row]; k < A.row_ptr[old_row+1]; k++) {
                size_t new_col_idx = current_col_perm[A.col_idx[k]];

                if (new_col_idx == best_move_col) {
                    size_t old_diag_idx = (new_col_idx + A.rows - new_row) % A.rows;
                    size_t new_diag_idx = (best_move_new_pos + A.rows - new_row) % A.rows;

                    num_nnz[old_diag_idx]--;
                    num_nnz[new_diag_idx]++;
                } else if (new_col_idx < best_move_col && new_col_idx >= best_move_new_pos) {
                    size_t old_diag_idx = (new_col_idx + A.rows - new_row) % A.rows;
                    size_t new_diag_idx = (old_diag_idx + 1) % A.rows;

                    num_nnz[old_diag_idx]--;
                    num_nnz[new_diag_idx]++;
                }
            }
        }
    } else {
        for (size_t old_row = 0; old_row < A.rows; old_row++) {
            size_t new_row = current_row_perm[old_row];

            for (size_t k = A.row_ptr[old_row]; k < A.row_ptr[old_row+1]; k++) {
                size_t new_col_idx = current_col_perm[A.col_idx[k]];

                if (new_col_idx == best_move_col) {
                    size_t old_diag_idx = (new_col_idx + A.rows - new_row) % A.rows;
                    size_t new_diag_idx = (best_move_new_pos + A.rows - new_row) % A.rows;

                    num_nnz[old_diag_idx]--;
                    num_nnz[new_diag_idx]++;
                } else if (new_col_idx > best_move_col && new_col_idx <= best_move_new_pos) {
                    size_t old_diag_idx = (new_col_idx + A.rows - new_row) % A.rows;
                    size_t new_diag_idx = (old_diag_idx + A.rows - 1) % A.rows;
                    
                    num_nnz[old_diag_idx]--;
                    num_nnz[new_diag_idx]++;
                }
            }
        }
    }
    
    update_perm_by_move(current_col_perm, best_move_col, best_move_new_pos);
}

/**
 * @brief Helper function for search algorithms to update variables accordingly.
 * Swap two rows
 * 
 * @param A Our CSRMatrix object
 * @param current_row_perm  vector of the row permutation (old->new) 
 * @param current_col_perm vector of the col permutation (old->new)
 * @param best_move_row0 the index of the row to swap
 * @param best_move_row1 the index of the other row to swap
 * @param num_nnz array holding the number of nonzeros in each HS diagonal
 */
void execute_row_swap(const CSRMatrix& A, std::vector<size_t>& current_row_perm, std::vector<size_t>& current_col_perm, 
                      size_t best_move_row0, size_t best_move_row1, unsigned int* num_nnz) {
    size_t idx_0, idx_1;
    for (size_t i = 0; i < A.rows; i++) {
        if (current_row_perm[i] == best_move_row0) idx_0 = i;
        if (current_row_perm[i] == best_move_row1) idx_1 = i;
    }
    std::swap(current_row_perm[idx_0], current_row_perm[idx_1]);

    for (size_t k = A.row_ptr[idx_0]; k < A.row_ptr[idx_0 + 1]; k++) {
        num_nnz[(current_col_perm[A.col_idx[k]] + A.rows - best_move_row0) % A.rows]--;
        num_nnz[(current_col_perm[A.col_idx[k]] + A.rows - best_move_row1) % A.rows]++;
    }
    for (size_t k = A.row_ptr[idx_1]; k < A.row_ptr[idx_1 + 1]; k++) {
        num_nnz[(current_col_perm[A.col_idx[k]] + A.rows - best_move_row1) % A.rows]--;
        num_nnz[(current_col_perm[A.col_idx[k]] + A.rows - best_move_row0) % A.rows]++;
    }
}

/**
 * @brief Helper function for search algorithms to update variables accordingly.
 * Swap two cols
 * 
 * @param A Our CSRMatrix object
 * @param current_row_perm  vector of the row permutation (old->new) 
 * @param current_col_perm vector of the col permutation (old->new)
 * @param best_move_col0 the index of the col to swap
 * @param best_move_col1 the index of the other col to swap
 * @param num_nnz array holding the number of nonzeros in each HS diagonal
 */
void execute_col_swap(const CSRMatrix& A, std::vector<size_t>& current_row_perm, std::vector<size_t>& current_col_perm, 
                      size_t best_move_col0, size_t best_move_col1, unsigned int* num_nnz) {
    size_t idx_0, idx_1;
    for (size_t i = 0; i < A.cols; i++) {
        if (current_col_perm[i] == best_move_col0) idx_0 = i;
        if (current_col_perm[i] == best_move_col1) idx_1 = i;
    }
    std::swap(current_col_perm[idx_0], current_col_perm[idx_1]);

    for (size_t row = 0; row < A.rows; row++) {
        size_t new_row0 = current_row_perm[row];
        size_t new_row1 = current_row_perm[row];

        for (size_t k = A.row_ptr[row]; k < A.row_ptr[row+1]; k++) {
            size_t idx = A.col_idx[k];
            if (idx == idx_0) {
                num_nnz[(best_move_col0 + A.rows - new_row0) % A.rows]--;
                num_nnz[(best_move_col1 + A.rows - new_row1) % A.rows]++;
            } else if (idx == idx_1) {
                num_nnz[(best_move_col1 + A.rows - new_row1) % A.rows]--;
                num_nnz[(best_move_col0 + A.rows - new_row0) % A.rows]++;
            }
        }
    }
}

/**
 * @brief Computes the potential gain of exeucting col0 <-> col1
 * 
 * @param A Our CSRMatrix object
 * @param current_row_perm vector of the row permutation (old->new) 
 * @param current_col_perm vector of the col permutation (old->new)
 * @param num_nnz array holding the number of nonzeros in each HS diagonal
 * @param col0 the index of the col to swap
 * @param col1 the index of the other col to swap
 * @return int - The potential gain
 */
int compute_gain_col_swap(const CSRMatrix& A, const std::vector<size_t>& current_row_perm, const std::vector<size_t>& current_col_perm, 
                          unsigned int* num_nnz, size_t col0, size_t col1) {
    int leave_gains = 0;
    int arrival_losses = 0;

    size_t idx_0, idx_1;
    for (size_t i = 0; i < A.cols; i++) {
        if (current_col_perm[i] == col0) idx_0 = i;
        if (current_col_perm[i] == col1) idx_1 = i;
    }
    
    int num_nnz_copy[A.rows];
    std::copy(num_nnz, num_nnz + A.rows, num_nnz_copy);

    for (size_t row = 0; row < A.rows; row++) {
        size_t new_row0 = current_row_perm[row];
        size_t new_row1 = current_row_perm[row];

        for (size_t k = A.row_ptr[row]; k < A.row_ptr[row+1]; k++) {
            size_t idx = A.col_idx[k];
            if (idx == idx_0) {
                if (num_nnz_copy[(col0 + A.rows - new_row0) % A.rows] == 1) leave_gains++;
                num_nnz_copy[(col0 + A.rows - new_row0) % A.rows]--;
                if (num_nnz_copy[(col1 + A.rows - new_row1) % A.rows] == 0) arrival_losses++;
                num_nnz_copy[(col1 + A.rows - new_row1) % A.rows]++;
            } else if (idx == idx_1) {
                if (num_nnz_copy[(col1 + A.rows - new_row1) % A.rows] == 1) leave_gains++;
                num_nnz_copy[(col1 + A.rows - new_row1) % A.rows]--;
                if (num_nnz_copy[(col0 + A.rows - new_row0) % A.rows] == 0) arrival_losses++;
                num_nnz_copy[(col0 + A.rows - new_row0) % A.rows]++;
            }
        }
    }

    return leave_gains - arrival_losses;
}

/**
 * @brief Computes the potential gain of exeucting row0 <-> row1
 * 
 * @param A Our CSRMatrix object
 * @param current_row_perm vector of the row permutation (old->new) 
 * @param current_col_perm vector of the col permutation (old->new)
 * @param num_nnz array holding the number of nonzeros in each HS diagonal
 * @param row0 the index of the row to swap
 * @param row1 the index of the other row to swap
 * @return int - The potential gain
 */
int compute_gain_row_swap(const CSRMatrix& A, const std::vector<size_t>& current_row_perm, const std::vector<size_t>& current_col_perm, 
                          unsigned int* num_nnz, size_t row0, size_t row1) {
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
        if (num_nnz_copy[(current_col_perm[A.col_idx[k]] + A.rows - row0) % A.rows] == 1) leave_gains++;
        num_nnz_copy[(current_col_perm[A.col_idx[k]] + A.rows - row0) % A.rows]--;
        if (num_nnz_copy[(current_col_perm[A.col_idx[k]] + A.rows - row1) % A.rows] == 0) arrival_losses++;
        num_nnz_copy[(current_col_perm[A.col_idx[k]] + A.rows - row1) % A.rows]++;
    }

    for (size_t k = A.row_ptr[idx_1]; k < A.row_ptr[idx_1+1]; k++) {
        if (num_nnz_copy[(current_col_perm[A.col_idx[k]] + A.rows - row1) % A.rows] == 1) leave_gains++;
        num_nnz_copy[(current_col_perm[A.col_idx[k]] + A.rows - row1) % A.rows]--;
        if (num_nnz_copy[(current_col_perm[A.col_idx[k]] + A.rows - row0) % A.rows] == 0) arrival_losses++;
        num_nnz_copy[(current_col_perm[A.col_idx[k]] + A.rows - row0) % A.rows]++;
    }

    return leave_gains - arrival_losses;
}

/**
 * @brief Computes the potential gain of exeucting col -> new_pos
 * 
 * @param A Our CSRMatrix object
 * @param current_row_perm vector of the row permutation (old->new) 
 * @param current_col_perm vector of the col permutation (old->new)
 * @param num_nnz array holding the number of nonzeros in each HS diagonal
 * @param col the index of the col to move
 * @param new_pos the index of the destination
 * @return int - The potential gain
 */
int compute_gain_displace_col(const CSRMatrix& A, const std::vector<size_t>& current_row_perm, const std::vector<size_t>& current_col_perm, 
                               unsigned int* num_nnz, size_t col, size_t new_pos) {
    int gain = 0;
    int loss = 0;

    int num_nnz_copy[A.rows];
    std::copy(num_nnz, num_nnz + A.rows, num_nnz_copy);

    if (col > new_pos) {
        for (size_t old_row = 0; old_row < A.rows; old_row++) {
            size_t new_row = current_row_perm[old_row];

            for (size_t k = A.row_ptr[old_row]; k < A.row_ptr[old_row+1]; k++) {
                size_t new_col_idx = current_col_perm[A.col_idx[k]];

                if (new_col_idx == col) {
                    size_t old_diag_idx = (new_col_idx + A.rows - new_row) % A.rows;
                    size_t new_diag_idx = (new_pos + A.rows - new_row) % A.rows;

                    num_nnz_copy[old_diag_idx]--;
                    num_nnz_copy[new_diag_idx]++;

                    if (num_nnz_copy[old_diag_idx] == 0) gain++;
                    if (num_nnz_copy[new_diag_idx] == 1) loss++;
                } else if (new_col_idx < col && new_col_idx >= new_pos) {
                    size_t old_diag_idx = (new_col_idx + A.rows - new_row) % A.rows;
                    size_t new_diag_idx = (old_diag_idx + 1) % A.rows;

                    num_nnz_copy[old_diag_idx]--;
                    num_nnz_copy[new_diag_idx]++;

                    if (num_nnz_copy[old_diag_idx] == 0) gain++;
                    if (num_nnz_copy[new_diag_idx] == 1) loss++;
                }
            }
        }
    } else {
        for (size_t old_row = 0; old_row < A.rows; old_row++) {
            size_t new_row = current_row_perm[old_row];

            for (size_t k = A.row_ptr[old_row]; k < A.row_ptr[old_row+1]; k++) {
                size_t new_col_idx = current_col_perm[A.col_idx[k]];

                if (new_col_idx == col) {
                    size_t old_diag_idx = (new_col_idx + A.rows - new_row) % A.rows;
                    size_t new_diag_idx = (new_pos + A.rows - new_row) % A.rows;

                    num_nnz_copy[old_diag_idx]--;
                    num_nnz_copy[new_diag_idx]++;

                    if (num_nnz_copy[old_diag_idx] == 0) gain++;
                    if (num_nnz_copy[new_diag_idx] == 1) loss++;
                } else if (new_col_idx > col && new_col_idx <= new_pos) {
                    size_t old_diag_idx = (new_col_idx + A.rows - new_row) % A.rows;
                    size_t new_diag_idx = (old_diag_idx + A.rows - 1) % A.rows;
                    
                    num_nnz_copy[old_diag_idx]--;
                    num_nnz_copy[new_diag_idx]++;

                    if (num_nnz_copy[old_diag_idx] == 0) gain++;
                    if (num_nnz_copy[new_diag_idx] == 1) loss++;
                }
            }
        }
    }
    return gain - loss;
}

/**
 * @brief Computes the potential gain of exeucting row -> new_pos
 * 
 * @param A Our CSRMatrix object
 * @param current_row_perm vector of the row permutation (old->new) 
 * @param current_col_perm vector of the col permutation (old->new)
 * @param num_nnz array holding the number of nonzeros in each HS diagonal
 * @param row the index of the row to move
 * @param new_pos the index of the destination
 * @return int - The potential gain
 */
int compute_gain_displace_row(const CSRMatrix& A, const std::vector<size_t>& current_row_perm, const std::vector<size_t>& current_col_perm, 
                               unsigned int* num_nnz, size_t row, size_t new_pos) {
    int gain = 0;
    int loss = 0;

    std::vector<size_t> inv_row_perm(A.rows);
    for (size_t i = 0; i < A.rows; i++) inv_row_perm[current_row_perm[i]] = i;

    int num_nnz_copy[A.rows];
    std::copy(num_nnz, num_nnz + A.rows, num_nnz_copy);
    
    for (size_t k = A.row_ptr[inv_row_perm[row]]; k < A.row_ptr[inv_row_perm[row] + 1]; k++) {
        size_t new_col_idx = current_col_perm[A.col_idx[k]];

        size_t new_diag_idx = (new_col_idx + A.rows - new_pos) % A.rows;
        size_t old_diag_idx = (new_col_idx + A.rows - row) % A.rows;

        num_nnz_copy[new_diag_idx]++;
        num_nnz_copy[old_diag_idx]--;

        if (num_nnz_copy[new_diag_idx] == 1) loss++;
        if (num_nnz_copy[old_diag_idx] == 0) gain++;
    }

    if (row > new_pos) {
        for (size_t current_row = new_pos; current_row < row; current_row++) {
            size_t old_row_idx = inv_row_perm[current_row];

            for (size_t k = A.row_ptr[old_row_idx]; k < A.row_ptr[old_row_idx+1]; k++) {
                size_t new_col_idx = current_col_perm[A.col_idx[k]];

                size_t old_diag_idx = (new_col_idx + A.rows - current_row) % A.rows;
                size_t new_diag_idx = (old_diag_idx + A.rows - 1) % A.rows;

                num_nnz_copy[new_diag_idx]++;
                num_nnz_copy[old_diag_idx]--;

                if (num_nnz_copy[new_diag_idx] == 1) loss++;
                if (num_nnz_copy[old_diag_idx] == 0) gain++;
            }
        }
    } else {
        for (size_t current_row = new_pos; current_row > row; current_row--) {
            size_t old_row_idx = inv_row_perm[current_row];

            for (size_t k = A.row_ptr[old_row_idx]; k < A.row_ptr[old_row_idx+1]; k++) {
                size_t new_col_idx = current_col_perm[A.col_idx[k]];

                size_t old_diag_idx = (new_col_idx + A.rows - current_row) % A.rows;
                size_t new_diag_idx = (old_diag_idx + 1) % A.rows;

                num_nnz_copy[new_diag_idx]++;
                num_nnz_copy[old_diag_idx]--;

                if (num_nnz_copy[new_diag_idx] == 1) loss++;
                if (num_nnz_copy[old_diag_idx] == 0) gain++;
            }
        }
    }

    return gain - loss;
}

/**
 * @brief Searches for the best possible displace move (row -> pos) or (col -> pos) and executes it on the current permutation of A
 * 
 * @param A Our CSRMatrix object
 * @param current_row_perm vector of the row permutation (old->new) 
 * @param current_col_perm vector of the col permutation (old->new)
 * @param num_nnz array holding the number of nonzeros in each HS diagonal
 * @param[out] best_move_gain the gain of the best move which is executed
 * @param current_diagonal_count the diag count of the current permutation
 * @param steps_without_gain used to backtrack after certain number of moves
 * @param row_turn bool to denote turn
 */
void find_and_execute_displace_move(const CSRMatrix& A, std::vector<size_t>& current_row_perm, std::vector<size_t>& current_col_perm,
                                unsigned int* num_nnz, int& best_move_gain, unsigned int& current_diagonal_count, 
                                unsigned int& steps_without_gain, bool &row_turn) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> row_dist(0, A.rows - 1);
    std::uniform_int_distribution<> col_dist(0, A.cols - 1);

    size_t best_move_col, best_move_row, best_move_new_pos;

    const int max_r_iter = std::min(static_cast<size_t>(500), std::max(A.rows, A.cols));
    if (row_turn) {
        for (size_t r_iter = 0; r_iter < max_r_iter; r_iter++) {
            size_t row = row_dist(gen);
            size_t new_pos = row_dist(gen);

            while (new_pos == row) new_pos = row_dist(gen);

            int move_gain = compute_gain_displace_row(A, current_row_perm, current_col_perm, num_nnz, row, new_pos);

            // only allow negative gain moves if we feel stuck
            if (r_iter < max_r_iter / 2 && move_gain < 0) continue;

            if (move_gain > best_move_gain) {
                best_move_gain = move_gain;
                best_move_row = row;
                best_move_new_pos = new_pos;
            }
        }
    } else {
        for (size_t r_iter = 0; r_iter < max_r_iter; r_iter++) {
            size_t col = col_dist(gen);
            size_t new_pos = col_dist(gen);

            while (col == new_pos) new_pos = col_dist(gen);

            int move_gain = compute_gain_displace_col(A, current_row_perm, current_col_perm, num_nnz, col, new_pos);

            // only allow negative gain moves if we feel stuck
            if (r_iter < max_r_iter / 2 && move_gain < 0) continue;

            if (move_gain > best_move_gain) {
                best_move_gain = move_gain;
                best_move_col = col;
                best_move_new_pos = new_pos;
            }
        }
    }

    if (best_move_gain == INT_MIN) {
        std::cout << "No beneficial moves found." << std::endl;
        return;
    } else {
        // if (row_turn) std::cout << "Best move gain: " << best_move_gain << " (Row " << best_move_row << " -> Row " << best_move_new_pos << ")" << std::endl;
        // else std::cout << "Best move gain: " << best_move_gain << " (Col " << best_move_col << " -> Col " << best_move_new_pos << ")" << std::endl;

        if (best_move_gain > 0) steps_without_gain = 0;
        else steps_without_gain++;

        current_diagonal_count -= best_move_gain;

        if (row_turn) {
            execute_move_row(A, current_row_perm, current_col_perm, best_move_row, best_move_new_pos, num_nnz);
        } else {
            execute_move_col(A, current_row_perm, current_col_perm, best_move_col, best_move_new_pos, num_nnz);
        }
    }
    
}

/**
 * @brief Searches for the best possible swap move (row0 <-> row1) or (col0 <-> col1) and executes it on the current permutation of A
 * 
 * @param A Our CSRMatrix object
 * @param current_row_perm vector of the row permutation (old->new) 
 * @param current_col_perm vector of the col permutation (old->new)
 * @param num_nnz array holding the number of nonzeros in each HS diagonal
 * @param[out] best_move_gain the gain of the best move which is executed
 * @param current_diagonal_count the diag count of the current permutation
 * @param steps_without_gain used to backtrack after certain number of moves
 * @param row_turn bool to denote turn
 */
void find_and_execute_swap_move(const CSRMatrix& A, std::vector<size_t>& current_row_perm, std::vector<size_t>& current_col_perm,
                                unsigned int* num_nnz, int& best_move_gain, unsigned int& current_diagonal_count, 
                                unsigned int& steps_without_gain, bool &row_turn) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> row_dist(0, A.rows - 1);
    std::uniform_int_distribution<> col_dist(0, A.cols - 1);
    size_t best_move_row0, best_move_row1, best_move_col0, best_move_col1;

    //TODO: This should be adaptive
    const int max_r_iter = std::min(static_cast<size_t>(500), A.values.size());
    if (row_turn) {
        for (size_t r_iter = 0; r_iter < max_r_iter; r_iter++) {
            // select 2 rows at random
            size_t row0 = row_dist(gen);
            size_t row1 = row_dist(gen);
            while (row1 == row0) row1 = row_dist(gen);

            int move_gain = compute_gain_row_swap(A, current_row_perm, current_col_perm, num_nnz, row0, row1);

            // only allow negative gain moves if we feel stuck
            if (r_iter < max_r_iter / 2 && move_gain < 0) continue;

            if (move_gain > best_move_gain) {
                best_move_gain = move_gain;
                best_move_row0 = row0;
                best_move_row1 = row1;
            }
        }
    } else {
        for (size_t r_iter = 0; r_iter < max_r_iter; r_iter++) {
            size_t col0 = col_dist(gen);
            size_t col1 = col_dist(gen);

            while (col1 == col0) col1 = col_dist(gen);

            int move_gain = compute_gain_col_swap(A, current_row_perm, current_col_perm, num_nnz, col0, col1);

            // only allow negative gain moves if we feel stuck
            if (r_iter < max_r_iter / 2 && move_gain < 0) continue;

            if (move_gain > best_move_gain) {
                best_move_gain = move_gain;
                best_move_col0 = col0;
                best_move_col1 = col1;
            }
        }
    }
    if (best_move_gain == INT_MIN) {
        std::cout << "No beneficial moves found." << std::endl;
        return;
    } else {
        // if (row_turn) std::cout << "Best move gain: " << best_move_gain << " (Row " << best_move_row0 << " <-> Row " << best_move_row1 << ")" << std::endl;
        // else std::cout << "Best move gain: " << best_move_gain << " (Col " << best_move_col0 << " <-> Col " << best_move_col1 << ")" << std::endl;
        
        if (best_move_gain > 0) steps_without_gain = 0;
        else steps_without_gain++;

        current_diagonal_count -= best_move_gain;

        if (row_turn) {
            execute_row_swap(A, current_row_perm, current_col_perm, best_move_row0, best_move_row1, num_nnz);
        } else {
            execute_col_swap(A, current_row_perm, current_col_perm, best_move_col0, best_move_col1, num_nnz);
        }
    }
}

/**
 * @brief Perform a search over the permutation space of A to place nonzeros on HS diagonals.
 * 
 * @param A CSRMatrix object
 * @param[out] row_perm result row permutation
 * @param[out] col_perm result col permutation
 * 
 * Alternatingly swap or displace rows and columns to fit nonzeros on HS diagonals. Search algorithm is greedy to 
 * always choose the best possible move. If it gets stuck in a local minimum, it allows negative moves and backtracks after
 * @p max_steps_without_gain moves are made.
 */
void hsorder(const CSRMatrix& A, std::vector<size_t>& row_perm, std::vector<size_t>& col_perm) {
    // number of nonzeros in each HS diagonal
    unsigned int num_nnz[A.rows] = {0};
    unsigned int best_num_nnz[A.rows] = {0};

    // initialize nonzero counts
    for (size_t row = 0; row < A.rows; row++) {
        for (size_t k = A.row_ptr[row]; k < A.row_ptr[row + 1]; k++) {
            num_nnz[(A.col_idx[k] + A.rows - row) % A.rows]++;
            best_num_nnz[(A.col_idx[k] + A.rows - row) % A.rows]++;
        }
    }

    std::vector<size_t> current_row_perm(A.rows);
    std::vector<size_t> current_col_perm(A.cols);
    for (size_t i = 0; i < A.rows; i++) current_row_perm[i] = i;
    for (size_t i = 0; i < A.cols; i++) current_col_perm[i] = i;

    const unsigned int max_steps_without_improvement = 100000 * ((A.rows + 255) / 256);
    const unsigned int max_steps_without_gain = max_steps_without_improvement / 1000;

    unsigned int steps_without_gain = 0;
    unsigned int steps_without_improvement = 0;

    std::vector<size_t> best_row_perm = current_row_perm;
    std::vector<size_t> best_col_perm = current_col_perm;
    unsigned int best_diagonal_count = count_nonempty_hs_diagonals(A);
    unsigned int current_diagonal_count = best_diagonal_count;

    // 0 -> swap
    // 1 -> move
    int move_type = 0;
    bool row_turn = true;

    std::mt19937 gen(std::random_device{}());          
    std::uniform_int_distribution<int> dist01(0, 1);
    while (steps_without_improvement < max_steps_without_improvement) {
        int best_move_gain = INT_MIN;
        
        if (move_type == 0) find_and_execute_swap_move(A, current_row_perm, current_col_perm, num_nnz, best_move_gain, current_diagonal_count, steps_without_gain, row_turn);
        else find_and_execute_displace_move(A, current_row_perm, current_col_perm, num_nnz, best_move_gain, current_diagonal_count, steps_without_gain, row_turn);
        
        if (best_move_gain == INT_MIN) {
            std::cout << "No beneficial moves found after " << move_type << ". Stopping..." << std::endl;
            break;
        }
        
        if (current_diagonal_count < best_diagonal_count) {
            best_diagonal_count = current_diagonal_count;
            std::copy(current_row_perm.begin(), current_row_perm.end(), best_row_perm.begin());
            std::copy(current_col_perm.begin(), current_col_perm.end(), best_col_perm.begin());
            std::copy(num_nnz, num_nnz + A.rows, best_num_nnz);
            steps_without_improvement = 0;
        } else if (steps_without_gain > max_steps_without_gain) {
            std::cout << "Backtracking..." << std::endl;
            std::copy(best_row_perm.begin(), best_row_perm.end(), current_row_perm.begin());
            std::copy(best_col_perm.begin(), best_col_perm.end(), current_col_perm.begin());
            std::copy(best_num_nnz, best_num_nnz + A.rows, num_nnz);
            steps_without_gain = 0;
        } else steps_without_improvement++;

        // // print debug info
        // std::cout << "Current diagonal count: " << current_diagonal_count << std::endl;
        // std::cout << "Best diagonal count: " << best_diagonal_count << std::endl;
        // std::cout << "Current row perm: ";
        // for (size_t i = 0; i < A.rows; i++) {
        //     std::cout << current_row_perm[i] << " ";
        // }
        // std::cout << "Current col perm: ";
        // for (size_t i = 0; i < A.cols; i++) {
        //     std::cout << current_col_perm[i] << " ";
        // }
        // std::cout << std::endl;
        // std::cout << "nnz array: ";
        // for (size_t i = 0; i < A.rows; i++) {
        //     std::cout << num_nnz[i] << " ";
        // }
        // std::cout << std::endl;

        move_type = dist01(gen);        // random 0 or 1
        row_turn  = static_cast<bool>(dist01(gen)); // random true or false

        #ifdef _WIN32
            if (_kbhit()) {
                stop_requested.store(true);
                break;
            }
        #else
            if (kbhit()) {
                stop_requested.store(true);
                break;
            }
        #endif
        
    }
    std::cout << "HSOrder complete. Active diagonals: " << best_diagonal_count << "/" << A.rows << std::endl;
    
    std::copy(best_row_perm.begin(), best_row_perm.end(), row_perm.begin());
    std::copy(best_col_perm.begin(), best_col_perm.end(), col_perm.begin());
}

/**
 * @brief Exhaustive search algorithm for fitting nonzeros on HS diagonals.
 * 
 * @param A CSRMatrix object
 * @param[out] row_perm result row permutation
 * @param[out] col_perm result col permutation
 * @param window_size 
 * @param debug_mode print debug info if set to true
 * 
 * At each step, pick a random row/col, do the best move available for that choice (swap or move) and execute it. Then move to
 * the next row/col. Rows/cols already swapped or displaced are not touched again. When all rows and cols are visited, this is 1 pass.
 * Do passes over the matrix until we are stuck. Row visit orders in each pass are random, column visit order is from least
 * populated to the most dense.
 * 
 * Returns the current best permutation if a key is pressed abruptly.
 */
void hsorder_long(const CSRMatrix& A, std::vector<size_t>& row_perm, std::vector<size_t>& col_perm, size_t window_size, bool debug_mode) {
    std::random_device rd;   // seed
    std::mt19937 g(rd()); 

    window_size = std::min(A.rows, window_size);

    const unsigned int initial_diag_count = count_nonempty_hs_diagonals(A);
    unsigned int current_diag_count = initial_diag_count;
    unsigned int best_diag_count = current_diag_count;

    std::vector<size_t> current_row_perm(A.rows), current_col_perm(A.cols);
    unsigned int num_nnz[A.rows] = {0};
    unsigned int best_num_nnz[A.rows] = {0};
    for (size_t row = 0; row < A.rows; row++) {
        for (size_t k = A.row_ptr[row]; k < A.row_ptr[row + 1]; k++) {
            num_nnz[(A.col_idx[k] + A.rows - row) % A.rows]++;
            best_num_nnz[(A.col_idx[k] + A.rows - row) % A.rows]++;
        }
    }
    for (size_t i = 0; i < A.rows; i++) {
        current_col_perm[i] = i;
        current_row_perm[i] = i;
    }

    bool allow_negative = false;
    
    size_t passes_without_improvement = 0;
    const size_t max_passes_without_improvement = 1000 * ((A.rows + 255 )/ 256);

    std::vector<size_t> best_row_perm(A.rows), best_col_perm(A.cols);

    for (size_t i = 0; i < A.rows; i++) best_row_perm[i] = i;
    for (size_t i = 0; i < A.cols; i++) best_col_perm[i] = i;

    size_t negative_passes = 0;
    const size_t backtracking_threshold = 100;

    // When a row or col is affected by a move, mark it to not move it again.
    std::vector<bool> is_row_placed(A.rows, false), is_col_placed(A.cols, false);

    std::vector<size_t> random_col_visit_order(A.cols), num_nnz_cols(A.cols, 0);
    for (size_t row = 0; row < A.rows; row++) {
        for (size_t k = A.row_ptr[row]; k < A.row_ptr[row+1]; k++) {
            num_nnz_cols[A.col_idx[k]]++;
        }
    }

    std::iota(random_col_visit_order.begin(),
          random_col_visit_order.end(),
          0);

    // sort indices based on num_nnz_cols values (descending)
    std::sort(random_col_visit_order.begin(),
        random_col_visit_order.end(),
        [&](size_t i, size_t j) {
            return num_nnz_cols[i] < num_nnz_cols[j];
        });

    while (passes_without_improvement < max_passes_without_improvement) {
        for (size_t i = 0; i < A.rows; i++) is_row_placed[i] = false;
        for (size_t i = 0; i < A.cols; i++) is_col_placed[i] = false;

        std::vector<size_t> random_row_visit_order(A.rows);
        int pass_gain = 0;
    
        std::shuffle(random_row_visit_order.begin(), random_row_visit_order.end(), g);
        // std::shuffle(random_col_visit_order.begin(), random_col_visit_order.end(), g);

        for (size_t i = 0; i < A.rows; i++) {
            size_t row0 = random_row_visit_order[i];
            size_t col0 = random_col_visit_order[i];

            if (!is_row_placed[row0]) {
                int best_move_gain = INT_MIN;
                size_t best_move_row1;
                bool is_best_swap = true;
                
                // choose random rows and columns for the window
                std::vector<size_t> trial_window(A.rows);
                for (size_t j = 0; j < A.rows; j++) trial_window[j] = j;
                std::shuffle(trial_window.begin(), trial_window.end(), g);
                trial_window.resize(window_size);
                
                // Find the best row move
                for (size_t w_i = 0; w_i < window_size; w_i++) {
                    if (trial_window[w_i] == row0) continue;
                    int swap_move_gain = INT_MIN, displace_move_gain = INT_MIN;
                    if (!is_row_placed[trial_window[w_i]]) swap_move_gain = compute_gain_row_swap(A, current_row_perm, current_col_perm, num_nnz, row0, trial_window[w_i]);
                    
                    bool allow_displace = true;
                    size_t high = max(row0, trial_window[w_i]), low = min(row0, trial_window[w_i]);
                    for (size_t r = low; r <= high; r++) if (is_row_placed[r]) {
                        allow_displace = false;
                        break;
                    }
                    if (allow_displace) displace_move_gain = compute_gain_displace_row(A, current_row_perm, current_col_perm, num_nnz, row0, trial_window[w_i]);

                    if (swap_move_gain > best_move_gain) {
                        best_move_gain = swap_move_gain;
                        best_move_row1 = trial_window[w_i];
                        is_best_swap = true;
                    } 
                    if (displace_move_gain > best_move_gain) {
                        best_move_gain = displace_move_gain;
                        best_move_row1 = trial_window[w_i];
                        is_best_swap = false;
                    }
                }

                // Execute the best row move
                if (best_move_gain == INT_MIN) continue;
                else if (best_move_gain < 0 && !allow_negative) continue;
                else {
                    current_diag_count -= best_move_gain;
                    pass_gain += best_move_gain;
                    if (is_best_swap) {
                        execute_row_swap(A, current_row_perm, current_col_perm, row0, best_move_row1, num_nnz);
                        is_row_placed[row0] = true;
                        is_row_placed[best_move_row1] = true;
                    } else {
                        execute_move_row(A, current_row_perm, current_col_perm, row0, best_move_row1, num_nnz);
                        size_t high = max(row0, best_move_row1), low = min(row0, best_move_row1);
                        for (size_t r = low; r < high; r++) is_row_placed[r] = true;
                    }

                    if (current_diag_count < best_diag_count) {
                        best_diag_count = current_diag_count;
                        std::copy(current_row_perm.begin(), current_row_perm.end(), best_row_perm.begin());
                        std::copy(current_col_perm.begin(), current_col_perm.end(), best_col_perm.begin());
                        std::copy(num_nnz, num_nnz + A.rows, best_num_nnz);
                        passes_without_improvement = 0;
                    }
                }
            }

            if (!is_col_placed[col0]) {
                int best_move_gain = INT_MIN;
                size_t best_move_col1;
                bool is_best_swap = true;

                std::vector<size_t> trial_window(A.cols);
                for (size_t j = 0; j < A.cols; j++) trial_window[j] = j;
                std::shuffle(trial_window.begin(), trial_window.end(), g);
                trial_window.resize(window_size);
                
                // find the best col move
                for (size_t w_i = 0; w_i < window_size; w_i++) {
                    if (trial_window[w_i] == col0) continue;
                    int swap_move_gain = INT_MIN, displace_move_gain = INT_MIN;
                    if (!is_col_placed[trial_window[w_i]]) swap_move_gain = compute_gain_col_swap(A, current_row_perm, current_col_perm, num_nnz, col0, trial_window[w_i]);
                    
                    bool allow_displace = true;
                    size_t high = max(col0, trial_window[w_i]), low = min(col0, trial_window[w_i]);
                    for (size_t c = low; c < high; c++) if (is_col_placed[c]) {
                        allow_displace = false;
                        break;
                    }
                    if (allow_displace) displace_move_gain = compute_gain_displace_col(A, current_row_perm, current_col_perm, num_nnz, col0, trial_window[w_i]);

                    if (swap_move_gain > best_move_gain) {
                        best_move_gain = swap_move_gain;
                        best_move_col1 = trial_window[w_i];
                        is_best_swap = true;
                    } 
                    if (displace_move_gain > best_move_gain) {
                        best_move_gain = displace_move_gain;
                        best_move_col1 = trial_window[w_i];
                        is_best_swap = false;
                    }
                }

                // execute the best col move
                if (best_move_gain == INT_MIN) continue;
                else if (best_move_gain < 0 && !allow_negative) continue;
                else {
                    current_diag_count -= best_move_gain;
                    pass_gain += best_move_gain;
                    if (is_best_swap) {
                        execute_col_swap(A, current_row_perm, current_col_perm, col0, best_move_col1, num_nnz);
                        is_col_placed[col0] = true;
                        is_col_placed[best_move_col1] = true;
                    } else {
                        execute_move_col(A, current_row_perm, current_col_perm, col0, best_move_col1, num_nnz);
                        size_t low = min(col0, best_move_col1), high = max(col0, best_move_col1);
                        for (size_t c = low; c < high; c++) is_col_placed[c] = true;
                    }

                    if (current_diag_count < best_diag_count) {
                        best_diag_count = current_diag_count;
                        std::copy(current_row_perm.begin(), current_row_perm.end(), best_row_perm.begin());
                        std::copy(current_col_perm.begin(), current_col_perm.end(), best_col_perm.begin());
                        std::copy(num_nnz, num_nnz + A.rows, best_num_nnz);
                        passes_without_improvement = 0;
                    }
                }
            }
        }

        if (passes_without_improvement >= 5) allow_negative = true;

        if (pass_gain <= 0) {
            negative_passes++;
            passes_without_improvement++;
        }
        else {
            allow_negative = false;
            negative_passes = 0;
        }

        if (negative_passes >= backtracking_threshold) {
            std::cout << "Backtracking...\n";
            std::copy(best_row_perm.begin(), best_row_perm.end(), current_row_perm.begin());
            std::copy(best_col_perm.begin(), best_col_perm.end(), current_col_perm.begin());
            std::copy(best_num_nnz, best_num_nnz + A.rows, num_nnz);
            negative_passes = 0;
        }

        if (debug_mode && passes_without_improvement % 50 == 0) {
            std::cout << "Finished a pass with gain: " << pass_gain << " current diag count: " << current_diag_count << " best diag count: " << best_diag_count;
            std::cout << " passes without impr: " << passes_without_improvement << " allow negative: " << allow_negative << " negative passes: " << negative_passes << std::endl;
        }

        #ifdef _WIN32
            if (_kbhit()) {
                stop_requested.store(true);
                break;
            }
        #else
            if (kbhit()) {
                stop_requested.store(true);
                break;
            }
        #endif
        
    }
    
    std::copy(best_row_perm.begin(), best_row_perm.end(), row_perm.begin());
    std::copy(best_col_perm.begin(), best_col_perm.end(), col_perm.begin());
}

void independent_set_order(const CSRMatrix& A, std::vector<size_t>& row_perm, std::vector<size_t>& col_perm, bool debug_mode) {    
    std::vector<size_t> current_row_perm(A.rows), current_col_perm(A.cols);
    for (size_t i = 0; i < A.rows; i++) current_row_perm[i] = i;
    for (size_t i = 0; i < A.rows; i++) current_col_perm[i] = i;

    std::vector<size_t> best_row_perm(current_row_perm), best_col_perm(current_col_perm);
    size_t best_diagonal_count = INT_MAX;

    // Sort by nnz descending
    sort(current_row_perm.begin(), current_row_perm.end(),
        [&A](size_t a, size_t b) {
            size_t nnz_a = A.row_ptr[a+1] - A.row_ptr[a];
            size_t nnz_b = A.row_ptr[b+1] - A.row_ptr[b];
            return nnz_a < nnz_b;
        });

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> dist(0, A.rows - 1);
    
    size_t trial_num = 0;
    const size_t max_trials_without_impr = 1000000 * ((A.rows + 255) / 256);

    while (trial_num < max_trials_without_impr) {
        std::shuffle(current_row_perm.begin() + 1, current_row_perm.end(), gen);
        std::set<size_t> occupied_diagonal_indices, fixed_cols;

        std::vector<size_t> inv_row_perm(A.rows);
        for (size_t i = 0; i < A.rows; i++) inv_row_perm[current_row_perm[i]] = i;
        size_t num_nnz = A.row_ptr[inv_row_perm[0]+1] - A.row_ptr[inv_row_perm[0]];

        while (occupied_diagonal_indices.size() < num_nnz) {
            occupied_diagonal_indices.insert(dist(gen));
        }

        bool made_a_move = true;
        size_t current_row = 0;
        for (current_row = 0; current_row < A.rows && made_a_move; current_row++) {
            size_t old_row_index = inv_row_perm[current_row];
            made_a_move = false;
            
            num_nnz = A.row_ptr[old_row_index+1] - A.row_ptr[old_row_index];
            std::vector<size_t> old_col_idxes(num_nnz);

            for (size_t k = A.row_ptr[old_row_index], j = 0; k < A.row_ptr[old_row_index+1]; k++, j++) old_col_idxes[j] = A.col_idx[k]; 
            std::shuffle(old_col_idxes.begin(), old_col_idxes.end(), gen);

            for (size_t k = 0; k < old_col_idxes.size(); k++) {
                size_t col_idx = current_col_perm[old_col_idxes[k]];

                if (fixed_cols.find(col_idx) != fixed_cols.end()) {
                    occupied_diagonal_indices.insert((col_idx + A.rows - current_row) % A.rows);
                } else {
                    bool moved_nnz = false;
                    for (size_t diag_to_place: occupied_diagonal_indices) {
                        size_t new_pos = (current_row + diag_to_place) % A.rows;
                        if (fixed_cols.find(new_pos) != fixed_cols.end()) continue;
                        if (new_pos == col_idx) {
                            fixed_cols.insert(col_idx);
                            made_a_move = true;
                            moved_nnz = true;
                            break;
                        } else {
                            size_t idx0, idx1;
                            for (size_t i = 0; i < A.cols; i++) {
                                if (current_col_perm[i] == col_idx) idx0 = i;
                                else if (current_col_perm[i] == new_pos) idx1 = i;
                            }
                            std::swap(current_col_perm[idx0], current_col_perm[idx1]);
                            fixed_cols.insert(new_pos);
                            made_a_move = true;
                            moved_nnz = true;
                            break;
                        }
                    }
                    if (!moved_nnz) {
                        occupied_diagonal_indices.insert((col_idx + A.rows - current_row) % A.rows);
                    }
                }
            }
        }

        CSRMatrix A_perm = permute(A, current_row_perm, current_col_perm);
        size_t current_diagonal_count = count_nonempty_hs_diagonals(A_perm);
        
        if (current_diagonal_count < best_diagonal_count) {
            std::copy(current_row_perm.begin(), current_row_perm.end(), best_row_perm.begin());
            std::copy(current_col_perm.begin(), current_col_perm.end(), best_col_perm.begin());
            best_diagonal_count = current_diagonal_count;
            trial_num = 0;
        }

        if (debug_mode && trial_num % 1000 == 0) std::cout << "Current best diagonal count: " << best_diagonal_count << std::endl;
        trial_num++;

        #ifdef _WIN32
            if (_kbhit()) {
                stop_requested.store(true);
                break;
            }
        #else
            if (kbhit()) {
                stop_requested.store(true);
                break;
            }
        #endif
    }
    
    std::cout << "ISOrder finished with " << best_diagonal_count << " diagonals\n";

    std::cout << "Row perm: ";
    for (auto p: best_row_perm) std::cout << p << " ";
    std::cout << std::endl << "Col perm: ";
    for (auto p: best_col_perm) std::cout << p << " ";
    std::cout << std::endl;

    std::copy(best_row_perm.begin(), best_row_perm.end(), row_perm.begin());
    std::copy(best_col_perm.begin(), best_col_perm.end(), col_perm.begin());
}