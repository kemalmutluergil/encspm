#include "csr.hpp"
#include "rcm.hpp"
#include "chunk.hpp"
#include "gorder.hpp"
#include "hsorder.hpp"
#include "homomorphic_spmv.hpp"

#include <map>
#include <unordered_set>
#include <assert.h>

#define NUM_TRIALS 5

void apply_GOrder(const CSRMatrix& A, int window_size, std::vector<size_t>& row_perm, std::vector<size_t>& col_perm) {
    CSRMatrix B = build_sym(A);
    CSRMatrix Bt = build_transpose(B);

    auto perm = gorder(B, Bt, window_size);

    for (auto p : perm) {
        row_perm.push_back(p);
        col_perm.push_back(p);
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " [-c] [-d] <input.mtx> <order_id>" << std::endl;
        std::cerr << "   -c    Check computations using CKKS" << std::endl;
        std::cerr << "   -d    Enable debug mode" << std::endl;
        return 1;
    }

    bool check_computations = false;
    bool debug_mode   = false;

    std::string filename;
    std::string order_id;

    // Index of first non-flag argument
    int argi = 1;

    // Parse flags
    while (argi < argc && argv[argi][0] == '-') {
        std::string flag(argv[argi]);
        if (flag == "-c") {
            check_computations = true;
        } else if (flag == "-d") {
            debug_mode = true;
        } else {
            std::cerr << "Unknown flag: " << flag << std::endl;
            return 1;
        }
        argi++;
    }

    // Now we must have two remaining arguments
    if (argc - argi < 2) {
        std::cerr << "Usage: " << argv[0] << " [-c] [-d] <input.mtx> <order_id>" << std::endl;
        std::cerr << "   -c    Check computations using CKKS" << std::endl;
        std::cerr << "   -d    Enable debug mode" << std::endl;
        return 1;
    }

    filename = argv[argi];
    order_id = argv[argi + 1];
    if (atoi(order_id.c_str()) == 2 && argc - argi < 3) {
        std::cerr << "Window size required for GOrder (id: 2)" << std::endl;
        return 1;
    }

    CSRMatrix A = read_mm(filename);
    std::cout << "Initial diagonal count: " << count_nonempty_hs_diagonals(A) << std::endl;

    if (debug_mode) {
        std::cout << "Initial matrix: " << std::endl;
        visualize_csr(A);
    }
    std::cout << "Read Matrix: Rows: " << A.rows << ", Cols: " << A.cols << ", Non-zeros: " << A.values.size() << std::endl;

    CSRMatrix A_perm;
    int order_id_int = atoi(order_id.c_str());

     std::vector<size_t> best_row_perm, best_col_perm;
    // initialize best_row_perm and best_col_perm
    best_row_perm.resize(A.rows);
    best_col_perm.resize(A.cols);
    for (size_t i = 0; i < A.rows; i++) best_row_perm[i] = i;
    for (size_t i = 0; i < A.cols; i++) best_col_perm[i] = i;

    const size_t m = A.rows, n = A.cols;

    if (order_id_int == 0) {
        A_perm = A;

    } if(order_id_int == 1) {   
        rcm::find_rcm_ordering(A, best_row_perm, best_col_perm);
        A_perm = permute_kemal(A, best_row_perm, best_col_perm);
    } else if (order_id_int == 2) {
        int window_size = atoi(argv[argc - 1]);

        if(m != n) {
            std::cerr << "Program does not support nonsquare matrices yet...\n";
            return 1;
            // CSRMatrix B = build_bipartite_sym(A);
            // CSRMatrix Bt = build_transpose(B);
            // auto perm = gorder(B, Bt, window_size);
        
            // for (auto p : perm) {
            //     if (p < m) {
            //         row_perm.push_back(p);
            //     } else {
            //         col_perm.push_back(p - m);
            //     }
            // }
        } else { 
            apply_GOrder(A, window_size, best_row_perm, best_col_perm);
            A_perm = permute_kemal(A, best_row_perm, best_col_perm);
        }
    }  else if (order_id_int == 3) {
        std::vector<size_t> rcm_row_perm(A.rows), rcm_col_perm(A.cols), hs_row_perm(A.rows), hs_col_perm(A.cols);
        for (size_t i = 0; i < A.rows; i++) {
            rcm_row_perm[i] = i;
            hs_row_perm[i] = i;
        }
        for (size_t i = 0; i < A.cols; i++) {
            rcm_col_perm[i] = i;
            hs_col_perm[i] = i;
        }
        rcm::find_rcm_ordering(A, rcm_row_perm, rcm_col_perm);
        A_perm = permute_kemal(A, rcm_row_perm, rcm_col_perm);
        visualize_csr(A_perm);
        std::cout << "Diagonal count: " << count_nonempty_hs_diagonals(A_perm) << std::endl;

        hsorder_kemal(A_perm, hs_row_perm, hs_col_perm);
        for (size_t i = 0; i < A.rows; i++) best_row_perm[i] = hs_row_perm[rcm_row_perm[i]];
        for (size_t i = 0; i < A.cols; i++) best_col_perm[i] = hs_col_perm[rcm_col_perm[i]];            
        
        A_perm = permute_kemal(A_perm, hs_row_perm, hs_col_perm);
    }

    std::cout << "Result Matrix: Rows: " << A_perm.rows << ", Cols: " << A_perm.cols << ", Non-zeros: " << A_perm.values.size() << std::endl;

    std::cout << "Diagonal count: " << count_nonempty_hs_diagonals(A_perm) << std::endl;
    if (debug_mode) {
        std::cout << "Permuted matrix: " << std::endl;
        visualize_csr(A_perm);
    }

    if (check_computations) {
        std::cout << "**********************************CKKS************************************\n";
        if (check_multiplications(A, A_perm, best_row_perm, best_col_perm)) {
            std::cout << "Success: The results match!" << std::endl;
        } else {
            std::cout << "Failed: The results do not match!" << std::endl;
        }
    }
    return 0;
}