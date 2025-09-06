#include "csr.hpp"
#include "rcm.hpp"
#include "chunk.hpp"
#include "gorder.hpp"
#include "hsorder.hpp"
#include "homomorphic_spmv.hpp"

#include <map>
#include <unordered_set>
#include <assert.h>
#include <chrono>

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
        std::cerr << "order_id: " << std::endl;
        std::cerr << "   0: Identity" << std::endl;
        std::cerr << "   1: RCM -> HSOrder" << std::endl;
        std::cerr << "   2: GOrder -> HSOrder" << std::endl;
        std::cerr << "   3: HSOrder" << std::endl;
        std::cerr << "   4: Greedy HSOrder" << std::endl;
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
        // RCM Order

        auto order_start = std::chrono::high_resolution_clock::now();
        rcm::find_rcm_ordering(A, best_row_perm, best_col_perm);     
        auto order_end = std::chrono::high_resolution_clock::now();
        auto order_dur = std::chrono::duration_cast<std::chrono::milliseconds>(order_end - order_start).count();

        std::cout << "RCM Ordering took " << order_dur << " ms" << std::endl;
        A_perm = permute(A, best_row_perm, best_col_perm);
    } else if (order_id_int == 2) {
        // GOrder

        int window_size = atoi(argv[argc - 1]);
        auto order_start = std::chrono::high_resolution_clock::now();
        gorder(A, best_row_perm, best_col_perm, window_size);
        auto order_end = std::chrono::high_resolution_clock::now();
        auto order_dur = std::chrono::duration_cast<std::chrono::milliseconds>(order_end - order_start).count();

        std::cout << "GOrder took " << order_dur << " ms" << std::endl;
        A_perm = permute(A, best_row_perm, best_col_perm);
    }  else if (order_id_int == 3) {
        // HSOrder

        auto order_start = std::chrono::high_resolution_clock::now();
        hsorder_long(A, best_row_perm, best_col_perm, 500, debug_mode);
        auto order_end = std::chrono::high_resolution_clock::now();
        auto order_dur = std::chrono::duration_cast<std::chrono::milliseconds>(order_end - order_start).count();

        std::cout << "HS Ordering took " << order_dur << " ms" << std::endl;
        A_perm = permute(A, best_row_perm, best_col_perm);
    } else if (order_id_int == 4) {
        // RCM then HSOrder

        std::vector<size_t> hs_row_perm(A.rows), rcm_row_perm(A.rows), hs_col_perm(A.cols), rcm_col_perm(A.cols);

        auto rcm_order_start = std::chrono::high_resolution_clock::now();
        rcm::find_rcm_ordering(A, rcm_row_perm, rcm_col_perm);
        auto rcm_order_end = std::chrono::high_resolution_clock::now();
        auto rcm_order_dur = std::chrono::duration_cast<std::chrono::milliseconds>(rcm_order_end - rcm_order_start).count();

        std::cout << "RCM Ordering took " << rcm_order_dur << " ms" << std::endl;
        A_perm = permute(A, rcm_row_perm, rcm_col_perm);
        std::cout << "Intermediate diagonal count after RCM: " << count_nonempty_hs_diagonals(A_perm) << std::endl;

        auto hs_order_start = std::chrono::high_resolution_clock::now();
        hsorder_long(A_perm, hs_row_perm, hs_col_perm, 500, debug_mode);
        auto hs_order_end = std::chrono::high_resolution_clock::now();
        auto hs_order_dur = std::chrono::duration_cast<std::chrono::milliseconds>(hs_order_end - hs_order_start).count();

        std::cout << "HS Ordering took " << hs_order_dur << " ms" << std::endl;
        A_perm = permute(A_perm, hs_row_perm, hs_col_perm);

        for (size_t i = 0; i < A.rows; i++) best_row_perm[i] = hs_row_perm[rcm_row_perm[i]];
        for (size_t i = 0; i < A.cols; i++) best_col_perm[i] = hs_col_perm[rcm_col_perm[i]];
    } else if (order_id_int == 5) {
        // GOrder then HSOrder

        int window_size = atoi(argv[argc - 1]);
        std::vector<size_t> hs_row_perm(A.rows), go_row_perm(A.rows), hs_col_perm(A.cols), go_col_perm(A.cols);

        auto go_order_start = std::chrono::high_resolution_clock::now();
        gorder(A, go_row_perm, go_col_perm, window_size);
        auto go_order_end = std::chrono::high_resolution_clock::now();
        auto go_order_dur = std::chrono::duration_cast<std::chrono::milliseconds>(go_order_end - go_order_start).count();

        std::cout << "GOrder took " << go_order_dur << " ms" << std::endl;
        A_perm = permute(A, go_row_perm, go_col_perm);
        std::cout << "Intermediate diagonal count after GOrder: " << count_nonempty_hs_diagonals(A_perm) << std::endl;

        auto hs_order_start = std::chrono::high_resolution_clock::now();
        hsorder_long(A_perm, hs_row_perm, hs_col_perm, 500, debug_mode);
        auto hs_order_end = std::chrono::high_resolution_clock::now();
        auto hs_order_dur = std::chrono::duration_cast<std::chrono::milliseconds>(hs_order_end - hs_order_start).count();

        std::cout << "HS Ordering took " << hs_order_dur << " ms" << std::endl;
        A_perm = permute(A_perm, hs_row_perm, hs_col_perm);

        for (size_t i = 0; i < A.rows; i++) best_row_perm[i] = hs_row_perm[go_row_perm[i]];
        for (size_t i = 0; i < A.cols; i++) best_col_perm[i] = hs_col_perm[go_col_perm[i]];
    } else if (order_id_int == 6) {
        // Independent Set Ordering

        std::vector<size_t> is_row_perm(A.rows), is_col_perm(A.cols), hs_row_perm(A.rows), hs_col_perm(A.cols);
        auto order_start = std::chrono::high_resolution_clock::now();
        independent_set_order(A, is_row_perm, is_col_perm, debug_mode);
        auto order_end = std::chrono::high_resolution_clock::now();
        auto order_dur = std::chrono::duration_cast<std::chrono::milliseconds>(order_end - order_start).count();
        
        std::cout << "Independent Set Ordering took " << order_dur << " ms" << std::endl;

        A_perm = permute(A, is_row_perm, is_col_perm);
        // std::cout << "Intermediate diagonal count after GOrder: " << count_nonempty_hs_diagonals(A_perm) << std::endl;

        // auto hs_order_start = std::chrono::high_resolution_clock::now();
        // hsorder_long(A_perm, hs_row_perm, hs_col_perm, 500, debug_mode);
        // auto hs_order_end = std::chrono::high_resolution_clock::now();
        // auto hs_order_dur = std::chrono::duration_cast<std::chrono::milliseconds>(hs_order_end - hs_order_start).count();
        
        // std::cout << "HS Ordering took " << hs_order_dur << " ms" << std::endl;
        // A_perm = permute(A_perm, hs_row_perm, hs_col_perm);

        // for (size_t i = 0; i < A.rows; i++) best_row_perm[i] = hs_row_perm[is_row_perm[i]];
        // for (size_t i = 0; i < A.cols; i++) best_col_perm[i] = hs_col_perm[is_col_perm[i]];
    } else if (order_id_int == 7) {
        // Our initial brute force HS search

        auto order_start = std::chrono::high_resolution_clock::now();
        hsorder(A, best_row_perm, best_col_perm);
        auto order_end = std::chrono::high_resolution_clock::now();
        auto order_dur = std::chrono::duration_cast<std::chrono::milliseconds>(order_end - order_start).count();

        std::cout << "HSBF Ordering took " << order_dur << " ms" << std::endl;
        A_perm = permute(A, best_row_perm, best_col_perm);
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