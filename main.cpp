#include "seal/seal.h"
#include "csr.hpp"
#include "rcm.hpp"
#include "chunk.hpp"

#include <map>

#define SLOT_COUNT 8192

using namespace seal;

bool approx_equal(const std::vector<double>& a, const std::vector<double>& b, double tol = 1e-12) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::fabs(a[i] - b[i]) > tol) return false;
    }
    return true;
}

vector<double> rotate_left(std::vector<double> v, size_t i) {
    if (v.empty()) return v;
    i = i % v.size(); // handle i >= v.size()
    std::rotate(v.begin(), v.begin() + i, v.end());
    return v;
}

// Helper: compute wraparound diagonal index
size_t diag_index(size_t r, size_t c, size_t n) {
    return (c - r + n) % n;
}

Ciphertext hs_csr_mult(const CSRMatrix& A, vector<double> x_plain, int slot_count,
                 Encryptor& encryptor, Evaluator& evaluator, CKKSEncoder& encoder,
                 GaloisKeys& gal_keys, double scale) {
    vector<Ciphertext> ct_diagonals(A.cols);
    vector<Plaintext> pt_diagonals(A.cols);
    vector<vector<double>> pt_diagonal_values(A.cols);
    int max_diag_idx = -1;
    
    for (size_t r = 0; r < A.rows; r++) {
        for (size_t idx = A.row_ptr[r]; idx < A.row_ptr[r+1]; idx++) {
            int diag_idx = diag_index(r, A.col_idx[idx], A.rows);
            // Pad with zeros until the vector's size is r
            while (pt_diagonal_values[diag_idx].size() < r) {
                pt_diagonal_values[diag_idx].push_back(0.0);
            }
            pt_diagonal_values[diag_idx].push_back(A.values[idx]);
            if (diag_idx > max_diag_idx) {
                max_diag_idx = diag_idx;
            }
        }
    }
    for (size_t i = 0; i < pt_diagonal_values.size(); i++) {
        if (pt_diagonal_values[i].size() < slot_count) {
            pt_diagonal_values[i].resize(slot_count, 0.0); // pad with zeros
        }
        encoder.encode(pt_diagonal_values[i], scale, pt_diagonals[i]);
        encryptor.encrypt(pt_diagonals[i], ct_diagonals[i]);

        // rotate x_plain vector i steps to the left cyclically
        Plaintext pt_x_rot;
        encoder.encode(rotate_left(x_plain, i), scale, pt_x_rot);
        evaluator.multiply_plain_inplace(ct_diagonals[i], pt_x_rot);
    }

    // sum all ct_diagonals together
    Ciphertext ct_result;
    evaluator.add_many(ct_diagonals, ct_result);

    return ct_result;

}

vector<Chunk> halevishoup_pretest(const CSRMatrix& A, int slot_count) {
    // Use a map to group non-zeros by their chunk identifier.
    // The key is {diagonal_index, row_chunk_index}.
    map<pair<size_t, size_t>, vector<pair<size_t, double>>> chunk_data;

    // 1. Iterate through the CSR matrix and populate the map.
    for (size_t r = 0; r < A.rows; ++r) {
        for (size_t ptr = A.row_ptr[r]; ptr < A.row_ptr[r+1]; ++ptr) {
            size_t c = A.col_idx[ptr];
            double val = A.values[ptr];
            
            // Calculate wraparound diagonal index: ((i + k) mod m, k) for 0 <= k < n
            size_t diag_idx = (((r - c) % A.rows) + A.rows) % A.rows;
            
            // Determine which row-chunk this non-zero belongs to.
            size_t chunk_row_index = c / slot_count;
            size_t row_in_chunk_slot = c % slot_count;

            // Group the non-zero value by its {diagonal, row_chunk} key.
            chunk_data[{chunk_row_index, diag_idx}].push_back({row_in_chunk_slot, val});
        }
    }

    vector<Chunk> chunks;
    chunks.reserve(chunk_data.size());

    // 2. Create a Chunk for each entry in the map.
    for (auto const& [key, nonzeros] : chunk_data) {
        size_t diag_idx = key.second;
        size_t chunk_row_idx = key.first;
        size_t start_row = chunk_row_idx * slot_count;

        // A diagonal chunk conceptually has width 1 and height slot_count.
        Chunk chunk(start_row, diag_idx, 1, slot_count, slot_count, DIAGONAL);

        for (auto const& nz : nonzeros) {
            size_t row_in_chunk_slot = nz.first;
            double val = nz.second;
            chunk.plaintext[row_in_chunk_slot] = val;

            // Store original coordinates for stats/debugging.
            size_t orig_row = start_row + row_in_chunk_slot;
            size_t orig_col = (diag_idx + orig_row) % A.cols;
            chunk.nonzeros.push_back({ {orig_row, orig_col}, val });
        }
        chunks.emplace_back(std::move(chunk));
    }
    return chunks;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input.mtx> <order_id>" << std::endl;
        return 1;
    }

    std::cout << "Reading matrix " << argv[1] << std::endl;
    
    CSRMatrix A = read_mm(argv[1]);

    // visualize_csr(A);

    std::cout << "Read matrix...\n";
    std::cout << "Rows: " << A.rows << ", Cols: " << A.cols << ", Non-zeros: " << A.values.size() << std::endl;
    std::vector<double> features_x(A.cols);
    for (size_t i = 0; i < A.cols; ++i) {
        // features_x[i] = static_cast<double>(rand()) / RAND_MAX;
        features_x[i] = 1.0 + i*0.01;
    }

    std::vector<double> y(A.rows);
    for (size_t r = 0; r < A.rows; ++r) {
        double sum = 0.0;
        for (size_t k = A.row_ptr[r]; k < A.row_ptr[r+1]; ++k) {
            sum += A.values[k] * features_x[A.col_idx[k]];
        }
        y[r] = sum;
    }
    
    // // Example: print the nonzeros of first 10 rows
    // size_t limit = std::min<size_t>(10, A.rows);
    // for (size_t r = 0; r < limit; ++r) {
    //     std::cout << "A(" << r << ", *): ";
    //     for(size_t ptr = A.row_ptr[r]; ptr < A.row_ptr[r+1]; ++ptr) {
    //         size_t col = A.col_idx[ptr];
    //         double v = A.values[ptr];
    //         std::cout << "(" << col << ", " << v << ") ";
    //     } 
    //     std::cout << std::endl;
    // }

    CSRMatrix A_perm;
    int order_id = atoi(argv[2]);

    if (order_id == 0) {
        A_perm = A;
    } else {
        size_t m = A.rows, n = A.cols;

        std::vector<size_t> best_row_perm, best_col_perm;
        size_t best_perm_chunk_count = m * n; 

        std::vector<size_t> row_perm, col_perm;

        if(order_id == 1) {   
            if(m != n) {
                CSRMatrix B = build_bipartite_sym(A);
                auto perm = rcm::rcm_ordering(B);
            
                for (auto p : perm) {
                    if (p < m) {
                        row_perm.push_back(p);
                    } else {
                        col_perm.push_back(p - m);
                    }
                }
            } else {
                CSRMatrix B = build_sym(A);
                auto perm = rcm::rcm_ordering(B);

                for (auto p : perm) {
                    row_perm.push_back(p);
                    col_perm.push_back(p);
                }
            }
        }

        // Suppose rcm_ordering gives P_inv: new -> old
        std::vector<size_t> row_perm_fwd(A.rows), col_perm_fwd(A.cols);

        for (size_t i = 0; i < A.rows; ++i)
            row_perm_fwd[row_perm[i]] = i;  // forward: old -> new

        for (size_t i = 0; i < A.cols; ++i)
            col_perm_fwd[col_perm[i]] = i;  // forward: old -> new

        A_perm = permute(A, row_perm_fwd, col_perm_fwd);

        // visualize_csr(A_perm);
        // print_csr_info(A_perm);

        // cout << "row perm: ";
        // for (const auto& val : row_perm_fwd) {
        //     cout << val << " ";
        // }
        // cout << std::endl;
        // cout << "col perm: ";
        // for (const auto& val : col_perm_fwd) {
        //     cout << val << " ";
        // }
        // cout << std::endl;

        std::vector<double> permuted_x(features_x.size());
        for (size_t i = 0; i < col_perm_fwd.size(); ++i) {
            permuted_x[i] = features_x[col_perm_fwd[i]];
        }

        // cout << "permuted x: ";
        // for (size_t i = 0; i < permuted_x.size(); ++i) {
        //     std::cout << permuted_x[i] << " ";
        // }
        // cout << std::endl;

        std::vector<double> result(A_perm.rows);
        for (size_t r = 0; r < A_perm.rows; ++r) {
            double sum = 0.0;
            for (size_t k = A_perm.row_ptr[r]; k < A_perm.row_ptr[r+1]; ++k) {
                sum += A_perm.values[k] * permuted_x[A_perm.col_idx[k]];
            }
            result[r] = sum;
        }

        std::vector<double> result_permuted(y.size());
        for (size_t i = 0; i < row_perm.size(); ++i) {
            result_permuted[i] = result[row_perm[i]];
        }

        if (approx_equal(result_permuted, y)) {
            std::cout << "Success: The results match!" << std::endl;
            for (size_t i = 0; i < result_permuted.size(); ++i) {
                std::cout << "Row " << i << ": " << result_permuted[i] << std::endl;
            }
        } else {
            std::cout << "Error: The results do not match!" << std::endl;
        }

        // auto chunksPermuted = halevishoup_pretest(A_perm, SLOT_COUNT);
        // std::cout << "Order trial: " << chunksPermuted.size() << std::endl;

        // if(chunksPermuted.size() < best_perm_chunk_count) {
        //     best_row_perm = row_perm;
        //     best_col_perm = col_perm;    
        //     best_perm_chunk_count = chunksPermuted.size();    
        // }
        // // visualize_csr(A_perm);

        cout << "**********************************CKKS************************************\n";

        size_t poly_modulus_degree = 8192;
        double scale = pow(2.0, 40);
        EncryptionParameters parms(scheme_type::ckks);
        parms.set_poly_modulus_degree(poly_modulus_degree);
        parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {60,40,40,60}));

        SEALContext context(parms);
        KeyGenerator keygen(context);
        auto secret_key = keygen.secret_key();
        PublicKey public_key;
        keygen.create_public_key(public_key);
        RelinKeys relin_keys;
        keygen.create_relin_keys(relin_keys);
        GaloisKeys gal_keys;
        keygen.create_galois_keys(gal_keys);
        Encryptor encryptor(context, public_key);
        Evaluator evaluator(context);
        Decryptor decryptor(context, secret_key);

        CKKSEncoder encoder(context);
        size_t slot_count = encoder.slot_count();

        Ciphertext ct_result = hs_csr_mult(A, features_x, 4, encryptor, evaluator, encoder, gal_keys, scale);

        // Decrypt
        Plaintext plain_result;
        decryptor.decrypt(ct_result, plain_result);
        vector<double> res_vec;
        encoder.decode(plain_result, res_vec);

        cout << "Result: ";
        for (size_t i = 0; i < A.cols; i++) cout << res_vec[i] << " ";
        cout << endl;

        }
    return 0;
}