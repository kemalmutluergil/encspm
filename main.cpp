#include "seal/seal.h"
#include "csr.hpp"
#include "rcm.hpp"
#include "chunk.hpp"

#include <map>
#include <unordered_set>

#define SLOT_COUNT 8192
#define NUM_TRIALS 1

using namespace seal;

bool approx_equal(const std::vector<double>& a, const std::vector<double>& b, double tol = 1e-6) {
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
    vector<Plaintext> pt_diagonals(A.cols);
    vector<vector<double>> pt_diagonal_values(A.cols);
    
    for (size_t r = 0; r < A.rows; r++) {
        for (size_t idx = A.row_ptr[r]; idx < A.row_ptr[r+1]; idx++) {
            int diag_idx = diag_index(r, A.col_idx[idx], A.rows);
            // Pad with zeros until the vector's size is r
            while (pt_diagonal_values[diag_idx].size() < r) {
                pt_diagonal_values[diag_idx].push_back(0.0);
            }
            pt_diagonal_values[diag_idx].push_back(A.values[idx]);
        }
    }

    vector<Ciphertext> ct_diagonals;
    vector<Plaintext> pt_rotated_xes;

    for (size_t i = 0; i < pt_diagonal_values.size(); i++) {
        if (pt_diagonal_values[i].size() == 0) continue;

        if (pt_diagonal_values[i].size() < slot_count) {
            pt_diagonal_values[i].resize(slot_count, 0.0); // pad with zeros
        }
        encoder.encode(pt_diagonal_values[i], scale, pt_diagonals[i]);
        Ciphertext ct;
        encryptor.encrypt(pt_diagonals[i], ct);
        ct_diagonals.push_back(ct);

        // rotate x_plain vector i steps to the left cyclically
        Plaintext pt_x_rot;
        encoder.encode(rotate_left(x_plain, i), scale, pt_x_rot);
        pt_rotated_xes.push_back(pt_x_rot);
    }

    for (size_t i = 0; i < ct_diagonals.size(); i++) {
        evaluator.multiply_plain_inplace(ct_diagonals[i], pt_rotated_xes[i]);
    }

    // sum all ct_diagonals together
    Ciphertext ct_result;
    evaluator.add_many(ct_diagonals, ct_result);

    return ct_result;

}

size_t count_nonempty_hs_diagonals(const CSRMatrix &A) {
    std::unordered_set<size_t> nonempty_diags;

    for (size_t r = 0; r < A.rows; ++r) {
        for (size_t k = A.row_ptr[r]; k < A.row_ptr[r + 1]; ++k) {
            size_t c = A.col_idx[k];
            size_t diag = (c + A.rows - r) % A.rows; // (j - i + n) mod n
            nonempty_diags.insert(diag);
        }
    }

    return nonempty_diags.size();
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

    CSRMatrix A_perm;
    int order_id = atoi(argv[2]);

    if (order_id == 0) {
        A_perm = A;
    } else {
        size_t m = A.rows, n = A.cols;

        std::vector<size_t> best_row_perm, best_col_perm;

        std::vector<size_t> row_perm, col_perm;

        if(order_id == 1) {   
            int min_diag_count = n + 1;
            for (size_t trial_num = 0; trial_num < NUM_TRIALS; trial_num++) {
                row_perm.clear();
                col_perm.clear();
                
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

                A_perm = permute(A, row_perm, col_perm);
                int diag_count = count_nonempty_hs_diagonals(A_perm);

                cout << "Trial " << trial_num << ": " << diag_count << " non-empty diagonals" << endl;               
                if (diag_count < min_diag_count) {
                    min_diag_count = diag_count;
                    best_row_perm = row_perm;
                    best_col_perm = col_perm;
                }
            }
        }
        A_perm = permute(A, best_row_perm, best_col_perm);
        
        // visualize_csr(A_perm);
        // print_csr_info(A_perm);

        // cout << "col perm: ";
        // for (const auto& val : best_col_perm) {
        //     cout << val << " ";
        // }
        // cout << std::endl;

        std::vector<double> permuted_x(features_x.size());
        for (size_t i = 0; i < best_col_perm.size(); ++i) {
            permuted_x[i] = features_x[best_col_perm[i]];
        }

        // cout << "permuted x: ";
        // for (size_t i = 0; i < permuted_x.size(); ++i) {
        //     std::cout << permuted_x[i] << " ";
        // }
        // cout << std::endl;

        // std::vector<double> result(A_perm.rows);
        // for (size_t r = 0; r < A_perm.rows; ++r) {
        //     double sum = 0.0;
        //     for (size_t k = A_perm.row_ptr[r]; k < A_perm.row_ptr[r+1]; ++k) {
        //         sum += A_perm.values[k] * permuted_x[A_perm.col_idx[k]];
        //     }
        //     result[r] = sum;
        // }

        // std::vector<double> result_permuted(y.size());
        // for (size_t i = 0; i < row_perm.size(); ++i) {
        //     result_permuted[i] = result[row_perm[i]];
        // }

        // if (approx_equal(result_permuted, y)) {
        //     std::cout << "Success: The results match!" << std::endl;
        // } else {
        //     std::cout << "Error: The results do not match!" << std::endl;
        // }

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

        auto hs_start = chrono::high_resolution_clock::now();
        Ciphertext ct_result = hs_csr_mult(A_perm, permuted_x, 4, encryptor, evaluator, encoder, gal_keys, scale);
        auto hs_end = chrono::high_resolution_clock::now();
        auto hs_dur = chrono::duration_cast<chrono::milliseconds>(hs_end - hs_start).count();
        cout << "Homomorphic mult duration: " << hs_dur << " milliseconds" << endl;
        // Decrypt
        Plaintext plain_result;
        decryptor.decrypt(ct_result, plain_result);
        vector<double> res_vec;
        encoder.decode(plain_result, res_vec);
        
        vector<double> permuted_result(res_vec.size());
        // Inverse row permutation
        std::vector<size_t> inv_row_perm(best_row_perm.size());
        for (size_t j = 0; j < best_row_perm.size(); j++) {
            inv_row_perm[best_row_perm[j]] = j;
        }
        for (size_t i = 0; i < best_row_perm.size(); ++i) permuted_result[i] = res_vec[inv_row_perm[i]];

        permuted_result.resize(A.cols);
        if (approx_equal(permuted_result, y)) {
            std::cout << "Success: The results match!" << std::endl;
        } else {
            std::cout << "Failed: The results do not match!" << std::endl;    
        }
    }
    return 0;
}