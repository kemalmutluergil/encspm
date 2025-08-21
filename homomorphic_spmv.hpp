#pragma once
#include "seal/seal.h"

#include "utils.hpp"
#include "csr.hpp"

using namespace std;
using namespace seal;

vector<Ciphertext> hs_csr_mult(const CSRMatrix& A, vector<double> x_plain, int slot_count,
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

    for (size_t i = 0; i < pt_diagonal_values.size(); i++) {
        if (pt_diagonal_values[i].size() == 0) continue;
        // Pad with zeros until the vector's size is slot_count
        pt_diagonal_values[i].resize(slot_count, 0.0);
    }

    if (A.rows > slot_count) {
        size_t num_chunks = (A.rows + slot_count - 1) / slot_count; // ceil(A.rows / slot_count)
        vector<Ciphertext> ct_results;

        for (size_t start = 0; start < A.rows; start += slot_count) {
            vector<Ciphertext> ct_diagonal_chunks;
            vector<Plaintext> pt_rotated_x_chunks;

            for (size_t i = 0; i < pt_diagonal_values.size(); i++) {
                if (pt_diagonal_values[i].size() == 0) continue;
                vector<double> diagonal_chunk = slice_and_pad(pt_diagonal_values[i], start, slot_count, slot_count);
                vector<double> rot_x = slice_and_pad(rotate_left(x_plain, i), start, slot_count, slot_count);
                Plaintext pt_x_rot, pt_diag;
                Ciphertext ct_diag;
                encoder.encode(rot_x, scale, pt_x_rot);
                encoder.encode(diagonal_chunk, scale, pt_diag);
                encryptor.encrypt(pt_diag, ct_diag);
                ct_diagonal_chunks.push_back(ct_diag);
                pt_rotated_x_chunks.push_back(pt_x_rot);
            }
            for (size_t i = 0; i < ct_diagonal_chunks.size(); i++) {
                evaluator.multiply_plain_inplace(ct_diagonal_chunks[i], pt_rotated_x_chunks[i]);
            }
            Ciphertext ct_result;
            evaluator.add_many(ct_diagonal_chunks, ct_result);
            ct_results.push_back(ct_result);
        }

        return ct_results;
        
    } else {
        vector<Ciphertext> ct_diagonals;
        vector<Plaintext> pt_rotated_xes;
        for (size_t i = 0; i < pt_diagonal_values.size(); i++) {
            if (pt_diagonal_values[i].size() == 0) continue;

            encoder.encode(pt_diagonal_values[i], scale, pt_diagonals[i]);
            Ciphertext ct;
            encryptor.encrypt(pt_diagonals[i], ct);
            ct_diagonals.push_back(ct);
            // rotate x_plain vector i steps to the left cyclically
            vector<double> rot_x = rotate_left(x_plain, i);
            Plaintext pt_x_rot;
            encoder.encode(rot_x, scale, pt_x_rot);
            pt_rotated_xes.push_back(pt_x_rot);
        }
        for (size_t i = 0; i < ct_diagonals.size(); i++) {
            evaluator.multiply_plain_inplace(ct_diagonals[i], pt_rotated_xes[i]);
        }
        // sum all ct_diagonals together
        Ciphertext ct_result;
        evaluator.add_many(ct_diagonals, ct_result);

        vector<Ciphertext> result(1, ct_result);
        return result;
    }

}

bool check_multiplications(const CSRMatrix& A, const CSRMatrix& A_perm, const std::vector<size_t>& best_row_perm, const std::vector<size_t>& best_col_perm) {
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

    std::vector<double> permuted_x(features_x.size());
    for (size_t i = 0; i < best_col_perm.size(); ++i) {
        permuted_x[i] = features_x[best_col_perm[i]];
    }

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
    vector<Ciphertext> ct_result = hs_csr_mult(A_perm, permuted_x, slot_count, encryptor, evaluator, encoder, gal_keys, scale);
    auto hs_end = chrono::high_resolution_clock::now();
    auto hs_dur = chrono::duration_cast<chrono::milliseconds>(hs_end - hs_start).count();
    std::cout << "Homomorphic mult duration: " << hs_dur << " milliseconds" << std::endl;
    
    vector<double> res_vec;
    size_t ct_idx = 0;
    for (size_t start = 0; start < A.rows; start += slot_count) {
        Plaintext pt_result;
        decryptor.decrypt(ct_result[ct_idx], pt_result);
        vector<double> temp;
        encoder.decode(pt_result, temp);
        res_vec.insert(res_vec.end(), temp.begin(), temp.end());
        ct_idx++;
    }        
    vector<double> permuted_result(res_vec.size());
    // Inverse row permutation
    std::vector<size_t> inv_row_perm(best_row_perm.size());
    for (size_t j = 0; j < best_row_perm.size(); j++) {
        inv_row_perm[best_row_perm[j]] = j;
    }
    for (size_t i = 0; i < best_row_perm.size(); ++i) permuted_result[i] = res_vec[inv_row_perm[i]];

    permuted_result.resize(A.cols);
    if (approx_equal(permuted_result, y)) {
        return true;
    } else {
        return false;
    }
}