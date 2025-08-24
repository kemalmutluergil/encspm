#pragma once
#include "seal/seal.h"

#include "utils.hpp"
#include "csr.hpp"

using namespace std;
using namespace seal;

struct chunk_kemal {
    vector<double> values;
    size_t index;
};

vector<chunk_kemal> chunkify_csr(const CSRMatrix& a) {
    vector<chunk_kemal> chunks(a.rows);
    for (size_t i = 0; i < a.rows; i++) {
        chunks[i].index = i;
    }
    for (size_t r = 0; r < a.rows; r++) {
        for (size_t k = a.row_ptr[r]; k < a.row_ptr[r+1]; k++) {
            const size_t diag_idx = (a.col_idx[k] + a.rows - r) % a.rows;
            if (chunks[diag_idx].values.size() < r) chunks[diag_idx].values.resize(r, 0.0);
            chunks[diag_idx].values.push_back(a.values[k]);
        }
    }

    // erase chunks with no values
    for (size_t i = 0; i < chunks.size(); i++) {
        if (chunks[i].values.size() == 0) {
            chunks.erase(chunks.begin() + i);
            i--;
        } else {
            // pad with zeros until the vector's size is a.cols
            chunks[i].values.resize(a.cols, 0.0);
        }
    }
    return chunks;
}

vector<Ciphertext> hs_chunks_mult(const vector<chunk_kemal>& chunks, vector<double> x_plain,
                                    size_t slot_count, Encryptor& encryptor, Evaluator& evaluator,
                                    CKKSEncoder& encoder, GaloisKeys& gal_keys, double scale, Decryptor& decryptor) {
    if (chunks.size() == 0) return vector<Ciphertext>();
    if (chunks[0].values.size() > slot_count) {

    } else if (chunks[0].values.size() <= slot_count / 2) {
        size_t pack_num = slot_count / chunks[0].values.size();
        size_t largest_pow2 = 1;
        while (largest_pow2 * 2 < pack_num) {
            largest_pow2 *= 2;
        }
        pack_num = largest_pow2;
        size_t ct_count = (chunks.size() + pack_num - 1) / pack_num; // ceil(chunks.size() / pack_num)

        std::vector<Ciphertext> ct_packed_chunks(ct_count);
        for (size_t i = 0; i < ct_count; i++) {
            size_t start = i * pack_num;
            size_t end = min(start + pack_num - 1, chunks.size() - 1);

            std::vector<double> packed_values;
            std::vector<double> packed_x;

            for (size_t j = start; j <= end; j++) {
                packed_values.insert(packed_values.end(), chunks[j].values.begin(), chunks[j].values.end());
                std::vector<double> rotated_x = rotate_left(x_plain, chunks[j].index);
                packed_x.insert(packed_x.end(), rotated_x.begin(), rotated_x.end());
            }
            packed_values.resize(slot_count, 0.0); // pad with zeros
            packed_x.resize(slot_count, 0.0); // pad with zeros

            Plaintext pt_chunk, pt_x;
            encoder.encode(packed_values, scale, pt_chunk);
            encryptor.encrypt(pt_chunk, ct_packed_chunks[i]);
            encoder.encode(packed_x, scale, pt_x);

            Plaintext temp_chunk;
            std::vector<double> temp_chunk_vec, temp_x_vec;
            decryptor.decrypt(ct_packed_chunks[i], temp_chunk);
            encoder.decode(temp_chunk, temp_chunk_vec);
            encoder.decode(pt_x, temp_x_vec);

            evaluator.multiply_plain_inplace(ct_packed_chunks[i], pt_x);
        }
        return ct_packed_chunks;
    } else {
        vector<Ciphertext> ct_chunks(chunks.size());
        vector<Plaintext> pt_chunks(chunks.size());

        for (size_t i = 0; i < chunks.size(); i++) {
            std::vector<double> vals = chunks[i].values;
            vals.resize(slot_count, 0.0);

            encoder.encode(vals, scale, pt_chunks[i]);
            encryptor.encrypt(pt_chunks[i], ct_chunks[i]);

            Plaintext rot_x;
            x_plain.resize(slot_count, 0.0);
            encoder.encode(rotate_left(x_plain, chunks[i].index), scale, rot_x);

            evaluator.multiply_plain_inplace(ct_chunks[i], rot_x);
        }

        return ct_chunks;
    }
}

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
    std::vector<double> inv_col_perm(best_col_perm.size());
    for (size_t i = 0; i < best_col_perm.size(); ++i) {
        inv_col_perm[best_col_perm[i]] = i;
    }
    for (size_t i = 0; i < inv_col_perm.size(); i++) permuted_x[i] = features_x[inv_col_perm[i]];

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
    std::vector<chunk_kemal> chunks = chunkify_csr(A_perm);
    vector<Ciphertext> ct_result = hs_chunks_mult(chunks, permuted_x, slot_count, encryptor, evaluator, encoder, gal_keys, scale, decryptor);
    auto hs_end = chrono::high_resolution_clock::now();
    auto hs_dur = chrono::duration_cast<chrono::milliseconds>(hs_end - hs_start).count();
    std::cout << "Homomorphic mult duration: " << hs_dur << " milliseconds" << std::endl;
    vector<double> res_vec;

    vector<double> temp_vec;
    Plaintext temp_pt;
    decryptor.decrypt(ct_result[0], temp_pt);
    encoder.decode(temp_pt, temp_vec);

    if (chunks[0].values.size() > slot_count) {

    } else if (chunks[0].values.size() <= slot_count / 2) {
        size_t pack_num = slot_count / chunks[0].values.size();
        size_t largest_pow2 = 1;
        while (largest_pow2 * 2 < pack_num) largest_pow2 *= 2;
        size_t pow2_2 = 1;
        while (pow2_2 < chunks.size()) pow2_2 *= 2;
        pack_num = min(largest_pow2, pow2_2);

        Ciphertext sum_res;
        evaluator.add_many(ct_result, sum_res);

        size_t j = 1;
        for (size_t i = 1; i < pack_num; i *= 2) {
            int offset = chunks[0].values.size() * j;
            Ciphertext rotated;
            evaluator.rotate_vector(sum_res, offset, gal_keys, rotated);
            evaluator.add_inplace(sum_res, rotated);
            j*=2;
        }

        Plaintext res;
        decryptor.decrypt(sum_res, res);
        encoder.decode(res, res_vec);
    } else {

    }
           


    vector<double> permuted_result(res_vec.size());
    // Inverse row permutation
    std::vector<size_t> inv_row_perm(best_row_perm.size());
    for (size_t j = 0; j < best_row_perm.size(); j++) {
        inv_row_perm[best_row_perm[j]] = j;
    }
    for (size_t old = 0; old < best_row_perm.size(); ++old) permuted_result[old] = res_vec[best_row_perm[old]];

    permuted_result.resize(y.size());
    if (approx_equal(permuted_result, y)) {
        return true;
    } else {
        return false;
    }
}