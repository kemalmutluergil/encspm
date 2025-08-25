#pragma once

#include <map>

static map<string, chrono::high_resolution_clock::time_point> __prof_start;
static map<string, chrono::nanoseconds> __prof_durations;
static map<string, size_t> __prof_counts;

#include <unordered_set>

using namespace std;

#define PROFILE_START(key) \
    __prof_start[key] = chrono::high_resolution_clock::now();
#define PROFILE_END(key) \
    do { \
        auto __prof_end = chrono::high_resolution_clock::now(); \
        auto __d = chrono::duration_cast<chrono::nanoseconds>(__prof_end - __prof_start[key]); \
        __prof_durations[key] += __d; \
        ++__prof_counts[key]; \
    } while(0)

inline void print_profiling() {
    cout << "\n--- Profiling Results ---\n";
    for (auto &entry : __prof_durations) {
        const string &key = entry.first;
        auto total_ns = entry.second.count();
        size_t calls = __prof_counts[key];
        double total_ms = total_ns / 1e6;
        double avg_us = (total_ns / 1000.0) / calls;
        cout << key
             << ": total " << total_ms << " ms"
             << ", calls " << calls
             << ", avg " << avg_us << " us\n";
    }
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

bool approx_equal(const std::vector<double>& a, const std::vector<double>& b, double tol = 1e-6) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::fabs(a[i] - b[i]) > tol) return false;
    }
    return true;
}

vector<double> rotate_left(const vector<double>& v, size_t i) {
    if (v.empty()) return {};
    i = i % v.size();
    vector<double> out(v.size());
    std::rotate_copy(v.begin(), v.begin() + i, v.end(), out.begin());
    return out;
}


// Helper: compute wraparound diagonal index
size_t diag_index(size_t r, size_t c, size_t n) {
    return (c - r + n) % n;
}

static inline vector<double> slice_and_pad(const vector<double>& v,
                                           size_t start, size_t len, size_t slot_count) {
    if (start >= v.size())
        return vector<double>(slot_count, 0.0);

    vector<double> out;
    out.reserve(slot_count);
    for (size_t k = 0; k < std::min(len, v.size() - start); ++k)
        out.push_back(v[start + k]);

    if (out.size() < slot_count)
        out.resize(slot_count, 0.0);

    return out;
}


/// Decrypts the packed SpMV result, compares to y_expected, and reports errors.
// void checkSpMVResult(
//     const vector<Ciphertext> &y_encrypted,
//     const vector<double>    &y_expected,
//     Decryptor               &decryptor,
//     CKKSEncoder             &encoder)
// {
//     const size_t slot_count = encoder.slot_count();
//     const size_t nrows      = y_expected.size();

//     double max_abs_err = 0.0;
//     double sum_abs_err = 0.0;
//     size_t count       = 0;

//     // loop over each ciphertext
//     for (size_t ct_idx = 0; ct_idx < y_encrypted.size(); ++ct_idx) {
//         // 1) Decrypt
//         Plaintext pt;
//         decryptor.decrypt(y_encrypted[ct_idx], pt);

//         // 2) Decode to vector<double>
//         vector<double> decoded;
//         encoder.decode(pt, decoded);

//         // 3) Extract each slot as a row result
//         for (size_t slot = 0; slot < slot_count; ++slot) {
//             size_t row = ct_idx * slot_count + slot;
//             if (row >= nrows) break;

//             double approx = decoded[slot];
//             double exact  = y_expected[row];
//             double err    = fabs(approx - exact);

//             // Track stats
//             max_abs_err = max(max_abs_err, err);
//             sum_abs_err += err;
//             ++count;

//             // Optionally: print any row with large error
//             if (err > 1e-3) {
//                 cout << "[Row " << row << "] "
//                      << "approx=" << approx
//                      << ", exact=" << exact
//                      << ", abs_err=" << err << "\n";
//             }
//         }
//     }

//     double avg_abs_err = (count > 0 ? sum_abs_err / double(count) : 0.0);
//     cout << fixed << setprecision(6)
//          << "\nSpMV Result Check:\n"
//          << "  Rows compared       : " << count << "\n"
//          << "  Max absolute error  : " << max_abs_err << "\n"
//          << "  Avg absolute error  : " << avg_abs_err << "\n";
// }