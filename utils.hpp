
static map<string, chrono::high_resolution_clock::time_point> __prof_start;
static map<string, chrono::nanoseconds> __prof_durations;
static map<string, size_t> __prof_counts;

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


/// Decrypts the packed SpMV result, compares to y_expected, and reports errors.
void checkSpMVResult(
    const vector<Ciphertext> &y_encrypted,
    const vector<double>    &y_expected,
    Decryptor               &decryptor,
    CKKSEncoder             &encoder)
{
    const size_t slot_count = encoder.slot_count();
    const size_t nrows      = y_expected.size();

    double max_abs_err = 0.0;
    double sum_abs_err = 0.0;
    size_t count       = 0;

    // loop over each ciphertext
    for (size_t ct_idx = 0; ct_idx < y_encrypted.size(); ++ct_idx) {
        // 1) Decrypt
        Plaintext pt;
        decryptor.decrypt(y_encrypted[ct_idx], pt);

        // 2) Decode to vector<double>
        vector<double> decoded;
        encoder.decode(pt, decoded);

        // 3) Extract each slot as a row result
        for (size_t slot = 0; slot < slot_count; ++slot) {
            size_t row = ct_idx * slot_count + slot;
            if (row >= nrows) break;

            double approx = decoded[slot];
            double exact  = y_expected[row];
            double err    = fabs(approx - exact);

            // Track stats
            max_abs_err = max(max_abs_err, err);
            sum_abs_err += err;
            ++count;

            // Optionally: print any row with large error
            if (err > 1e-3) {
                cout << "[Row " << row << "] "
                     << "approx=" << approx
                     << ", exact=" << exact
                     << ", abs_err=" << err << "\n";
            }
        }
    }

    double avg_abs_err = (count > 0 ? sum_abs_err / double(count) : 0.0);
    cout << fixed << setprecision(6)
         << "\nSpMV Result Check:\n"
         << "  Rows compared       : " << count << "\n"
         << "  Max absolute error  : " << max_abs_err << "\n"
         << "  Avg absolute error  : " << avg_abs_err << "\n";
}