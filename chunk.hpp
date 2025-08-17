#pragma once 

enum Orientation {
  ROWWISE,
  COLWISE,
  DIAGONAL,
  BICYCLIC,
  ELLPACK,
  HALEVISHOUP,
  ARBITRARY
};

struct Chunk {
    public: 
        size_t start_i;
        size_t start_j;
        size_t width;
        size_t height; 
        size_t inner_sz;   
        enum Orientation orientation;
        
        //normally these are not stored - for experimental purposes
        vector< std::pair<std::pair<size_t, size_t>, double>> nonzeros;
        vector<double> plaintext; 

        // Ciphertext ciphertext;

        Chunk(int i, int j, int width, int height, int inner_sz, Orientation orientation) : start_i(i), start_j(j), 
                                                                                            width(width), height(height), inner_sz(inner_sz),
                                                                                            orientation(orientation) {
            plaintext.resize(width * height, 0);
        }
};


void printChunkStats(const std::vector<Chunk> &chunks, const CSRMatrix& A, Orientation opt) {
    // 1) Compute number of chunks per row
    std::vector<size_t> chunks_per_row(A.rows, 0);
    std::vector<size_t> chunks_per_col(A.cols, 0);
    std::vector<size_t> chunks_per_diag(A.cols, 0);

    for (const auto &chunk : chunks) {
        if(chunk.orientation == ROWWISE) { 
            chunks_per_row[chunk.start_i]++;
        } else if(chunk.orientation == COLWISE) { 
            chunks_per_col[chunk.start_j]++;
        } else if(chunk.orientation == DIAGONAL) { 
            chunks_per_diag[chunk.start_j]++;
        }
    }

    // 2) Compute nonzeros per chunk ("density")
    std::vector<size_t> nnz_per_chunk;
    nnz_per_chunk.reserve(chunks.size());
    for (const auto &chunk : chunks) {
        nnz_per_chunk.push_back(chunk.nonzeros.size());
    }

    auto compute_stats = [](std::vector<size_t> &data) {
        size_t n = data.size();
        std::sort(data.begin(), data.end());
        size_t min_v = data.front();
        size_t max_v = data.back();
        double sum = std::accumulate(data.begin(), data.end(), 0ull);
        double avg = sum / double(n);
        double median = (n % 2 == 0)
            ? (data[n/2 - 1] + data[n/2]) / 2.0
            : data[n/2];
        double sq_sum = 0.0;
        for (auto v : data) {
            double d = double(v) - avg;
            sq_sum += d * d;
        }
        double stdev = std::sqrt(sq_sum / double(n));
        return std::tuple<size_t,size_t,double,double,double>(min_v, max_v, avg, median, stdev);
    };


    // Print
    std::cout << "Number of chunks: " << chunks.size() << "\n";

    // Stats for chunks-per-row
    if(opt == ROWWISE || opt == ARBITRARY) {
        auto [min_cpr, max_cpr, avg_cpr, med_cpr, sd_cpr] = compute_stats(chunks_per_row);
        std::cout << "Min chunks in a row:    " << min_cpr << "\n";
        std::cout << "Max chunks in a row:    " << max_cpr << "\n";
        std::cout << "Avg chunks per row:     " << avg_cpr << "\n";
        std::cout << "Median chunks per row:  " << med_cpr << "\n";
        std::cout << "Stdev chunks per row:   " << sd_cpr << "\n\n";
    }
    // Stats for chunks-per-col
    if(opt == COLWISE || opt == ARBITRARY) {
        auto [min_cpr_c, max_cpr_c, avg_cpr_c, med_cpr_c, sd_cpr_c] = compute_stats(chunks_per_col);
        std::cout << "Min chunks in a col:    " << min_cpr_c << "\n";
        std::cout << "Max chunks in a col:    " << max_cpr_c << "\n";
        std::cout << "Avg chunks per col:     " << avg_cpr_c << "\n";
        std::cout << "Median chunks per col:  " << med_cpr_c << "\n";
        std::cout << "Stdev chunks per col:   " << sd_cpr_c << "\n\n";
    }
    // Stats for chunks-per-col
    if(opt == DIAGONAL || opt == ARBITRARY) {
        auto [min_cpr_d, max_cpr_d, avg_cpr_d, med_cpr_d, sd_cpr_d] = compute_stats(chunks_per_diag);
        std::cout << "Min chunks in a diag:  " << min_cpr_d << "\n";
        std::cout << "Max chunks in a diag:  " << max_cpr_d << "\n";
        std::cout << "Avg chunks in a diag:  " << avg_cpr_d << "\n";
        std::cout << "Median chunks in a diag:  " << med_cpr_d << "\n";
        std::cout << "Stdev chunks in a diag:   " << sd_cpr_d << "\n\n";
    }

    // Stats for nnz-per-chunk
    std::cout << "Slot count: " << chunks[0].width * chunks[0].height << "\n";
    auto [min_nnz, max_nnz, avg_nnz, med_nnz, sd_nnz] = compute_stats(nnz_per_chunk);
        std::cout << "Min nonzeros in chunk:  " << min_nnz << "\n";
        std::cout << "Max nonzeros in chunk:  " << max_nnz << "\n";
        std::cout << "Avg nonzeros per chunk: " << avg_nnz << " - must be " << (A.row_ptr[A.rows] + 0.0f) / chunks.size() << "\n";
        std::cout << "Median nonzeros:        " << med_nnz << "\n";
        std::cout << "Stdev nonzeros:         " << sd_nnz << "\n";

    // After `vector<Chunk> chunks = rowChunksNaive(...);`
    size_t total_bytes = 0;
    for (auto const &chunk : chunks) {
        std::ostringstream oss;
        // chunk.ciphertext.save(oss);       // serialize ct into the stream
        total_bytes += oss.tellp();       // add the number of bytes written
    }

    std::cout << "Total ciphertext size: " << total_bytes / (1024.0 * 1024.0) << " MB\n";
}