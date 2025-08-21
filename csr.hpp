#pragma once 

#include <vector>
#include <iostream>
#include <random>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <tuple>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <set>

using namespace std;
// Simple CSR sparse matrix representation
struct CSRMatrix {
    public:
        size_t rows;
        size_t cols;
        vector<size_t> row_ptr;
        vector<size_t> col_idx;
        vector<double> values;
};

void print_csr_info(const CSRMatrix& A) {
    for (size_t r = 0; r < A.rows; ++r) {
        cout << "Row " << r << ": ";
        for (size_t ptr = A.row_ptr[r]; ptr < A.row_ptr[r+1]; ++ptr) {
            cout << "(" << A.col_idx[ptr] << ", " << A.values[ptr] << ") ";
        }
        cout << endl;
    }
}

void visualize_csr(const CSRMatrix& A, size_t max_height = 40, size_t max_width = 100) {
    if (A.rows == 0 || A.cols == 0) {
        cout << "[Empty matrix]" << endl;
        return;
    }
    
    cout << "*********************************************************************************\n";
    
    // Scale rows/cols into a fixed grid
    size_t out_rows = min(max_height, A.rows);
    size_t out_cols = min(max_width, A.cols);

    vector<string> grid(out_rows, string(out_cols, '.'));

    for (size_t r = 0; r < A.rows; ++r) {
        size_t start = A.row_ptr[r];
        size_t end   = A.row_ptr[r+1];
        for (size_t k = start; k < end; ++k) {
            size_t c = A.col_idx[k];
            // Scale down to output grid
            size_t rr = r * out_rows / A.rows;
            size_t cc = c * out_cols / A.cols;
            grid[rr][cc] = 'x';
        }
    }

    // Print result
    for (auto& line : grid) {
        cout << line << "\n";
    }
    cout << "*********************************************************************************\n";
}

// Reads a Matrix Market (.mtx) file and returns the matrix in CSR format
CSRMatrix read_mm(const std::string &filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open file: " + filename);

    bool symmetric = false;
    bool pattern_only = false;
    std::string line;

    // Read header
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '%') {
            if (line.find("symmetric") != std::string::npos) symmetric = true;
            if (line.find("pattern") != std::string::npos) pattern_only = true;
            continue;
        } else {
            break;  // first non-comment line is the size line
        }
    }

    std::istringstream header(line);
    size_t M, N, NNZ;
    header >> M >> N >> NNZ;

    std::vector<std::tuple<size_t, size_t, double>> entries;
    entries.reserve(symmetric ? 2 * NNZ : NNZ);

    size_t i, j;
    double val;

    for (size_t k = 0; k < NNZ; ++k) {
        if (pattern_only) {
            in >> i >> j;
            val = 1.0; // default value for pattern
        } else {
            in >> i >> j >> val;
        }

        // Convert to 0-based
        size_t row = i - 1;
        size_t col = j - 1;
        entries.emplace_back(row, col, val);

        if (symmetric && row != col) {
            entries.emplace_back(col, row, val);
        }
    }
    in.close();

    // Sort entries by row then column
    std::sort(entries.begin(), entries.end(), [](auto &a, auto &b) {
        if (std::get<0>(a) != std::get<0>(b))
            return std::get<0>(a) < std::get<0>(b);
        return std::get<1>(a) < std::get<1>(b);
    });

    // Remove duplicates (keep first value)
    entries.erase(std::unique(entries.begin(), entries.end(),
        [](auto &a, auto &b) { return std::get<0>(a) == std::get<0>(b) && std::get<1>(a) == std::get<1>(b); }),
        entries.end());

    size_t NNZ_final = entries.size();

    // Build CSR
    CSRMatrix A;
    A.rows = M;
    A.cols = N;
    A.row_ptr.assign(M + 1, 0);
    A.col_idx.resize(NNZ_final);
    A.values.resize(NNZ_final);

    // Count entries per row
    for (auto &t : entries) {
        A.row_ptr[std::get<0>(t) + 1]++;
    }

    // Cumulative sum for row_ptr
    for (size_t r = 1; r <= M; ++r) {
        A.row_ptr[r] += A.row_ptr[r - 1];
    }

    // Fill col_idx and values
    std::vector<size_t> offset = A.row_ptr;
    for (auto &t : entries) {
        size_t row = std::get<0>(t);
        size_t dest = offset[row]++;
        A.col_idx[dest] = std::get<1>(t);
        A.values[dest] = std::get<2>(t);
    }

    return A;
}

/**
 * @brief For an m x n matrix A, builds the symmetric (m+n) x (m+n) adjacency matrix.
 */
CSRMatrix build_bipartite_sym(const CSRMatrix &A) {
    size_t m = A.rows, n = A.cols;
    size_t dim = m + n;
    std::vector<std::vector<size_t>> adj(dim);

    for (size_t i = 0; i < m; ++i) {
        for (size_t k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
            size_t j = A.col_idx[k];
            adj[i].push_back(m + j); 
            adj[m + j].push_back(i);
        }
    }

    CSRMatrix B;
    B.rows = B.cols = dim;
    B.row_ptr.resize(dim + 1, 0);
    size_t nnz = 0;
    for (size_t i = 0; i < dim; ++i) {
        std::sort(adj[i].begin(), adj[i].end());
        B.row_ptr[i + 1] = B.row_ptr[i] + adj[i].size();
        nnz += adj[i].size();
    }
    B.col_idx.reserve(nnz);
    for (const auto &lst : adj) {
        B.col_idx.insert(B.col_idx.end(), lst.begin(), lst.end());
    }
    B.values.assign(nnz, 1.0);
    return B;
}

/**
 * @brief For an n x n matrix A, builds the symmetric (n) x (n) adjacency matrix.
 */
CSRMatrix build_sym(const CSRMatrix &A) {
    size_t n = A.rows;
    std::vector<std::vector<size_t>> adj(n);

    for (size_t i = 0; i < n; ++i) {
        for (size_t k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
            size_t j = A.col_idx[k];
            adj[i].push_back(j); 
            adj[j].push_back(i);
        }
    }
    
    CSRMatrix B;
    B.rows = B.cols = n;
    B.row_ptr.resize(n + 1, 0);
    size_t nnz = 0;
    for (size_t i = 0; i < n; ++i) {
        std::sort(adj[i].begin(), adj[i].end());
        adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
        B.row_ptr[i + 1] = B.row_ptr[i] + adj[i].size();
        nnz += adj[i].size();
    }
    B.col_idx.reserve(nnz);
    for (const auto &lst : adj) {
        B.col_idx.insert(B.col_idx.end(), lst.begin(), lst.end());
    }
    B.values.assign(nnz, 1.0);
    return B;
}


/**
 * @brief For an n x n matrix A, builds the transposed (n) x (n) adjacency matrix.
 */
CSRMatrix build_transpose(const CSRMatrix &A) {
    size_t n = A.rows;
    std::vector<std::vector<size_t>> adj(n);

    for (size_t i = 0; i < n; ++i) {
        for (size_t k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
            size_t j = A.col_idx[k];
            adj[j].push_back(i);
        }
    }
    
    CSRMatrix B;
    B.rows = B.cols = n;
    B.row_ptr.resize(n + 1, 0);
    size_t nnz = 0;
    for (size_t i = 0; i < n; ++i) {
        B.row_ptr[i + 1] = B.row_ptr[i] + adj[i].size();
        nnz += adj[i].size();
    }
    B.col_idx.reserve(nnz);
    for (const auto &lst : adj) {
        B.col_idx.insert(B.col_idx.end(), lst.begin(), lst.end());
    }
    B.values.assign(nnz, 1.0);
    return B;
}


/**
 * @brief Finds a pseudo-peripheral node in the graph, a good starting point for RCM.
 * This is done by repeatedly performing a Breadth-First Search (BFS) to find a node
 * in the last level of the traversal, which approximates a node at the "edge" of the graph.
 * @param A The input sparse matrix (graph).
 * @param start_node An initial node to begin the search.
 * @return The index of a pseudo-peripheral node.
 */
size_t find_pseudo_peripheral_node(const CSRMatrix &A, size_t start_node) {
    int trial = 0;
    size_t u = start_node;
    std::vector<bool> visited(A.rows, false);
    size_t prev_clevel = 0;
    // Loop until the width of the BFS levels no longer increases
    while (true) {
        std::fill(visited.begin(), visited.end(), false);
        std::vector<size_t> current_level = {u};
        visited[u] = true;
        size_t clevel = 1;
        while (true) {
            std::vector<size_t> next_level;
            for (auto v : current_level) {
                for (size_t k = A.row_ptr[v]; k < A.row_ptr[v + 1]; ++k) {
                    size_t w = A.col_idx[k];
                    if (!visited[w]) {
                        visited[w] = true;
                        next_level.push_back(w);
                    }
                }
            }
            if (next_level.empty()) break;
            current_level = std::move(next_level);
            clevel++;
        }
        //std::cout << "u: " << u << "\tclevel: " << clevel << std::endl;
        if (clevel < prev_clevel || (clevel == prev_clevel && trial == 10)) break;

        if(clevel == prev_clevel) trial++;
        if(clevel > prev_clevel) trial = 0;
        prev_clevel = clevel;

        u = current_level[0];
        std::vector<size_t> out;

        std::sample(current_level.begin(), current_level.end(), std::back_inserter(out), 1, std::mt19937{std::random_device{}()});
        u = out[0];
    }
    return u;
}


void printRowNnzStats(const CSRMatrix &A) {
    size_t n = A.rows;
    // 1) Compute nonzeros per row
    std::vector<size_t> nnz_per_row(n);
    for (size_t i = 0; i < n; ++i) {
        nnz_per_row[i] = A.row_ptr[i+1] - A.row_ptr[i];
    }

    // 2) Total nonzeros
    size_t total_nnz = A.values.size();

    // 3) Min & Max
    auto [min_it, max_it] = std::minmax_element(nnz_per_row.begin(), nnz_per_row.end());
    size_t min_nnz = *min_it;
    size_t max_nnz = *max_it;

    // 4) Average
    double sum = std::accumulate(nnz_per_row.begin(), nnz_per_row.end(), 0ull);
    double avg = sum / static_cast<double>(n);

    // 5) Median
    std::sort(nnz_per_row.begin(), nnz_per_row.end());
    double median;
    if (n % 2 == 0) {
        median = (nnz_per_row[n/2 - 1] + nnz_per_row[n/2]) / 2.0;
    } else {
        median = nnz_per_row[n/2];
    }

    // 6) Standard Deviation (population)
    double sq_sum = 0.0;
    for (auto v : nnz_per_row) {
        double diff = static_cast<double>(v) - avg;
        sq_sum += diff * diff;
    }
    double stdev = std::sqrt(sq_sum / static_cast<double>(n));

    // 7) Output results
    std::cout << "Total nonzeros:          " << total_nnz << "\n";
    std::cout << "Min nonzeros per row:    " << min_nnz  << "\n";
    std::cout << "Max nonzeros per row:    " << max_nnz  << "\n";
    std::cout << "Avg nonzeros per row:    " << avg      << "\n";
    std::cout << "Median nonzeros per row: " << median   << "\n";
    std::cout << "Stdev nonzeros per row:  " << stdev    << "\n";
}

/**
 * @brief Permutes a matrix A given forward permutations, producing B = P_r A P_c^T.
 * @param row_P_fwd The forward row permutation vector (old->new).
 * @param col_P_fwd The forward column permutation vector (old->new).
 * @return The permuted matrix B.
 */
CSRMatrix permute(const CSRMatrix &A, const std::vector<size_t> &row_perm, const std::vector<size_t> &col_perm) {
    size_t m = A.rows;
    size_t n = A.cols;
    CSRMatrix B;
    B.rows = m; 
    B.cols = n;

    // Temporary storage for each row
    std::vector<std::vector<size_t>> new_cols(m);
    std::vector<std::vector<double>> new_vals(m);
    B.row_ptr.assign(m + 1, 0);

    // Inverse column permutation
    std::vector<size_t> inv_col_perm(n);
    for (size_t j = 0; j < n; j++) {
        inv_col_perm[col_perm[j]] = j;
    }

    // Place non-zero elements into their new positions
    for (size_t i_old = 0; i_old < m; ++i_old) {
        size_t i_new = row_perm[i_old];
        for (size_t k = A.row_ptr[i_new]; k < A.row_ptr[i_new + 1]; ++k) {
            size_t j_old = A.col_idx[k];
            size_t j_new = inv_col_perm[j_old];
            new_cols[i_old].push_back(j_new);
            new_vals[i_old].push_back(A.values[k]);  // preserve original value
        }
    }

    // Compute row_ptr
    size_t nnz = 0;
    for (size_t i = 0; i < m; ++i) {
        // Sort columns and values together
        std::vector<size_t>& cols = new_cols[i];
        std::vector<double>& vals = new_vals[i];
        std::vector<size_t> idx(cols.size());
        for (size_t k = 0; k < idx.size(); ++k) idx[k] = k;
        std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return cols[a] < cols[b]; });

        std::vector<size_t> sorted_cols(cols.size());
        std::vector<double> sorted_vals(vals.size());
        for (size_t k = 0; k < idx.size(); ++k) {
            sorted_cols[k] = cols[idx[k]];
            sorted_vals[k] = vals[idx[k]];
        }
        cols.swap(sorted_cols);
        vals.swap(sorted_vals);

        B.row_ptr[i + 1] = B.row_ptr[i] + cols.size();
        nnz += cols.size();
    }

    // Flatten columns and values
    B.col_idx.reserve(nnz);
    B.values.reserve(nnz);
    for (size_t i = 0; i < m; ++i) {
        B.col_idx.insert(B.col_idx.end(), new_cols[i].begin(), new_cols[i].end());
        B.values.insert(B.values.end(), new_vals[i].begin(), new_vals[i].end());
    }

    return B;
}


void printColNnzStats(const CSRMatrix &A) {
    size_t m = A.rows;
    size_t n = A.cols;

    // 1) Compute nonzeros per column
    std::vector<size_t> nnz_per_col(n, 0);
    for (size_t row = 0; row < m; ++row) {
        for (size_t ptr = A.row_ptr[row]; ptr < A.row_ptr[row+1]; ++ptr) {
            ++nnz_per_col[A.col_idx[ptr]];
        }
    }

    // 2) Total nonzeros
    size_t total_nnz = A.values.size();

    // 3) Min & Max
    auto [min_it, max_it] = std::minmax_element(nnz_per_col.begin(), nnz_per_col.end());
    size_t min_nnz = *min_it;
    size_t max_nnz = *max_it;

    // 4) Average
    double sum = std::accumulate(nnz_per_col.begin(), nnz_per_col.end(), 0ull);
    double avg = sum / static_cast<double>(n);

    // 5) Median
    std::sort(nnz_per_col.begin(), nnz_per_col.end());
    double median;
    if (n % 2 == 0) {
        median = (nnz_per_col[n/2 - 1] + nnz_per_col[n/2]) / 2.0;
    } else {
        median = nnz_per_col[n/2];
    }

    // 6) Standard Deviation (population)
    double sq_sum = 0.0;
    for (auto v : nnz_per_col) {
        double diff = static_cast<double>(v) - avg;
        sq_sum += diff * diff;
    }
    double stdev = std::sqrt(sq_sum / static_cast<double>(n));

    // 7) Output results
    std::cout << "Total nonzeros:              " << total_nnz << "\n";
    std::cout << "Min nonzeros per column:     " << min_nnz   << "\n";
    std::cout << "Max nonzeros per column:     " << max_nnz   << "\n";
    std::cout << "Avg nonzeros per column:     " << avg       << "\n";
    std::cout << "Median nonzeros per column:  " << median    << "\n";
    std::cout << "Stdev nonzeros per column:   " << stdev     << "\n";
}

struct PartitionResult {
    CSRMatrix row_matrix;    // NNZs from high-degree rows
    CSRMatrix col_matrix;    // NNZs from high-degree columns (excluding those already in row_matrix)
    CSRMatrix remaining_matrix; // All other NNZs
};

PartitionResult partitionCSR(const CSRMatrix& original, int mr, int mc) {
    // Calculate row degrees (number of non-zeros per row)
    vector<pair<size_t, size_t>> row_degrees; // (degree, row_index)
    for (size_t i = 0; i < original.rows; i++) {
        size_t degree = original.row_ptr[i + 1] - original.row_ptr[i];
        row_degrees.push_back({degree, i});
    }
    
    // Calculate column degrees
    vector<size_t> col_degree_count(original.cols, 0);
    for (size_t col : original.col_idx) {
        col_degree_count[col]++;
    }
    
    vector<pair<size_t, size_t>> col_degrees; // (degree, col_index)
    for (size_t j = 0; j < original.cols; j++) {
        col_degrees.push_back({col_degree_count[j], j});
    }
    
    // Sort by degree (descending)
    sort(row_degrees.rbegin(), row_degrees.rend());
    sort(col_degrees.rbegin(), col_degrees.rend());
    
    // Get top mr rows and top mc columns
    set<size_t> high_deg_rows, high_deg_cols;
    
    for (int i = 0; i < min(mr, (int)original.rows); i++) {
        high_deg_rows.insert(row_degrees[i].second);
    }
    
    for (int i = 0; i < min(mc, (int)original.cols); i++) {
        high_deg_cols.insert(col_degrees[i].second);
    }
    
    // Partition the non-zeros
    vector<vector<pair<size_t, double>>> row_matrix_data(original.rows);
    vector<vector<pair<size_t, double>>> col_matrix_data(original.rows);
    vector<vector<pair<size_t, double>>> remaining_matrix_data(original.rows);
    
    for (size_t i = 0; i < original.rows; i++) {
        for (size_t idx = original.row_ptr[i]; idx < original.row_ptr[i + 1]; idx++) {
            size_t j = original.col_idx[idx];
            double val = original.values[idx];
            
            if (high_deg_rows.count(i)) {
                // NNZ belongs to high-degree row
                row_matrix_data[i].push_back({j, val});
            } else if (high_deg_cols.count(j)) {
                // NNZ belongs to high-degree column (but not high-degree row)
                col_matrix_data[i].push_back({j, val});
            } else {
                // NNZ belongs to remaining matrix
                remaining_matrix_data[i].push_back({j, val});
            }
        }
    }
    
    // Helper function to build CSR matrix from row data
    auto buildCSRMatrix = [&](const vector<vector<pair<size_t, double>>>& data) -> CSRMatrix {
        CSRMatrix matrix;
        matrix.rows = original.rows;
        matrix.cols = original.cols;
        matrix.row_ptr.resize(original.rows + 1);
        
        size_t nnz_count = 0;
        matrix.row_ptr[0] = 0;
        
        for (size_t i = 0; i < original.rows; i++) {
            for (const auto& entry : data[i]) {
                matrix.col_idx.push_back(entry.first);
                matrix.values.push_back(entry.second);
                nnz_count++;
            }
            matrix.row_ptr[i + 1] = nnz_count;
        }
        
        return matrix;
    };
    
    PartitionResult result;
    result.row_matrix = buildCSRMatrix(row_matrix_data);
    result.col_matrix = buildCSRMatrix(col_matrix_data);
    result.remaining_matrix = buildCSRMatrix(remaining_matrix_data);
    
    return result;
}

// Example usage and verification function
void printMatrixInfo(const CSRMatrix& matrix, const string& name) {
    size_t nnz = matrix.values.size();
    cout << name << ": " << matrix.rows << "x" << matrix.cols 
         << " with " << nnz << " non-zeros" << endl;
}

// Verification function to ensure partition is valid
bool verifyPartition(const CSRMatrix& original, const PartitionResult& result) {
    size_t original_nnz = original.values.size();
    size_t total_partitioned = result.row_matrix.values.size() + 
                               result.col_matrix.values.size() + 
                               result.remaining_matrix.values.size();
    
    return original_nnz == total_partitioned;
}