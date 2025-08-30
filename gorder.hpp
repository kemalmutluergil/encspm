#pragma once
#include <queue>

using namespace std;

#include "csr.hpp"

/**
 * @brief GOrder reordering (Shun et al., ASPLOS'15 simplified version)
 * 
 * @param A CSRMatrix (assumed square, adjacency-like)
 * @param row_perm Output: row forward permutation (old -> new)
 * @param col_perm Output: col forward permutation (old -> new)
 * @param window_size Sliding window size (default = 50)
 */
void gorder(const CSRMatrix& A,
            vector<size_t>& row_perm,
            vector<size_t>& col_perm,
            size_t window_size = 50) {
    
    size_t n = A.rows;
    row_perm.assign(n, n);
    col_perm.assign(n, n);

    // --- Step 1: build adjacency ---
    vector<vector<size_t>> adj(n);
    for (size_t i = 0; i < n; i++) {
        for (size_t k = A.row_ptr[i]; k < A.row_ptr[i+1]; k++) {
            size_t j = A.col_idx[k];
            if (i != j) {
                adj[i].push_back(j);
                adj[j].push_back(i); // assume undirected
            }
        }
    }

    // --- Step 2: initial ordering by degree (descending) ---
    vector<size_t> vertices(n);
    iota(vertices.begin(), vertices.end(), 0);
    sort(vertices.begin(), vertices.end(),
         [&](size_t a, size_t b) { return adj[a].size() > adj[b].size(); });

    // --- Step 3: sliding window placement ---
    vector<size_t> ordering;
    ordering.reserve(n);
    vector<bool> placed(n, false);

    for (size_t v : vertices) {
        if (placed[v]) continue;

        // place this vertex
        ordering.push_back(v);
        placed[v] = true;

        // maintain a queue as a sliding window
        deque<size_t> window;
        window.push_back(v);

        while (!window.empty()) {
            size_t u = window.front();
            window.pop_front();

            // sort neighbors by degree, descending
            vector<size_t> neigh;
            for (auto w : adj[u]) {
                if (!placed[w]) neigh.push_back(w);
            }
            sort(neigh.begin(), neigh.end(),
                 [&](size_t a, size_t b){ return adj[a].size() > adj[b].size(); });

            for (auto w : neigh) {
                if (!placed[w]) {
                    ordering.push_back(w);
                    placed[w] = true;
                    window.push_back(w);

                    if (window.size() > window_size)
                        window.pop_front(); // keep window size fixed
                }
            }
        }
    }

    // --- Step 4: map old -> new ---
    for (size_t new_idx = 0; new_idx < n; new_idx++) {
        size_t old_idx = ordering[new_idx];
        row_perm[old_idx] = new_idx;
        col_perm[old_idx] = new_idx; // same if symmetric
    }
}
