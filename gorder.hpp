#include <queue>

using namespace std;

#include "csr.hpp"

vector<size_t> gorder(const CSRMatrix& A, const CSRMatrix& At, size_t window_size) {

    size_t start_vertex = find_pseudo_peripheral_node(A, 0);

    size_t n = A.rows;
    size_t min_unpermed = 0;
    
    vector<size_t> scores;          // scores (for reinsertion - no pq updates will be performed/can be more efficient?)
    scores.reserve(n);
    for(size_t i = 0; i < n; i++) scores[i] = 0;

    vector<size_t> perm;           // resulting ordering
    perm.reserve(n);
    
    vector<bool> permed;          // inverse permutation
    permed.reserve(n);
    for(size_t i = 0; i < n; i++) permed[i] = false;

    using P = pair<size_t, size_t>;    // (score, vertex)
    priority_queue<P> pq;              // max‐heap, lazy updates

    perm.push_back(start_vertex);
    permed[start_vertex] = true;

    //initialize the priority queue with incoming/outgoing neighbors of the start vertex
    for (size_t idx = A.row_ptr[start_vertex]; idx < A.row_ptr[start_vertex + 1]; ++idx) {
        size_t u = A.col_idx[idx]; //Matrix has the entry A[start_vertex, u]
        scores[u]++;
        pq.emplace(scores[u], u);
    }

    for (size_t idx = At.row_ptr[start_vertex]; idx < At.row_ptr[start_vertex + 1]; ++idx) {
        size_t u = At.col_idx[idx]; //Matrix has the entry [u, start_vertex]
        scores[u]++;
        pq.emplace(scores[u], u);
    }

    while (perm.size() < n) {
        size_t v = n;  
        while (!pq.empty()) {
            auto [s, u] = pq.top();
            pq.pop();
            //std::cout << u << " " << s << " " << scores[u] << " " << permed[u] << std::endl;

            if (!permed[u] && scores[u] == s) {
                v = u;
                //std::cout << "found a new vertex " << u << std::endl;
                break;
            }
        }
        //std::cout << v << std::endl;
        // if none in pq, pick the smallest-index unplaced vertex
        if (v == n) {
            for (; min_unpermed < n; ++min_unpermed) {
                if (!permed[min_unpermed]) {
                    v = min_unpermed;
                    break;
                }
            }
        } 
        //std::cout << "Disconnected: - continue from: " << v << std::endl;

        // place v
        perm.push_back(v);
        permed[v] = true;

        // add v’s contribution
        for (size_t idx = A.row_ptr[v]; idx < A.row_ptr[v + 1]; ++idx) {
            size_t u = A.col_idx[idx]; //Matrix has the entry A[start_vertex, u]
            if(!permed[u]) { 
                scores[u]++;
                pq.emplace(scores[u], u);
            }
         }

        for (size_t idx = At.row_ptr[v]; idx < At.row_ptr[v + 1]; ++idx) {
            size_t u = At.col_idx[idx]; //Matrix has the entry [u, start_vertex]
             if(!permed[u]) { 
                scores[u]++;
                pq.emplace(scores[u], u);
            }
        }
        
       if(perm.size() > window_size) {
            size_t r = perm[perm.size() - window_size];
            // remove r's contribution
            for (size_t idx = A.row_ptr[r]; idx < A.row_ptr[r + 1]; ++idx) {
                size_t u = A.col_idx[idx]; //Matrix has the entry A[start_vertex, u]
                if(!permed[u]) { 
                    scores[u]--;
                    pq.emplace(scores[u], u);
                }
            }

            for (size_t idx = At.row_ptr[r]; idx < At.row_ptr[r + 1]; ++idx) {
                size_t u = At.col_idx[idx]; //Matrix has the entry [u, start_vertex]
                if(!permed[u]) { 
                    scores[u]--;
                    pq.emplace(scores[u], u);
                }
            }
       }
    }
   /* for(auto x : perm) {
        std::cout << x << " ";
    }
    std::cout  << std::endl;*/
    return perm;
}