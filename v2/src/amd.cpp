#include "../include/matrix.hpp"
#include "../include/amd.hpp"
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>


AMD_Graph::AMD_Graph(const Matrix& matrix)
    : n(matrix.n), eliminated(0), min_degree(0), stamp(0),
      status(matrix.n, Status::VARIABLE), nv(matrix.n, 1), degrees(matrix.n, 0),
      var_adj(matrix.n), elem_adj(matrix.n), pattern(matrix.n), pattern_weight(matrix.n, 0),
      bucket_head(matrix.n, NONE), bucket_next(matrix.n, NONE), bucket_prev(matrix.n, NONE),
      member_next(matrix.n, NONE), member_tail(matrix.n, NONE),
      mark(matrix.n, 0), elem_stamp(matrix.n, 0), elem_external(matrix.n, 0), outside_weight(matrix.n, 0) {
    for (std::size_t c = 0; c < n; c++) {
        for (auto idx = matrix.col_ptr[c]; idx < matrix.col_ptr[c + 1]; idx++) {
            auto r = matrix.row_indices[idx];
            if (r != c) {
                var_adj[c].push_back(r);
                var_adj[r].push_back(c);
            }
        }
    }

    for (std::size_t i = 0; i < n; i++) {
        std::sort(var_adj[i].begin(), var_adj[i].end());
        var_adj[i].erase(std::unique(var_adj[i].begin(), var_adj[i].end()), var_adj[i].end());
        degrees[i] = var_adj[i].size();
        member_tail[i] = i;
        bucketInsert(i);
    }
}


void AMD_Graph::bucketInsert(std::size_t vertex) {
    auto d = degrees[vertex];
    bucket_prev[vertex] = NONE;
    bucket_next[vertex] = bucket_head[d];
    if (bucket_head[d] != NONE) {
        bucket_prev[bucket_head[d]] = vertex;
    }
    bucket_head[d] = vertex;
}


void AMD_Graph::bucketRemove(std::size_t vertex) {
    auto prev = bucket_prev[vertex];
    auto next = bucket_next[vertex];
    if (prev != NONE) {
        bucket_next[prev] = next;
    } else {
        bucket_head[degrees[vertex]] = next;
    }
    if (next != NONE) {
        bucket_prev[next] = prev;
    }
    bucket_prev[vertex] = NONE;
    bucket_next[vertex] = NONE;
}


std::size_t AMD_Graph::selectPivot() {
    while (min_degree < n && bucket_head[min_degree] == NONE) {
        min_degree++;
    }
    if (min_degree >= n) {
        throw std::runtime_error("No active vertices found in the graph.");
    }

    return bucket_head[min_degree];
}


bool AMD_Graph::isDone() const {
    return eliminated >= n;
}


void AMD_Graph::mergeSupervariables(std::vector<std::size_t>& lp) {
    std::unordered_map<std::size_t, std::vector<std::size_t>> buckets;
    for (auto i : lp) {
        std::size_t hash = 0;
        for (auto j : var_adj[i]) {
            hash += j;
        }
        for (auto e : elem_adj[i]) {
            hash += e;
        }
        buckets[hash].push_back(i);
    }

    for (auto& entry : buckets) {
        auto& group = entry.second;
        for (std::size_t a = 0; a + 1 < group.size(); a++) {
            auto i = group[a];
            if (status[i] != Status::VARIABLE) {
                continue;
            }

            auto i_stamp = ++stamp;
            for (auto j : var_adj[i]) {
                mark[j] = i_stamp;
            }
            for (auto e : elem_adj[i]) {
                mark[e] = i_stamp;
            }

            for (std::size_t b = a + 1; b < group.size(); b++) {
                auto j = group[b];
                if (status[j] != Status::VARIABLE) {
                    continue;
                }
                if (var_adj[j].size() != var_adj[i].size() || elem_adj[j].size() != elem_adj[i].size()) {
                    continue;
                }

                auto same = true;
                for (auto k : var_adj[j]) {
                    if (mark[k] != i_stamp) {
                        same = false;
                        break;
                    }
                }
                for (auto e : elem_adj[j]) {
                    if (mark[e] != i_stamp) {
                        same = false;
                        break;
                    }
                }
                if (!same) {
                    continue;
                }

                nv[i] += nv[j];
                nv[j] = 0;
                status[j] = Status::MERGED;
                member_next[member_tail[i]] = j;
                member_tail[i] = member_tail[j];
                std::vector<std::size_t>().swap(var_adj[j]);
                std::vector<std::size_t>().swap(elem_adj[j]);
            }
        }
    }

    lp.erase(std::remove_if(lp.begin(), lp.end(), [this](std::size_t v) {
        return status[v] != Status::VARIABLE;
    }), lp.end());
}


void AMD_Graph::eliminatePivot(std::size_t pivot, std::vector<std::size_t>& order) {
    bucketRemove(pivot);

    // Lp = (A_pivot U union of L_e for e in E_pivot) \ {pivot}; the elements in E_pivot are absorbed.
    auto lp_stamp = ++stamp;
    mark[pivot] = lp_stamp;
    std::vector<std::size_t> lp;
    for (auto j : var_adj[pivot]) {
        if (status[j] == Status::VARIABLE && mark[j] != lp_stamp) {
            mark[j] = lp_stamp;
            lp.push_back(j);
        }
    }
    for (auto e : elem_adj[pivot]) {
        if (status[e] != Status::ELEMENT) {
            continue;
        }
        for (auto j : pattern[e]) {
            if (status[j] == Status::VARIABLE && mark[j] != lp_stamp) {
                mark[j] = lp_stamp;
                lp.push_back(j);
            }
        }
        status[e] = Status::DEAD_ELEMENT;
        std::vector<std::size_t>().swap(pattern[e]);
    }
    std::vector<std::size_t>().swap(var_adj[pivot]);
    std::vector<std::size_t>().swap(elem_adj[pivot]);

    status[pivot] = Status::ELEMENT;
    eliminated += nv[pivot];
    nv[pivot] = 0;
    for (auto v = pivot; v != NONE; v = member_next[v]) {
        order.push_back(v);
    }

    for (auto i : lp) {
        bucketRemove(i);
    }

    // |L_e \ Lp| for every element adjacent to Lp, stored in elem_external.
    for (auto i : lp) {
        for (auto e : elem_adj[i]) {
            if (status[e] != Status::ELEMENT) {
                continue;
            }
            if (elem_stamp[e] != lp_stamp) {
                elem_stamp[e] = lp_stamp;
                elem_external[e] = pattern_weight[e];
            }
            elem_external[e] -= nv[i];
        }
    }

    // Clean the lists of every variable in Lp, absorb elements contained in Lp (aggressive
    // absorption) and eliminate the variables that became indistinguishable from the pivot.
    std::vector<std::size_t> survivors;
    survivors.reserve(lp.size());
    for (auto i : lp) {
        std::size_t outside = 0;

        auto& elements = elem_adj[i];
        std::size_t kept = 0;
        for (auto e : elements) {
            if (status[e] != Status::ELEMENT) {
                continue;
            }
            if (elem_external[e] == 0) {
                status[e] = Status::DEAD_ELEMENT;
                std::vector<std::size_t>().swap(pattern[e]);
                continue;
            }
            elements[kept++] = e;
            outside += elem_external[e];
        }
        elements.resize(kept);

        auto& neighbors = var_adj[i];
        kept = 0;
        for (auto j : neighbors) {
            if (status[j] != Status::VARIABLE || mark[j] == lp_stamp) {
                continue;
            }
            neighbors[kept++] = j;
            outside += nv[j];
        }
        neighbors.resize(kept);

        if (elements.empty() && neighbors.empty()) {
            status[i] = Status::ELIMINATED;
            eliminated += nv[i];
            nv[i] = 0;
            for (auto v = i; v != NONE; v = member_next[v]) {
                order.push_back(v);
            }
            continue;
        }

        elements.push_back(pivot);
        outside_weight[i] = outside;
        survivors.push_back(i);
    }

    mergeSupervariables(survivors);

    std::size_t pivot_weight = 0;
    for (auto i : survivors) {
        pivot_weight += nv[i];
    }

    // Approximate external degree of every variable in the new element.
    auto remaining = n - eliminated;
    for (auto i : survivors) {
        auto external = pivot_weight - nv[i];
        auto d = std::min({remaining - nv[i], degrees[i] + external, outside_weight[i] + external});
        degrees[i] = d;
        bucketInsert(i);
        min_degree = std::min(min_degree, d);
    }

    if (survivors.empty()) {
        status[pivot] = Status::DEAD_ELEMENT;
    } else {
        pattern[pivot] = std::move(survivors);
        pattern_weight[pivot] = pivot_weight;
    }
}


std::vector<std::size_t> AMD(const Matrix& matrix) {
    AMD_Graph graph = AMD_Graph(matrix);

    std::vector<std::size_t> permutation;
    permutation.reserve(matrix.n);
    while (!graph.isDone()) {
        auto pivot = graph.selectPivot();
        graph.eliminatePivot(pivot, permutation);
    }

    return permutation;
}
