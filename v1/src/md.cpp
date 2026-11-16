#include "../include/matrix.hpp"
#include "../include/md.hpp"
#include <vector>
#include <limits>
#include <stdexcept>


MD_Graph::MD_Graph(const Matrix& matrix) : n(matrix.n), degrees(matrix.n, 0), edges(matrix.n), active(matrix.n, true) {
    for (std::size_t i = 0; i < matrix.n; i++) {
        for (std::size_t j = matrix.col_ptr[i]; j < matrix.col_ptr[i + 1]; j++) {
            if (matrix.row_indices[j] != i) {
                edges[i].insert(matrix.row_indices[j]);
                degrees[i]++;
            }
        }
    }
}


std::size_t MD_Graph::findMinimumDegree() const {
    std::size_t min_degree = std::numeric_limits<std::size_t>::max();
    std::size_t min_index = std::numeric_limits<std::size_t>::max();

    for (std::size_t i = 0; i < n; i++) {
        if (active[i] && degrees[i] < min_degree) {
            min_index = i;
            min_degree = degrees[i];
        }
    }

    return min_index;
}


void MD_Graph::eliminateVertex(std::size_t vertex) {
    active[vertex] = false;

    for (auto p = edges[vertex].begin(); p != edges[vertex].end(); p++) {
        edges[*p].erase(vertex);
        degrees[*p]--;

        for (auto q = std::next(p); q != edges[vertex].end(); q++) {
            if (edges[*p].insert(*q).second) {
                edges[*q].insert(*p);
                degrees[*p]++;
                degrees[*q]++;
            }
        }
    }

    edges[vertex].clear();
}


bool MD_Graph::isActive(std::size_t vertex) const {
    return active[vertex];
}


std::vector<std::size_t> MD(const Matrix& matrix) {
    MD_Graph graph = MD_Graph(matrix);

    std::vector<std::size_t> permutation(matrix.n);
    for (std::size_t i = 0; i < matrix.n; i++) {
        auto vertex = graph.findMinimumDegree();
        if (vertex == std::numeric_limits<std::size_t>::max()) {
            throw std::runtime_error("No active vertices found in the graph.");
        }

        graph.eliminateVertex(vertex);
        permutation[i] = vertex;
    }

    return permutation;
}