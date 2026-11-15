#pragma once

#include <vector>
#include <cstddef>
#include <limits>
#include "matrix.hpp"

// Approximate Minimum Degree ordering (Amestoy, Davis, Duff) on a quotient graph.
class AMD_Graph {
private:
    enum class Status { VARIABLE, ELEMENT, MERGED, ELIMINATED, DEAD_ELEMENT };

    static constexpr std::size_t NONE = std::numeric_limits<std::size_t>::max();

    std::size_t n;
    std::size_t eliminated;
    std::size_t min_degree;
    std::size_t stamp;

    std::vector<Status> status;
    std::vector<std::size_t> nv;
    std::vector<std::size_t> degrees;

    // Quotient graph: variable i keeps adjacent variables in var_adj[i] and adjacent
    // elements in elem_adj[i]; element e keeps its variables in pattern[e].
    std::vector<std::vector<std::size_t>> var_adj;
    std::vector<std::vector<std::size_t>> elem_adj;
    std::vector<std::vector<std::size_t>> pattern;
    std::vector<std::size_t> pattern_weight;

    // Degree buckets (doubly linked lists indexed by approximate degree).
    std::vector<std::size_t> bucket_head;
    std::vector<std::size_t> bucket_next;
    std::vector<std::size_t> bucket_prev;

    // Members of a supervariable: linked list starting at the principal variable.
    std::vector<std::size_t> member_next;
    std::vector<std::size_t> member_tail;

    // Work arrays that are valid only for the current stamp.
    std::vector<std::size_t> mark;
    std::vector<std::size_t> elem_stamp;
    std::vector<std::size_t> elem_external;
    std::vector<std::size_t> outside_weight;

    void bucketInsert(std::size_t vertex);
    void bucketRemove(std::size_t vertex);
    void mergeSupervariables(std::vector<std::size_t>& lp);

public:
    AMD_Graph(const Matrix& matrix);
    std::size_t selectPivot();
    void eliminatePivot(std::size_t pivot, std::vector<std::size_t>& order);
    bool isDone() const;
};

std::vector<std::size_t> AMD(const Matrix& matrix);
