#include "../include/matrix.hpp"
#include "../include/md.hpp"
#include <iostream>

int main() {
    Matrix matrix = Matrix(4);
    matrix.values = 
    {4, 1,
     1, 4, 1,
        1, 4, 1,
           1, 4};
    matrix.row_indices = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3};
    matrix.col_ptr = {0, 2, 5, 8, 10};

    auto permutation = MD(matrix);
    for(auto& e : permutation) {
        std::cout << e << " ";
    }

    return 0;
}