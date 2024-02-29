#ifndef LINEAR_ALGEBRA_HEADER
#define LINEAR_ALGEBRA_HEADER
#include <iostream>
#include <cassert>

#include <linear_algebra/vector.hpp>
#include <linear_algebra/matrix.hpp>
#include <linear_algebra/vector_support.hpp>
#include <tools/permutation.hpp>

template <typename T_TYPE_BLOCK, size_t T_ROWS, size_t T_COLS>
bool partial_echelon_form(t_matrix_row<T_TYPE_BLOCK,T_ROWS,T_COLS>& A,size_t nb_cols, size_t& stop){
    assert(nb_cols <= A.rows() && nb_cols <= A.cols());
    t_vector<T_TYPE_BLOCK,T_COLS> tmp;
    size_t current_col = 0;
    while (current_col < nb_cols){
        size_t j = current_col;
        stop = current_col;
        while((j < A.rows()) && !A[j][current_col]){
            j += 1;
        }
        if (j == A.rows()){
            return false;
        }
        if(j != current_col){

            A[current_col] ^= A[j];
        }
        for(size_t z = 0; z < A.rows(); ++z){
            if((z != current_col) && A[z][current_col]){
                    A[z] ^= A[current_col];
            }
        }
        current_col +=1;
        stop = current_col;
    }
    return true;
}

#endif