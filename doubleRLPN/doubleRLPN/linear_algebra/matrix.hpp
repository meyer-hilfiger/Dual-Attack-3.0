#ifndef MATRIX_HEADER
#define MATRIX_HEADER

#include <iostream>
#include <cassert>
#include <random>

#include <linear_algebra/vector.hpp>
#include <linear_algebra/vector_support.hpp>
#include <tools/permutation.hpp>


#include <cmath>


template <typename T_TYPE_BLOCK, size_t T_ROWS, size_t T_COLS, bool T_COL_REPR>
struct t_matrix
{

    using TYPE_VECTOR = t_vector<T_TYPE_BLOCK,T_ROWS>;
    static constexpr size_t COLS = T_COLS;
    static constexpr size_t ROWS = T_ROWS;

    TYPE_VECTOR& operator[](size_t c){
         return mat[c];
    }

    t_matrix() : mat(){

    }
    template <size_t T_START, size_t T_V_PADDING, typename T_TYPE_BLOCK_1,size_t T_ROWS_1, size_t T_COLS_1>
    void extract_last(t_matrix<T_TYPE_BLOCK_1,T_ROWS_1,T_COLS_1,T_COL_REPR>& m){
        for (int i = T_START; i < T_COLS_1 + T_START ; ++i){
            mat[i].template extract_last<T_V_PADDING>(m[i-T_START]);
        }
    }

    void init(){
        for(size_t i = 0; i < T_COLS; ++i){
            mat[i].init();
        }
    }

    void print(){
        bool res;
        for(int i= 0; i < T_ROWS; ++i){
            for(int j = 0; j < T_COLS; ++j){
                res = mat[j][i];
                std::cout << (uint16_t) res;
            }
            std::cout << std::endl;
        }
    }

    // Set bit to value v in row i and column j
    void set(int i, int j, bool v){
        if constexpr(T_COL_REPR){
            mat[j].set(i,v);
        }else{
            mat[i].set(j,v);
        }
    }
    // Get value in row i and column j
    bool get(int i, int j){
        if constexpr(T_COL_REPR){
            return mat[j].get(i);
        }else{
            return mat[i].get(j);
        }
    }

    // Switch from column to row representation (or vis-versa)
    t_matrix<T_TYPE_BLOCK,T_COLS,T_ROWS,!T_COL_REPR> inverse_representation(){
        t_matrix<T_TYPE_BLOCK,T_COLS,T_ROWS,!T_COL_REPR> M;
        for(size_t i = 0; i < T_ROWS; ++i){
            for(size_t j = 0; j < T_COLS; ++j){
                M[i].set(j,mat[j][i]);
            }
        }
        return M;
    }

    size_t rows(){
        if constexpr (T_COL_REPR){
            return T_ROWS;
        }else{
            return T_COLS;
        }
    }

    size_t cols(){
        if constexpr (!T_COL_REPR){
            return T_ROWS;
        }else{
            return T_COLS;
        }
    }

    void rand(){
        for(size_t i = 0; i < T_COLS; ++i){
            mat[i].rand();
            
        }
    }

    t_matrix permute_columns(permutation<T_COL_REPR ? T_COLS : T_ROWS >& p){
        if constexpr(T_COL_REPR){
            t_matrix mat_permuted;
            for(size_t i= 0; i < T_COLS; ++i){
                mat_permuted.mat[p[i]] = mat[i];
            }
            return mat_permuted;      
        }else{
            t_matrix<T_TYPE_BLOCK,T_COLS,T_ROWS,!T_COL_REPR> A_col = this->inverse_representation();
            t_matrix<T_TYPE_BLOCK,T_COLS,T_ROWS,!T_COL_REPR> A_col_permuted = A_col.permute_columns(p);
            return A_col_permuted.inverse_representation();

        }

    }

    template <typename TYPE_BLOCK_V_MUL>
    TYPE_VECTOR operator*(t_vector<TYPE_BLOCK_V_MUL,T_COLS>& v){
        TYPE_VECTOR mul;
        for(size_t i = 0; i < T_COLS; ++i){
            if(v[i]){
                mul ^= mat[i];
            }
        }
        return mul;
    }



    template <size_t T_P>
    TYPE_VECTOR operator*(vector_support<T_P>& p){
        TYPE_VECTOR mul;
        for(size_t i = 0; i < T_P; ++i){
                mul ^= mat[p[i]];
        }
        return mul;
    }

    TYPE_VECTOR mat[T_COLS];

};

#endif






