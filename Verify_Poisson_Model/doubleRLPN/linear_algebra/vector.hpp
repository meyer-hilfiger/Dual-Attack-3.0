#ifndef VECTOR_HEADER
#define VECTOR_HEADER

#include <iostream>
#include <cassert>
#include <random>


#include <linear_algebra/vector_support.hpp>
#include <tools/permutation.hpp>

#include <cmath>
//using namespace std;

template <typename T_TYPE_BLOCK, size_t T_ROWS, size_t T_COLS, bool T_COL_REPR>
struct t_matrix;

template <typename T_TYPE_BLOCK, size_t T_ROWS, size_t T_COLS>
using t_matrix_row = t_matrix<T_TYPE_BLOCK,T_COLS,T_ROWS,false>;

template <typename T_TYPE_BLOCK, size_t T_ROWS, size_t T_COLS>
using t_matrix_col = t_matrix<T_TYPE_BLOCK,T_ROWS,T_COLS,true>;

inline int popcount(uint8_t x){
  return __builtin_popcount(x);
}
inline int popcount(uint16_t x){
  return __builtin_popcount(x);
}
inline int popcount(uint64_t x){
  return __builtin_popcountll(x);
}
inline int popcount(uint32_t x){
  return __builtin_popcountl(x);
}

// Compute the number of required blocks of type T_TYPE_BLOCK to store a binary vector of length n
template<typename T_TYPE_BLOCK>
constexpr size_t nb_Blocks(size_t n){
    if (n % (sizeof(T_TYPE_BLOCK) * 8) == 0){
        return n / (sizeof(T_TYPE_BLOCK) * 8);
    }
    else{
        return 1 + n / (sizeof(T_TYPE_BLOCK) * 8);
    }
}

template<typename T_TYPE_BLOCK>
constexpr size_t offset(size_t n){
    return ((n-1) % (sizeof(T_TYPE_BLOCK)*8));
}



template <typename T_TYPE_BLOCK, size_t T_LEN>
struct t_vector{
    using TYPE_BLOCK = T_TYPE_BLOCK;
    static constexpr size_t nb_block = nb_Blocks<T_TYPE_BLOCK>(T_LEN);

    public:
    size_t get_k(size_t i) const {
        const int off = sizeof(T_TYPE_BLOCK)*8 - (T_LEN % (sizeof(T_TYPE_BLOCK)*8));
        if (0 == (T_LEN % (sizeof(T_TYPE_BLOCK)*8))){
            return i/(sizeof(T_TYPE_BLOCK)*8);
        }else{
            return (off + i)/(sizeof(T_TYPE_BLOCK)*8);
        }
    }

    size_t get_mod(size_t i) const {
        const int off = (sizeof(T_TYPE_BLOCK)*8) - (T_LEN % (sizeof(T_TYPE_BLOCK)*8));
        return ((off+i) % (sizeof(T_TYPE_BLOCK)*8));
    }

    void init(){
        for(size_t i = 0; i < nb_block; ++i){
            vect[i] = 0;
        }
    }
    t_vector() : vect(){
        init();
    }
    t_vector(const t_vector& v) : vect() {
        for(size_t i = 0; i < nb_block; ++i){
            vect[i] = 0;
            vect[i] ^= v.vect[i];
        }
    }

    // Create vector of hamming weight T_P with support given by p
    template <size_t T_P>
    t_vector(const vector_support<T_P>& p) : vect(){
        init();
        for(size_t i = 0; i < T_P; ++i){
            set(p[i],true);
        }
    }
    size_t len(){
        return T_LEN;
    }
    

    t_vector& operator^=(const t_vector& right){
        for(size_t i = 0; i < nb_block; ++i){
            this->vect[i] ^= right.vect[i];
        }
        return *this;
    }

    t_vector operator^(const t_vector& right){
        t_vector v;
        for(size_t i = 0; i < nb_block; ++i){
            v.vect[i] = this->vect[i] ^ right.vect[i];
        }
        return v;
    }

    t_vector operator&(const t_vector& right){
        t_vector v;
        for(size_t i = 0; i < nb_block; ++i){
            v.vect[i] = this->vect[i] & right.vect[i];
        }
        return v;
    }


    t_vector& operator=(const t_vector& right){
        for(size_t i = 0; i < nb_block; ++i){
            this->vect[i] = right.vect[i];
        }
        return *this;
    }

    template <size_t T_P>
    t_vector& operator=(const vector_support<T_P>& p) {
        init();
        for(size_t i = 0; i < T_P; ++i){
            set(p[i],true);
        }
        return (*this);
    }

    T_TYPE_BLOCK get_block (size_t i) const{
        return vect[i];
    }

    void set(size_t i, bool value){
        int k = get_k(i);
        int pos = get_mod(i);

        vect[k] &= ~((T_TYPE_BLOCK)1 << (sizeof(T_TYPE_BLOCK)*8 - 1 - pos));
        if (value){
            vect[k] |= ((T_TYPE_BLOCK)1 << (sizeof(T_TYPE_BLOCK)*8 - 1 - pos));
        }
    }

    void invert(size_t i){
        int k = get_k(i);
        int pos = get_mod(i);
        vect[k] ^= ((T_TYPE_BLOCK)1 << (sizeof(T_TYPE_BLOCK)*8 - 1 - pos));
    }
    
    size_t hamming(){
        size_t res = 0;
        for(size_t i = 0; i < nb_block; ++i){
            res += popcount(vect[i]);
        }
        return res;
    }

    bool dot_product(const t_vector& right){
        size_t res = ((*this) & right).hamming();
        return res % 2;
    }
    bool operator[](size_t i) const {
        size_t k = get_k(i);
        size_t m = get_mod(i);
        T_TYPE_BLOCK mask = ((T_TYPE_BLOCK) 1 << (sizeof(T_TYPE_BLOCK)*8 - 1 - m));
        T_TYPE_BLOCK c = get_block(k);
        return (( c & mask) >> (sizeof(T_TYPE_BLOCK)*8 - 1 - m));
    }


    bool get(size_t i) const{
        return (*this)[i];
    }
    bool operator==(const t_vector<T_TYPE_BLOCK,T_LEN> &other) const {
    for(int i = 0; i < nb_block; ++i){
        if(vect[i] != other.vect[i]){
            return false;
        }
    }
        return true;
    }
    
    bool operator!=(const t_vector<T_TYPE_BLOCK,T_LEN>  &other) const {
    return !(*this == other);
    }
    
    bool operator<(const t_vector<T_TYPE_BLOCK,T_LEN>  &other) const {
    for(int i = 0; i < nb_block; ++i){
        if(vect[i] < other.vect[i]){
            return true;
        }else if (vect[i] > other.vect[i]){
            return false;
        }
    }
        return false;
    }
    
    bool operator>(const t_vector<T_TYPE_BLOCK,T_LEN> &other) const {
    for(int i = 0; i < nb_block; ++i){
        if(vect[i] > other.vect[i]){
            return true;
        }else if (vect[i] < other.vect[i]){
            return false;
        }
    }
        return false;
    }
    bool operator<=(const t_vector<T_TYPE_BLOCK,T_LEN>  &other) const {
    return !(*this > other);
    }
    bool operator>=(const t_vector<T_TYPE_BLOCK,T_LEN>  &other) const {
    return !(*this < other);
    }

    template <size_t T_ROWS>
    t_vector<T_TYPE_BLOCK,T_ROWS> operator*(const t_matrix_row<T_TYPE_BLOCK,T_LEN,T_ROWS>& m){
        t_vector<T_TYPE_BLOCK,T_ROWS> res;
        for(size_t i = 0; i < T_LEN; ++i){
            if((*this)[i]){
                res ^= m.mat[i];
            }
        }
        return res;
    }
    
    uint64_t to_integer(){
        if constexpr (nb_block != 1){
            std::cout << "Vector as integer is not implemented for vectors with more than one block" << std::endl;
            exit(0);
        }else{
            return vect[0];
        }
        
    }
    void integer_to_vector(T_TYPE_BLOCK b){
        if constexpr (nb_block != 1){
            std::cout << "Integer to vector is not implemented for vectors with more than one block" << std::endl;
            exit(0);
        }else{
            vect[0] = b;
        }
    }
   
    // SLOW TO -> to be bitsliced
    // Extract last L bits of this into v
    template <size_t T_PADDING, typename T_TYPE_BLOCK_1,size_t L>
    void extract_last(t_vector<T_TYPE_BLOCK_1,L>& v){
        v.init();
        size_t start = T_LEN - L - T_PADDING;
        for(size_t i = 0; i < L; ++i){
            v.set(i,(*this)[start+i]);
        }
    }
    template <typename T_TYPE_BLOCK_1,size_t L>
    void set_first_coordinates(t_vector<T_TYPE_BLOCK_1,L>& v){
        for(size_t i = 0; i < L; ++i){
            this->set(i,v[i]);
        }
    }
    void print(){
        for (int i = 0; i < T_LEN ; ++i){
            std::cout << (uint64_t) (*this) [i]; //<< endl;
        }
        std::cout << std::endl;
    }

    void rand(){
        size_t j = this->get_mod(0);
        for(size_t i = 0; i < T_LEN; ++i){
            set(i,std::rand()%2);
        }
        vect[0] = vect[0] << j;
        vect[0] = vect[0] >> j;

    }

    void rand_weight(size_t t){
        init();
        for(size_t i = 0; i < t; ++i){
            set(i,1);
        }
        permutation<T_LEN> per;
        per.rand();
        permute(per);
    }

    t_vector permute(permutation<T_LEN>& p){
        t_vector tmp;
        for(size_t i= 0; i < T_LEN; ++i){
            tmp.set(p[i],(*this)[i]);
        }
        return tmp;
    }

    T_TYPE_BLOCK vect[nb_Blocks<T_TYPE_BLOCK>(T_LEN)];


};

#endif






