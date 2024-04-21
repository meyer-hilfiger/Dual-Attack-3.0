#ifndef HEADER_PERMUTATION
#define HEADER_PERMUTATION
#include <iostream>

template <int T_N>
struct permutation{
    
  uint16_t permu[T_N];
  static constexpr int P = T_N;

    void init(){
        for(int i = 0; i < T_N; ++i){
            permu[i] = i;
        }
    }
  permutation(){
      init();
  };

  uint16_t& operator[](size_t i) {
    return permu[i];
  }

  uint16_t operator[](size_t i) const {
    return permu[i];
  }
  void rand(){
    init();
    int index;
    int temp;
    for(int i = 0; i< T_N; ++i){
        int tmp = std::rand();
        index = tmp % (T_N-i);
        temp = permu[T_N-i-1];
        permu[T_N-i-1] = permu[index];
        permu[index] = temp;
    }
  }
  permutation inverse(){
    permutation p_inv;
    for(size_t i = 0; i < T_N; ++i){
      p_inv[permu[i]] = i;
    }
    return p_inv;
  }

  permutation composition(permutation& p){
    permutation ret;
    for(size_t i = 0; i < T_N; ++i){
      ret[i] = (*this)[p[i]];
    }
    return ret;
  }

    void print(){
    for(size_t i = 0; i < T_N; ++i){
      std::cout << permu[i] << "   ";
    }
    std::cout << std::endl;
  }




};



#endif