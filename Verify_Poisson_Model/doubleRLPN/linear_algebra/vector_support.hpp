

#ifndef HEADER_VECTOR_SUPPORT
#define HEADER_VECTOR_SUPPORT
#include <algorithm>

// Set of position of the support of a vector (of low weight)
template <size_t T_P>
struct vector_support;

template <>
struct vector_support<0>{
  static constexpr size_t T_P = 0;
  void init(){
    
  }

  vector_support(){

  }

  template <size_t T_P1>
  vector_support<T_P + T_P1> operator+(vector_support<T_P1>& p){
    vector_support<T_P + T_P1> new_p;
    for(size_t i = 0; i < T_P; ++i){
      new_p[i] = permu[i];
    }
    for(size_t i = 0; i < T_P1; ++i){
      new_p[i+T_P] = p.permu[i];
    }
    return new_p;
  }

  void print(){
   
  }

  void sort(){
   
  }
  uint16_t& operator[](size_t i) {
    return permu[i];
  }

  uint16_t operator[](size_t i) const {
    return permu[i];
  }


  uint16_t permu[T_P];

  static constexpr size_t P = T_P;
};




template <size_t T_P>
struct vector_support{

  void init(){
    for(int i = 0; i < T_P; ++i){
      permu[i] = 0;
    }
  }

  vector_support(){
    init();
  }

  template <size_t T_P1>
  vector_support<T_P + T_P1> operator+(vector_support<T_P1>& p){
    vector_support<T_P + T_P1> new_p;
    for(size_t i = 0; i < T_P; ++i){
      new_p[i] = permu[i];
    }
    for(size_t i = 0; i < T_P1; ++i){
      new_p[i+T_P] = p.permu[i];
    }
    return new_p;
  }

  void print(){
    for(size_t i = 0; i < T_P; ++i){
      std::cout << permu[i] << "   ";
    }
    std::cout << std::endl;
  }

  void sort(){
    std::sort(permu,permu+T_P);
  }
  uint16_t& operator[](size_t i) {
    return permu[i];
  }

  uint16_t operator[](size_t i) const {
    return permu[i];
  }


  uint16_t permu[T_P];

  static constexpr size_t P = T_P;
};

#endif