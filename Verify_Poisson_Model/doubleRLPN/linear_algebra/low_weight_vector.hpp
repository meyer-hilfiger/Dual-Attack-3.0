#ifndef HEADER_LOWWEIGHTVECTOR
#define HEADER_LOWWEIGHTVECTOR

#include <linear_algebra/linear_algebra.hpp>

// low_weight_vector is the representation of a word v of F_2^n of low weight T_w on the part N.

// v_P (v restricted to P) is represented by a vector of F_2^{T_s} of type t_vector
// v_N (v restricted to N) is represented by its support (of size T_w)
template <size_t T_n, size_t T_s, size_t T_w>
struct low_weight_vector{
	using type_word_P = t_vector<uint64_t,T_s>;
	using type_word_N = vector_support<T_w>;
	using type_word = t_vector<uint64_t,T_n>;

	type_word_P word_P;
	type_word_N word_N;

	low_weight_vector(type_word_P _word_P,type_word_N _word_N ) : word_P(_word_P), word_N(_word_N){
	}

	// Convert this representation into the usual vector representation
	type_word to_vector(){
		type_word ret;
		ret.set_first_coordinates(word_P);
		for(size_t i = 0; i < T_w; ++i){
			ret.set(word_N[i] + T_s,1);
		}
		return ret;
	}
};

#endif