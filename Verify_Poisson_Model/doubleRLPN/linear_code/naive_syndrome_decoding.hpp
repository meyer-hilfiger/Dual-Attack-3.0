#ifndef HEADER_NAIVESYNDROMEDECODING
#define HEADER_NAIVESYNDROMEDECODING

#include <linear_algebra/linear_algebra.hpp>
#include <tools/enumerate.hpp>

// Solve syndrome decoding problem naively
template <size_t T_n, size_t T_k, size_t T_w>
struct naive_syndrome_decoding{
	using type_syndrome = t_vector<uint64_t,T_n - T_k>;
	using type_word = vector_support<T_w>;
	using type_parity_check_matrix = t_matrix_col<uint64_t,T_n - T_k,T_n>;

	type_parity_check_matrix& parity_check_matrix;
	type_syndrome& syndrome;

	naive_syndrome_decoding(type_parity_check_matrix& _parity_check_matrix,type_syndrome& _syndrome) : parity_check_matrix(_parity_check_matrix), syndrome(_syndrome){
	}

	void execute(auto&& callback){
		auto callback2 = [&callback,this](type_syndrome& sum_columns,type_word& error){
			if(sum_columns == syndrome){
				callback(error);
			}
		};
		enumerate3<T_w>(parity_check_matrix,0,T_n,callback2);
	}
};

#endif