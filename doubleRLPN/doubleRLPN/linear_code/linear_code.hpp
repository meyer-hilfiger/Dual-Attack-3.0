#ifndef HEADER_LINEARCODE
#define HEADER_LINEARCODE

#include <linear_algebra/linear_algebra.hpp>
#include <tools/permutation.hpp>

// RANDOM LINEAR CODE
template <size_t T_n, size_t T_k>
struct linear_code{
	using type_message = t_vector<uint64_t,T_k>;
	using type_syndrome = t_vector<uint64_t,T_n-T_k>;
	using type_codeword = t_vector<uint64_t,T_n>;
	using type_parity_check_matrix = t_matrix_col<uint64_t,T_n-T_k,T_n>;
	using type_generator_matrix = t_matrix_row<uint64_t,T_k,T_n>;
	
	permutation<T_n> permutation_positions;
	type_parity_check_matrix parity_check_matrix;
	type_generator_matrix generator_matrix;

	public:
	linear_code(){
		// Choose support of weight k uniformly at random in [1,n]
		parity_check_matrix.init();
		generator_matrix.init();
		for(size_t i = 0; i < T_n - T_k; ++i){
			parity_check_matrix.set(i,i + T_k,1);
		}
		for(size_t i = 0; i < T_k; ++i){
			generator_matrix.set(i,i,1);
		}
		for(size_t i = 0; i < T_k; ++i){
			for(size_t j = 0; j < T_n - T_k; ++j){
				bool r = rand()%2;
				generator_matrix.set(i,T_k + j,r);
				parity_check_matrix.set(j,i,r);
			}
		}

		permutation_positions.rand();

		parity_check_matrix = parity_check_matrix.permute_columns(permutation_positions);
		generator_matrix = generator_matrix.permute_columns(permutation_positions);

	}

	// c must be a codeword of the code for decode to return the right message
	type_message decode(type_codeword c){
		type_message m;
		for(int i = 0; i < T_k; ++i){
			m.set(i,c[permutation_positions[i]]);
		}
		return m;
	}
	type_codeword encode(type_message m){
		return m*generator_matrix;
	}

	// Verify that the code is of full rank (dimension s) when restricted to its s first coordinates
	// s <= k
	bool is_full_rank(size_t s){
		size_t stop=0;
		type_generator_matrix copy_generator_matrix(generator_matrix);
		return partial_echelon_form(copy_generator_matrix,s, stop);
	}

	linear_code permute(permutation<T_n>& perm){
		linear_code code_permuted;
		code_permuted.generator_matrix = generator_matrix.permute_columns(perm);
		code_permuted.parity_check_matrix = parity_check_matrix.permute_columns(perm);
		code_permuted.permutation_positions = perm.composition(permutation_positions);
		return code_permuted;
	}
};

#endif