#ifndef HEADER_PARITY_CHECKS_EQUATION
#define HEADER_PARITY_CHECKS_EQUATION

#include <linear_algebra/linear_algebra.hpp>
#include <linear_algebra/low_weight_vector.hpp>
#include <linear_code/naive_syndrome_decoding.hpp>
#include <tools/math.hpp>

// Implements the procedure PARITY-CHECKS-EQUATIONS
template <size_t T_n, size_t T_k, size_t T_s, size_t T_w>
struct parity_checks_equations{
 
	using type_generator_matrix_row = t_matrix_row<uint64_t,T_k,T_n>;
	using type_generator_matrix_col = t_matrix_col<uint64_t,T_k,T_n>;

	using type_generator_matrix_shortened_code_col = t_matrix_col<uint64_t,T_k - T_s,T_n - T_s>;
	using type_matrix_R_col = t_matrix_col<uint64_t, T_s, T_n - T_s>;

	using type_low_weight_vector = low_weight_vector<T_n,T_s,T_w>;
	using type_decoder = naive_syndrome_decoding<T_n-T_s,T_n-T_k,T_w>;



	type_generator_matrix_row& generator_matrix;

	type_generator_matrix_shortened_code_col generator_matrix_shortened_code_col;
	type_matrix_R_col matrix_R;
	

	type_decoder syndrome_decoder;

	t_vector<uint64_t, T_k-T_s> null_vector;


	parity_checks_equations(type_generator_matrix_row& _generator_matrix) : generator_matrix_shortened_code_col(), matrix_R(), null_vector(), generator_matrix(_generator_matrix), syndrome_decoder(generator_matrix_shortened_code_col,null_vector) {
	}
	void iterate(auto&& callback){
		size_t nbSupress = 0;
		type_generator_matrix_row copy_generator_matrix(generator_matrix);
		bool success = partial_echelon_form(copy_generator_matrix, T_s,nbSupress);

		type_generator_matrix_col generator_matrix_col = copy_generator_matrix.inverse_representation();

		generator_matrix_col.template extract_last<T_s,0>(generator_matrix_shortened_code_col);
		generator_matrix_col.template extract_last<T_s,T_k - T_s>(matrix_R);

		auto callback2 = [&callback,this](typename type_low_weight_vector::type_word_N& h_N){
			typename type_low_weight_vector::type_word_P h_P(matrix_R*h_N);
			type_low_weight_vector _low_weight_vector(h_P,h_N);
			callback(_low_weight_vector);
		};

		syndrome_decoder.execute(callback2);
	}
};
#endif