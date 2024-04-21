#ifndef DECODING_PROBLEM_HEADER
#define DECODING_PROBLEM_HEADER
template <size_t T_n,size_t T_k,size_t T_t>
struct decoding_problem{
	using type_code = linear_code<T_n,T_k>;
	using type_word = t_vector<uint64_t,T_n>;
	using type_message = t_vector<uint64_t,T_k>;

	
	// Secret
	type_message message;
	type_word codeword;
	type_word error;

	// Public
	type_code code;
	type_word y;

	decoding_problem() : code() {
		message.rand();
		codeword = message*code.generator_matrix;
		error.rand_weight(T_t);
		y = codeword^error;
	}
};

#endif