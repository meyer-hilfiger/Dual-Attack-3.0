#ifndef HEADER_doubleRLPN_iteration
#define HEADER_doubleRLPN_iteration

#include <doubleRLPN/LPN_samples.hpp>
#include <doubleRLPN/parity_checks_equations.hpp>
#include <doubleRLPN/fft_decode.hpp>

#include <linear_code/linear_code.hpp>
#include <linear_algebra/linear_algebra.hpp>

#include <tools/permutation.hpp>



// Run one iteration of Algorithm 4.4 but does not finish the decoding: it stops at the "Recover-e" procedure.
template <size_t T_n,size_t T_k,size_t T_kaux, size_t T_s,size_t T_w, size_t T_taux>
struct doubleRLPN_iteration{
	using type_linear_code = linear_code<T_n,T_k>;
	using type_word = t_vector<uint64_t,T_n>;
	using type_LPN_samples = LPN_samples<T_n,T_k,T_kaux,T_s,T_w,T_taux>;
	using type_fft_decode = fft_decode<T_n, T_k, T_kaux, T_s, T_w,T_taux>;
	using type_permutation = permutation<T_n>;

	type_permutation permutation_positions;
	type_linear_code& code;
	type_word& y;

	type_linear_code code_permuted;
	type_word y_permuted;

	type_LPN_samples _LPN_samples;
	type_fft_decode _fft_decode;

	doubleRLPN_iteration(type_linear_code& _code,type_word& _y) : code(_code), y(_y), y_permuted(), code_permuted(code), _LPN_samples(code_permuted), _fft_decode(_LPN_samples,y_permuted){
	}

	bool execute(){
		permutation_positions.rand();

		code_permuted = code.permute(permutation_positions);
		y_permuted = y.permute(permutation_positions);

		// Fail is the permuted code is not of full rank dimension s
		if(!code_permuted.is_full_rank(T_s)){
			return false;
		}
		_fft_decode.execute();
		return true;
	}
};
#endif