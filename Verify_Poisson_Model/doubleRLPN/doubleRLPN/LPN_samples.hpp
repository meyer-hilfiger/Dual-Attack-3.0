#ifndef HEADER_LPN_samples
#define HEADER_LPN_samples

#include <iostream>
#include <linear_algebra/linear_algebra.hpp>
#include <linear_code/linear_code.hpp>
#include <linear_code/naive_list_decoder.hpp>
#include <linear_algebra/low_weight_vector.hpp>
#include <doubleRLPN/parity_checks_equations.hpp>

// Implements Algorithm 1 of the paper
template <size_t T_n, size_t T_k, size_t T_kaux, size_t T_s,size_t T_w, size_t T_taux>
struct LPN_samples{
	using type_code = linear_code<T_n,T_k>;
	using type_auxilary_code = linear_code<T_s,T_kaux>;
	using type_auxilary_decoder = naive_list_decoder<T_s,T_kaux,T_taux>;

	using type_parity_check = low_weight_vector<T_n, T_s, T_w>;
	using type_parity_checks_equations = parity_checks_equations<T_n,T_k,T_s,T_w>;

	type_code& code;

	type_auxilary_code auxilary_code;
	type_auxilary_decoder *auxilary_decoder;
	type_parity_checks_equations function_parity_checks_equations;

	size_t number_LPN_samples;



	LPN_samples(type_code& _code) : code(_code), auxilary_code(), function_parity_checks_equations(code.generator_matrix){
		number_LPN_samples = 0;
		auxilary_decoder = new type_auxilary_decoder(auxilary_code);
	}
	~LPN_samples(){
		delete auxilary_decoder;
	}
	void execute(auto&& callback){
		number_LPN_samples = 0;
		auto callback2 = [&callback,this](type_parity_check& h){
			auto callback3 = [&h,&callback,this](typename type_auxilary_code::type_codeword& auxilary_codeword){
				number_LPN_samples += 1;
				callback(h,auxilary_codeword);
			};
			auxilary_decoder->decode_codeword_callback(h.word_P,callback3);
		};
		function_parity_checks_equations.iterate(callback2);
	}
};

#endif