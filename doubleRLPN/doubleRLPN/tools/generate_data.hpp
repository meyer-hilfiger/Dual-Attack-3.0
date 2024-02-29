#ifndef HEADER_GENERATEDATA
#define HEADER_GENERATEDATA

#include <linear_algebra/linear_algebra.hpp>
#include <linear_code/linear_code.hpp>
#include <linear_code/decoding_problem.hpp>
#include <doubleRLPN/doubleRLPN_iteration.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
template <typename type_int>
void write_to_file(auto& file, type_int *fft_f, size_t T_kaux, size_t secretLPN){
	type_int *copy_fft_f = new type_int[power(2,T_kaux)];
	copy (fft_f, fft_f + power(2,T_kaux), copy_fft_f); 
	std::sort(copy_fft_f, copy_fft_f +power(2,T_kaux) );
	file << fft_f[secretLPN] << ",";

	type_int base = copy_fft_f[0];
	size_t counter = 1;

	for(size_t i = 1; i < power(2,T_kaux); ++i){
		if(copy_fft_f[i] == base){
			counter += 1;
		}else{
			file << base << "," << counter;
			file << ",";
			counter = 1;
			base = copy_fft_f[i];
		}
	}
	file << base << "," << counter;
}

template <size_t T_n, size_t T_k,size_t T_kaux, size_t T_s,size_t T_w, size_t T_taux,size_t T_t>
void generate_data(size_t number_instances){
	using type_decoding_problem = decoding_problem<T_n,T_k,T_t>;
	using type_doubleRLPN_iteration = doubleRLPN_iteration<T_n, T_k, T_kaux, T_s, T_w, T_taux>;
	using type_word = t_vector<uint64_t,T_n>;
	using type_word_P = t_vector<uint64_t,T_s>;
	using type_word_LPN = t_vector<uint64_t,T_kaux>;

	using type_linear_code = linear_code<T_n,T_k>;
	using type_auxilary_code = linear_code<T_s,T_kaux>;
	

	string filename = "experimental_data/doubleRLPN_" + to_string(T_w) + "_" + to_string(T_taux) + "_" + to_string(T_kaux) + "_" + to_string(T_s) + "_" +to_string(T_k) + "_" +to_string(T_n) + "_" +to_string(T_t) + ".csv";
	
	std::ofstream file;
	file.open(filename);

	size_t i = 0;
	while(i <  number_instances){
		type_decoding_problem problem;
		type_doubleRLPN_iteration doubleRLPN_iter(problem.code, problem.y);


		bool sucess = doubleRLPN_iter.execute();
		if(sucess){
			// Get LPN secret : e_P*(G_aux^{\top})
			type_word e_permuted = problem.error.permute(doubleRLPN_iter.permutation_positions);
			type_word_P e_P;
			e_permuted.template extract_last<T_n-T_s>(e_P);
			type_word_LPN LPN_secret = doubleRLPN_iter._LPN_samples.auxilary_code.generator_matrix.inverse_representation()*e_P;
			
			write_to_file(file, doubleRLPN_iter._fft_decode.fft_f,T_kaux, LPN_secret.to_integer()) ;
			file << std::endl;
			i+=1;
		}
		
	}
	file.close();
}

#endif