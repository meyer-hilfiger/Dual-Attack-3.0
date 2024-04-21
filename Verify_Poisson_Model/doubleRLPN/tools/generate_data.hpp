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
#include <vector>
#include <sstream>

template <typename type_v>
void add_vector(std::vector<type_v>& v, type_v value){
	if(value >= 0){
		size_t current_capacity = 2*(v.size()-1);
		size_t necessary_capacity = ((uint64_t)value) - (((uint64_t)value)%2) + 2;
		if (current_capacity < necessary_capacity){
			size_t offset = (((necessary_capacity - current_capacity)/2));
			v.resize(v.size() + offset, 0);
		}
		v[(uint64_t)(value/2)] += 1;
	}
}

template <typename type_int,typename type_survival>
void add_to_survival(type_int *fft_f, size_t T_kaux, size_t secretLPN,std::vector<type_survival>& survival){
	for(size_t i = 0; i < power(2,T_kaux); ++i){
		if((i != secretLPN) && (fft_f[i] >= 0)){
			add_vector(survival,(long double)fft_f[i]);
		}
	}
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

	std::vector<long double> survival_function;
	survival_function.push_back(0);

	size_t i = 0;
	while(i < number_instances){
		type_decoding_problem problem;
		type_doubleRLPN_iteration doubleRLPN_iter(problem.code, problem.y);


		bool sucess = doubleRLPN_iter.execute();
		if(sucess){
			std::cout << "Generate data doubleRLPN: " + std::to_string(i) + " / " + std::to_string(number_instances) << std::endl;
			// Get LPN secret : e_P*(G_aux^{\top})
			type_word e_permuted = problem.error.permute(doubleRLPN_iter.permutation_positions);
			type_word_P e_P;
			e_permuted.template extract_last<T_n-T_s>(e_P);
			type_word_LPN LPN_secret = doubleRLPN_iter._LPN_samples.auxilary_code.generator_matrix.inverse_representation()*e_P;
			add_to_survival( doubleRLPN_iter._fft_decode.fft_f,T_kaux, LPN_secret.to_integer(),survival_function);
			i+=1;
		}
		
	}
	for (long long int j = survival_function.size() - 2; j >= 0; j--){
		survival_function[j] += survival_function[j+1];
	}
	
	string filename = "data/doubleRLPN_" + to_string(T_w) + "_" + to_string(T_taux) + "_" + to_string(T_kaux) + "_" + to_string(T_s) + "_" +to_string(T_k) + "_" +to_string(T_n) + "_" +to_string(T_t) + "_" +to_string(number_instances) + ".csv";
	std::ofstream file;
	file.open(filename);
	int precision = 15;
	file.precision(precision);
	file << std::scientific;

	for(size_t j = 0; j < survival_function.size(); ++j){
		survival_function[j] /= (( long double ) number_instances);

		file << 2*j << ","<< survival_function[j] << std::endl;

	}
	file.close();
}

#endif