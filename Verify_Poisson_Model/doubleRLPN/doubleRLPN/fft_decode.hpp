#ifndef HEADER_FFT_DECODE
#define HEADER_FFT_DECODE

#include <iostream>

#include <linear_algebra/linear_algebra.hpp>
#include <linear_code/linear_code.hpp>
#include <doubleRLPN/LPN_samples.hpp>
#include <tools/fft.hpp>
#include <tools/math.hpp>

template <size_t T_n, size_t T_k, size_t T_kaux, size_t T_s, size_t T_w,size_t T_taux>
struct fft_decode{
	using type_auxilary_code = linear_code<T_s,T_kaux>;
	using type_LPN_samples = LPN_samples<T_n,T_k,T_kaux,T_s,T_w,T_taux>;
	using type_word =  t_vector<uint64_t,T_n>;
	using type_parity_check =  typename type_LPN_samples::type_parity_check;

	using type_auxilary_codeword = typename type_auxilary_code::type_codeword;
	using type_auxilary_message = typename type_auxilary_code::type_message;

	type_LPN_samples& _LPN_samples;
	type_word& y;


	long int *f;

	long int *fft_f;

	fft_decode(type_LPN_samples& v_LPN_samples, type_word& _y) : _LPN_samples(v_LPN_samples), y(_y){
		if constexpr(T_kaux >= 63){
			std::cout << "Not yet implemented for kaux >= 63. Change long int into long long int in fft_decode" << std::endl;
			exit(0);
		}
		f = new long int[power(2,T_kaux)];
		fft_f = new long int[power(2,T_kaux)];
		for(size_t i = 0; i < power(2,T_kaux); ++i){
			f[i] = 0;
		}
	}
	~fft_decode(){
		delete[] f;
		delete[] fft_f;
	}
	void execute(){
		auto callback1 = [this](type_parity_check& h, type_auxilary_codeword& c_aux){
	
			uint8_t sample = y.dot_product(h.to_vector());
			type_auxilary_message m_aux = _LPN_samples.auxilary_decoder->code.decode(c_aux);
			if(sample == 0){
				f[m_aux.to_integer()] += 1;
			}else{
				f[m_aux.to_integer()] -= 1;
			}
		};
		_LPN_samples.execute(callback1);
		fft<T_kaux>(f,fft_f);
	}
};

#endif