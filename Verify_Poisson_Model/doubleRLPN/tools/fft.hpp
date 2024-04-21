#ifndef HEADER_FFT
#define HEADER_FFT
#include <cassert>
#include <cstring>
#include <iostream>

#include <tools/math.hpp>
using namespace std;
template <size_t n, typename eval>
void fft(eval *f, eval* fft_f){
	assert(n < 64);
	memcpy(fft_f,f,sizeof(eval)*((size_t)  power(2,n)));
	eval *ep = new eval[(size_t) power(2,n)];
	for(size_t k = 1; k <= n; ++k){
		size_t kpow = power(2,k);
		for(size_t i = 0; i < power(2,n-k) ; ++i){
			size_t kpowi = kpow*i;
			for(size_t j = 0; j < power(2,k-1); ++j){
				ep[kpowi + j] = fft_f[kpowi+j] + fft_f[kpowi+j + kpow/2];
				ep[kpowi + kpow/2 + j] = fft_f[kpowi+j] - fft_f[kpowi+j + kpow/2];
			}
		}
		memcpy(fft_f,ep,sizeof(eval)*((size_t)  power(2,n)));
	}
	delete[] ep;
	
}

	
#endif