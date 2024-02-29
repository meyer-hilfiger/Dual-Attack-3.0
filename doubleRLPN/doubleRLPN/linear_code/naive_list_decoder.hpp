#ifndef HEADER_NAIVELISTDECODER
#define HEADER_NAIVELISTDECODER

#include <linear_algebra/linear_algebra.hpp>
#include <tools/enumerate.hpp>
#include <tools/math.hpp>


template <size_t T_n, size_t T_k, size_t T_w>
struct naive_list_decoder;
// Naive list decoder for an unstructured code
template <size_t T_n, size_t T_k, size_t T_w>
struct naive_list_decoder{
	using type_linear_code = linear_code<T_n, T_k>;

	using type_word = typename type_linear_code::type_codeword;
	using type_message = typename type_linear_code::type_message;
	using type_syndrome = typename type_linear_code::type_syndrome;

	using type_error = vector_support<T_w>;

	type_linear_code& code;

	std::vector<type_word> *syndrome_errors;

	private:
		void compute_syndromes(){
			auto add_error_to_syndrome = [this](type_syndrome& syndrome, type_error& error){
				syndrome_errors[syndrome.to_integer()].push_back(error);
			};
			if constexpr (T_w >= 1){
				enumerate3<T_w>(code.parity_check_matrix, 0, T_n, add_error_to_syndrome);
			}else{
				type_word v;
				syndrome_errors[0].push_back(v);
			}
		}
	public:
		naive_list_decoder(type_linear_code& _code) : code(_code){
			syndrome_errors = new std::vector<type_word>[power(2,T_n-T_k)]();
			compute_syndromes();
		}
		~naive_list_decoder(){
			delete[] syndrome_errors;
		}
		void decode_error_callback(type_word& w,auto&& callback){
			type_syndrome s = code.parity_check_matrix*w;
			size_t size = syndrome_errors[s.to_integer()].size();
			for(size_t i = 0; i < size; ++i){
				type_word err(syndrome_errors[s.to_integer()][i]);
				callback(err);
			}
		}
		void decode_codeword_callback(type_word& w, auto&& callback){
			auto callback2 = [&w,&callback](type_word& e){
				type_word v_r(e);
				v_r ^= w;
				callback(v_r);
			};
			decode_error_callback(w,callback2);
		}

};

#endif