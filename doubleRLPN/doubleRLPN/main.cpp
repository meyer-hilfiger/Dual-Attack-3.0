#include <iostream>
#include <cmath>

#include <vector>

#include <tools/generate_data.hpp>

using namespace std;

int main(int argc, char **argv){
	srand(time(NULL));
	int number_instance = atoi(argv[1]);
	generate_data<PARAM_N, PARAM_K,PARAM_KAUX, PARAM_S,PARAM_W, PARAM_TAUX,PARAM_T>(number_instance);

}