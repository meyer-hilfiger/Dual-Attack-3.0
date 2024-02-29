# Authors: Charles Meyer-Hilfiger

name = "RLPN_decoder.out"
import re
import subprocess
import sys

def compile(N,K,S,W,TAUX,KAUX,T,name):
	prec = " -D PARAM_N=" + str(N) + " -D PARAM_K=" + str(K) + " -D PARAM_S=" + str(S) + " -D PARAM_W=" + str(W) + " -D PARAM_TAUX=" + str(TAUX) + " -D PARAM_KAUX=" + str(KAUX)+ " -D PARAM_T=" + str(T)
	optim = " -O3 -march=native -I ."
	fil = " main.cpp -std=c++20"
	compiler = "g++"
	command = compiler + fil + optim + prec + " -o " + name
	comp = subprocess.Popen(command,shell=True)
	comp.wait()



def launch(name,number_instances):
	command = "./" +name  + " " + str(number_instances)
	comp = subprocess.Popen(command,shell=True)
	comp.wait()


def main():
		compile(n,k,s,w,taux,kaux,t,name)
		launch(name,number_instances)

if __name__ == "__main__":
    # execute only if run as a script
    main()