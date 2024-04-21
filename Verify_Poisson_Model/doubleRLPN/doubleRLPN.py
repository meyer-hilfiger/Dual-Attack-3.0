# Maintainer: Charles Meyer-Hilfiger

import subprocess
import sys
import pathlib


def make(N,K,S,W,TAUX,KAUX,T,name):
	currentDir = str(pathlib.Path(__file__).parent.resolve())
	prec = " -D PARAM_N=" + str(N) + " -D PARAM_K=" + str(K) + " -D PARAM_S=" + str(S) + " -D PARAM_W=" + str(W) + " -D PARAM_TAUX=" + str(TAUX) + " -D PARAM_KAUX=" + str(KAUX)+ " -D PARAM_T=" + str(T)
	optim = " -O3 -march=native -I ."
	fil = " main.cpp -std=c++20"
	compiler = "g++"
	command = compiler + fil + optim + prec + " -o " + name
	comp = subprocess.Popen(command,shell=True,cwd=currentDir)
	comp.wait()


def launch(name,number_instances):
	currentDir = str(pathlib.Path(__file__).parent.resolve())
	command = "./" +name  + " " + str(number_instances)
	comp = subprocess.Popen(command,shell=True,cwd=currentDir)
	comp.wait()


def make_launch(w,taux,kaux,s,k,n,t,Niter):
	name = "RLPN_decoder.out"
	make(n,k,s,w,taux,kaux,t,name)
	launch(name,Niter)

def main():
	assert(len(sys.argv) == 9)
	#[w,taux,kaux,s,k,n,t] = [5,2,20,28,30,60,8]
	[w,taux,kaux,s,k,n,t,Niter] = [int(sys.argv[i]) for i in range(1,9)]
	#print([w,taux,kaux,s,k,n,t,Niter])
	make_launch(w,taux,kaux,s,k,n,t,Niter)

if __name__ == "__main__":
    main()