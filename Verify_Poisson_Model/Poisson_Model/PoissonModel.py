# Maintainer: Charles Meyer-Hilfiger

from tools import *
import math
import numpy
import sys
import pathlib


def generate_data(param_doubleRLPN,Niter = None,batch = None):
	[w,taux,kaux,s,k,n,t] = param_doubleRLPN
	if(Niter== None):
		Niter = 10
	if(batch == None):
		batch = min(Niter*(2**kaux),2**10)
	nbIt = math.ceil(Niter*(2**kaux)/batch)
	list_fourier = doubleRLPN_Poisson(n,k,s,w,kaux,taux,batch,nbIt)
	L = [[2*i, list_fourier[i]] for i in range(len(list_fourier))]
	currentDir = str(pathlib.Path(__file__).parent.resolve())
	name = currentDir  + "/data/PoissonModel_"  + str(w) + "_" + str(taux)+ "_" + str(kaux) + "_" + str(s)+ "_" + str(k) + "_"+ str(n)+ "_" + str(t) + "_" + str(Niter) + ".csv"
	with open(name, 'w') as data:
		numpy.savetxt(data, L, delimiter=",")

def main():
	assert(len(sys.argv) == 9)
	[w,taux,kaux,s,k,n,t,Niter] = [int(sys.argv[i]) for i in range(1,9)]
	generate_data([w,taux,kaux,s,k,n,t],Niter)

if __name__ == "__main__":
    main()