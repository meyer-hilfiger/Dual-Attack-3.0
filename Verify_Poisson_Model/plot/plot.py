# Maintainer: Charles Meyer-Hilfiger
import numpy as np
import matplotlib.pyplot as plt
import csv
from math import log2
import sys
import os

import importlib
def read_file(filename):
	x = []
	y = []
	with open(filename, 'r') as f:
		reader = csv.reader(f,delimiter = ',')
		for row in reader:
			x.append(float(row[0]))
			y.append(float(row[1]))
	return x,y

def get_limit(y,minimum):
	for i in range(len(y)):
		if y[i] < minimum:
			return i
	return len(y) - 1

def make_plot(L_parameters,Niter):
	[w,taux,kaux,s,k,n,t] = L_parameters
	filename_doubleRLPN = "../doubleRLPN/data/doubleRLPN_" + str(w) + "_" + str(taux)+ "_" + str(kaux) + "_" + str(s)+ "_" + str(k) + "_"+ str(n)+ "_" + str(t) + "_" + str(Niter) + ".csv"
	filename_PoissonModel = "../Poisson_Model/data/PoissonModel_" + str(w) + "_" + str(taux)+ "_" + str(kaux) + "_" + str(s)+ "_" + str(k) + "_"+ str(n)+ "_" + str(t) + "_" + str(Niter) + ".csv"

	x_d,y_d = read_file(filename_doubleRLPN)
	x_p,y_p = read_file(filename_PoissonModel)


	minimum = 300/Niter
	x_max_d = x_d[get_limit(y_d,minimum)]
	x_max_p = x_p[get_limit(y_p,minimum)]

	x_max = min(x_max_d,x_max_p) - 1

	fig, axs = plt.subplots()

	axs.set_yscale('log', base=2)
	axs.set_xlim(right=x_max)
	

	axs.plot(x_d,y_d,color='red',linewidth=1.5,label='doubleRLPN',alpha=0.6)
	axs.plot(x_p,y_p,color='blue',linewidth=1.5,label='Poisson Model',alpha=0.6)

	axs.legend()

	axs.grid(visible=True)

	axs.set_xlabel(r'T')
	axs.set_ylabel(r'$\mathbb{E}\left( \left| \{ \mathbf{x} \in \mathbb{F}_2^{k_{aux}} \backslash  \{ \mathbf{e}_{P} \mathbf{G}_{aux} \} \: : \: \widehat{f} \left( \mathbf{x} \right) \geq T \} \right|\right)$')
	outputname = "plot_" + str(w) + "_" + str(taux)+ "_" + str(kaux) + "_" + str(s)+ "_" + str(k) + "_"+ str(n)+ "_" + str(t) + "_" + str(Niter) + ".pdf"
	plt.savefig(outputname,dpi=1000,bbox_inches='tight',transparent=True)


def main():
	assert(len(sys.argv) == 9)
	[w,taux,kaux,s,k,n,t,Niter] = [int(sys.argv[i]) for i in range(1,9)]
	L_parameters = [w,taux,kaux,s,k,n,t]

	filename_doubleRLPN = "../doubleRLPN/data/doubleRLPN_" + str(w) + "_" + str(taux)+ "_" + str(kaux) + "_" + str(s)+ "_" + str(k) + "_"+ str(n)+ "_" + str(t) + "_" + str(Niter) + ".csv"
	filename_PoissonModel = "../Poisson_Model/data/PoissonModel_" + str(w) + "_" + str(taux)+ "_" + str(kaux) + "_" + str(s)+ "_" + str(k) + "_"+ str(n)+ "_" + str(t) + "_" + str(Niter) + ".csv"

	exists_doubleRLPN = os.path.exists(filename_doubleRLPN)
	exists_PoissonModel = os.path.exists(filename_PoissonModel)
	if(not exists_doubleRLPN):
		sys.path.append('../doubleRLPN/')
		doubleRLPN = importlib.import_module('doubleRLPN')
		doubleRLPN.make_launch(w,taux,kaux,s,k,n,t,Niter)
	if(not exists_PoissonModel):
		sys.path.append('../Poisson_Model/')
		PoissonModel = importlib.import_module('PoissonModel')
		PoissonModel.generate_data(L_parameters,Niter)
	make_plot(L_parameters,Niter)

if __name__ == "__main__":
	main()