# Maintainer: Charles Meyer-Hilfiger
import numpy as np
import matplotlib.pyplot as plt
import csv
from math import log2
import sys
import os

import importlib
import argparse
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
	filename_doubleRLPN = "doubleRLPN/data/doubleRLPN_" + str(w) + "_" + str(taux)+ "_" + str(kaux) + "_" + str(s)+ "_" + str(k) + "_"+ str(n)+ "_" + str(t) + "_" + str(Niter) + ".csv"
	filename_PoissonModel = "Poisson_Model/data/PoissonModel_" + str(w) + "_" + str(taux)+ "_" + str(kaux) + "_" + str(s)+ "_" + str(k) + "_"+ str(n)+ "_" + str(t) + "_" + str(Niter) + ".csv"

	x_d,y_d = read_file(filename_doubleRLPN)
	x_p,y_p = read_file(filename_PoissonModel)


	minimum = 500/Niter
	x_max_d = x_d[get_limit(y_d,minimum)]
	x_max_p = x_p[get_limit(y_p,minimum)]

	i_min = min(get_limit(y_d,minimum),get_limit(y_p,minimum))

	x_max = min(x_max_d,x_max_p)

	fig, axs = plt.subplots()
	

	axs.set_yscale('log', base=2)
	axs.set_xlim(right=x_max)
	x_dd = x_d[:i_min]
	y_dd = y_d[:i_min]

	y_pp = y_p[:i_min]
	x_pp = x_p[:i_min]
	

	axs.plot(x_dd,y_dd,color='red',linewidth=1.5,label='doubleRLPN',alpha=0.6)
	axs.plot(x_pp,y_pp,color='blue',linewidth=1.5,label='Poisson Model',alpha=0.6)

	axs.legend()

	axs.grid(visible=True)

	axs.set_xlabel(r'T')
	axs.set_ylabel(r'$\mathbb{E}\left( \left| \{ \mathbf{x} \in \mathbb{F}_2^{k_{aux}} \backslash  \{ \mathbf{e}_{P} \mathbf{G}_{aux} \} \: : \: \widehat{f} \left( \mathbf{x} \right) \geq T \} \right|\right)$')
	outputname = "plot/plot_" + str(w) + "_" + str(taux)+ "_" + str(kaux) + "_" + str(s)+ "_" + str(k) + "_"+ str(n)+ "_" + str(t) + "_" + str(Niter) + ".pdf"
	plt.savefig(outputname,dpi=1000,bbox_inches='tight',transparent=True)


def main():
	parser = argparse.ArgumentParser(description='Verify the Poisson Model. Read README.pdf at root for more information.')
	parser.add_argument('--d1', dest='d1', action='store_true',
						default=False,
						help='Create a dataset containing the expected number of false candidates given experimentally by doubleRLPN.')
	parser.add_argument('--d2', dest='d2', action='store_true',
						default=False,
						help='Create a dataset containing the expected number of false candidates given by the Poisson Model.')
	parser.add_argument('--plot', dest='plot', action='store_true',
						default=False,
						help='Combine the two previous datasets into a plot in the folder plot/. This option must be either combined with option d1 and d2 if the corresponding datasets do not already exist, or can be used alone if the datasets already exist.')

	parser.add_argument('w', metavar='w', type=int, help='Parameter of doubleRLPN')
	parser.add_argument('taux', metavar='taux', type=int, help='Parameter of doubleRLPN')
	parser.add_argument('kaux', metavar='kaux', type=int, help='Parameter of doubleRLPN')
	parser.add_argument('s', metavar='s', type=int, help='Parameter of doubleRLPN')
	parser.add_argument('k', metavar='k', type=int, help='Parameter of doubleRLPN')
	parser.add_argument('n', metavar='n', type=int, help='Parameter of doubleRLPN')
	parser.add_argument('t', metavar='t', type=int, help='Parameter of doubleRLPN')
	parser.add_argument('Niter', metavar='Niter', type=int, help='Number of times the algorithm runs')

	args = parser.parse_args()

	w = args.w
	taux = args.taux
	kaux = args.kaux
	s = args.s
	k = args.k
	n = args.n
	t = args.t
	Niter = args.Niter
	L_parameters = [w,taux,kaux,s,k,n,t]
	if(not args.d1 and not args.d2 and not args.plot):
		args.d1 = True
		args.d2 = True
		args.plot = True
	filename_doubleRLPN = "doubleRLPN/data/doubleRLPN_" + str(w) + "_" + str(taux)+ "_" + str(kaux) + "_" + str(s)+ "_" + str(k) + "_"+ str(n)+ "_" + str(t) + "_" + str(Niter) + ".csv"
	filename_PoissonModel = "Poisson_Model/data/PoissonModel_" + str(w) + "_" + str(taux)+ "_" + str(kaux) + "_" + str(s)+ "_" + str(k) + "_"+ str(n)+ "_" + str(t) + "_" + str(Niter) + ".csv"


	if(args.d1):
		sys.path.append('doubleRLPN/')
		doubleRLPN = importlib.import_module('doubleRLPN')
		doubleRLPN.make_launch(w,taux,kaux,s,k,n,t,Niter)
	if(args.d2):
		sys.path.append('Poisson_Model/')
		PoissonModel = importlib.import_module('PoissonModel')
		PoissonModel.generate_data(L_parameters,Niter)
	if(args.plot):
		exists_doubleRLPN = os.path.exists(filename_doubleRLPN)
		exists_PoissonModel = os.path.exists(filename_PoissonModel)
		ok = 1
		if(not exists_doubleRLPN):
			print("DoubleRLPN dataset missing. You must use option --d1 at least once with these parameters.")
			ok = 0
		if(not exists_PoissonModel):
			print("Poisson Model dataset missing. You must use option --d2 at least once with these parameters.")
			ok = 0
		if(ok == 1):
			make_plot(L_parameters,Niter)
if __name__ == "__main__":
	main()