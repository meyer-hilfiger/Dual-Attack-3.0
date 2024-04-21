import numpy
import math
from math import log2,comb
#from gmpy2 import bincoef, log2,log
#from gmpy2 import mpz,mpfr
import scipy.stats as stat

from scipy.optimize import brentq

def binary_search(arr,low, high, x):
	mid = (high + low) // 2
	if(high == low + 1):
		if arr[high] <= x:
			return high
		else:
			return low
	if arr[mid] == x:
		return mid
	elif arr[mid] > x:
		return binary_search(arr, low, mid, x)
	else:
		return binary_search(arr, mid, high, x)
def launch_binary_search(arr,x):
	if(x < arr[0]):
		return -1
	else:
		return binary_search(arr,0, len(arr) - 1, x)

def inverse_sum_partial(L):
	t = 0
	tmp = 0
	for i in range(len(L)-1, -1, -1) :
		tmp = L[i]
		L[i] += t
		t += tmp
		tmp = 0
	return L
def H2_G(x,a):
	if(x == 0 or x == 1):
		return -a
	return -x*log2(x) - (1-x)*log2(1-x) - a
	
def H2(x):
	assert x >= 0
	assert x <= 1
	return H2_G(x,0)

def H2_I(x):
	if(x == 0):
		return 0

	if(x == 1):
		return 0.5
	return brentq(H2_G,a=0,b=0.5,args=(x))
"""
def kraw(N,W,X):
	res = 0
	for j in range(0,W+1):
		res += mpz(((-1)**j))*bincoef(mpz(X),mpz(j))*bincoef(mpz(N-X),mpz(W-j))
	return res

def eps(N,W,X):
	return kraw(N,W,X)/bincoef(N,W)
"""

def kraw(N,W,X):
	res = 0
	for j in range(0,W+1):
		res += ((-1)**j)*comb(X,j)*comb(N-X,W-j)
	return res

def eps(N,W,X):
	return kraw(N,W,X)/comb(N,W)

def binom(n,k):
	if(n == k or k == 0):
		return 1
	return math.sqrt(1/(2*math.pi*n*(k/n)*(1-k/n)))*(2**(H2(k/n)*n))



def compute_list_eps_kraw(n,w):
	L = [numpy.zeros(n+1), numpy.zeros(n+1)]
	for i in range(n+1):
		L[1][i] = kraw(n,w,i)
		L[0][i] = L[1][i]/comb(n,w)
	return L




def compute_superior_tresholds(G,tresholds):
	X = numpy.zeros(len(tresholds))
	for g in G:
			idx = launch_binary_search(tresholds,g)
			if idx >= 0:
				X[idx] += 1
	return X
"""
def compute_superior_tresholds(G,tresholds):
	X = numpy.zeros(len(tresholds))
	for g in G:
			idx = launch_binary_search(tresholds,g)
			if idx >= 0:
				X[idx] += 1
	return X
"""

def compute_survival(L, G):
	max_treshold = 2*(len(G) - 1)
	max_tresh_current = int(max(L)) - (int(max(L))%2) + 2
	if(max_tresh_current > max_treshold):
		G.extend([0] * (int((max_tresh_current - max_treshold)/2)))
		max_treshold = max_tresh_current
	for l in L:
		if(l >= 0):
			G[int((int(l) - (int(l)%2))/2)] += 1
def sum_progressive(G):
	for i in range(len(G) - 2, -1, -1):
		G[i] += G[i+1]


def monteCarlo_bias(n,k,w,batch,optims,multiplier_in):
	L = numpy.zeros(batch)
	list_eps = optims[0]
	list_kraw = optims[1]
	weight_enumerator= numpy.zeros((n+1,batch))
	L_poisson = []
	for i in range(0,n+1,1):
		kr = list_kraw[i]
		ep = list_eps[i]
		E_point = comb(n,i)/(2**(n-k))
		P = stat.poisson.rvs(multiplier_in*float(E_point),loc=0, size=batch,random_state=None)
		weight_enumerator[i] = numpy.copy(P)
		L_poisson.append(kr*P)
	for LP in L_poisson:
		L = L + LP
	return L,weight_enumerator



def doubleRLPN_Poisson(n,k,s,w,kaux,taux,batch,number_It):
	survival_function = [0]
	list_eps_par = compute_list_eps_kraw(n-s,w)
	list_eps_aux = compute_list_eps_kraw(s,taux)
	eps_aux = list_eps_aux[0]
	kr_aux = list_eps_aux[1]
	j = 0
	for i in range(number_It):
		if(i*batch >= j):
			print("Generate data Poisson Model: " +str(int(j/(2**kaux))) + " / " + str(math.ceil(number_It*batch/(2**kaux))))
			j += 2**kaux
		list_bias = numpy.zeros(batch)
		list_bias_aux,list_enum_aux =  monteCarlo_bias(s,s-kaux,taux,batch,list_eps_aux,1.0)
		for z in range(s+1):
			assert(len(list_enum_aux[z]) == batch)
			list_bias_par,list_enum_par = monteCarlo_bias(n-s,k-s,w,batch,list_eps_par,list_enum_aux[z])  
			list_bias += (kr_aux[z])*list_bias_par
		list_bias *= 2**(kaux-k)
		compute_survival(list_bias,survival_function)
	sum_progressive(survival_function)
	return [survival_function[i]*((2**kaux - 1)/(number_It*batch)) for i in range(len(survival_function))]