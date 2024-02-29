
import numpy as np
import matplotlib.pyplot as plt
import csv
from math import log2

# Plot the average experimental size of the set of candidates in function of the treshold
# First execute  ../doubleRLPN/make.py on your desired parameter to create a dataset
# Then plug your parameters here
L_parameters = [5,2,20,28,30,60,8]
# The x-axis of graph being plotted
tresholds = np.arange(0,4000,20)


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


def associated_treshold(Value,Number,tresholds):
	X = np.zeros(len(tresholds))
	for i in range(len(Value)):
			n = Number[i]
			v = Value[i]
			idx = launch_binary_search(tresholds,v)
			if idx >= 0:
				X[idx] += n
	return X

def survival_function_from_experimental_data(L_parameters, tresholds):
	[w,taux,kaux,s,k,n,t] = L_parameters
	filename = "../doubleRLPN/experimental_data/doubleRLPN_" + str(w) + "_" + str(taux)+ "_" + str(kaux) + "_" + str(s)+ "_" + str(k) + "_"+ str(n)+ "_" + str(t) + ".csv"
	
	S = np.zeros(len(tresholds))
	with open(filename, 'r') as f:
		reader = csv.reader(f,delimiter = ',')
		number_instances = 0
		for row in reader:
			number_instances += 1
			data = list(row)
			data = np.array(data)
			nb_disinct_values = int((len(data)-1)/2)
			fft_value_LPN_secret = int(data[0])
			fft_value = np.zeros(nb_disinct_values)
			fft_number = np.zeros(nb_disinct_values)
			for i in range(nb_disinct_values):
				fft_value[i] = data[2*i+1]
				fft_number[i] = data[2*i + 2]
				if(int(fft_value[i]) == fft_value_LPN_secret):
					fft_number[i] -= 1
			S += associated_treshold(fft_value,fft_number,tresholds)
	return inverse_sum_partial(S)/number_instances

survival = survival_function_from_experimental_data(L_parameters,tresholds)

fig, axs = plt.subplots()
axs.set_yscale('log', base=2)

axs.plot(tresholds,[x for x in survival],color='red',linewidth=2.2,label='Experiments',alpha=0.6)

[w,taux,kaux,s,k,n,t] = L_parameters
outputname = "experiment_" + str(w) + "_" + str(taux)+ "_" + str(kaux) + "_" + str(s)+ "_" + str(k) + "_"+ str(n)+ "_" + str(t) + ".pdf"
plt.savefig(outputname,dpi=1000,bbox_inches='tight',transparent=True)
	
