# Maintainer: Charles Meyer-Hilfiger
import sys
from function_utilitaries import *
from function_doubleRLPN_baseComplexity import *
import matplotlib.pyplot as plt

############# BEGIN BJMM12 FUNCTIONS

# complexity_BJMM12(R, w, l_1,l_2,p_1,p_2)
# Asymptotic exponent for computing parity-checks of relative weight w of a code of rate R with BJMM12 technique, Proposition 11

# Inputs: R,w,l_1,l_2,p_1,p_2 (Float)

# --------------
def complexity_asym_BJMM12(R, w, l_1,l_2,p_1,p_2):
		# Size of the lists of BJMM12 at each step
		S0 = asym_binomial(0.5,p_1/2.0)
		S1 = asym_binomial(1.0,p_1) - l_1
		S2 = asym_binomial(1.0,p_2) - l_2
		S3 = asym_binomial(1.0,w) - R
		

		# Time complexity of BJMM12 at each step
		T0 = S0
		T1 = max(S0, 2*S0 - l_1)
		T2 = max(S1, 2*S1 - (l_2 - l_1))
		T3 = max(S2, 2*S2 - (R-l_2))

		complexity = max(T0, max(T1,max(T2,max(T3,S3))))
		memory = max(S0,S1,S2)

		return np.array([complexity,memory])


# --------------
# constraints_BJMM12(R,w,l_1,l_2,p_1,p_2, tolerance)
# BJM12 parameters (l_1,l_2,p_1,p_2) must verify constraints (See Proposition 11) in order for BJMM12 method to compute all parity-checks of relative weight w of a code of rate R

# Output: True if constraints are verified (Up to tolerance)
#		  False else

# Inputs: R,w,l_1,l_2,p_1,p_2,tolerance (Float)

# --------------
def constraints_BJMM12(R,w,l_1,l_2,p_1,p_2, tolerance):

	assert(l_1 >= 0 and l_1 <= 1)
	assert(l_2 >= 0 and l_2 <= 1)
	assert(p_1 >= 0 and p_1 <= 1)
	assert(p_2 >= 0 and p_2 <= 1)

	assert(l_1 <= l_2)
	assert(l_2 <= R)

	assert(((p_1 - p_2/2) >= 0) and ((p_1 - p_2/2) <= 1.0-p_2))
	assert(((p_2 - w/2) >= 0) and ((p_2 - w/2) <= 1.0-w))

	representation_1 = ((asym_binomial(1.0-p_2,p_1 - p_2/2) + p_2))
	assert(l_1 <= (representation_1 + tolerance)) # (38)

	representation_2 = ((asym_binomial(1.0-w,p_2 - w/2) + w))
	assert(l_2 <= (representation_2 + tolerance)) # (38)


	return True

############# END BJMM12 FUNCTIONS

def plot_complexity(R,C):
	
	f, ax = plt.subplots()
	plt.plot([0] + R + [1],[0] + C + [0],linestyle="solid",color="black",label="doubleRLPN")
	ax.yaxis.grid(color='gray', linestyle='dashed',linewidth=0.2,alpha=0.3)
	ax.xaxis.grid(color='gray', linestyle='dashed',linewidth=0.2,alpha=0.3)
	plt.xlim(0,1)
	plt.ylim(0,0.13)

	ax.legend(fontsize=25)
	plt.xlabel(r'Rate',fontsize=25)
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.ylabel(r'$\frac{\log_2 Complexity}{n}$',fontsize=25)
	plt.gcf().set_size_inches(13, 11)
	plt.savefig("complexity_doubleRLPN_BJMM12.pdf")



import csv
filename = "doubleRLPN_BJMM12.csv"
file1 = open(filename)
csv1 = csv.reader(file1,delimiter=";")
headers = next(csv1)
headersString = [str(k).strip() for k in headers]
rows = []
for row in csv1:
		rowFloat = [float(k.strip()) for k in row]
		rows.append(rowFloat)
file1.close()

# Index of each parameter in the parameter file
def main():
	idx_R = 0
	idx_s = 3
	idx_k_aux = 4
	idx_u = 5
	idx_w = 6
	idx_l_1 = 7
	idx_l_2 = 8
	idx_p_1 = 9
	idx_p_2 = 10
	idx_t = 11
	tolerance = 1e-7

	print("Complexity exponents of the Double RLPN decoder when decoding at Gilbert-Varshamov distance and using BJMM12 to compute low-weight parity-checks")
	list_R = []
	list_C = []
	for i in range(len(rows)):
		#Rate
		R = rows[i][idx_R]

		#Decoding distance
		tau = rows[i][idx_t]



		# Core double-RLPN parameters
		sigma = rows[i][idx_s]
		R_aux = rows[i][idx_k_aux]
		mu = rows[i][idx_u]
		omega = rows[i][idx_w]

		#BJMM12 related parameters to compute all parity checks of weight omega
		l_1 = rows[i][idx_l_1]
		l_2 = rows[i][idx_l_2]
		p_1 = rows[i][idx_p_1]
		p_2 = rows[i][idx_p_2]

		#We decode the auxilary code at Gilbert-Varshamov distance
		t_aux = sigma*Hinv(1-(R_aux/sigma))
		# N_aux= 1 works for every parameter
		N_aux = 1

		len_punctured_code = 1.0 - sigma
		R_punctured_code = R - sigma
		rel_BJMM12_R = R_punctured_code/len_punctured_code
		rel_BJMM12_w = omega/len_punctured_code
		rel_BJMM12_l_1 = l_1/len_punctured_code
		rel_BJMM12_l_2 = l_2/len_punctured_code
		rel_BJMM12_p_1 = p_1/len_punctured_code
		rel_BJMM12_p_2 = p_2/len_punctured_code
		rel_BJMM12_tolerance = tolerance
		

		constraints_BJMM12(rel_BJMM12_R,rel_BJMM12_w,rel_BJMM12_l_1,rel_BJMM12_l_2,rel_BJMM12_p_1,rel_BJMM12_p_2, rel_BJMM12_tolerance)
		complexity_parity_check_BJMM12 = len_punctured_code*complexity_asym_BJMM12(rel_BJMM12_R,rel_BJMM12_w,rel_BJMM12_l_1,rel_BJMM12_l_2,rel_BJMM12_p_1,rel_BJMM12_p_2)

		Comp,Memory = doubleRLPN_baseComplexity(R, sigma, R_aux, t_aux, mu, omega, tau, N_aux, complexity_parity_check_BJMM12, tolerance)
		print("Rate: " + "{:.5f}".format(R) + ";	" + "Complexity: " + "{:.5f}".format(Comp))# + ";	" + "Memory: " + "{:.5f}".format(Memory) )
if __name__ == '__main__':
    sys.exit(main())
