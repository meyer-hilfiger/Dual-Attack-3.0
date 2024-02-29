
# This script computes (using the formula in Section 3.3, Proposition 5) the complexity exponent of the Double RLPN algorithm (Section 3.2, Algorithm 1) when using BJMM12 technique to compute low-weight parity checks.
# This script also verifies that the parameters meet the constraints of Proposition 5.

# The parameters of the algorithm are taken from the file doubleRLPN_BJMM12.csv

# The Sections, Algorithms, Propositions mentioned here all refers to the eprint version of the article "Reduction from sparse to plain LPN and application to RLPN decoding, statistical Decoding 3.0": https://eprint.iacr.org/2023/1852.pdf.

# $R$ = $\frac{k}{n}$ the rate of the code to decode?
# $\alpha_{\text{DoubleRLPN-BJMM12}}$ = $\frac{log2(Time complexity of Double RLPN with BJMM12 decoder)}{n}$ when n -> Infinity. This value is only here for reference, this script will recompute it from the parameters of the algorithm.
# $\beta_{\text{DoubleRLPN-BJMM12}}$ = $\frac{log2(Memory complexity of Double RLPN with BJMM12 decoder)}{n}$ when n -> Infinity. This value is only here for reference, this script will recompute it from the parameters of the algorithm.
# $\sigma$ = $\frac{s}{n}$
# $R_{\text{aux}}$ = $\frac{k_aux}{n}$
# $\nu$ = $\frac{u}{n}$
# $\omega$ = $\frac{w}{n}$

# Internal parameters of the BJMM12 method to compute low-weight parity-checks. They are expressed relatively to n:
# $\ell_1$ = $\frac{l_1}{n}$
# $\ell_2$ = $\frac{l_2}{n}$
# $\pi_1$ = $\frac{p_1}{n}$
# $\pi_2$ = $\frac{p_2}{n}$
# $\tau$ = $\frac{t}{n}$ (Here it is equal to the relative Gilbert-Varshamov distance of a code of rate R, namely, $Hinv(1-R)$ where Hinv is the inverse of the binary entropy function)

# We take $\rho$ = $\frac{r}{n}$ as being equal to the relative Gilbert-Varshamov distance of an [s, k_aux] linear code.

############# BEGIN UTILITIES FUNCTIONS
from math import log2,sqrt
from scipy.optimize import brentq
import numpy as np
#from doubleRLPN_BJMM12_parameters_ISD_Dumer import parameters_Decode_Dumer
# --------------
# H2(x)
# Binary entropy function
# Output: -x*log2(x) - (1-x)*log2(1-x)

# Inputs: x (Float)

# Constraints:
# 0 <= x <= 1

# --------------	
def H2(x):
	assert x >= 0
	assert x <= 1
	return H2_G(x,0)


# --------------
# Hinv(x)
# Inverse of binary entropy function in (0,0.5)

# Inputs: x (Float)

# Constraints:
# 0 <= x <= 1

# --------------	
def H2_G(x,a):
	if(x == 0 or x == 1):
		return -a
	return -x*log2(x) - (1-x)*log2(1-x) - a
def Hinv(x):
	if(x == 0):
		return 0

	if(x == 1):
		return 0.5
	return brentq(H2_G,a=0,b=0.5,args=(x),rtol=0.00000000001)



# --------------
# asym_binomial(n,k)
# Asymptotic exponent of the binomial coefficient
# Output: n*H2(k/n) ~ \frac{log2(\bincoef{n*p}{k*p})}{p} when p -> Infinity

# Inputs: n,k (Float)

# Constraints:
# 0 <= n <= 1
# 0 <= k <= n

# --------------
def asym_binomial(n,k):
	return n*H2(k/n)


# --------------
# asym_bias(w,u)
# Asymptotic exponent of bias(<e,h>) where e is a word of relative weight u and h is a random word of relative weight w
# Output: See Proposition 2

# Inputs: w,u (Float)

# ---- 

# Constraints:
# 0 <= w <= 1
# 0 <= u <= 1

# --------------
def asym_bias(w,u):
	nuP = u
	omegaP = w
	if (nuP >= 1/2):
		nuP = 1-nuP
	if ( nuP < 1/2 - sqrt( omegaP*(1-omegaP) ) ):
		Delta = pow(1-2*nuP,2) - 4*omegaP*(1-omegaP)
		z = ( 1-2*nuP - sqrt(Delta) )/( 2*(1-omegaP) )
		kraw = (nuP*log2(1-z) + (1-nuP)*log2(1+z) - omegaP*log2(z))
	else:
		kraw = (1/2)*(1+ H2(omegaP) - H2(nuP))
	return kraw - asym_binomial(1.0, omegaP)


def expoAsymptKrawNormalisation(omega,nu):
    nuP = nu
    omegaP = omega
    if (nuP >= 1/2):
        nuP = 1-nuP
    if ( nuP < 1/2 - sqrt( omegaP*(1-omegaP) ) ):
        Delta = pow(1-2*nuP,2) - 4*omegaP*(1-omegaP)
        z = ( 1-2*nuP - sqrt(Delta) )/( 2*(1-omegaP) )
        return (nuP*log2(1-z) + (1-nuP)*log2(1+z) - omegaP*log2(z))
    else:
        return (1/2)*(1+ H2(omegaP) - H2(nuP))

def range_kraw(norm,omega):
	return norm*0.5*H2(omega/norm),norm*H2(omega/norm)


def kraw_inverse(norm,omega,y):
	mi,ma = range_kraw(norm,omega)
	if (y < mi or y > ma):
		return -1
	return brentq(lambda x : norm*expoAsymptKrawNormalisation(omega/norm,x/norm) - y,a=0,b=norm/2,maxiter=1000,rtol=1e-15,xtol=1e-15)


############# END UTILITIES FUNCTIONS


############# BEGIN BJMM12 FUNCTIONS

# complexity_BJMM12(R, w, l_1,l_2,p_1,p_2)
# Asymptotic exponent for computing parity-checks of relative weight w of a code of rate R with BJMM12 technique (Equation (5.4) in Section 5 of https://eprint.iacr.org/2022/1000.pdf)

# Inputs: R,w,l_1,l_2,p_1,p_2 (Float)


# Constraints:
# See section 4.2 and constraints_BJMM12 function

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
# BJM12 parameters (l_1,l_2,p_1,p_2) must verify constraints (specified in Proposition A.8) in order for BJMM12 method to compute all parity-checks of relative weight w of a code of rate R

# Output: True if constraints are verified (Up to tolerance)
#		  False else

# Inputs: R,w,l_1,l_2,p_1,p_2,tolerance (Float)

# --------------
def constraints_BJMM12(R,w,l_1,l_2,p_1,p_2, tolerance):
	assert(0 <= l_1)
	assert(l_1 <= l_2)
	assert(l_2 <= R)

	assert(0 <= p_1)
	assert(0 <= p_2)


	assert(((p_1 - p_2/2) >= 0) and ((p_1 - p_2/2) <= 1.0-p_2))
	assert(((p_2 - w/2) >= 0) and ((p_2 - w/2) <= 1.0-w))


	
	representation_1 = ((asym_binomial(1.0-p_2,p_1 - p_2/2) + p_2))
	assert(l_1 <= (representation_1 + tolerance)) # (38)

	representation_2 = ((asym_binomial(1.0-w,p_2 - w/2) + w))
	assert(l_2 <= (representation_2 + tolerance)) # (38)


	return True

############# END BJMM12 FUNCTIONS

# Asymptotic complexity exponent of Dumer decoder to return *all* solutions to the decoding problem at relative distance t in a code of rate R as given by Proposition C.6
def complexity_asym_ISD_Dumer(R, t, l, w):
		assert(l <= (1-R))
		assert(w <= (R+l)/2.0)
		bet = (R+l)*H2(w/(R+l)) + (1.0-R-l)*H2((t-w)/(1.0-R-l)) - H2(t)

		S0 = ((R+l)/2)*H2(w/(R+l))
		S1 = 2*S0 - l
		T_sub = max(S0,S1)
		#Number_solution = max(H2(t) - (1-R),0.0)
		memory = S0
		return np.array([T_sub - bet,memory])


############# BEGIN DOUBLE RLPN FUNCTIONS


# Compute the expected number of candidates as given by $\nu_{candidate}$ of Proposition C.4 (And Proposition 5.2)
def number_candidates(R,sigma, mu, omega, taux,kaux, tau,v_eps):
	v_max = -1
	i_max = -1
	j_max = -1

	limit_max = (1-sigma)*expoAsymptKrawNormalisation(omega/(1-sigma),mu/(1-sigma)) + (sigma)*expoAsymptKrawNormalisation(taux/(sigma),(tau-mu)/(sigma))
	abs_eps = limit_max*v_eps

	nb_step_aux = 10000

	step_aux = (sigma/2)/nb_step_aux

	for j in np.arange(0,sigma/2,step_aux):
		bias_j = (sigma)*expoAsymptKrawNormalisation(taux/(sigma),j/(sigma))
		target = limit_max - bias_j - abs_eps
		i = kraw_inverse((1-sigma),omega,target)
		if(i != -1):
			N_ij = (1-sigma)*H2(i/(1-sigma)) - (1-R) + sigma*H2(j/sigma)
			if( N_ij > v_max):
				v_max = N_ij
				i_max = i
				j_max = j
	return (v_max,i_max,j_max)


# --------------
# constraints_Double_RLPN(R, s, k_aux, u, w, t, tolerance)
# Double RLPN parameters must verify constraints

# Output: True if constraints are verified (Up to tolerance)
#		  False else

# Inputs: R, s, k_aux, u, w, t, tolerance (Float)

# --------------
def constraints_Double_RLPN(R, s, k_aux, t_aux, u, w, t, tolerance):

	assert(0 < k_aux and k_aux <= s and s <= R)

	assert(u <= 1.0-s)
	assert(t-u >= 0 and t-u <= s)
	assert(w >= 0 and w <= (1.0-u)/2.0)

	#tau = GV(s,k_aux)
	#rel_GV = s*Hinv(1-(k_aux/s))
	# Constraints bias 

	assert(asym_binomial(1.0-s,w) + asym_binomial(s,t_aux) <= R + tolerance) # (Constraint C.10)
	assert(asym_binomial(s,t_aux) <= s - k_aux + tolerance) # (Constraint C.11)

	len_punctured_code = 1.0 - s
	R_punctured_code = R - s

	eps_1 = len_punctured_code*asym_bias(u/len_punctured_code,w/len_punctured_code)

	
	eps_2 = s*asym_bias(t_aux/s,(t-u)/s)

	eps = eps_1 + eps_2

	number_avaible_parity_checks =asym_binomial(len_punctured_code,w) - R_punctured_code

	assert(-2*eps <= number_avaible_parity_checks + tolerance) # (Constraint C.9)

	return True

# --------------
# complexity_Double_RLPN(R, s, k_aux, u, w, t, tolerance)
# Asymptotic complexity exponent of the Double RLPN decoder (see Section 3.2, Algorithm 1) when using BJMM12 technique to compute low-weight parity checks

# Inputs: R, s, k_aux, u, w, t, tolerance (Float)

# Constraints:
# See Section 3.3, Proposition 4 and constraints_Double_RLPN function

# --------------
def complexity_Double_RLPN_BJMM12(R, s, k_aux,t_aux, u, w, l_1, l_2, p_1, p_2, t, N_aux, Decode_dumer_l,Decode_dumer_w,SolveSubProblem_Dumer_l,SolveSubProblem_Dumer_w, tolerance):
	constraints_Double_RLPN(R, s, k_aux, t_aux, u, w, t, tolerance)
	
	# Exponent of the number of expected needed iterations for the bet on the weight of e to be verified (e of weight u on N)
	_asym_number_iteration = asym_binomial(1.0,t) - (asym_binomial(s,t-u) + asym_binomial(1.0-s,u))

	# Complexity exponent of Fast Fourier Transform
	_complexity_asym_FFT = k_aux

	# Complexity exponent of BJMM12 technique to compute all parity-checks of weight w
	# Assert that BJMM12 parameters verifies the constraints

	len_punctured_code = 1.0 - s
	R_punctured_code = R - s

	assert(len_punctured_code > 0 and len_punctured_code < 1)

	rel_BJMM12_R = R_punctured_code/len_punctured_code
	assert(rel_BJMM12_R > 0 and rel_BJMM12_R < 1)

	rel_BJMM12_w = w/len_punctured_code
	assert(rel_BJMM12_w > 0 and rel_BJMM12_w < 1)

	rel_BJMM12_l_1 = l_1/len_punctured_code
	assert(rel_BJMM12_l_1 > 0 and rel_BJMM12_l_1 < 1)

	rel_BJMM12_l_2 = l_2/len_punctured_code
	assert(rel_BJMM12_l_2 > 0 and rel_BJMM12_l_2 < 1)

	rel_BJMM12_p_1 = p_1/len_punctured_code
	assert(rel_BJMM12_p_1 > 0 and rel_BJMM12_p_1 < 1)

	rel_BJMM12_p_2 = p_2/len_punctured_code
	assert(rel_BJMM12_p_2 > 0 and rel_BJMM12_p_2 < 1)

	rel_BJMM12_tolerance = tolerance
	
	constraints_BJMM12(rel_BJMM12_R,rel_BJMM12_w,rel_BJMM12_l_1,rel_BJMM12_l_2,rel_BJMM12_p_1,rel_BJMM12_p_2, rel_BJMM12_tolerance)

	_complexity_asym_BJMM12, _memory_asym_BJMM12 = len_punctured_code*complexity_asym_BJMM12(rel_BJMM12_R,rel_BJMM12_w,rel_BJMM12_l_1,rel_BJMM12_l_2,rel_BJMM12_p_1,rel_BJMM12_p_2)

	_number_parity_checks= asym_binomial(len_punctured_code,w) - R_punctured_code

	_number_samples = _number_parity_checks + asym_binomial(s,t_aux) - (s-k_aux)

	_number_candidate = number_candidates(R,s, u, w, t_aux,k_aux, t,0.0)[0]

	_number_solution_Decode_Dumer = max(asym_binomial(s,t-u) - N_aux*k_aux,0)
	#Complexity of the two ISD dumer routines

	# Time and space complexity of the call to ISD dumer
	_alpha_DecodeDumer, _beta_DecodeDumer = s*complexity_asym_ISD_Dumer(1-(N_aux * k_aux)/s, (t-u)/s,Decode_dumer_l/s,Decode_dumer_w/s)
	_alpha_SolveSubProblem, _beta_SolveSubProblem = (1.0-s)*complexity_asym_ISD_Dumer((R-s)/(1-s), u/(1-s),SolveSubProblem_Dumer_l/(1.0-s),SolveSubProblem_Dumer_w/(1.0-s))

	_alpha_ISD  = max(_alpha_DecodeDumer , _number_solution_Decode_Dumer + _alpha_SolveSubProblem)
	_beta_ISD  = max(_beta_DecodeDumer ,  _beta_SolveSubProblem)

	_complexity_asym_Double_RLPN_BJMM12 = _asym_number_iteration +  max(_complexity_asym_BJMM12, _complexity_asym_FFT , _number_samples, _number_candidate*N_aux + _alpha_ISD)
	_memory_asym_Double_RLPN_BJMM12 = max(_memory_asym_BJMM12,k_aux,_beta_ISD)

	return np.array([_complexity_asym_Double_RLPN_BJMM12, _memory_asym_Double_RLPN_BJMM12])


############# END DOUBLE RLPN FUNCTIONS


############# LOAD CSV FILE
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
Tab = []
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
	lambda_1 = rows[i][idx_l_1]
	lambda_2 = rows[i][idx_l_2]
	pi_1 = rows[i][idx_p_1]
	pi_2 = rows[i][idx_p_2]

	#We decode the auxilary code at Gilbert-Varshamov distance
	t_aux = sigma*Hinv(1-(R_aux/sigma))
	# N_aux= 1 works for every parameter
	N_aux = 1


	
	parameters_Decode_Dumer = [[0.0018644849, 0.0006728533], [0.0037552195, 0.0013313797], [0.0049190036, 0.0017252408], [0.0062609503, 0.0021727683], [0.0077161297, 0.0026497793], [0.0099420333, 0.003384119], [0.0113651531, 0.0038449976], [0.0160265599, 0.0053794], [0.0200636221, 0.0066709164], [0.024082797, 0.0079699156], [0.0262912985, 0.0086634656], [0.0281614211, 0.0092318443], [0.0313795416, 0.0102177289], [0.0347020738, 0.0112264855], [0.0367416333, 0.0118893625], [0.0398680174, 0.012875869], [0.0403934835, 0.0129977195], [0.0429431973, 0.0137006519], [0.0446051422, 0.0142209134], [0.0489945799, 0.0154037109], [0.0478411498, 0.0151500738], [0.0471415351, 0.0148234542], [0.0500542333, 0.0157820617], [0.0517908392, 0.0162376337], [0.0528462258, 0.0165301157], [0.053238913, 0.016584109], [0.054990691, 0.0171227944], [0.0550108602, 0.0169928407], [0.0548827609, 0.0168335265], [0.0564036585, 0.0173555673], [0.0549532097, 0.0165591457], [0.0522434819, 0.0154727072], [0.0518684322, 0.0152893187], [0.0511365423, 0.0148434428], [0.0526142387, 0.0153923117], [0.0498311218, 0.0142896911], [0.0490169348, 0.0140442867], [0.0491869436, 0.0139926977], [0.0497309277, 0.0140956116], [0.0506441316, 0.0144119106], [0.047434143, 0.0132062802], [0.047604162, 0.0131948032], [0.0470497293, 0.0129770461], [0.0454139744, 0.0123798292], [0.0443707875, 0.0119983588], [0.0437863882, 0.0117486411], [0.0428254314, 0.0114339635], [0.0427974132, 0.0113824557], [0.0417024248, 0.0110001173], [0.0405699016, 0.0106102781], [0.0407816577, 0.0106670761], [0.0398189654, 0.0103212538], [0.0393405516, 0.0101543249], [0.0380335976, 0.009722741], [0.0366941395, 0.0092870953], [0.0367605996, 0.0092875827], [0.0349617414, 0.0087179775], [0.0347198556, 0.0086351403], [0.0337307634, 0.0083169273], [0.032822591, 0.008041232], [0.0322907621, 0.0078681999], [0.0323227375, 0.0078767104], [0.0309820087, 0.0074662309], [0.0297278535, 0.0070883891], [0.0295709612, 0.0070374589], [0.0285025377, 0.0067271619], [0.0277363138, 0.0064970532], [0.0272426194, 0.0063579877], [0.026206634, 0.006059416], [0.0254210269, 0.0058365334], [0.0244879205, 0.005573547], [0.023688297, 0.0053502829], [0.0232556784, 0.0052314047], [0.0224334839, 0.0050072247], [0.0216345902, 0.0047894283], [0.0206254904, 0.0045212467], [0.0200517179, 0.0043726321], [0.0193054172, 0.004174936], [0.0184065164, 0.003941645], [0.0178245408, 0.0037934883], [0.0173793596, 0.0036801658], [0.0160975122, 0.003358933], [0.0154941093, 0.0032104385], [0.0147260974, 0.0030221393], [0.0139787823, 0.0028418541], [0.013094396, 0.0026305268], [0.0122276618, 0.0024264119], [0.0112023105, 0.0021895992], [0.0106746244, 0.0020695676], [0.0098372537, 0.0018812306], [0.0090332994, 0.0017035866], [0.008060961, 0.0014926663], [0.0072717288, 0.0013250006], [0.0062569485, 0.0011142542], [0.0054147026, 0.0009438458], [0.0042244787, 0.0007107243], [0.003408338, 0.0005567886], [0.0025087145, 0.0003936422], [0.0013257206, 0.0001923705]]
	parameters_SolveSubProblem = [[0.0008496007, 0.0003062593], [0.0017058591, 0.0006031176], [0.0025480648, 0.0008875807], [0.0033641223, 0.0011572553], [0.0041519133, 0.0014126065], [0.0048914617, 0.0016476907], [0.0056253369, 0.0018774265], [0.0062069259, 0.002053836], [0.0067509866, 0.0022160943], [0.0072458526, 0.0023601763], [0.0077873194, 0.0025177576], [0.0083116749, 0.0026682969], [0.0087093241, 0.0027772975], [0.0090461837, 0.0028662072], [0.0094705811, 0.0029805199], [0.0097570463, 0.0030511433], [0.0102303017, 0.0031794151], [0.010439076, 0.0032265579], [0.0107688227, 0.0033081993], [0.0105483712, 0.0032263549], [0.0112921962, 0.0034298447], [0.0117603243, 0.0035532723], [0.0118652986, 0.0035634044], [0.011959655, 0.0035738878], [0.0121825589, 0.0036208285], [0.0124541929, 0.003682496], [0.0125498061, 0.0036904866], [0.0126522532, 0.0037008108], [0.0128168071, 0.0037300632], [0.0131289683, 0.0038036111], [0.0128795907, 0.0037154707], [0.0129486039, 0.0037135313], [0.0131542463, 0.0037524992], [0.0128850486, 0.003663358], [0.0133320268, 0.0037704649], [0.0131742887, 0.0037066175], [0.0136187041, 0.0038044808], [0.0134956767, 0.0037578552], [0.0134583887, 0.003735321], [0.0137846104, 0.0038061485], [0.013478843, 0.0037030104], [0.013420055, 0.0036751411], [0.0135053116, 0.0036797498], [0.0133940039, 0.003631108], [0.0133651798, 0.0036050801], [0.0132017973, 0.0035499717], [0.0133643582, 0.0035687173], [0.0133112391, 0.0035434936], [0.013233847, 0.0035052905], [0.0131322309, 0.0034610172], [0.0132981806, 0.0034869987], [0.0130797065, 0.0034192655], [0.0130983301, 0.0034069051], [0.0129258862, 0.0033452987], [0.0127259227, 0.0032771396], [0.0127526943, 0.003273756], [0.0124214147, 0.0031729199], [0.01247287, 0.0031698303], [0.0122473315, 0.0031029539], [0.0121808299, 0.0030644075], [0.0120791635, 0.0030292888], [0.0122489542, 0.0030497912], [0.0119227034, 0.0029594123], [0.0116098907, 0.0028728083], [0.0116819633, 0.0028754275], [0.0114802908, 0.0028055629], [0.0113006322, 0.0027581972], [0.0112860184, 0.0027290022], [0.0110295463, 0.0026582935], [0.0108676871, 0.0026052731], [0.0106430854, 0.0025428751], [0.0104749669, 0.002494122], [0.010442969, 0.0024676887], [0.0102504104, 0.0024087761], [0.0101073446, 0.0023713926], [0.0098166882, 0.0022903111], [0.0096114703, 0.0022113287], [0.0095347963, 0.0021942967], [0.0093555216, 0.0021490726], [0.0091910079, 0.0020899689], [0.009140407, 0.0020612064], [0.0087363235, 0.0019662735], [0.0084866799, 0.0018861029], [0.0083459851, 0.0018465357], [0.0080413961, 0.0017598654], [0.0078806216, 0.0017200742], [0.0078127591, 0.0017037126], [0.0074549076, 0.0016173788], [0.0072170832, 0.0015408114], [0.007054184, 0.0014971745], [0.0067037401, 0.0014050338], [0.0064687133, 0.0013496376], [0.0061609238, 0.0012678949], [0.0058724559, 0.0012014818], [0.0055466748, 0.0011173313], [0.005185586, 0.0010448555], [0.0046859331, 0.0009186056], [0.0040964332, 0.0007744445], [0.0034916107, 0.0006440705]]
	
	#Parameters for the Dumer decoder in the Decode-Dumer procedure (Line 5 of Algorithm 4.3))
	Decode_dumer_l = sigma*parameters_Decode_Dumer[i][0]
	Decode_dumer_w = sigma*parameters_Decode_Dumer[i][1]

	#Parameters for the Dumer decoder in the SolveSubProblem procedure (Line 6 of Algorithm 4.3))
	SolveSubProblem_Dumer_l = (1.0-sigma)*parameters_SolveSubProblem[i][0]
	SolveSubProblem_Dumer_w = (1.0-sigma)*parameters_SolveSubProblem[i][1]
	
	Comp,Memory = complexity_Double_RLPN_BJMM12(R, sigma, R_aux, t_aux, mu, omega, lambda_1, lambda_2, pi_1, pi_2, tau, N_aux,Decode_dumer_l ,Decode_dumer_w,SolveSubProblem_Dumer_l,SolveSubProblem_Dumer_w, tolerance)
	print("Rate: " + "{:.5f}".format(R) + ";	" + "Complexity: " + "{:.5f}".format(Comp)  + ";	" + "Memory: " + "{:.5f}".format(Memory) )