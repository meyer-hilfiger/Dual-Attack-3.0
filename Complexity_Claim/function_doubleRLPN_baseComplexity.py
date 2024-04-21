# Proposition number refers to https://eprint.iacr.org/archive/2023/1852/1701452846.pdf

# These functions computes using Proposition 9 the complexity exponent of the Double RLPN algorithm when using BJMM12 technique to compute low-weight parity checks.
# Thsese functions verifies that the parameters meet the constraints of Proposition 9.

# The parameters of the algorithm are taken from the file doubleRLPN_BJMM12.csv


# $R$ = $\frac{k}{n}$ the rate of the code to decode?
# $TimeComplexity = $\frac{log2(Time complexity of Double RLPN with BJMM12 decoder)}{n}$ when n -> Infinity. This value is only here for reference, this script will recompute it from the parameters of the algorithm.
# $\sigma$ = $\frac{s}{n}$
# $R_{\text{aux}}$ = $\frac{k_aux}{n}$
# $\nu$ = $\frac{u}{n}$
# $\omega$ = $\frac{w}{n}$

# Internal parameters of the BJMM12 method to compute low-weight parity-checks. They are expressed relatively to n:
# $\lambda_1$ = $\frac{l_1}{n}$
# $\lambda_2$ = $\frac{l_2}{n}$
# $\pi_1$ = $\frac{p_1}{n}$
# $\pi_2$ = $\frac{p_2}{n}$
# $\tau$ = $\frac{t}{n}$ (Here it is equal to the relative Gilbert-Varshamov distance of a code of rate R, namely, $Hinv(1-R)$ where Hinv is the inverse of the binary entropy function)

from function_utilitaries import *
# Asymptotic complexity exponent of Dumer decoder to return *all* solutions to the decoding problem at relative distance t in a code of rate R as given by Proposition 10
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


# Compute the expected number of candidates as given by Proposition 5 (and asymptotically by $\nu_{\text{candidates}}$ in Proposition 9)
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
# Double RLPN parameters must verify constraints (Proposition 9)

# Output: True if constraints are verified (Up to tolerance)
#		  False else

# Inputs: R, s, k_aux, u, w, t, tolerance (Float)

# --------------
def constraints_Double_RLPN(R, s, k_aux, t_aux, u, w, t, tolerance):

	assert(0 < k_aux and k_aux <= s and s <= R)

	assert(u <= 1.0-s)
	assert(t-u >= 0 and t-u <= s)
	assert(w >= 0 and w <= (1.0-u)/2.0)


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
# Asymptotic complexity exponent of the Double RLPN decoder when using BJMM12 technique to compute low-weight parity checks

# Inputs: R, s, k_aux, u, w, t, tolerance (Float)

# Constraints:
# See Proposition 9 and constraints_Double_RLPN function

# --------------
def doubleRLPN_baseComplexity(R, s, k_aux,t_aux, u, w, t, N_aux, complexity_compute_parity_check, tolerance):
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

	rel_punctured_code_R = R_punctured_code/len_punctured_code
	assert(rel_punctured_code_R > 0 and rel_punctured_code_R < 1)

	rel_punctured_code_w = w/len_punctured_code
	assert(rel_punctured_code_w > 0 and rel_punctured_code_w < 1)

	_alpha_complexity_compute_parity_check, _beta_complexity_compute_parity_check = complexity_compute_parity_check

	_number_parity_checks= asym_binomial(len_punctured_code,w) - R_punctured_code

	_number_samples = _number_parity_checks + asym_binomial(s,t_aux) - (s-k_aux)

	_number_candidate = number_candidates(R,s, u, w, t_aux,k_aux, t,0.0)[0]

	_number_solution_Decode_Dumer = max(asym_binomial(s,t-u) - N_aux*k_aux,0)
	#Complexity of the two ISD dumer routines

	_,rel_Decode_dumer_w,rel_Decode_dumer_l = optimize_dumer(1-(N_aux * k_aux)/s, (t-u)/s)
	_,rel_SolveSubProblem_Dumer_w,rel_SolveSubProblem_Dumer_l = optimize_dumer((R-s)/(1-s), u/(1-s))


	# Time and space complexity of the call to ISD dumer
	_alpha_DecodeDumer, _beta_DecodeDumer = s*complexity_asym_ISD_Dumer(1-(N_aux * k_aux)/s, (t-u)/s,rel_Decode_dumer_l,rel_Decode_dumer_w)
	_alpha_SolveSubProblem, _beta_SolveSubProblem = (1.0-s)*complexity_asym_ISD_Dumer((R-s)/(1-s), u/(1-s),rel_SolveSubProblem_Dumer_l,rel_SolveSubProblem_Dumer_w)
	
	"""
	_alpha_DecodeDumer, _beta_DecodeDumer = s*isd(1-(N_aux * k_aux)/s, (t-u)/s,"B"),0.0
	_alpha_SolveSubProblem, _beta_SolveSubProblem = (1.0-s)*isd((R-s)/(1-s), u/(1-s),"B"),0.0
	"""
	_alpha_ISD  = max(_alpha_DecodeDumer , _number_solution_Decode_Dumer + _alpha_SolveSubProblem)
	_beta_ISD  = max(_beta_DecodeDumer ,  _beta_SolveSubProblem)

	_alpha_process_set_candidates = _number_candidate*N_aux + _alpha_ISD
	_beta_process_set_candidates = _beta_ISD

	_complexity_iteration = max(_alpha_complexity_compute_parity_check, _complexity_asym_FFT , _number_samples,_alpha_process_set_candidates)
	_memory_iteration = max(_beta_complexity_compute_parity_check,k_aux,_beta_process_set_candidates)
	
	if(_alpha_process_set_candidates == _complexity_iteration):
		print("Warning processing the set of false candidates dominate the time complexity... Change the subroutine to get a better time complexity?")
		print("Cost of checking false candidates : " + str(_alpha_process_set_candidates))
		print("Cost of core of one iteration of doubleRLPN" + str(max(_alpha_complexity_compute_parity_check, _complexity_asym_FFT , _number_samples)))
	if(_beta_process_set_candidates == _memory_iteration):
		print("Warning processing the set dominate the memory complexity... Change the subroutine to get a better memory complexity?")

	_complexity_asym_Double_RLPN = _asym_number_iteration + _complexity_iteration
	_memory_asym_Double_RLPN = _memory_iteration
	

	return np.array([_complexity_asym_Double_RLPN, _memory_asym_Double_RLPN])
############# END DOUBLE RLPN FUNCTIONS