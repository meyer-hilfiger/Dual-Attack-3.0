############# BEGIN UTILITIES FUNCTIONS
from math import log2,sqrt
from scipy.optimize import brentq
import numpy as np

from scipy.optimize import minimize_scalar

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
### OPTIMIZE ISD DUMER

def dumer_w(R,t,w):
	nb_solution = max(asym_binomial(1.0,t) - (1-R),0.0)
	f = lambda l : - (asym_binomial(R+l,w) + asym_binomial(1.0-(R+l),t-w) - asym_binomial(1.0,t)) + max(0.5*((R+l)*H2(w/(R+l))),((R+l)*H2(w/(R+l))) - l) - nb_solution
	l_max = min((1.0-R) + (w-t),1-R)
	l_min = max(0.0,w-R)
	res = minimize_scalar(f,bounds=(l_min,l_max),method='bounded')
	return [res.fun, res.x]
def optimize_dumer(R,t):
	f = lambda w : dumer_w(R,t,w)[0]
	w_min = max(t+R-1,0.0)
	w_max = t
	res = minimize_scalar(f,bounds=(w_min,w_max),method='bounded')
	complexity,l = dumer_w(R,t,res.x)
	return [complexity,res.x,l]