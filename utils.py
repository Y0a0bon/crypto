############################
#         utils.py         #
############################

# Contains the following functions :
# test_utils()
# euclide(a, b, verbose)
# calcpgcd(a, b, verbose)
# mod_inv(b, n, verbose)
# ind_euler(N, verbose)
# fast_exp(m, e, n, verbose)
# phi(p, q)
# rabin_miller(n)
# is_prime(n)
# iroot(k, n)
# int_list(a, b)
# baby_giant_step(p, g, verbose)
# pollard_rho(n, x1, f, verbose)
# pollard_rho_lambda(x)


import os, sys
import math, random


#
# Test function for utils.py
#
def test_utils():
	print("********************************************************************"
              "\n Test of functions in utils.py:\n\n")
	if(euclide(141, 255, 0) == [38, -21, 3]):
		print(" Test of euclide(a, b, verbose) ...\t\t\t\tOK")
	else:
		print(" Test of euclide(a, b, verbose) ...\t\t\t\tFAILED")
	if(calcpgcd(45,89, 0) == 1 and calcpgcd(78, 1547, 0) == 13):
		print(" Test of calcpgcd(a, b, verbose) ...\t\t\t\tOK")
	else:
		print(" Test of calcpgcd(a, b, verbose) ...\t\t\t\tFAILED")
	if(mod_inv(547, 696964823868233, 0) == 44595555457748):		  
		print(" Test of mod_inv(b, n verbose) ...\t\t\t\tOK")
	else:
		print(" Test of mod_inv(b, n, verbose) ...\t\t\t\tFAILED")
	if(fast_exp(1254894561, 44595555457748, 696964823868233, 0) == 269024829525699 and fast_exp(11, 13, 19, 0) == pow(11, 13)%19):
		print(" Test of fast_exp(m, e, n, verbose) ...\t\t\t\tOK")
	else:
		print(" Test of fast_exp(m, e, n, verbose) ...\t\t\t\tFAILED")
	if(phi(70457971, 9891923) == 696964743518340):
		print(" Test of phi(p, q) ...\t\t\t\t\t\tOK")
	else:
		print(" Test of phi(p, q) ...\t\t\t\tFAILED")
	if(is_prime(999982) == False and is_prime(999979)== True and is_prime(516484451321516846844515313513848486461351351468498435159) == True):
		print(" Test of is_prime(n) ...\t\t\t\t\tOK")
	else:
		print(" Test of is_prime(n) ...\t\t\t\t\tFAILED")
	if(baby_giant_step(31, 3, 0) == 25):
		print(" Test of baby_giant_step(p, g, verbose) ...\t\t\tOK")
	else:
		print(" Test of baby_giant_step(p, g, verbose) ...\t\t\tFAILED")
	if(pollard_rho(15770708441, 1, pollard_rho_lambda, 0) == 135979):
		print(" Test of pollard_rho(n, x1, f, verbose) ...\t\t\tOK")
	else:
		print(" Test of pollard_rho(n, x1, f, verbose) ...\t\t\tFAILED")
	print("\n\n End of tests\n"
              "********************************************************************")


#
# Apply Euclide algorithm
#
# a 1st number
# b 2nd number
# v verbose if not 0 (doesn't work ?)
#
# return [u, v, pgcd]
#
def euclide(a, b, verbose):
	r=-1
	pgcd=0
	q=0
	a1=a
	b1=b
	u=[]
	v=[]
	u.append(0)
	u.append(1)
	v.append(1)
	v.append(0)
	i=1
	while(r!=0):
		u.append(u[i-1]-q*u[i])
		v.append(v[i-1]-q*v[i])
		i+=1
		q=int(a1/b1)
		r=a1%b1
		a1=b1
		b1=r
		if(b1!=0):
			pgcd=b1
	if(verbose!=0):
		print(str(a)+' * '+str(u[i])+' - '+str(b)+' * '+str(v[i])+' = '+str(pgcd)) 
	return [u[i],v[i],pgcd]


#
# Compute the gcd of a and b
#
def calcpgcd(a, b, verbose):
	pgcd = euclide(a, b, verbose)
	if(verbose!=0):
		print('gcd('+str(a)+', '+str(b)+') = '+str(pgcd[2]))
	return pgcd[2]

#
# Compute inverse of b mod n
#
def mod_inv(b, n, verbose):
	n0 = n
	b0 = b
	t0 = 0
	t = 1 
	q  = int(math.floor(n0 / b0)) # floor useless here ?
	#print(str(n0)+' / '+str(b0)+' = '+str(q))
	r = n0 - q * b0
	while(r > 0):
		temp = t0 - q * t
		if(temp >= 0):
			temp =  temp % n
		else:
			temp = n - ((-temp) % n)
		t0 = t
		t = temp
		n0 = b0
		b0 = r
		q = int(math.floor(n0 / b0))
		r = n0 -q * b0
	if(b0 != 1):
		print('Pas de solution')
		return None
	else:
		if(verbose != 0):
			print(str(b)+'^-1 mod '+str(n)+' = '+str(t))
		return t


#
# Euler indicative
#
def ind_euler(N, verbose):
	nbr=0
	for i in range(2,N):
		if (N%i==0):
			if(verbose!=0):
				print(i, "n'est pas premier avec",N)
		nbr+=1
	if(verbose!=0):
		print("l'indicateur d'Euler est :",nbr)
	return nbr


#
# Fast exponentiation
#
def fast_exp(m, e, n, verbose):
	c , e_sv , m_tmp = 1 , e, m
	while(e_sv>0):
		if(e_sv % 2 == 0):
			e_sv//=2
		else:
			c = (c * m_tmp) % n
			e_sv = (e_sv - 1)//2
		#print("e = " + str(e_sv) + " and m = " + str(m_tmp))
		m_tmp = (m_tmp * m_tmp)
		#print("m = " + str(m_tmp))
		m_tmp = m_tmp % n
	if(verbose != 0):
		print(str(m) + " ^ " + str(e) + " % " + str(n) + " = " + str(c))
	return c


#
# Compute phi(n) = phi(pq)
#
def phi(p, q):
	return(p-1)*(q-1)


#
# Compute prime with Rabin-Miller primality test
#
def rabin_miller(n):
		# Returns True if n is a prime number.

		s = n - 1
		t = 0
		while s % 2 == 0:
				# keep halving s while it is even (and use t
				# to count how many times we halve s)
				s = s // 2
				t += 1

		for trials in range(5): # try to falsify n's primality 5 times
				a = random.randrange(2, n - 1)
				v = pow(a, s, n)
				if v != 1: # this test does not apply if v is 1.
					i = 0
					while v != (n - 1):
						if i == t - 1:
							return False
						else:
							i = i + 1
							v = (v ** 2) % n
		return True

#
# Return True if n is a prime number. This function does a quicker
# prime number check before calling rabinMiller().
def is_prime(n):
		if (n < 2):
				return False # 0, 1, and negative numbers are not prime

		# About 1/3 of the time we can quickly determine if n is not prime
		# by dividing by the first few dozen prime numbers. This is quicker
		# than rabinMiller(), but unlike rabinMiller() is not guaranteed to
		# prove that a number is prime.
		low_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997]

		if n in low_primes:
				return True

		# See if any of the low prime numbers can divide n
		for prime in low_primes:
				if (n % prime == 0):
					return False

		# If all else fails, call rabinMiller() to determine if n is a prime.
		return rabin_miller(n)


#
# Compute k nth root
#
def iroot(k, n):
	u, s = n, n+1
	while u < s:
		s = u
		t = (k-1) * s + n / pow(s, k-1)
		u = t / k
	return u

#print(iroot(2, 16))


#
# List of integers from a to b
#
def int_list(a, b):
	return [i for i in range(a, b+1)]


#
# Function applying Baby Step Giant Step algorithm
#
def baby_giant_step(p, g, verbose):
	m = isqrt(p)+1
	print(m)
	g_t = [0]*m
	h_t = [0]*m
	g_inv = 0
	for i in range(0, m):
		g_t[i] = fast_exp(g, i, p, verbose) % p
	if(verbose != 0):
		print("g ={0}".format(g_t))
	# g^(-m) mod p = (g ^ (-1)) ^ 6 mod p :
	g_inv_m = fast_exp(mod_inv(g, p, verbose), m, p, verbose)
	if(verbose != 0):
		print("g ^ (-m) mod p = {0}^(-{1}) mod {2} = {3}".format(g, m, p, g_inv_m))
	j = ret = cont = 0
	while(j < m and cont == 0):
		h_t[j] = m * fast_exp(g_inv_m, j, p, verbose) % p
		if(verbose != 0):
			print("h[{0}] = {1}".format(j, h_t[j]))
			print("h is {0} and g is {1}".format(h_t, g_t))
		if(h_t[j] in g_t):
			ret = j * m + g_t.index(h_t[j]) # NON OPIMISE
			if(verbose != 0):
				print("Found it in g_t[{0}] which is {1} ".format(g_t.index(h_t[j]), h_t[j]))
			cont = 1
		j = j+1
	return ret

#
# BSGS gen a, res b, order n
#
def BSGS(a, b, n, verbose):
    m = isqrt(n)+1
    print(m)
    g_t = [0]*m
    h_t = [0]*m
    for i in range(0,m):
        g_t[i] = fast_exp(a, i, n, verbose)
    if(verbose != 0):
        print("g ={0}".format(g_t))
    # g^(-m) mod n = (g ^ (-1)) ^ m mod n :
    g_inv_m = fast_exp(mod_inv(a, n, verbose), m, n, verbose)
    if(verbose != 0):
        print("g ^ (-m) mod n = {0}^(-{1}) mod {2} = {3}".format(a, m, n, g_inv_m))
    j = ret = cont = 0
    while(j < m and cont == 0):
        h_t[j] = m * fast_exp(g_inv_m, j, n, verbose) % n
        if(verbose != 0):
            print("h[{0}] = {1}".format(j, h_t[j]))
            print("h is {0} and g is {1}".format(h_t, g_t))
        if(h_t[j] in g_t):
            ret = j * m + g_t.index(h_t[j]) # NON OPIMISE
            if(verbose != 0):
                print("Found it in g_t[{0}] which is {1} ".format(g_t.index(h_t[j]), h_t[j]))
            cont = 1
        j = j+1
    return ret

        
#
# Pollard's Rho algorithm
#
def pollard_rho(n, x1, f, verbose):
        x = x1
        y = f(x) % n
        p = calcpgcd(y - x, n, verbose)
        while(p == 1):
                x = f(x) % n
                y = f(f(y)) % n
                p = calcpgcd(y - x, n, verbose)
        if(p == n):
                if(verbose != 0):
                        print("No solution found !")
                return None
        return p

#
# Function to use in Pollard's Rho algorithm
#
def pollard_rho_lambda(x):
        return (x**2 + 1)


#
# Function to generate a random prime number
#
def gen_prime(n_max=3513514684985159):
        random.seed()
        mod = 2
        while(mod == 0 and is_prime(mod) == False):
                mod = random.randint(100000, n_max)



#
# Mod_inv for ECDSA
#
def mod_inv_ecdsa(k, p):
        if k == 0:
            raise ZeroDivisionError('division by zero')
        if k < 0:
            # k ** -1 = p - (-k) ** -1  (mod p)
            return p - self.mod_inv_ecdsa(-k, p)

        # Extended Euclidean algorithm.
        s, old_s = 0, 1
        t, old_t = 1, 0
        r, old_r = p, k

        while r != 0:
            quotient = old_r // r
            old_r, r = r, old_r - quotient * r
            old_s, s = s, old_s - quotient * s
            old_t, t = t, old_t - quotient * t

        gcd, x, y = old_r, old_s, old_t

        assert gcd == 1
        assert (k * x) % p == 1

        return x % p
                
        
#
# CRT
#
def crt(a, moda, b, modb):
        r = euclide(moda, modb)
        u = r[0]
        v = r[1]
        x =(a*v*modb) + (b*u*moda)
        return [x, moda*modb]

#
# multiple crt
#
def crt_list(values, mods):
        lg = len(values)
        while(len >1):
                i = 0
                while((i*2) < length):
                        if((i*2)+1 == length):
                                values[i] = values[i*2]
                                mods[i] = mods[i*2]
                        else:
                                res = crt(values[i*2], mods[i*2], values[i*2+1], mods[i*2+1])
                        i = i+1
                lg =(lg >> 1) + lg%2
        return [values[0], mods[0]]


def cont_frac(x : int, y : int) -> int:
	"""
	Convert a rational x/y fraction into
	a list of partial quotients [a0, ... , an]
	"""

	L = []
	a = x
	b = y
	r = log2(y) #Maximum partial quotients

	while(a%b != 0 and r > 0 ):
		q = a // b
		a_t = a
		a = b
		b = a_t % b
		L.append(q)
		r = r-1

	L.append(a // b)
	if(r>0) :
		L.append(b)

	return L


def frac_seq(L : list) -> list:
	"""
	Give the fractional sequence from a list of partial quotient
	"""
	
	seq = []

	for i in range(len(L)) :
		l = L[0:i]
		if len(l) == 0 :
			continue
		elif len(l) == 1 :
			seq.insert(i,(L[0],1))
		else :
			seq.insert(i,compute_frac(l))

	return seq


def compute_frac(l : list) -> tuple:
	"""
	compute symbolically a sequence of fractions
	"""	
	i = len(l)-2

	c = 1
	a = l[i]
	b = l[i+1]
	

	while ( i > 0 ) :
		a = l[i]
		a = a * b + c
		c = b
		b = a
		i = i - 1
	
	a = l[0]
	a = a * b + c

	return (a,b)


def isqrt(n):
	'''
	Calculates the integer square root
	for arbitrary large nonnegative integers
	'''
	if n < 0:
		raise ValueError('square root not defined for negative numbers')
    
	if n == 0:
		return 0

	a = int(math.log2(n)) + 1 // 2
	b = int(math.log2(n)) + 1 % 2

	x = 2**(a+b)
	while True:
		y = (x + n//x)//2
		if y >= x:
			return x
		x = y        




