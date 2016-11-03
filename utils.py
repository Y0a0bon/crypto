############################
#         utils.py         #
############################

# Contains the following functions :
# euclide(a, b, verbose)
# calcpgcd(a, b, verbose)
# mod_inv(b, n, verbose)
# ind_euler(N, verbose)
# fast_exp(m, e, n, verbose)
# phi(p, q)
# rabin_miller(n)
# is_prime(n)
# test_utils()
# iroot(k, n)
# eratho(a, b, ch)
# int_list(a, b)
# baby_giant_step(p, g, verbose)
# rec_eratho(x)


import os, sys
import math, random


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
	while(e>0):
		if(e % 2 == 0):
			e//=2
		else:
			c = (c * m_tmp) % n
			e = (e - 1)//2
		m_tmp = (m_tmp * m_tmp) % n
	if(verbose != 0):
		print(str(m) + " ^ " + str(e_sv) + " % " + str(n) + " = " + str(c))
	return c


#
# Compute phi(n) = phi(pq)
#
def phi(p, q):
	return(p-1)*(q-1)


#
# Compute prime with Rabin-Miller algorithm
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
# Test function for utils.py
#
def test_utils():
	if(calcpgcd(45,89, 0) == 1 and calcpgcd(78, 1547, 0) == 13):
		print("Test de calcpgcd(a, b, verbose) ...\t\tOK")
	else:
		print("Test de calcpgcd(a, b, verbose) ...\t\tFAILED")
	if(mod_inv(547, 696964823868233, 0) == 44595555457748):		  
		print("Test de mod_inv(b, n verbose) ...\t\tOK")
	else:
		print("Test de mod_inv(b, n, verbose) ...\t\tFAILED")
	if(fast_exp(1254894561, 44595555457748, 696964823868233, 0) == 269024829525699 and fast_exp(11, 13, 19, 0) == pow(11, 13)%19):
		print("Test de fast_exp(m, e, n, verbose) ...\t\tOK")
	else:
		print("Test de fast_exp(m, e, n, verbose) ...\t\tFAILED")
	if(phi(70457971, 9891923) == 696964743518340):
		print("Test de phi(p, q) ...\t\t\t\tOK")
	else:
		print("Test de phi(p, q) ...\t\t\t\tFAILED")
	if(is_prime(999982) == False and is_prime(999979)== True and is_prime(516484451321516846844515313513848486461351351468498435159) == True):
		print("Test de is_prime(n) ...\t\t\t\tOK")
	else:
		print("Test de is_prime(n) ...\t\t\t\tFAILED")


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
	m = math.ceil(math.sqrt(p))
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



## ----------------------------
##      NOT TO BE USED
## ----------------------------

#
# List all prime numbers under n if 0, else the greatest
#
##def primes(n, ch): 
##	if n==2: return [2]
##	elif n<2: return []
##
##	s=list(range(3,n+1,2))
##	mroot = n ** 0.5
##	half=(n+1)/2-1
##	i=0
##	m=3
##	while(m <= mroot):
##	    if s[i]:
##		j=int((m*m-3)/2)
##		s[j]=0
##		while(j<half):
##		    s[j]=0
##		    j+=m
##		i=i+1
##		m=2*i+3
##	pr = [2]+[x for x in s if x]
##	if(ch==0):
##		return pr
##	else :
##		return pr[-1]


#
# Apply the Sieve of Eratosthenes
# Find prime numbers between a and b (1 is forbidden)
# 0 for list, else return the greatest one
# Less effective than primes()
#
def eratho(a, b, ch):
	if(ch==0):
		return rec_eratho(int_list(a,b))
	else:
		return rec_eratho(int_list(a,b))[-1] # wtf?

# print(primes(100000, 1))
# print(eratho(2, 100000, 1))


# print(int_list(1,100))


#
# Recursive function of the Sieve of Eratosthenes
#
def rec_eratho(x):
	if(x[0]*x[0] > x[-1]):
		return x
	else:
		l=[]
		for i,elt in enumerate(x):
			if(elt%x[0] != 0):
				l = l + [elt]
		return [x[0]] + rec_eratho(l)

# print("Crible d'Erathostene : {}".format(rec_eratho(int_list(2,100))))



		
