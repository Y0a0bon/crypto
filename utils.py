############################
#         utils.py         #
############################

# Contains the following functions :
# euclide(a, b, verbose)
# calcpgcd(a, b, verbose)
# mod_inv(b, n, verbose)
# ind_euler(N, verbose)
# fast_exp(m, e, n verbose)
# phi(p, q)
# primes(n)
# iroot(k, n)
# eratho(a, b, ch)
# int_list(a, b)
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

# calcpgcd(45,89, 1)
# print(calcpgcd(696964823868233, 547, 0))


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

#mod_inv(5, 77, 1)
# Compute 4.45955554577e+13
# e_1 = mod_inv(547, 696964823868233, 0)
# print(e_1)


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
    c , e_sv = 1 , e
    while(e>0):
        if(e % 2 == 0):
            e//=2
        else:
            c = (c * m) % n
            e = (e - 1)//2
        m = (m * m) % n
    if(verbose != 0):
        print(str(m) + " ^ " + str(e_sv) + " % " + str(n) + " = " + str(c))
    return c

# fast_exp(11, 13, 19, 1)
# print(pow(11,13)%19)
# print(fast_exp(1254894561, 44595555457748, 696964823868233, 0))


#
# Compute phi(n) = phi(pq)
#
def phi(p, q):
    return(p-1)*(q-1)

# print(phi(70457971, 9891923))


#
# List all prime numbers under n if 0, else the greatest
#
def primes(n, ch): 
	if n==2: return [2]
	elif n<2: return []

	s=list(range(3,n+1,2))
	mroot = n ** 0.5
	half=(n+1)/2-1
	i=0
	m=3
	while m <= mroot:
		if s[i]:
			j=int((m*m-3)/2)
			s[j]=0
			while j<half:
				s[j]=0
				j+=m
		i=i+1
		m=2*i+3
	pr = [2]+[x for x in s if x]
        if(ch==0):
            return pr
        else :
            return pr[-1]

#print(primes(70457984, 1))
#print(primes(9891945, 1))

# to try
#for i in range(0, int(mroot/2)) :
#    m=2*i+3
#    if s[i]:
#        for j in range((m*m-3)/2, half, m) :
#            s[j]=0


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


#
# List of integers from a to b
#
def int_list(a, b):
    return [i for i in range(a, b+1)]

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
        
