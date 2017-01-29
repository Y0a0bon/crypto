############################
#          RSA.py          #
############################

import os, sys
import math
import utils


#
# Test function for rsa.py
#
def test_rsa():
    if(e_rsa(1254894561, 547, 696964823868233, 0) == 467534266279873):
        # p = 70457971 and q = 9891923 and n = 696964823868233
        print("Test de e_rsa(m, e, n, verbose) ...\t\tOK")
    else:
        print("Test de e_rsa(m, e, n, verbose) ...\t\tFAILED")
    if(d_rsa(436 ,1079, 1073, 0) == 726):
        print("Test de d_rsa(c, d, n, verbose) ...\t\tOK")
    else:
        print("Test de d_rsa(c, d, n, verbose) ...\t\tFAILED")
    if(rsa_crt(8363, 11787, 137, 131, 0) == 513):
        print("Test de rsa_crt(c, d, p, q, verbose) ...\tOK")
    else:
        print("Test de rsa_crt(c, d, p, q, verbose) ...\tFAILED")


#
# RSA key
#
class RSA_key:
    """ RSA_key class
         p & q prime (numbers distincts)
         n = p*q
         phi = phi(p, q)
         e prime with phi and < phi
         d = mod_inv(e, phi, 0)
    """

    def __init__(self, p, q, e, d):
        """ Constructor """
        self.p = p
        self.q= q
        self.e =e
        self._d= d
        self.n = p *q

    """ Methods """
    def get_n(self):
        return self.n

    def get_e(self):
        return self.e

    def _get_d(self):
        return self._d

    d = property(_get_d)

    def encrypt_message(self, m, verbose):
        """ Apply RSA to m """
        return e_rsa(m, self.e, self.n, verbose)

    def decrypt_message(self, c, verbose):
        """ Apply RSA to c """
        return d_rsa(c, self.d, self.n, verbose)

    def __del__(self):
        print("RSA key deleted")

    


#
# Apply RSA
#
def e_rsa(m, e, n, verbose):
    if(utils.calcpgcd(n, e, 0)!=1):
        print(str(e)+" and "+str(n)+"should not have common factors.")
        return 0
    c = utils.fast_exp(m, e, n, verbose)
    if(verbose!=0):
        print('Message is now '+str(c))
    return c


#
# Decode RSA
#
def d_rsa(c, d, n, verbose):
    m = utils.fast_exp(c, d, n, verbose)
    if(verbose!=0):
        print('Message was '+str(m))
    return m


#
# Decode RSA with CRT
#
def rsa_crt(c, d, p, q, verbose):
    #d_p = mod_inv(e, p-1, verbose)
    #d_q = mod_inv(e, q-1, verbose)
    d_p = d % (p-1) 
    d_q = d % (q-1)
    q_inv = utils.mod_inv(q, p, verbose)
    
    m_p = utils.fast_exp(c, d_p, p, 0)
    m_q = utils.fast_exp(c, d_q, q, 0)
    h = q_inv * (m_p - m_q) % p
    m = m_q + h * q

    if(verbose != 0):
        print('Message was '+str(m))
    return m


#
# Bellcore attack on RSA
#
def bellcore(c, d, p, p_fault, q, verbose):
    d_p = d % (p-1)
    d_q = d % (q-1)
    d_p_fault = d % (p_fault-1)
    
    q_inv = mod_inv(q, p, verbose)
    q_inv_fault = mod_inv(q, p_fault, verbose)
    
    m_p = utils.fast_exp(c, d_p, p, 0)
    m_q = utils.fast_exp(c, d_q, q, 0)
    m_p_fault = utils.fast_exp(c, d_p_fault, p, 0)
    
    h = q_inv * (m_p - m_q) % p
    h_fault = q_inv_fault * (m_p_fault - m_q) % p_fault
    
    m = m_q + h * q
    m_fault = m_q + h_fault * q

    q_r = calcpgcd(m, m-m_fault)

    if(verbose != 0):
        print('q is '+str(q_r))
    return q_r


#
# Broadcast attack on RSA
# C[] is m message encrypted with the same exponent e and a different modulo n
#
def broadcast_attack(C, N, e, verbose):
    V = []
    M = []
    for i in range(1,e):
        V.append(C[i])
        M.append(N[i])
    c = utils.crt_list(V, M)
    return utils.iroot(c[0], e)


#
# Broadast attack on RSA (m,e,n) e_rsa(m, e, n, verbose)
# NOT WORKING YET
#
def broadcast_attack_old(m, e, n, verbose):
    c=[]
    M = 1
    random.seed()
    for i in range(1,e):
        m_tmp = gen_prime(n)
        c[i] = e_rsa(m, e, gen_prime(n), verbose)
        M = M * m_tmp
    m_r = 0
    if(verbose != 0):
        a=0


#
# Wiener attack
# Source : ndiab
#
def Wiener(N : int, e : int) -> int :
        """
        Wiener Attack
        return the private key if it's vulnerable
        if not, return 0
        """
        
        #1 - initiate the list of partial quotients
        quot = utils.cont_frac(e,N)
        
        #2 - initiate a list of continued fractions
        seq = utils.frac_seq(quot)

        #check if d is actually the key
        for (k,d) in seq :
                # Now we are going to check if the couple (k,d)
                # produce a factorization of N (p * q)
                
                if k == 0:
                        continue

                # compute phi
                # phi(N) = (ed - 1)/ k
                phi = (e * d - 1) // k
                
                # Now from phi and N we can deduce the eventual factorization
                # by solving the equation : x^2 - ((N-phi)+1)x + N
                a = 1
                b = (N - phi) + 1
                c = N
                
                discr = b*b - 4 * a * c
        
                # in that case we got solutions to the equation
                if discr > 0 :
                        # to works, the solutions must be integers
                        # first, the root of the discriminent must be integer 
                        rdiscr = utils.isqrt(discr)

                        if (rdiscr * rdiscr) == discr :
                                #the root of the discr is integer, we compute the solutions 
                                p = -((-b - rdiscr) // (2 * a))
                                q = -((-b + rdiscr) // (2 * a))
                                if N == (p * q) :
                                        print("Hacked")
                                        return d

        print("Not vulnerable to Wiener")

        return 0


#
# Attack on exo2
# Source : ndiab
#
def attack_exo_2(N : int, e : int, N_bis : int, e_bis : int) -> int :
	d = Wiener (N,e)
	nblb = int(log2(d))  # number of least significant bits k known

	offset = 1 << (nblb +1)	 # Offset for each iteration
	max_iter = 1 << (276 - nblb)

	print(max_iter)
	print(offset)

	d_bis = d + offset

	m_ref = randint(10,1000)

	m = [m_ref,N_bis]
	
	# brute force the high order bits
	for i in range (max_iter):
		# to check if the key works, we compare the encryption and the decryption	
		
		if i%100000 == 0 :
			print (i)
	 
		m[0] = fast_exp(m[0], e_bis, m[1], 0)
		m[0] = fast_exp(m[0], d_bis, m[1], 0)		
		if (m[0] == m_ref) :
			print("Hacked")
			print(d_bis)
			return d_bis


		d_bis = d_bis + offset
		m[0] = m_ref

	print("doesn't work")
	return 0
    





        

