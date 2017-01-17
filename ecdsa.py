############################
#         ecdsa.py         #
############################

import os, sys
import math, random
import utils


class ECDSA_key:
    """ ECDSA_key class
    """

    def __init__(p, a, b, h, n, G):
        """ Constructor """
        if(is_prime(p) ==  False):
            print("Error : You have chosen a non prime number for p")
        if(is_prime(n) ==  False):
            print("Error : You have chosen a non prime number for n")
        if(n<pow(2,160)):
            print("Warning : n should be greater than 2^160")
        if(s == 0 or s>(n-1)):
            print("Error : s should be between 1 and n-1")
        # Initialize variables #
        self.G = G # (coordinates (g_x, g_y)) #
        self._s = 0
        self.Q = 0 # (coordinates (q_x, q_y)) #
        self.n = n
        self.p = p
        self.a = a
        self.b = b
        self.h = h
        # Generate key #
        generate_key()
        

    """ Methods """
    def get_pkey():
        return self.Q

    def get_gen():
        return self.G

    def get_n():
        return n

    def _get_skey():
        return self._s

    def _set_skey(s):
        self._s = s

    def set_pkey(Q):
        self.Q = Q

    s = property(_get_skey)
    s = property(_set_skey)

    #
    # Generate key
    #
    def generate_key():
        random.seed()
        _set_skey(random.randint(1, get_n()-1))
        set_pkey(_get_skey()*get_gen())

    #
    # Sign m with ECDSA
    #
    def sign(m):
        r , t = 0 # t is s in the paper #
        # calculated truncated hash of m #
        z = hashed(m, int(sys.getsizeof(get_n())))
        # r != 0 #
        while(r == 0 or t == 0):
            # Take random int between 1 and n-1 #
            k = random.randint(1, get_n()-1)
            # P = k*G #
            P = [i * k for i in get_gen()]
            # r = x_P mod n  #
            r = P[1] % get_n()
            # t = k^-1(z+r*d_A) mod n #
            t = mod_inv(k, get_n(), 0)*(z+r*_get_skey()) % get_n()
        return (r, t)
            
    #
    # Verify signature Sgn
    #
    def verify(Sgn, Q, z):
        # u1 = t^-1*z mod n #
        # u2 = t^-1*r mod n #
        t_inv = mod_inv(Sgn[1], get_n(), 0)
        u1 = t_inv * z % get_n()
        u2 = t_inv * Sgn[0] % get_n()
        # P = u1 * G + u2 * Q
        P_x = u1 * get_gen()[0] + u2 * get_pkey()[0]
        P_y = u1 * get_gen()[1] + u2 * get_pkey()[1]
        P = [P_x, P_y] # useless #
        if(Sgn[0] != P[0] % get_n() ):
            print("WRONG")
            return False
        else:
            print("OK")
            return True

    #
    # Hash function
    #
    def hashed(m, size):
        return hash(m)

"""
def verify(c, sig):
        if(self.Q == O): # /!\ verify that it's not the infinite point (not 0...)
            print("Doesn't match")
            return False
        if(self.n*self.Q != O): # /!\ Same
            print("Dosen't match")
            return False
        if(sig[1] == 0 or y == 0 or sig[1]>n-1 or sig[2]>n-1):
            print("Dosen't match")
            return False
        (i,j)=(hashed(c)/sig[2] % self.n)*self.G + (sig[1]/sig[2] % self.n)*self.Q
        if(sig[1] != i % self.n):
            print("Doesn't match")
            return False
        print("OK")
        return True
"""






            
    


