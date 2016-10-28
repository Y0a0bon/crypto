############################
#          RSA.py          #
############################

import os, sys
import math
import utils


p = 70457971
q = 9891923
# n = p * q = 696964823868233
n = 696964823868233
e = 547
# d = e^-1 mod phi(n)
#   = e^-1 mod phi(p, q)
#   = 44595555457748 mod 696964743518340
d = 269024829525699
#
# (e,n) = (547,696964823868233)
# (d,n = (2.32972487969e+14,696964823868233)
c = 0


#
# Apply RSA
#
def e_rsa(m, e, p, q, verbose):
    n = p*q
    if(utils.calcpgcd(n, e, 0)!=1):
        print(str(e)+" and "+str(n)+"should not have common factors.")
        return 0
    #c = pow(m, e)%n
    c = utils.fast_exp(m, e, n, verbose)
    if(verbose!=0):
        print('Message is now '+str(c))
    return c
# m = "HELLO" in ASCII
# m = 7269767679
# Blocks have to be sliced in size not equal or greater than the number that compose n (1073 : 4 numbers)
# Answer should be 436
#print("Should be 436 : ")
#e_rsa(726, 71, 29, 37, 1)

print("Big numbers woohoo :")
m = 1254894561
print("m is "+str(m)+" and c is :")
c = e_rsa(m, e, p, q, 1)
print(c)


#
# Decode RSA
#
def d_rsa(c, d, n, verbose):
    #m = pow(c, d)%n
    m = utils.fast_exp(c, d, n, verbose)
    if(verbose!=0):
        print('Message was '+str(m))
    return m

# Answer should be 726
# print("Should be 726 : ")
# d_rsa(436 ,1079, 1073, 1)
print("c is "+str(c)+" and m should be "+str(m))
print(d_rsa(c, d, n, 1))


#
# Decode RSA with CRT
#
def rsa_crt(c, d, p, q, verbose):
    #d_p = mod_inv(e, p-1, verbose)
    #d_q = mod_inv(e, q-1, verbose)
    d_p = d % (p-1) 
    d_q = d % (q-1)
    q_inv = utils.mod_inv(q, p, verbose)
    
    #m_p = pow(c, d_p) % p
    m_p = utils.fast_exp(c, d_p, p, 0)
    #m_q = pow(c, d_q) % q
    m_q = utils.fast_exp(c, d_q, q, 0)
    h = q_inv * (m_p - m_q) % p
    m = m_q + h * q

    if(verbose != 0):
        print('Message was '+str(m))
    return m

#print(rsa_crt(436, 1079, 29, 37, 0))
# Should return 513
#rsa_crt(8363, 11787, 137, 131, 1)
# Error :
print("c is "+str(c)+" and m should be "+str(m))
print(rsa_crt(c, d, p, q, 1))


#
# Bellcore attack on RSA
#
def bellcore(c, d, p, p_fault, q, verbose):
    d_p = d % (p-1)
    d_q = d % (q-1)
    d_p_fault = d % (p_fault-1)
    
    q_inv = mod_inv(q, p, verbose)
    q_inv_fault = mod_inv(q, p_fault, verbose)
    
    #m_p = pow(c, d_p) % p
    m_p = utils.fast_exp(c, d_p, p, 0)
    #m_q = pow(c, d_q) % q
    m_q = utils.fast_exp(c, d_q, q, 0)
    #m_p_fault = pow(c, d_p_fault) % p_fault
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
# Broadast attack on RSA (m,e,n) e_rsa(m, e, n, verbose)
#
def broadcast_attack(m, e, n, verbose):
    c=[]
    M = 1
    random.seed
    for i in range(1,e):
        m_tmp = gen_prime(n)
        c[i] = e_rsa(m, e, gen_prime(n))
        M = M * m_tmp
    m_r = 0
    if(verbose != 0):
        a=0

