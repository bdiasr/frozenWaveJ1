import scipy 
import numpy as np

import math

from numpy import pi

import scipy.special as s
import scipy.integrate as integrate

import cmath
import matplotlib.pyplot as plt

import specialFunc
import cte

from cte import epslon2, m

lamb = 1064*10**(-9)
k = 2*pi/lamb
L = 400*10**(-6)

def F_z(z, L):

  if (z >= L/12 and z <= 11*L/12):
    expoente = -5*((z - 0.5*L)**2*1/L**2)
    value = np.exp(expoente)*np.cos(6*pi*z/L)
  else:
    value = 0
  
  return(value)

'''
print("Fz z = 0", F_z(0, L))
print("Fz z = L/10", F_z(L/10, L))
print("Fz z = 10*L/11", F_z(10*L/11, L))
print("Fz z = 10*L/12", F_z(10*L/12, L))
'''


def complex_integrate(L, q, a, b):

    real_func = lambda z: F_z(z, L)*(np.cos((((2*pi*q)/L)*z)))
    imag_func = lambda z: F_z(z, L)*(np.sin((((2*pi*q)/L)*z)))

    real_integral = integrate.quad(real_func, a, b)
    imag_integral = integrate.quad(imag_func, a, b)

    if imag_integral[0] == 0.0:
        return real_integral[0]

    return real_integral[0] + 1j*imag_integral[0]

def Aq(q, L):

    frac = 1/L
    integ = complex_integrate(L, q, L/12, (11*L)/12)

    return frac*integ

#print("Aq para q = 75 -> ",Aq(75, L))



def psi(rho, z):

  Q = 0.8*k

  soma = []
  total = 0
  N = 75
  qs = np.arange(-75, 76)

  for q in qs:

      k_zq = Q + 2*np.pi*q/L
      k_pq = math.sqrt(k**2 - k_zq**2)
      a = Aq(q, L)
      j0 = scipy.special.j0(k_pq *rho) 
      exponencial = np.exp(-1j * k_zq * z)
      soma = a * j0 * exponencial
      
      total += soma
  

  return total

def pi_n(m, n, x):
    
    fator = cmath.sqrt(1 - x**2)
    
    
    frac =  1/fator

    aux = s.lpmn(m, n, x)[0]

    pi_val = aux[m][n]*frac
    
    return pi_val

def tau_n(m, n, x):

    if x == 1:
      x = 0.99

    pmn = s.lpmn(m, n, x)[0][m][n]
    pmn1 = s.lpmn(m, n+1, x)[0][m][n+1]

    fator = cmath.sqrt(1 - x**2)

    termo1 = -(n + 1)*x*pmn
    termo2 = (n - m + 1)*pmn1

    
    return (termo1 + termo2)/fator

#print("tau", tau_n(1, 30, 0.9999))
#print("pi", pi_n(1, 30, 0.9999))

def kzq(q, L, Q):
  return Q + ((2*np.pi * q)/L)

def gn(n, z0):

  soma = []
  total = 0
  Q = 0.8*k

  primeiroTermoDem = (-2)/(n*(n + 1))

  qs = np.arange(-75, (75 + 1))

  for q in qs:

    k_zq = kzq(q, L, Q)

    k_termo = k_zq/k

    primeiroTermoSum = Aq(q, L)/(1 + k_termo)

    primeiroTermoMul = pi_n(1, n, k_termo) + tau_n(1, n, k_termo)

    exponencial = np.exp(1j * k_zq * z0)

    soma = primeiroTermoSum * primeiroTermoMul * exponencial
    total += soma

  return primeiroTermoDem*total

print('z0 = 0, n = 10', gn(10, 0))
print('z0 = 0, n = 50', gn(50, 0))

ceillingX = (lambda x: math.ceil(x + 4.05*x**(1/3)+2))

def j_any(epslon, m, x, z0):
    
  a = (3*epslon)/((abs(m)**2)*(x**3))
  soma = []
  total = 0
  nMax = ceillingX(x)

  for i in range(1, (nMax+1)):

    b = (i*(i+2)/m)
    c1 = (gn(i, z0))*(np.conjugate(gn(i, z0)))*(specialFunc.cs_n(m, cte.mi, i+1, x))*(np.conjugate(specialFunc.cs_n(m, cte.mi, i, x)))*(specialFunc.r_n(m, i+1, x))
    c2 = (gn(i+1, z0))*(np.conjugate(gn(i, z0)))*(abs(1/m)**2)*(specialFunc.ds_n(m, cte.mi, i+1, x)*(np.conjugate(specialFunc.ds_n(m, cte.mi, i, x)))*specialFunc.r_n(m, i, x))
        
    d = ((i*(i+2)/(i+1)))
    d1 = (gn(i, z0))*(np.conjugate(gn(i+1, z0)))*(specialFunc.cs_n(m, cte.mi, i, x))*(np.conjugate(specialFunc.cs_n(m, cte.mi, i+1, x)))
    d2 = (gn(i+1, z0))*(np.conjugate(gn(i, z0)))*((abs(1/m)**2))*(specialFunc.ds_n(m, cte.mi, i+1, x))*(np.conjugate(specialFunc.ds_n(m,  cte.mi, i, x)))
       
    f = (((2*i+1)/(i*(i+1)))*((gn(i, z0))**2))*(specialFunc.cs_n(m, cte.mi, i, x))*np.conjugate(specialFunc.ds_n(m, cte.mi, i, x))

    soma = b*(c1+c2) - (((d*(d1 + d2))+ f*(1/m))*specialFunc.s_n(m, i, x))
    total += soma

  return total.imag * a 


#print("z0 = 0", j_any(epslon2, m, 0.1, 0))
#print("z0 = 10", j_any(epslon2, m, 0.1, 10*10**(-6)))
#print("z0 = 200", j_any(epslon2, m, 0.1, 200*10**(-6)))
#print("z0 = 400", j_any(epslon2, m, 0.1, 400*10**(-6)))

def graph():

    Z0 = np.linspace(0, L, 50)

    jn1_gauss = []
    jn2_gauss = []
    jn3_gauss = []
   
    #gauss bean M = 1.57 - 0.038j
    for z0 in Z0:
        jn1_gauss.append((j_any(epslon2, m, 0.1,z0))*2500)
      
    #print(jn1_gauss)
   
    #gauss bean M = 1.57 - 0.19j
    for z0 in Z0:
        jn2_gauss.append(j_any(epslon2, m, 3, z0))

    #gauss bean M = 1.57 - 0.95j
    for z0 in Z0:
        jn3_gauss.append(j_any(epslon2, m, 8, z0))


    plt.plot(Z0, jn1_gauss,"b", label= "x = 0.1" )
    plt.plot(Z0, jn2_gauss,"r", label= "x = 3" )    
    plt.plot(Z0, jn3_gauss,"g", label= "x = 8" )

    plt.xlabel('Z0')
    plt.ylabel('J1')
    plt.title("Assymetry Factor J1 x Z0")
    plt.legend()
    plt.grid()
    plt.show()

graph()