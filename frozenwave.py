from turtle import color
import numpy
import scipy 
import numpy as np

import math

from scipy import fabs
from numpy import pi

import scipy.special as s
import scipy.integrate as integrate

import cmath


import matplotlib.pyplot as plt

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

def figura1():
  F_values = []

  Z = np.linspace(0, L, 200)
  for z in Z: 
    F_values.append(abs(F_z(z, L))**2)

  plt.plot(Z, F_values, 'r')

  plt.xlabel('z')
  plt.grid()
  plt.show()

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

print("Aq para q = 75 -> ",Aq(75, L))

def figura2a():
    values = []
    Q = np.arange(-75, 76)

    for q in Q:
        values.append(np.log(abs(Aq(q, L))))

    plt.figure(figsize=[7,5])
    plt.plot(Q, values, 'r.-')

    plt.ylim(-14, 0)
    plt.xlabel('q')
    plt.grid()
    plt.show()

def figura2b():

    values = []
    Q = np.arange(-75, 75)
    L = 400*10**(-6)

    for q in Q:
        values.append(np.log(abs(Aq(q, L))))
    
    plt.xlabel('q')
    plt.grid()

    valores = []
    for q in Q:
        valores.append(numpy.angle(Aq(q, L)))
    
    print(Aq(-30, L))

def psi(rho, z):

  Q = 0.8*k

  soma = []
  total = 0
  N = 75
  qs = np.arange(-75, 75)

  for q in qs:

      k_zq = Q + 2*np.pi*q/L
      k_pq = math.sqrt(k**2 - k_zq**2)
      a = Aq(q, L)
      j0 = scipy.special.j0(k_pq *rho) 
      exponencial = np.exp(-1j * k_zq * z)
      soma = a * j0 * exponencial
      
      total += soma
  

  return total

def figura3():

  Psi_values = []

  Z = np.linspace(0, L, 100)

  for z in Z: 
    Psi_values.append(abs(psi(0, z))**2)

  plt.plot(Z, Psi_values, 'b')
  
  plt.xlabel('z')
  plt.grid()
  plt.show()

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

print("tau", tau_n(1, 30, 0.9999))
print("pi", pi_n(1, 30, 0.9999))

def figura4a():

  n1 = []
  n2 = []
  n3 = []

  X = np.linspace(0.9, 1, 200)

  for x in X:
    n1.append(tau_n(1, 10, x))

  for x in X:
    n2.append(tau_n(1, 20, x))
  
  for x in X:
    n3.append(tau_n(1, 30, x))

  plt.figure(figsize=[7,5])
  plt.plot(X, n1, 'b', label='10')
  plt.plot(X, n2, 'r', label='20')
  plt.plot(X, n3, 'g', label='30')

  plt.ylim(-400, 400)
  plt.xlabel('x')
  plt.ylabel('tau')
  plt.legend()
  plt.grid()
  plt.show()

def figura4b():

  n1 = []
  n2 = []
  n3 = []

  X = np.linspace(0.9, 0.999999, 200)

  for x in X:
    n1.append(pi_n(1, 10, x))

  for x in X:
    n2.append(pi_n(1, 20, x))
  
  for x in X:
    n3.append(pi_n(1, 30, x))

  plt.plot(X, n1, 'b', label ='10')
  plt.plot(X, n2, 'r', label='20')
  plt.plot(X, n3, 'g', label='30')


  plt.ylim(-400, 100)
  plt.xlabel('x')
  plt.ylabel("pi")

  plt.legend()
  plt.grid()
  plt.show()

figura4a()
figura4b()

def gn(n, z0):

  soma = []
  total = 0
  Q = 0.8*k

  primeiroTermoDem = -2/(n*(n + 1))

  qs = np.arange(-75, 75)

  for q in qs:

    k_zq = Q + 2*np.pi*q/L

    k_termo = k_zq/k

    primeiroTermoSum = Aq(q, L)/(1 + k_termo)

    primeiroTermoMul = (pi_n(1, n, k_termo)) + (tau_n(1, n, k_termo))

    exponencial = np.exp(-1j * k_zq * z0)

    soma = primeiroTermoSum * primeiroTermoMul * exponencial
    total += soma

  return primeiroTermoDem*(total)

##print('z0 = 0, n = 10', gn(10, 0))
#print('z0 = 0, n = 50', gn(50, 0))

def figura5a():
    n1 = []
    n2 = []
    n3 = []

    Z0 = np.linspace(0, 400*10**(-6), 100)

    for z0 in Z0:
      n1.append(abs(gn(10, z0)))

    for z0 in Z0:
      n2.append(abs(gn(20, z0)))
    
    for z0 in Z0:
      n3.append(abs(gn(30, z0)))

    plt.plot(Z0, n1, 'b', label ='n1')
    plt.plot(Z0, n2, 'r', label='n2')
    plt.plot(Z0, n3, 'g', label='n3')

    plt.ylim(0, 0.4)
    plt.xlabel('z')
    plt.ylabel('|$g_n$|')
    plt.legend()
    plt.grid()
    plt.show()

#figura5a()

def figura5b():
    n1 = []
    n2 = []
    n3 = []

    Z0 = np.linspace(0, 400*10**(-6), 200)

    for z0 in Z0:
      n1.append(cmath.phase(gn(1, z0)))


    for z0 in Z0:
      n2.append(cmath.phase(gn(2, z0)))
    
    for z0 in Z0:
      n3.append(cmath.phase(gn(3, z0)))

    plt.plot(Z0, n1, 'b', label ='n1')
    plt.plot(Z0, n2, 'r', label='n2')
    plt.plot(Z0, n3, 'g', label='n3')

    plt.xlabel('z')
    plt.ylabel('|$g_n$|')
    plt.legend()
    plt.grid()
    plt.show()


mi = 1
n = np.linspace(1, 2, 100)
z = np.linspace(0, L, 50)


m = (1.57 - 0.038j)
m2 = (1.57 - 0.19j)
m3 = (1.57 - 0.95j)


epslon1 = (m**2).real
epslon2 = -(m**2).imag


epslon22 = -(m2**2).imag
epslon23 = -(m3**2).imag

#epsilon24 = -(m4**2).imag

#spherical bessel function 1 type 
# j_n = spherical_jn(n, z, derivative = bool) = j_n
def js_n(n):
   return s.spherical_jn(n, z, False)

#spherical bessel function 2 type
#y_n = spherical_jn(n, z, derivative = bool)
def ys_n(n):
    return s.spherical_yn(n, z, False)

#Ricatti Bessel 1type given by = x*(j_n(x))
def psiBessel(n, x):
    a = x*(s.spherical_jn(n, x, False))
    return a

#Ricatti Bessel 1type derivative
derivativePsiBesseln = (lambda n, x: ((1 + n)*s.spherical_jn(n, x, False) - (x*s.spherical_jn((n + 1), x, False))))

#Spherical Hankel H2 function 
def sphericalHankel_n(n, x):
    js_n = s.spherical_jn(n, x, False)
    ys_n = s.spherical_yn(n, x, False)
    return js_n - (ys_n*1j)

#Ricatti Bessel 2 type 
def RiccatiBessel_2Type(n, x):
    return x*(sphericalHankel_n(n, x))

#spherical hankel H2 derivative
def derivativeSphericalHankel_n(n, x):
    return (1+n)*sphericalHankel_n(n, x) - x*sphericalHankel_n((n+1), x)

#coefficient C_n
def cs_n(M, mi, n, x):
    m = M
    dem = (m*mi*((RiccatiBessel_2Type(n, x)*derivativePsiBesseln(n, x)) - derivativeSphericalHankel_n(n, x)*psiBessel(n, x)))
    num = (mi*RiccatiBessel_2Type(n, x)*derivativePsiBesseln(n, x*m) - m*derivativeSphericalHankel_n(n, x)*psiBessel(n, x*m))
    return dem/num

#list of cs values 
cs = []

#coefficient D_n
def ds_n(M, mi, n, x):
    m = M
    dem = (m*m)*(RiccatiBessel_2Type(n, x)*derivativePsiBesseln(n, x) - derivativeSphericalHankel_n(n, x)*psiBessel(n, x))
    num = (m*RiccatiBessel_2Type(n, x)*derivativePsiBesseln(n, m*x)) - (mi*derivativeSphericalHankel_n(n, x)*psiBessel(n, m*x))
    return dem/num

#list of ds values 
ds = []

#coefficient r_n
def r_n(M, n, x):
    m = M
    psiConj = np.conjugate(psiBessel(n, m*x))
    dem = (m*psiBessel(n+1, m*x)*psiConj)
    num = m*m
    return (dem.imag)/(num.imag)

#coefficient s_n
def s_n(m, n, x):
  
    a = ((1j)/(2*((m**2).imag)))
    b = (m*((abs(psiBessel(n, m*x))**2)))
    c = (np.conjugate(m)) * (abs(psiBessel(n+1, m*x)))**2
    d = (m +((((2*(n + 1)*((m*m).real)/m)))))*r_n(m, n, x)
    e = (2*n + 1)*(np.conjugate(m)*r_n(m, n+1, x))
    
    return -a*(x*(b+c) - d + e)
  

ceillingX = (lambda x: math.ceil(x + 4.05*x**(1/3)+2))

def j_any(epslon, m, x, z0):
    
  a = (3*epslon)/((abs(m)**2)*(x**3))
  soma = []
  total = 0
  nMax = ceillingX(x)


  for i in range(1, (nMax+1)):

    b = (i*(i+2)/m)
    c1 = (gn(i, z0))*(np.conjugate(gn(i, z0)))*(cs_n(m, mi, i+1, x))*(np.conjugate(cs_n(m, mi, i, x)))*(r_n(m, i+1, x))
    c2 = (gn(i+1, z0))*(np.conjugate(gn(i, z0)))*(abs(1/m)**2)*(ds_n(m, mi, i+1, x)*(np.conjugate(ds_n(m, mi, i, x)))*r_n(m, i, x))
        
    d = ((i*(i+2)/(i+1)))
    d1 = (gn(i, z0))*(np.conjugate(gn(i+1, z0)))*(cs_n(m, mi, i, x))*(np.conjugate(cs_n(m, mi, i+1, x)))
    d2 = (gn(i+1, z0))*(np.conjugate(gn(i, z0)))*((abs(1/m)**2))*(ds_n(m, mi, i+1, x))*(np.conjugate(ds_n(m,  mi, i, x)))
       
    f = (((2*i+1)/(i*(i+1)))*((gn(i, z0))**2))*(cs_n(m, mi, i, x))*np.conjugate(ds_n(m, mi, i, x))

    soma = b*(c1+c2) - (((d*(d1 + d2))+ f*(1/m))*s_n(m, i, x))
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
        jn1_gauss.append(j_any(epslon2, m, 0.1,z0))
      
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