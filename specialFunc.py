import numpy as np


import scipy.special as s


import cte

mi = cte.mi
z = cte.z

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
  


