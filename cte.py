import numpy as np


mi = 1
n = np.linspace(1, 2, 100)
z = np.linspace(0, 400*10**(-6), 50)


m = (1.57 - 0.038j)
m2 = (1.57 - 0.19j)
m3 = (1.57 - 0.95j)


epslon1 = (m**2).real
epslon2 = -(m**2).imag


epslon22 = -(m2**2).imag
epslon23 = -(m3**2).imag