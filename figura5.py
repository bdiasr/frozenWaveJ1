import numpy as np
import cmath

from frozenwave import gn

import matplotlib.pyplot as plt


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
