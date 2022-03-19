import numpy as np
from frozenwave import L, Aq
import matplotlib.pyplot as plt

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
        valores.append(np.angle(Aq(q, L)))
    
    print(Aq(-30, L))