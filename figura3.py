import numpy as np
from frozenwave import L, Aq, psi
import matplotlib.pyplot as plt


def figura3():

  Psi_values = []

  Z = np.linspace(0, L, 100)

  for z in Z: 
    Psi_values.append(abs(psi(0, z))**2)

  plt.plot(Z, Psi_values, 'b')
  
  plt.xlabel('z')
  plt.grid()
  plt.show()