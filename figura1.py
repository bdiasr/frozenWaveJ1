
import numpy as np
from frozenwave import F_z, L
import matplotlib.pyplot as plt


def figura1():
  F_values = []

  Z = np.linspace(0, L, 200)
  for z in Z: 
    F_values.append(abs(F_z(z, L))**2)

  plt.plot(Z, F_values, 'r')

  plt.xlabel('z')
  plt.grid()
  plt.show()