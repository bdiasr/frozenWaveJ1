import numpy as np
from frozenwave import tau_n, pi_n
import matplotlib.pyplot as plt

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