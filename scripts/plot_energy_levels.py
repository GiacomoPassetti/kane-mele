import numpy as np
import matplotlib.pyplot as plt
from lattice_utils import Hamiltonian_generator, I_J_conv, atom_positions_generator
from numpy.linalg import eigh
import matplotlib


a = 1
v1 = a*np.array([np.sqrt(3)/2, 3/2])
v2 = a*np.array([np.sqrt(3)/2, -3/2])

rA = a*[0, 0.5]
rB = a*[0, -0.5]

n1 = 18
n2 = n1

t = 1
t_so = 0.1
b_0 = 0.2

H = Hamiltonian_generator(n1, n2, t, t_so, b_0)
u, v = eigh(H)

plt.plot(u, ls = "", marker = "o")
plt.show()