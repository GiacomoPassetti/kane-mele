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

n1 = 16
n2 = n1

t = 1
t_so = 0.2
b_0 = 0.1



H = Hamiltonian_generator(n1, n2, t, t_so, b_0)
u, v = eigh(H)

ground_state = v[:, 5]


atom_positions = atom_positions_generator(v1, v2, n1, n2, rA, rB)

Labels = np.ones(shape=(n1, n2))
indeces_labels = np.column_stack(np.where(Labels==1))
#print(atom_positions)

xs = atom_positions[:, 0]
ys = atom_positions[:, 1]
fig, ax = plt.subplots()
cmap = matplotlib.cm.get_cmap('viridis')

print(ground_state.shape[0], xs.shape[0])
for i in range(xs.shape[0]):
    print(np.abs(ground_state[i]))
    ax.plot(xs[i], ys[i], color = cmap(np.abs(ground_state[2*i])), marker = "o", ls = "")

plt.show()