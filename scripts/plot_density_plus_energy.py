import numpy as np
import matplotlib.pyplot as plt
from lattice_utils import Hamiltonian_generator, I_J_conv, atom_positions_generator
from numpy.linalg import eigh
import matplotlib
from matplotlib import gridspec


a = 1
v1 = a*np.array([np.sqrt(3)/2, 3/2])
v2 = a*np.array([np.sqrt(3)/2, -3/2])

rA = a*[0, 0.5]
rB = a*[0, -0.5]

n1 = 16
n2 = n1

t = 1
t_so = 0.1
b_0 = 0.2



H = Hamiltonian_generator(n1, n2, t, t_so, b_0)
u, v = eigh(H)

index_of_state = int(u.shape[0]/2) - 10

selected_state = v[:,index_of_state]


atom_positions = atom_positions_generator(v1, v2, n1, n2, rA, rB)

Labels = np.ones(shape=(n1, n2))
indeces_labels = np.column_stack(np.where(Labels==1))
#print(atom_positions)

xs = atom_positions[:, 0]
ys = atom_positions[:, 1]
nrow = 1
ncol = 3
fig = plt.figure(figsize=(5.2, 7), dpi = 800) 

gs = gridspec.GridSpec(nrow, ncol, width_ratios = [0.2, 0.1, 1],
         wspace=0.0, hspace=0.0, top=0.97, bottom=0.08, left=0.18, right=0.98) 
ax = plt.subplot(gs[0, 2])
ax_energy = plt.subplot(gs[0, 0])



cmap = matplotlib.cm.get_cmap('OrRd')


print(selected_state.shape[0], xs.shape[0])
for i in range(xs.shape[0]):
    #print("local site occupation: ", (np.abs(selected_state[2*i])**2)+(np.abs(selected_state[2*i + 1])**2))
    ax.plot(xs[i], ys[i], color = cmap(((np.abs(selected_state[2*i])**2)+(np.abs(selected_state[2*i + 1])**2))*10), marker = "o", ls = "")
ax.set_aspect("equal")
ax_energy.plot(np.arange((index_of_state - 10), (index_of_state + 10)), u[(index_of_state - 10):(index_of_state + 10)], marker = "o", ls = "")
ax_energy.plot(index_of_state, u[index_of_state], color = "red", marker = "o", ls = "")
plt.savefig("plot_energy_plu_density.png")