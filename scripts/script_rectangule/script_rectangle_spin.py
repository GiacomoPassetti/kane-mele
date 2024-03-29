import numpy as np
import matplotlib.pyplot as plt
from lattice_utils import Hamiltonian_generator, I_J_conv, atom_positions_generator, Hamiltonian_boundaries_generator, Circle_label_generator,rectangular_label_generator , atom_position_IJ
from numpy.linalg import eigh
import matplotlib
from matplotlib import gridspec
import matplotlib as mpl




a = 1
v1 = a*np.array([np.sqrt(3)/2, 3/2])
v2 = a*np.array([np.sqrt(3)/2, -3/2])

rA = a*[0, 0.5]
rB = a*[0, -0.5]

n1 = 36
n2 = n1

t = 1
t_so = 0.1
b_0 = 0.2
l1 = 48
l2 = l1



atom_positions = atom_positions_generator(v1, v2, n1, n2, rA, rB)
Atom_Labels = np.ones(shape=(n1, n2))
indeces_labels = np.column_stack(np.where(Atom_Labels==1))


Labels = rectangular_label_generator(v1, v2, n1, n2, l1, l2)
print(Labels)
H, indeces_labels = Hamiltonian_boundaries_generator(n1, n2, t, t_so, b_0, Labels)

u, v = eigh(H)
index_of_state = int(u.shape[0]/2) - 4
selected_state = v[:,index_of_state]

fig = plt.figure(figsize=(7, 7), dpi = 800) 
nrow = 2
ncol = 3
gs = gridspec.GridSpec(nrow, ncol, width_ratios = [0.3, 0.1, 1], height_ratios=[1, 0.05],

         wspace=0.0, hspace=0.0, top=0.97, bottom=0.08, left=0.18, right=0.98) 
ax = plt.subplot(gs[0, 2])
ax_energy = plt.subplot(gs[0, 0])
ax_colorbar = plt.subplot(gs[1, 2])

index_of_state = int(u.shape[0]/2) 

cmap = matplotlib.cm.get_cmap('OrRd')
norm = mpl.colors.Normalize(vmin=0, vmax=1)

cb1 = mpl.colorbar.ColorbarBase(ax_colorbar, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
cb1.set_label('P_occupation*10')


for i in range(indeces_labels.shape[0]):
    
    coordA, coordB = atom_position_IJ(v1, v2, indeces_labels[i, 0], indeces_labels[i, 1], rA, rB)
    if Labels[tuple(indeces_labels[i])] == 1:
        ax.plot(coordA[0], coordA[1], color = cmap(((np.abs(selected_state[4*i])**2)+(np.abs(selected_state[4*i + 1])**2))*10), marker = "o", ls = "")
        ax.plot(coordB[0], coordB[1], color = cmap(((np.abs(selected_state[4*i+2])**2)+(np.abs(selected_state[4*i + 3])**2))*10), marker = "o", ls = "")

ax.set_aspect("equal")
ax_energy.plot(np.arange((index_of_state - 10), (index_of_state + 10)), u[(index_of_state - 10):(index_of_state + 10)], marker = "o", ls = "")
ax_energy.plot(index_of_state, u[index_of_state], color = "red", marker = "o", ls = "")
ax_energy.set_xlabel("Number of States")
ax_energy.set_ylabel("E[t]")
ax.set_xticks([])
ax.set_yticks([])
import os
plt.savefig(os.path.join("plots", "Second_attempt_energy_occupation_n1_"+str(n1)+"_n2_"+str(n2)+"l1_"+str(l1)+"l2_"+str(l2)+".png"))