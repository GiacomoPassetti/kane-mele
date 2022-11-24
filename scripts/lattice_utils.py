import numpy as np
import matplotlib.pyplot as plt
a = 1
v1 = a*np.array([np.sqrt(3)/2, 3/2])
v2 = a*np.array([np.sqrt(3)/2, -3/2])
n1 = 10
n2 = 10
rA = a*[0, 0.5]
rB = a*[0, -0.5]

def atom_positions_generator(a1, a2, n1, n2, rA, rB):
    positions_array = np.zeros(shape = (n1,n2, 2))
    for i in range(n1):
        for j in range(n2):
            positions_array[i, j, :] = ((i*a1)+(j*a2))
    positions_array = positions_array.reshape(n1*n2, 2)
    final_positions = np.zeros(shape =(2*n1*n2, 2))
    for i in range(n1*n2):
        final_positions[2*i, :] = positions_array[i, :] + rA
        final_positions[(2*i) + 1, :] = positions_array[i, :] + rB
    return final_positions

def unitary_cells_generator(a1, a2, n1, n2, rA, rB):
    positions_array = np.zeros(shape = (n1,n2, 2))
    for i in range(n1):
        for j in range(n2):
            positions_array[i, j, :] = ((i*a1)+(j*a2))
    return positions_array

def Hamiltonian_generator(n1, n2):
    number_of_states = n1*n2*4
    Labels = np.ones(shape=(n1, n2))
    indeces_labels = np.column_stack(np.where(Labels==1))
    H = np.zeros(shape=(number_of_states, number_of_states))
Hamiltonian_generator(4, 4,)




#pos = atom_positions_generator(v1, v2, n1, n2, rA, rB)
#fig, ax = plt.subplots()
#ax.plot(pos[:, 0], pos[:, 1], ls = "", marker = "o")
#ax.set_xlim(-2, 14)
#ax.set_ylim(-8, 8)
#ax.set_aspect('equal')
#plt.show()

 

      
                 