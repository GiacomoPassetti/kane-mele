import numpy as np
import matplotlib.pyplot as plt
a = 1
v1 = a*np.array([np.sqrt(3)/2, 3/2])
v2 = a*np.array([np.sqrt(3)/2, -3/2])

rA = a*[0, 0.5]
rB = a*[0, -0.5]



def I_J_conv(i, j, n):
    return int((n*i)+j)

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

def Hamiltonian_generator(n1, n2, t, t_so, b_0):
    number_of_states = n1*n2*4
    Labels = np.ones(shape=(n1, n2))
    indeces_labels = np.column_stack(np.where(Labels==1))
    H = np.zeros(shape=(number_of_states, number_of_states), dtype=complex)
    
    for i in range(indeces_labels.shape[0]):
        # Intra-cell Hopping terms
        H[4*i, 4*i+2] = -t 
        H[4*i+ 1, 4*i+3] = -t

        # Magnetic Field Term
        H[4*i, 4*i + 1] = -1j*b_0 
        H[4*i+2, 4*i + 3] = -1j*b_0 
        
        try:  # Top Right Next Cell
          if Labels[tuple(indeces_labels[i]+ np.array([1, 0]))]==1:
             new_lab = indeces_labels[i]+ np.array([1, 0])
             j = I_J_conv(new_lab[0], new_lab[1], n1)
             # Nearest Neighbour hoppings
             H[4*i, 4*j+2] = -t   # Particle A_up hops to B_up
             H[4*i + 1, 4*j+3] = -t  # Particle A_down hops to B_down
             
             # Next-nearest Neighbour hoppings
             H[4*i, 4*j+1] = 1j*t_so   # Particle A_up hops to A_down
             H[4*i + 1, 4*j] = 1j*t_so  # Particle A_down hops to A_up
             H[4*i+2, 4*j+3] = -1j*t_so   # Particle B_up hops to B_down
             H[4*i + 3, 4*j + 2] = -1j*t_so  # Particle B_down hops to B_up


        except:
            _=0
        try:  # Bottom right Next Cell
          if Labels[tuple(indeces_labels[i]+ np.array([0, 1]))]==1:
             new_lab = indeces_labels[i]+ np.array([0, 1])
             j = I_J_conv(new_lab[0], new_lab[1], n1)
             H[4*i+2, 4*j] = -t   # Particle B_up hops to A_up
             H[4*i + 3, 4*j+1] = -t  # Particle B_down hops to A_down

             # Next-nearest Neighbour hoppings
             H[4*i, 4*j+1] = 1j*t_so   # Particle A_up hops to A_down
             H[4*i + 1, 4*j] = 1j*t_so  # Particle A_down hops to A_up
             H[4*i+2, 4*j+3] = -1j*t_so   # Particle B_up hops to B_down
             H[4*i + 3, 4*j + 2] = -1j*t_so  # Particle B_down hops to B_up

        except:
            _=0
        try: # Top Left Next cell
          if (indeces_labels[i]- np.array([0, 1])[0] >= 0) and (indeces_labels[i] - np.array([0, 1])[1] >= 0):
           if Labels[tuple(indeces_labels[i]- np.array([0, 1]))]==1:
             new_lab = indeces_labels[i]- np.array([0, 1])
             j = I_J_conv(new_lab[0], new_lab[1], n1)
             H[4*i, 4*j+2] = -t   # Particle A_up hops to B_up
             H[4*i + 1, 4*j+3] = -t  # Particle A_down hops to B_down

             # Next-nearest Neighbour hoppings
             H[4*i, 4*j+1] = -1j*t_so   # Particle A_up hops to A_down
             H[4*i + 1, 4*j] = -1j*t_so  # Particle A_down hops to A_up
             H[4*i+2, 4*j+3] = 1j*t_so   # Particle B_up hops to B_down
             H[4*i + 3, 4*j + 2] = 1j*t_so  # Particle B_down hops to B_up
        except:
            _=0
        try: # Bottom Left Next cell 
          if (indeces_labels[i]- np.array([1, 0])[0] >= 0) and (indeces_labels[i] - np.array([1, 0])[1] >= 0):
           if Labels[tuple(indeces_labels[i]- np.array([1, 0]))]==1:
             new_lab = indeces_labels[i]- np.array([1, 0])
             j = I_J_conv(new_lab[0], new_lab[1], n1)
             H[4*i+2, 4*j] = -t   # Particle B_up hops to A_up
             H[4*i + 3, 4*j+1] = -t  # Particle B_down hops to A_down

             # Next-nearest Neighbour hoppings
             H[4*i, 4*j+1] = -1j*t_so   # Particle A_up hops to A_down
             H[4*i + 1, 4*j] = -1j*t_so  # Particle A_down hops to A_up
             H[4*i+2, 4*j+3] = 1j*t_so   # Particle B_up hops to B_down
             H[4*i + 3, 4*j + 2] = 1j*t_so  # Particle B_down hops to B_up

        except:
            _=0

    H = H + H.transpose().conj()
    assert np.allclose(H, H.conj().T)
    return H 

def Hamiltonian_boundaries_generator(n1, n2, t, t_so, b_0, Labels):
    number_of_states = n1*n2*4
    
    indeces_labels = np.column_stack(np.where(Labels==1))
    H = np.zeros(shape=(number_of_states, number_of_states), dtype=complex)
    
    for i in range(indeces_labels.shape[0]):
        # Intra-cell Hopping terms
        H[4*i, 4*i+2] = -t 
        H[4*i+ 1, 4*i+3] = -t

        # Magnetic Field Term
        H[4*i, 4*i + 1] = -1j*b_0 
        H[4*i+2, 4*i + 3] = -1j*b_0 
        
        try:  # Top Right Next Cell
          if Labels[tuple(indeces_labels[i]+ np.array([1, 0]))]==1:
             new_lab = indeces_labels[i]+ np.array([1, 0])
             j = I_J_conv(new_lab[0], new_lab[1], n1)
             # Nearest Neighbour hoppings
             H[4*i, 4*j+2] = -t   # Particle A_up hops to B_up
             H[4*i + 1, 4*j+3] = -t  # Particle A_down hops to B_down
             
             # Next-nearest Neighbour hoppings
             H[4*i, 4*j+1] = 1j*t_so   # Particle A_up hops to A_down
             H[4*i + 1, 4*j] = 1j*t_so  # Particle A_down hops to A_up
             H[4*i+2, 4*j+3] = -1j*t_so   # Particle B_up hops to B_down
             H[4*i + 3, 4*j + 2] = -1j*t_so  # Particle B_down hops to B_up


        except:
            _=0
        try:  # Bottom right Next Cell
          if Labels[tuple(indeces_labels[i]+ np.array([0, 1]))]==1:
             new_lab = indeces_labels[i]+ np.array([0, 1])
             j = I_J_conv(new_lab[0], new_lab[1], n1)
             H[4*i+2, 4*j] = -t   # Particle B_up hops to A_up
             H[4*i + 3, 4*j+1] = -t  # Particle B_down hops to A_down

             # Next-nearest Neighbour hoppings
             H[4*i, 4*j+1] = 1j*t_so   # Particle A_up hops to A_down
             H[4*i + 1, 4*j] = 1j*t_so  # Particle A_down hops to A_up
             H[4*i+2, 4*j+3] = -1j*t_so   # Particle B_up hops to B_down
             H[4*i + 3, 4*j + 2] = -1j*t_so  # Particle B_down hops to B_up

        except:
            _=0
        try: # Top Left Next cell
          if (indeces_labels[i]- np.array([0, 1])[0] >= 0) and (indeces_labels[i] - np.array([0, 1])[1] >= 0):
           if Labels[tuple(indeces_labels[i]- np.array([0, 1]))]==1:
             new_lab = indeces_labels[i]- np.array([0, 1])
             j = I_J_conv(new_lab[0], new_lab[1], n1)
             H[4*i, 4*j+2] = -t   # Particle A_up hops to B_up
             H[4*i + 1, 4*j+3] = -t  # Particle A_down hops to B_down

             # Next-nearest Neighbour hoppings
             H[4*i, 4*j+1] = -1j*t_so   # Particle A_up hops to A_down
             H[4*i + 1, 4*j] = -1j*t_so  # Particle A_down hops to A_up
             H[4*i+2, 4*j+3] = 1j*t_so   # Particle B_up hops to B_down
             H[4*i + 3, 4*j + 2] = 1j*t_so  # Particle B_down hops to B_up
        except:
            _=0
        try: # Bottom Left Next cell 
          if (indeces_labels[i]- np.array([1, 0])[0] >= 0) and (indeces_labels[i] - np.array([1, 0])[1] >= 0):
           if Labels[tuple(indeces_labels[i]- np.array([1, 0]))]==1:
             new_lab = indeces_labels[i]- np.array([1, 0])
             j = I_J_conv(new_lab[0], new_lab[1], n1)
             H[4*i+2, 4*j] = -t   # Particle B_up hops to A_up
             H[4*i + 3, 4*j+1] = -t  # Particle B_down hops to A_down

             # Next-nearest Neighbour hoppings
             H[4*i, 4*j+1] = -1j*t_so   # Particle A_up hops to A_down
             H[4*i + 1, 4*j] = -1j*t_so  # Particle A_down hops to A_up
             H[4*i+2, 4*j+3] = 1j*t_so   # Particle B_up hops to B_down
             H[4*i + 3, 4*j + 2] = 1j*t_so  # Particle B_down hops to B_up

        except:
            _=0

    H = H + H.transpose().conj()
    assert np.allclose(H, H.conj().T)
    return H, indeces_labels






 

      
                 