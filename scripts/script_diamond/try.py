import numpy as np


sigma_x = np.array([[0, 1],[1, 0]])
sigma_y = np.array([[0, -1j],[1j, 0]])
sigma_z = np.array([[1, 0],[0, -1]])


print(np.tensordot(sigma_x, sigma_z, 0))
print(np.tensordot(sigma_z, sigma_x, 0))
