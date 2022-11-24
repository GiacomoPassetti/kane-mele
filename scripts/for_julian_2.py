import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as LA

# should include magnetic flux
def winkel(x, y):
	return np.arctan2(y,x)%(2*np.pi)


def Hamiltonian(f, alpha, a, m, t, magmom, Nx, Ny):
	sig0=np.array([[1, 0],[0, 1]])
	sigx=np.array([[0, 1],[1, 0]])
	sigy=np.array([[0, -1j],[1j, 0]])
	sigz=np.array([[1, 0],[0, -1]])
	
	x = np.arange(Nx) - int(Nx/2) + 0.5
	y = np.arange(Ny) - int(Ny/2) + 0.5
	X = np.kron(x, np.ones(Ny))
	Y = np.kron(np.ones(Ny), y)
	args = winkel(X, Y)
	args = args.reshape((Nx, Ny))
	hoppingx = np.zeros((Nx, Ny), dtype = complex)
	hoppingy = np.zeros((Nx, Ny), dtype = complex)
	deltax = args[1:, :] - args[:-1, :]
	deltax = deltax % (2*np.pi)
	deltax = deltax + np.pi
	deltax = deltax % (2*np.pi)
	deltax = deltax - np.pi
	print(deltax)
	deltay = args[:, 1:] - args[:, :-1]
	deltay = deltay % (2*np.pi)
	deltay = deltay + np.pi
	deltay = deltay % (2*np.pi)
	deltay = deltay - np.pi
	print(deltay)
	hoppingx[:-1, :]  = np.exp(f*1j*deltax)
	hoppingx = hoppingx.flatten()
	print(hoppingx)
	hoppingx1 = np.diag(hoppingx[:-Ny], k=Ny)
	hoppingx2 = np.diag(np.conjugate(hoppingx)[:-Ny], k=-Ny)
	hoppingy[:, :-1] = np.exp(f*1j*deltay)
	hoppingy = hoppingy.flatten()
	hoppingy1 = np.diag(hoppingy[:-1], k=1)
	hoppingy2 = np.diag(np.conjugate(hoppingy[:-1]), k=-1)
		
	H0=1j*alpha/(2*a)*np.kron(hoppingx1-hoppingx2, np.kron(sigx, sigx))
	plt.matshow(np.imag(H0))
	H0+=1j*alpha/(2*a)*np.kron(hoppingy1-hoppingy2, np.kron(sigx, sigy))
	H0+=(4*t-1)*np.kron(np.eye(Nx), np.kron(np.eye(Ny), np.kron(sigz, sig0)))
	H0+=-t*np.kron(hoppingx1+hoppingx2, np.kron(sigz, sig0))
	H0+=-t*np.kron(hoppingy1+hoppingy2, np.kron(sigz, sig0))
	Hz=-2/4/m*np.kron(np.eye(Nx), np.kron(np.eye(Ny), np.kron(sig0, magmom[0]*sigx+magmom[1]*sigy+magmom[2]*sigz)))
	H=H0+Hz #zum testen
	return H
	
import sys
	
alpha=1
a=1
m=1/2
t= 1 / (2*m)
magmom=1/np.sqrt(2)*np.array([0.01,0.01,0]) #for the different directions of the Zeeman field

Nx=10
Ny=10
for f in [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]:

	H=Hamiltonian(f, alpha, a, m, t, magmom, Nx, Ny)
	w, v = LA.eigh(H)
	np.save('Data/'+str(f), w)
	
'''
import matplotlib.pyplot as plt
plt.figure()
plt.plot(w, 'o')
plt.show()
'''

#erstmal mit for_julian.py bei f=0 
