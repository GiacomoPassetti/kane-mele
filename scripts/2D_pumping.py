import numpy as np
import cmath
from numpy import linalg as LA
#import matplotlib.pyplot as plt
#Code for 2D SSH model with different boundary conditions
#First m, then sigma

def Boundary_fluctuation(corr, f):
    Q_B_2=np.dot(f,np.dot(corr,f))
    return Q_B_2

def Boundary_charge(dens, f):
    Q_B=np.dot(f,dens)
    return Q_B


def env(r, r1, r2):
    if r<=r1:
        return 1.
    elif r>=r2:
        return 0.
    else:
        return 1-(r-r1)/(r2-r1) 

def envelope_function_r(Nx, Ny, r1, r2):
    r=np.zeros(Nx*Ny)
    rx=np.kron(np.arange(1, Nx+1), np.ones(Ny))
    ry=np.kron(np.ones(Nx), np.arange(1, Ny+1))
    r=np.sqrt(rx**2+ry**2)
    f=np.zeros(Nx*Ny)
    for i in range(Nx*Ny):
        f[i]=env(r[i], r1, r2)
    return f

def Correlation(H):
	w, v = LA.eigh(H)   
	N=len(w)
	n_diag=np.zeros(N)
	n_diag[:int(N/2)]=1.0
	n_diag=np.diag(n_diag)
	vdag=np.conjugate(np.transpose(v))
	mat1=np.multiply(np.diag(n_diag)[:,None], vdag)
	n_mat=np.dot(v,mat1)
	nn_mat=np.abs(n_mat)**2
	nn_mat=-nn_mat+np.diag(np.diag(n_mat))
	return nn_mat

def density(H):
    w, v = LA.eigh(H)
    vec = v[:, :int(N/2)]**2
    dens = np.sum(vec, axis=1)
    return dens, w, v


####### to adjust ######
import sys
########################

L=30
Nx=L
Ny=L
N=Nx*Ny
Z=2
F=[]
phi_list=np.pi*np.arange(0.0, 2.0, 0.02).round(3)# p.linspace(0, 2*np.pi, 100)

for phi in phi_list:
    tx=1.0
    ty=1.0
    dtx=0.1*np.cos(phi)
    dty=0.1
    phi_0=np.pi/4
    dv=0.2*np.cos(phi+phi_0)
    # We need to pump in the x direction towards the y axis.
    
    r1=10
    r2=20
    f=envelope_function_r(Nx, Ny, r1, r2)
    
    tm_y=np.array([ty-dty, ty+dty]*int(Ny/2))[:-1]
    tm_x=np.array([tx-dtx, ty+dtx]*int(Nx/2))[:-1]
    
    H_x=-np.kron(np.diag(tm_x, k=1), np.diag([-1, +1]*int(Ny/2)))-np.kron(np.diag(tm_x, k=-1), np.diag([-1, +1]*int(Ny/2)))
    H_y=-np.kron(np.eye(Nx), np.diag(tm_y, k=1))-np.kron(np.eye(Nx), np.diag(tm_y, k= -1))
    Pot=dv*np.kron(np.diag([+1, -1]*int(Nx/2)),np.diag([+1, -1]*int(Ny/2)))
    
    H=H_x+H_y+Pot
    
    # Calculate Boundary charge fluctuations for open boundary conditions
    nn_mat, w, v=density(H)
    fluc=Boundary_charge(nn_mat-1/2, f)
    print(fluc)
    F+=[fluc]
    np.save('Daten_pump/QB_L'+str(L)+'_r1.'+str(r1)+'_r2.'+str(r2)+'_phi'+str(phi)+'.npy', fluc)
    np.save('Daten_pump/w_L'+str(L)+'_r1.'+str(r1)+'_r2.'+str(r2)+'_phi'+str(phi)+'.npy', w)
    np.save('Daten_pump/v_L'+str(L)+'_r1.'+str(r1)+'_r2.'+str(r2)+'_phi'+str(phi)+'.npy', v)

#plt.figure()
#plt.plot(phi_list, F, '-o')
