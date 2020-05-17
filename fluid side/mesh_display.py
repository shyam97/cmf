import numpy as np
import matplotlib.pyplot as plt

'''
Objective:
1. Display the mesh


'''

from cellcenter import *
from init_mesh import *
from facecenter import *

#Mesh parameters
L = 1e-2
W = 5e-3
h_x = 1e-5
h_y = 1e-5
dt=0.0004
nstep=100
maxiter=1000
maxError=0.001
beta=1.2

n_x = int(L/h_x)
n_y = int(W/h_y)
r=np.ones((nx+2,ny+2))*rho1 #Residual matrix for density, not necessary if incompressible

#tangential velocities
unorth=2.
usouth=0.
veast=0.
vwest=0.

grid = init_mesh(L,W,h_x,h_y)
pgrid = cellcenter(grid)
ugrid = facecenter(grid,0)
vgrid = facecenter(grid,1)

# print(grid)
# print(pgrid)


#Displaying mesh
for i in range(n_x+1):
    for j in range(n_y+1):
        plt.figure(num=1)
        plt.scatter(grid[0,i,j],grid[1,i,j],c='k',s=10)

        if i<n_x and j<n_y:
            plt.scatter(pgrid[0,i,j],pgrid[1,i,j],c='r',s=10)

        if i<n_x:
            plt.scatter(vgrid[0,i,j],vgrid[1,i,j],c='b',s=10,marker='x')

        if j<n_y:
            plt.scatter(ugrid[0,i,j],ugrid[1,i,j],c='g',s=10,marker='x')


#Defining velocities and pressures

u=np.zeros((n_x+1,n_y+2))
v=np.zeros((n_x+2,n_y+1))
p=np.zeros((n_x+2,n_y+2))

ut=np.zeros((n_x+1,n_y+2))
vt=np.zeros((n_x+2,n_y+1))

#temp variables for pressure calculations
tmp1=np.zeros((n_x+2,n_y+2))
tmp2=np.zeros((n_x+2,n_y+2))

#velocities at center of grid for plotting
uu=np.zeros((n_x+1,n_y+1))
vv=np.zeros((n_x+1,n_y+1))
x = np.linspace(0.,Lx,n_x)
y = np.linspace(0.,Ly,n_y)


time=0.
for steps in range(nstep):
    
    #tangential velocity at boundaries
    u[:,0]=2.*usouth-u[:,1]
    u[:,-1]=2.*unorth-u[:,-2]
    v[0,:]=2.*vwest-v[1,:]
    v[-1,:]=2.*veast-v[-2,:]

    
    # TEMPORARY u-velocity                               
    for i in range(1,nx):
        for j in range(1,ny+1):
            ut[i,j]=u[i,j]+dt*(-0.25*(((u[i+1,j]+u[i,j])**2-(u[i,j]+\
                u[i-1,j])**2)/dx+((u[i,j+1]+u[i,j])*(v[i+1,j]+\
                v[i,j])-(u[i,j]+u[i,j-1])*(v[i+1,j-1]+v[i,j-1]))/dy)+\
                m0/(0.5*(r[i+1,j]+r[i,j]))*((u[i+1,j]-2*u[i,j]+u[i-1,j])/dx**2+\
                (u[i,j+1]-2*u[i,j]+u[i,j-1])/dy**2 )+gx)

    # TEMPORARY v-velocity                               
    for i in range(1,nx+1):
        for j in range(1,ny):
            vt[i,j]=v[i,j]+dt*(-0.25*(((u[i,j+1]+u[i,j])*(v[i+1,j]+\
                v[i,j])-(u[i-1,j+1]+u[i-1,j])*(v[i,j]+v[i-1,j]))/dx+\
                ((v[i,j+1]+v[i,j])**2-(v[i,j]+v[i,j-1])**2)/dy)+\
                m0/(0.5*(r[i,j+1]+r[i,j]))*((v[i+1,j]-2*v[i,j]+v[i-1,j])/dx**2+\
                (v[i,j+1]-2*v[i,j]+v[i,j-1])/dy**2 )+gy)






plt.show()
