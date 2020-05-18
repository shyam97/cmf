import numpy as np
import matplotlib.pyplot as plt
from mesh import *

l = 0.2
w = 0.1
nx = 20
ny = 20
md = 0
h_x = l/(nx)
h_y = w/(ny)

rho_s = 1
p_a = 101325
mu_s = 0.05

dt = 1e-4
tsteps = 1
beta = 0.8 #Over-relaxation parameter
maxiter = 1000
maxError = 1e-3

grid,pgrid,ugrid,vgrid = create_mesh(l,w,h_x,h_y,md)

p = matrix(pgrid) + p_a
p1 = p
p2 = p
u = matrix(ugrid) + 0.5
ustar = matrix(ugrid)
vstar = matrix(vgrid)
v = matrix(vgrid)
rho = matrix(pgrid) + rho_s
mu = matrix(pgrid) + mu_s
gx = 0
gy = 0

print(nx)
print(ny)
print(p.shape)
print(u.shape)
print(v.shape)

v[:,0] = 0 #Vertical component of velocity at bottom wall
v[:,-1] = 0 #Vertical component of velocity at top wall
u[0,:] = 1 #Horizontal component of velocity defined at the inlet, to be changed to another profile later
p[:,0] = 0 #Pressure at the bottom wall ghost cells
p[:,-1] = 0 #Pressure at the top wall ghost cells
# rho[:,0] = 1e10 #Density at the bottom wall ghost cells
# rho[:,-1] = 1e10 #Density at the top wall ghost cells
time = 0
for t in range(0,tsteps):
    time = time + dt
    print(time)

    u[:,0] = -u[:,1] #Horizontal component of velocity at bottom wall ghost
    u[:,-1] = -u[:,-2] #Horizontal component of velocity at top wall ghost
    v[0,:] = -v[1,:] #Vertical component of velocity at ghost cell at inlet
    v[-1,:] = v[-2,:] #Vertical component of velocity at ghost cell at outlet
    p[-1,:] = 2*p_a - p[-2,:]  #Pressure at the ghost cells at the outlet
    p[0,:] = 2*p_a #Pressure at the inlet (to be changed later)

    for i in range(1,len(u)-1):
        for j in range(1,len(u[0])-1):

            #Advection terms - x component
            Ax = ((u[i+1,j] + u[i,j])**2 - (u[i,j]+u[i-1,j])**2)/(4*h_x) +\
                 ((u[i,j+1] + u[i,j])*(v[i+1,j] + v[i,j]) - (u[i,j] + u[i,j-1])*\
                   (v[i+1,j-1] + v[i,j-1]))/(4*h_y)

            #Diffusion terms - x component
            Dx = 0.5*(mu[i+1,j] + mu[i,j])*\
                ((u[i+1,j] -2*u[i,j] + u[i-1,j])/(h_x**2) +\
                 (u[i,j+1] -2*u[i,j] + u[i,j-1])/(h_y**2))

            #Predictor step - x component
            ustar[i,j] = u[i,j] + dt*(-Ax + (2/(rho[i,j] + rho[i+1,j]))*Dx + gx)

    for i in range(1,len(v)-1):
        for j in range(1,len(v[0])-1):

            #Advection terms - y component
            Ay = ((u[i,j+1] + u[i,j])*(v[i,j] + v[i+1,j]) - (u[i-1,j+1] + u[i-1,j])*\
                  (v[i,j] + v[i-1,j]))/(4*h_x) + ((v[i,j+1] + v[i,j])**2 -\
                  (v[i,j] + v[i,j-1])**2)/(4*h_y)

            #Diffusion terms - y component
            Dy = 0.5*(mu[i,j+1] + mu[i,j])*\
                ((v[i+1,j] - 2*v[i,j] + v[i-1,j])/(h_x**2) +\
                 (v[i,j+1] - 2*v[i,j] + v[i,j-1])/(h_y**2))

            #Predictor step - y component
            vstar[i,j] = v[i,j] + dt*(-Ay + gy + 2*Dy/((rho[i,j+1] + rho[i,j])))

    r = rho.copy()
    r[:,0] = 1e3
    r[:,-1] = 1e3

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    plt.figure()
    ax = plt.gca()
    im = ax.imshow(ustar.T, origin = 'lower', extent=[-h_x,l+h_x,-h_y,w+h_y])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    plt.show()

#Defining the constants in the corrector step
    for i in range(1,len(p)-1):
        for j in range(1,len(p[0])-1):
            #Constant on the LHS
            p1[i,j]= (0.5/dt)*((ustar[i,j]-ustar[i-1,j])/h_x +\
                                (vstar[i,j]-vstar[i,j-1])/h_y)
            #Constant on the RHS
            p2[i,j]=1.0/((1./h_x)*(1./(h_x*(r[i+1,j]+r[i,j]))+\
                1./(h_x*(r[i-1,j]+r[i,j])))+\
                (1./h_y)*(1./(h_y*(r[i,j+1]+r[i,j]))+\
                1./(h_y*(r[i,j-1]+r[i,j]))))

    #Pressure equation - to be solved implicitly

    iter=0
    while True:
        pn=p.copy()
        iter=iter+1
        for i in range(1,len(p)-1):
            for j in range(1,len(p[0])-1):
                p[i,j]=(1.0-beta)*p[i,j]+beta*p2[i,j]*(\
                    (1./h_x)*( p[i+1,j]/(h_x*(r[i+1,j]+r[i,j]))+\
                    p[i-1,j]/(h_x*(r[i-1,j]+r[i,j])))+\
                    (1./h_y)*( p[i,j+1]/(h_y*(r[i,j+1]+r[i,j]))+\
                    p[i,j-1]/(h_y*(r[i,j-1]+r[i,j])))- p1[i,j])
        print(np.abs(pn-p).max())
        if np.abs(pn-p).max()<maxError:
            break
        if iter>maxiter:
            break

    #Corrector step - to be done implicitly

    #CORRECT THE u-velocity
    for i in range(1,len(u)-1):
        for j in range(1,len(u[0])-1):
            u[i,j]=ustar[i,j]-dt*(2.0/h_x)*(p[i+1,j]-p[i,j])/(rho[i+1,j]+rho[i,j])

    #CORRECT THE v-velocity
    for i in range(1,len(v)-1):
        for j in range(1,len(v[0])-1):
            v[i,j]=vstar[i,j]-dt*(2.0/h_y)*(p[i,j+1]-p[i,j])/(rho[i,j+1]+rho[i,j])

from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.figure()
ax = plt.gca()
im = ax.imshow(u.T, origin = 'lower', extent=[-h_x,l+h_x,-h_y,w+h_y])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.show()
