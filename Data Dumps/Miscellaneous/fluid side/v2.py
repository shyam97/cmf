import numpy as np
import matplotlib.pyplot as plt
import time
start = time.time()
from mpl_toolkits.axes_grid1 import make_axes_locatable

# GRID SETTINGS
l = 0.10
w = 0.010
nx = 50
ny = 100
h_x = l/nx
h_y = w/ny
md = 0

x = np.linspace(0.,l,nx)
y = np.linspace(0.,w,ny)

# TIME STEP SETTINGS
dt = 1e-6
tsteps = 100
flow = dt*tsteps
clock = 0
switch = 0

# ITERATION SETTINGS
maxiter = 1000
maxerror = 1e-6
beta = 1.2 #Over-relaxation parameter
iter = 0

# STATIC PARAMETERS
rho_s = 1
rho_w = 1e10
mu_s = 1e-5
p_s = 1
gx = 0
gy = 0

# INITIALIZING MATRICES, defined using length of the interval and not number of grid points
u = np.zeros((nx+2,ny+2)) #Having an extra set of rows at the outlet for having ghost cells to implement boundary conditions
v = np.zeros((nx+2,ny+1))
p = np.zeros((nx+2,ny+2)) + p_s

# INTERMEDIATE VELOCITIES
ustar = u.copy()
vstar = v.copy()

# COEFFICIENTS FOR PRESSURE
p1 = p.copy()
p2 = p.copy()

# DENSITY AND VISCOSITY MATRIX
rho = np.ones((nx+2,ny+2))*rho_s
mu = np.ones((nx+2,ny+2))*mu_s

# INLET STATIC BOUNDARY
u[0,:] = 0
p[0,:] = 2

# OUTLET STATIC BOUNDARY
p[-1,:] = 1

# BOTTOM WALL STATIC BOUNDARY
v[:,0] = 0
p[:,0] = 0

# TOP WALL STATIC BOUNDARY
v[:,-1] = 0
p[:,-1] = 0

for t in range(0,tsteps):
    clock = clock + dt
    start = time.time()

    # INLET DYNAMIC BOUNDARY
    v[0,:] = -1*v[1,:]

    # OUTLET DYNAMIC BOUNDARY
    v[-1,:] = 1*v[-2,:]
    u[-1,:] = u[-2,:]

    # BOTTOM WALL DYNAMIC BOUNDARY
    u[:,0] = -1*u[:,1]

    # TOP WALL DYNAMIC BOUNDARY
    u[:,-1] = -1*u[:,-2]

    # U STAR
    for i in range(1,nx+1):
        for j in range(1,ny+1):
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


    #V STAR
    for i in range(1,nx+1):
        for j in range(1,ny):
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

    rt=rho.copy()
    rt[:,0] = rho_w
    rt[:,-1] = rho_w

    for i in range(1,nx+1):
        for j in range(1,ny+1):
            p1[i,j]= (0.5/dt)*( (ustar[i,j]-ustar[i-1,j])/h_x+(vstar[i,j]-vstar[i,j-1])/h_y)
            p2[i,j]=1.0/( (1./h_x)*( 1./(h_x*(rt[i+1,j]+rt[i,j]))+\
                1./(h_x*(rt[i-1,j]+rt[i,j])))+\
                (1./h_y)*(1./(h_y*(rt[i,j+1]+rt[i,j]))+\
                1./(h_y*(rt[i,j-1]+rt[i,j]))))

    iter=0
    while True:
        # p[0,:] = p[1,:]
        pn=p.copy()
        iter=iter+1
        for i in range(1,nx+1):
            for j in range(1,ny+1):
                p[i,j]=(1.0-beta)*p[i,j]+beta*p2[i,j]*(\
                    (1./h_x)*( p[i+1,j]/(h_x*(rt[i+1,j]+rt[i,j]))+\
                    p[i-1,j]/(h_x*(rt[i-1,j]+rt[i,j])))+\
                    (1./h_y)*( p[i,j+1]/(h_y*(rt[i,j+1]+rt[i,j]))+\
                    p[i,j-1]/(h_y*(rt[i,j-1]+rt[i,j])))-p1[i,j])

        # print(np.abs(pn-p).max())
        if np.abs(pn-p).max()<maxerror:
            break
        if iter>maxiter:
            break

    #CORRECT THE u-velocity
    for i in range(1,nx+1):
        for j in range(1,ny+1):
            u[i,j]=ustar[i,j]-dt*(2.0/h_x)*(p[i+1,j]-p[i,j])/(rho[i+1,j]+rho[i,j])

    #CORRECT THE v-velocity
    for i in range(1,nx+1):
        for j in range(1,ny):
            v[i,j]=vstar[i,j]-dt*(2.0/h_y)*(p[i,j+1]-p[i,j])/(rho[i,j+1]+rho[i,j])


    u[0,:] = u[1,:]

    print(clock, time.time() - start, iter, np.abs(pn-p).max())
    # uu=0.5*(u[0:nx,1:ny+1]+u[0:nx,0:ny])
    # vv=0.5*(v[1:nx+1,0:ny]+v[0:nx,0:ny])
    # yy,xx=np.mgrid[0:(nx-1)*h_x:nx*1j,0:(ny-1)*h_y:ny*1j]

    # plt.figure(num=1, figsize = (10,5))
    # plt.clf()
    # plt.quiver(xx,yy,uu.T,vv.T)
    # ax = plt.gca()
    # im = ax.imshow(u.T, origin = 'lower', extent=[-h_x,l+h_x,-h_y,w+h_y])
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("right", size="5%", pad=0.05)
    # plt.colorbar(im, cax=cax)
    # plt.streamplot(x,y,uu.T,vv.T);
    # plt.pause(0.001)

plt.figure(num=1, figsize = (10,5))
plt.clf()
plt.quiver(xx,yy,uu.T,vv.T)
plt.show()

plt.figure(num=3, figsize = (10,5))
ax = plt.gca()
im = ax.imshow(p[1:nx+1,1:ny+1].T, origin = 'lower', extent=[-h_x,l+h_x,-h_y,w+h_y])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.show()
