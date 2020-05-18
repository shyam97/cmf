import numpy as np
import matplotlib.pyplot as plt

# GRID SETTINGS
l = 10
w = 5
nx = 100
ny = 50
h_x = l/nx
h_y = w/ny

# TIME STEP SETTINGS
dt = 1e-4
tsteps = 3
time = 0

# ITERATION SETTINGS
maxiter = 100
maxerror = 1e-3
beta = 1.0 #Over-relaxation parameter

# STATIC PARAMETERS
rho_s = 1000
rho_w = 1e10
mu_s = 1e-5
p_s = 1e5
gx = 0
gy = 0

# INITIALIZING MATRICES, defined using length of the interval and not number of grid points
u = np.zeros((nx+1,ny+2))
v = np.zeros((nx+2,ny+1))
p = np.zeros((nx+2,ny+2))

# INTERMEDIATE VELOCITIES
ustar = u.copy()
vstar = v.copy()

# COEFFICIENTS FOR PRESSURE
p1 = p.copy()
p2 = p.copy()

# DENSITY AND VISCOSITY MATRIX
rho = np.ones((nx+2,ny+2))*rho_s
mu = np.ones((nx+2,ny+2))*mu_s

# STATIC BOUNDARIES
u[0,:] = 1

# GRIDS FOR PLOTTING
uu=np.zeros((nx+1,ny+1))
vv=np.zeros((nx+1,ny+1))
uu=0.5*(u[0:nx,1:ny+1]+u[0:nx,0:ny])
vv=0.5*(v[1:nx+1,0:ny]+v[0:nx,0:ny])
x = np.linspace(0.,l,nx)
y = np.linspace(0.,w,ny)

for t in range(0,tsteps):
    time = time + dt
    print(time)

    # TANGENTIAL VELOCITY AT BOUNDARIES
    u[:,0] = -u[:,1]
    u[:,-1] = -u[:,-2]

    v[0,:] = -v[1,:]
    v[-1,:] = v[-2,:]

    p[-1,:] = 2*p_s - p[-2,:]

    # U STAR
    for i in range(1,nx):
        for j in range(1,ny+1):
            ustar[i,j]=u[i,j]+dt*(-0.25*(((u[i+1,j]+u[i,j])**2-(u[i,j]+\
                u[i-1,j])**2)/h_x +((u[i,j+1]+u[i,j])*(v[i+1,j]+\
                v[i,j])-(u[i,j]+u[i,j-1])*(v[i+1,j-1]+v[i,j-1]))/h_y)+\
                0.5*(mu[i+1,j] + mu[i,j])/(0.5*(rho[i+1,j]+rho[i,j]))*\
                ((u[i+1,j]-2*u[i,j]+u[i-1,j])/h_x**2+\
                (u[i,j+1]-2*u[i,j]+u[i,j-1])/h_y**2 )+gx)
                #To find the u velocity component at the right end of the pipe,
                # use the discretized form of the continuity equation

    #V STAR
    for i in range(1,nx+1):
        for j in range(1,ny):
            vstar[i,j]=v[i,j]+dt*(-0.25*(((u[i,j+1]+u[i,j])*(v[i+1,j]+\
                v[i,j])-(u[i-1,j+1]+u[i-1,j])*(v[i,j]+v[i-1,j]))/h_x+\
                ((v[i,j+1]+v[i,j])**2-(v[i,j]+v[i,j-1])**2)/h_y)+\
                0.5*(mu[i,j+1] + mu[i,j])/(0.5*(rho[i,j+1]+rho[i,j]))*\
                ((v[i+1,j]-2*v[i,j]+v[i-1,j])/h_x**2+\
                (v[i,j+1]-2*v[i,j]+v[i,j-1])/h_y**2 )+gy)

    rt=rho.copy()
    rt[:,0]=rho_w
    rt[:,-1]=rho_w

    for i in range(1,nx+1):
        for j in range(1,ny+1):
            p1[i,j]= (0.5/dt)*( (ustar[i,j]-ustar[i-1,j])/dx+(vt[i,j]-vt[i,j-1])/dy)
            p2[i,j]=1.0/( (1./dx)*( 1./(dx*(rt[i+1,j]+rt[i,j]))+\
                1./(dx*(rt[i-1,j]+rt[i,j])))+\
                (1./dy)*(1./(dy*(rt[i,j+1]+rt[i,j]))+\
                1./(dy*(rt[i,j-1]+rt[i,j]))))
