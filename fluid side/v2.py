import numpy as np
import matplotlib.pyplot as plt
import time
start = time.time()

# GRID SETTINGS
l = 1
w = 0.5
nx = 20
ny = 20
h_x = l/nx
h_y = w/ny

# TIME STEP SETTINGS
dt = 2e-3
tsteps = 1000
clock = 0

# ITERATION SETTINGS
maxiter = 1000
maxerror = 1e-5
beta = 1.2 #Over-relaxation parameter
iter = 0

# STATIC PARAMETERS
rho_s = 1000
rho_w = 1e10
mu_s = 1e-5
p_s = 1
gx = 0
gy = 0

# INITIALIZING MATRICES, defined using length of the interval and not number of grid points
u = np.zeros((nx+2,ny+2)) #Having an extra set of rows at the outlet for having ghost cells to implement boundary conditions
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
u[0,:] = 0.1 # Hortizontal velocity at the inlet
p[-1,:] = p_s # Pressure boundary condition at the outlet

# GRIDS FOR PLOTTING
uu=np.zeros((nx+1,ny+1))
vv=np.zeros((nx+1,ny+1))
uu=0.5*(u[0:nx,1:ny+1]+u[0:nx,0:ny])
vv=0.5*(v[1:nx+1,0:ny]+v[0:nx,0:ny])
x = np.linspace(0.,l,nx)
y = np.linspace(0.,w,ny)

for t in range(0,tsteps):
    clock = clock + dt
    start = time.time()

    #DYNAMIC BOUNDARIES
    # Tangential velocities at the wall
    u[:,0] = -u[:,1]
    u[:,-1] = -u[:,-2]

    # Horizontal velocities at the outlet
    u[-1,:] = u[-2,:]

    # Vertical velocities at the inlet and outlet
    v[0,:] = v[1,:]
    v[-1,:] = v[-2,:]

    # Pressure boundary condition at the inlet
    p[0,:] = p[1,:]

    # U STAR
    for i in range(1,nx+1):
        for j in range(1,ny+1):
            ustar[i,j]=u[i,j]+dt*(-0.25*(((u[i+1,j]+u[i,j])**2-(u[i,j]+\
                u[i-1,j])**2)/h_x +((u[i,j+1]+u[i,j])*(v[i+1,j]+\
                v[i,j])-(u[i,j]+u[i,j-1])*(v[i+1,j-1]+v[i,j-1]))/h_y)+\
                0.5*(mu[i+1,j] + mu[i,j])/(0.5*(rho[i+1,j]+rho[i,j]))*\
                ((u[i+1,j]-2*u[i,j]+u[i-1,j])/h_x**2+\
                (u[i,j+1]-2*u[i,j]+u[i,j-1])/h_y**2 )+gx)


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
            p1[i,j]= (0.5/dt)*( (ustar[i,j]-ustar[i-1,j])/h_x+(vstar[i,j]-vstar[i,j-1])/h_y)
            p2[i,j]=1.0/( (1./h_x)*( 1./(h_x*(rt[i+1,j]+rt[i,j]))+\
                1./(h_x*(rt[i-1,j]+rt[i,j])))+\
                (1./h_y)*(1./(h_y*(rt[i,j+1]+rt[i,j]))+\
                1./(h_y*(rt[i,j-1]+rt[i,j]))))

    iter=0
    while True:
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
            

    print(clock, time.time() - start, iter, np.abs(pn-p).max())

    uu=0.5*(u[0:nx,1:ny+1]+u[0:nx,0:ny])
    vv=0.5*(v[1:nx+1,0:ny]+v[0:nx,0:ny])
    yy,xx=np.mgrid[0:(nx-1)*h_x:nx*1j,0:(ny-1)*h_y:ny*1j]
    plt.figure(num = 1)
    plt.quiver(xx,yy,uu.T,vv.T)
    plt.figure(num = 2)
    plt.streamplot(x,y,uu.T,vv.T);
    plt.figure()
    plt.imshow(p.T, origin='lower')
