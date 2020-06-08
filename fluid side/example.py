from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
# from IPython import display
# %matplotlib inline

#domain size and physical variables
Lx=5.0
Ly=5.0
gx=0.0
gy=0.0
rho1=1.0

#dynamic viscosity
m0=0.05

#tangential velocities
unorth=0.
usouth=0.
veast=0.
vwest=0.

#Numerical variables
nx=100
ny=100
dt=0.0004
nstep=100
maxiter=1000
maxError=0.001
beta=1.2

u=np.zeros((nx+1,ny+2))
v=np.zeros((nx+2,ny+1))
p=np.zeros((nx+2,ny+2))

ut=np.zeros((nx+1,ny+2))
vt=np.zeros((nx+2,ny+1))

#temp variables for pressure calculations
tmp1=np.zeros((nx+2,ny+2))
tmp2=np.zeros((nx+2,ny+2))

#velocities at center of grid for plotting
uu=np.zeros((nx+1,ny+1))
vv=np.zeros((nx+1,ny+1))
x = np.linspace(0.,Lx,nx)
y = np.linspace(0.,Ly,ny)

#Define the grid
dx=Lx/(nx)
dy=Ly/(ny)


r=np.ones((nx+2,ny+2))*rho1

time=0.
for steps in range(nstep):

    #tangential velocity at boundaries
    u[:,0]=2.*usouth-u[:,1]
    u[:,-1]=2.*unorth-u[:,-2]
    v[0,:]=2.*vwest-v[1,:]
    v[-1,:]=2.*veast-v[-2,:]
    u[0,:] = 1
    u[-1,:] = 1

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


    ## Smagorinsky model 
    
    τ_0, ρ, ω, ω_total, ω_0 = sp.symbols("tau_0 rho omega omega_total omega_0", positive=True, real=True)
    ν_0, C_S, S, Π = sp.symbols("nu_0, C_S, |S|, Pi", positive=True, real=True)
    
    Seq = sp.Eq(S, 3 * ω / 2 * Π)
    print(Seq)
    
    def relaxation_rate_from_lattice_viscosity(ν_0, C_S,S):
        omega = 2/(6*(ν_0 + C_S ** 2 * S) + 1)
        return omega
    
    ω = relaxation_rate_from_lattice_viscosity(ν_0, C_S,S)
    
    Seq2 = Seq.subs(ω, relaxation_rate_from_lattice_viscosity(ν_0, C_S,S ))
    print(Seq2)
    
    solveRes = sp.solve(Seq2, S)
    print(solveRes)
    assert len(solveRes) == 1
    SVal = solveRes[0]
    
    def lattice_viscosity_from_relaxation_rate(τ_0):
        ω_0 = 1/τ_0
        return ω_0
    
    SVal = SVal.subs(ν_0, lattice_viscosity_from_relaxation_rate(τ_0)).expand()
    
    def second_order_moment_tensor(function_values, stencil):
        assert len(function_values) == len(stencil)
        dim = len(stencil[0])
        return sp.Matrix(dim, dim, lambda i, j: sum(c[i] * c[j] * f for f, c in zip(function_values, stencil)))
    
    
    def frobenius_norm(matrix, factor=1):
        return sp.sqrt(sum(i*i for i in matrix) * factor)


    #Compute source term and the coefficient for p(i,j)
    rt=r.copy()
    lrg=1000.
    rt[:,0]=lrg
    rt[:,-1]=lrg
    rt[0,:]=lrg
    rt[-1,:]=lrg

    for i in range(1,nx+1):
        for j in range(1,ny+1):
            tmp1[i,j]= (0.5/dt)*( (ut[i,j]-ut[i-1,j])/dx+(vt[i,j]-vt[i,j-1])/dy)
            tmp2[i,j]=1.0/( (1./dx)*( 1./(dx*(rt[i+1,j]+rt[i,j]))+\
                1./(dx*(rt[i-1,j]+rt[i,j])))+\
                (1./dy)*(1./(dy*(rt[i,j+1]+rt[i,j]))+\
                1./(dy*(rt[i,j-1]+rt[i,j]))))

    iter=0
    while True:
        pn=p.copy()
        iter=iter+1
        for i in range(1,nx+1):
            for j in range(1,ny+1):
                p[i,j]=(1.0-beta)*p[i,j]+beta*tmp2[i,j]*(\
                    (1./dx)*( p[i+1,j]/(dx*(rt[i+1,j]+rt[i,j]))+\
                    p[i-1,j]/(dx*(rt[i-1,j]+rt[i,j])))+\
                    (1./dy)*( p[i,j+1]/(dy*(rt[i,j+1]+rt[i,j]))+\
                    p[i,j-1]/(dy*(rt[i,j-1]+rt[i,j])))-tmp1[i,j])
        if np.abs(pn-p).max()<maxError:
            break
        if iter>maxiter:
            break

    #CORRECT THE u-velocity
    for i in range(1,nx):
        for j in range(1,ny+1):
            u[i,j]=ut[i,j]-dt*(2.0/dx)*(p[i+1,j]-p[i,j])/(r[i+1,j]+r[i,j])

    #CORRECT THE v-velocity
    for i in range(1,nx+1):
        for j in range(1,ny):
            v[i,j]=vt[i,j]-dt*(2.0/dy)*(p[i,j+1]-p[i,j])/(r[i,j+1]+r[i,j])

    time=time+dt

    #Plot the velocity field
    uu=0.5*(u[0:nx,1:ny+1]+u[0:nx,0:ny])
    vv=0.5*(v[1:nx+1,0:ny]+v[0:nx,0:ny])
    yy,xx=np.mgrid[0:(nx-1)*dx:nx*1j,0:(ny-1)*dx:ny*1j]
    plt.clf()
    plt.quiver(xx,yy,uu.T,vv.T)
    plt.show()
    # display.clear_output(wait=True)
    # display.display(plt.gcf())
    print(steps)
    
    

plt.clf()
plt.streamplot(x,y,uu.T,vv.T);
# plt.show()
