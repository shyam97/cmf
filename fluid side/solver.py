import numpy as np
import matplotlib.pyplot as plt
from mesh import *

l = 0.2
w = 0.1
h_x = 0.01
h_y = 0.01
md = 0

rho_s = 1.25
p_a = 101325
nu_s = 1.48e-5

dt = 1e-2
tsteps = 1e3

grid,pgrid,ugrid,vgrid = create_mesh(l,w,h_x,h_y,md)

p = matrix(pgrid) + p_a
u = matrix(ugrid)
ustar = matrix(ugrid)
vstar = matrix(vgrid)
v = matrix(vgrid)
rho = matrix(pgrid) + rho_s
nu = matrix(pgrid) + nu_s

v[:,0] = 0
v[:,-1] = 0
u[0,:] = 1
p[:,0] = 0
p[:,-1] = 0
rho[:,0] = 1e10 #
rho[:,-1] = 1e10

for t in range(1,tsteps):
    time = time + dt

    u[:,0] = -u[:,1]
    u[:,-1] = -u[:,-2]
    p[-1,:] = 2*p_a - p[-2,:]
