import numpy as np
import matplotlib.pyplot as plt
from mesh import *

l = 0.2
w = 0.1
h_x = 0.01
h_y = 0.001
md = 0

rho_s = 1.25
p_a = 101325

grid,pgrid,ugrid,vgrid = create_mesh(l,w,h_x,h_y,md)

p = matrix(pgrid) + p_a
u = matrix(ugrid)
v = matrix(vgrid)
rho = matrix(pgrid) + rho_s
visc = matrix(pgrid)

v[:,0] = 0
v[:,-1] = 0
u[0,:] = 1
p[:,0] = 0
p[:,-1] = 0
rho[:,0] = 1e10
rho[:,-1] = 1e10
