import numpy as np
import matplotlib.pyplot as plt
from mesh import *

l = 2
w = 1
h_x = 0.1
h_y = 0.1
md = 0

rho_s = 1.25
p_a = 101325

grid,pgrid,ugrid,vgrid = create_mesh(l,w,h_x,h_y,md)

p = matrix(pgrid) + p_a
u = matrix(ugrid)
v = matrix(vgrid)
rho = matrix(pgrid) + rho_s
visc = matrix(pgrid)

plt.figure()
plt.imshow(p.T, origin = 'lower', extent=[0,l,0,w])
plt.colorbar()
plt.show()
