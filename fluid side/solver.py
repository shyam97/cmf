import numpy as np
import matplotlib.pyplot as plt
from mesh import *

l = 2
w = 1
h_x = 0.1
h_y = 0.001

grid,pgrid,ugrid,vgrid = create_mesh(l,w,h_x,h_y)
