import numpy as np
import matplotlib.pyplot as plt

'''
Objective:
1. Create the mesh coordinates


'''

# Uniform mesh for now
def init_mesh(L,W,n_x,n_y):

    step = L/n_x

    grid = step*np.mgrid[0:n_x+1,0:n_y+1]

    return grid
