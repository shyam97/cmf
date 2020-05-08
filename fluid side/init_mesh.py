import numpy as np
import matplotlib.pyplot as plt

'''
Objective:
1. Create the mesh coordinates


'''

# Uniform mesh for now
def init_mesh(L,W,h_x,h_y):
    grid = np.mgrid[-h_x:L+h_x:h_x,-h_y:W+h_y:h_y]
    return grid
