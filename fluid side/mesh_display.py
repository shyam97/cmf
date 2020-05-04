import numpy as np
import matplotlib.pyplot as plt

'''
Objective:
1. Display the mesh


'''

from cellcenter import *
from init_mesh import *

L = 10
W = 5
n_x = 10
n_y = 10

grid = init_mesh(L,W,n_x,n_y)
pgrid = cellcenter(grid)

# print(grid)
# print(pgrid)

for i in range(n_x+1):
    for j in range(n_y+1):
        plt.figure(num=1)
        plt.scatter(grid[0,i,j],grid[1,i,j],c='k')

        if i<n_x and j<n_y:
            plt.scatter(pgrid[0,i,j],pgrid[1,i,j],c='r')

plt.show()
