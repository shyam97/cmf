import numpy as np
import matplotlib.pyplot as pyplot

'''
Objective:
1. Create the mesh coordinates


'''

def facecenter(grid,axis):

    if axis==0:
        n_x = len(grid[0])-1
        n_y = len(grid[0][0])-1

    else:
        n_x = len(grid[0])-1
        n_y = len(grid[0][0])-1

    fgrid = np.zeros((2,n_x,n_y))

    for i in range(0,n_x):
        for j in range(0,n_y):

            if axis==0:

                fgrid[0][i][j] = (grid[0][i+1][j]+grid[0][i+1][j+1])/2
                fgrid[1][i][j] = (grid[1][i+1][j]+grid[1][i+1][j+1])/2

            else:
                fgrid[0][i][j] = (grid[0][i][j+1]+grid[0][i+1][j+1])/2
                fgrid[1][i][j] = (grid[1][i][j+1]+grid[1][i+1][j+1])/2

    return fgrid
