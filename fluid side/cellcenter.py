import numpy as np
import matplotlib.pyplot as pyplot

'''
Objective:
1. Interpolate cell centers


'''

def cellcenter(grid):
    n_x = len(grid[0])-1
    n_y = len(grid[0][0])-1

    pgrid = np.zeros((2,n_x,n_y))

    for i in range(0,n_x):
        for j in range(0,n_y):
            pgrid[0][i][j] = (grid[0][i][j]+grid[0][i+1][j]+grid[0][i][j+1]\
                                +grid[0][i+1][j+1])/4
                                
            pgrid[1][i][j] = (grid[1][i][j]+grid[1][i+1][j]+grid[1][i][j+1]\
                                +grid[1][i+1][j+1])/4

    return pgrid
