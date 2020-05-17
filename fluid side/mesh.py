import numpy as np
import matplotlib.pyplot as plt

'''
Objective:
1. Take inputs of L,W,H
2. Return grids for pressure and velocities
3. Comment out displaying the grids

'''

def init_mesh(L,W,h_x,h_y):
    grid = np.mgrid[-h_x:L+h_x:h_x,-h_y:W+h_y:h_y]
    return grid

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

def create_mesh(l,w,h_x,h_y,viewgrid):

    L = l+h_x
    W = w+h_y

    n_x = int(L/h_x)
    n_y = int(W/h_y)

    grid = init_mesh(L,W,h_x,h_y)
    pgrid = cellcenter(grid)
    ugrid = facecenter(grid,0)
    vgrid = facecenter(grid,1)

    # print(grid)
    # print(grid)
    # print(pgrid)

    if viewgrid==1:
        for i in range(n_x+2):
            for j in range(n_y+2):
                plt.figure(num=1)
                plt.scatter(grid[0,i,j],grid[1,i,j],c='k',s=10)

                if i<n_x+1 and j<n_y+1:
                    plt.scatter(pgrid[0,i,j],pgrid[1,i,j],c='r',s=10)

                if i<n_x+1 and j<n_y:
                    plt.scatter(vgrid[0,i,j],vgrid[1,i,j],c='b',s=10,marker='x')

                if j<n_y+1 and i<n_x:
                    plt.scatter(ugrid[0,i,j],ugrid[1,i,j],c='g',s=10,marker='x')

        plt.plot([0,l],[0,0],'k')
        plt.plot([l,l],[0,w],'k')
        plt.plot([l,0],[w,w],'k')
        plt.plot([0,0],[w,0],'k')
        plt.show()

    return grid,pgrid,ugrid,vgrid

def matrix(grid):
    n_x = len(grid[0])-1
    n_y = len(grid[0][0])-1

    m = np.zeros((n_x,n_y))

    return m
