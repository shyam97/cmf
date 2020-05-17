import numpy as np
import matplotlib.pyplot as plt

'''
Objective:
1. Display the mesh


'''

def meshes:
    from cellcenter import *
    from init_mesh import *
    from facecenter import *

    l = 2
    w = 1
    h_x = 0.1
    h_y = 0.1

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

    # for i in range(n_x+2):
    #     for j in range(n_y+2):
    #         plt.figure(num=1)
    #         plt.scatter(grid[0,i,j],grid[1,i,j],c='k',s=10)
    #
    #         if i<n_x+1 and j<n_y+1:
    #             plt.scatter(pgrid[0,i,j],pgrid[1,i,j],c='r',s=10)
    #
    #         if i<n_x+1 and j<n_y:
    #             plt.scatter(vgrid[0,i,j],vgrid[1,i,j],c='b',s=10,marker='x')
    #
    #         if j<n_y+1 and i<n_x:
    #             plt.scatter(ugrid[0,i,j],ugrid[1,i,j],c='g',s=10,marker='x')
    #
    # plt.plot([0,l],[0,0],'k')
    # plt.plot([l,l],[0,w],'k')
    # plt.plot([l,0],[w,w],'k')
    # plt.plot([0,0],[w,0],'k')
    # plt.show()
    #
    # plt.figure()
    # plt.plot([0,l],[0,0],'k')
    # plt.plot([l,l],[0,w],'k')
    # plt.plot([l,0],[w,w],'k')
    # plt.plot([0,0],[w,0],'k')
    # plt.scatter(ugrid[0,0,1],ugrid[1,0,1])
    # plt.show()

    return grid,pgrid,ugrid,vgrid
