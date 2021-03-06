## IMSHOW MODULE ---------------------------------------------------------------
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.figure()
ax = plt.gca()
im = ax.imshow(u.T, origin = 'lower', extent=[-h_x,l+h_x,-h_y,w+h_y])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.show()

## P,U,V SHAPE CHECK -----------------------------------------------------------
print(p.shape)
print(u.shape)
print(v.shape)
