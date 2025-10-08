## Import critical libraries/modules
import pencil as pc             # Local module for handling Pencil Code outputs
import numpy as np              # Arrays
import matplotlib.pyplot as plt # Plotting

## Read in VAR files from Pencil Code output
ff   = pc.read.var(trimall=True) # Pulling full VAR file output
rhop = ff.potself[0,:,:]         # Particle density at z=0 (z,y,x)
#amax = np.log10(rhop.max())
#amin = amax - 4
#rmin = 10**amin
#i = np.where(rhop < rmin)
#rhop[i] = rmin

## Check the grid size to apply to the plot
grid   = pc.read.grid()
x1, y1 = grid.x[3], grid.y[3]   # first corner coordinates
x2, y2 = grid.x[-4], grid.y[-4] # second corner coordinates

## Plot potential at midplane
plt.imshow(rhop,origin='lower',cmap='inferno',interpolation='none',
		extent=[x1,x2,y1,2])                                  # Plotting as image
plt.xlabel(r'$x$', fontsize=15)                               # x-axis label
plt.ylabel(r'$y$', fontsize=15)                               # y-axis label
bar = plt.colorbar()   # Colorbar
bar.set_label(r"potself", fontsize=15) # Colorbar label
plt.title("")                      # Plot title
plt.savefig("potential.pdf", dpi=300)                      # Save figure as a PDF file with set quality
plt.close()                                                   # Close figure after saving
