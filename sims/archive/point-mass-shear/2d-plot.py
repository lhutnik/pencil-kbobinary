## Import critical libraries/modules
import pencil as pc             # Local module for handling Pencil Code outputs
import numpy as np              # Arrays
import matplotlib.pyplot as plt # Plotting

## Read in VAR files from Pencil Code output
ff   = pc.read.var(trimall=True) # Pulling full VAR file output
rho  = ff.rho[0,:,:]            # Particle density at z=0 (z,y,x)
#amax = np.log10(rhop.max())
#amin = amax - 4
#rmin = 10**amin
#i = np.where(rhop < rmin)
#rhop[i] = rmin

## Check the grid size to apply to the plot
grid = pc.read.grid() # Read relevant information from grid data
lx = grid.Lx          # Length of x-axis side
ly = grid.Ly          # Length of y-axis side

## Plot 2d density field at midplane
plt.imshow(rho,origin='lower',cmap='inferno',interpolation='none',
		extent=[-lx/2,lx/2,-ly/2,ly/2])               # Plotting as image
plt.xlabel(r'$x$', fontsize=15)                               # x-axis label
plt.ylabel(r'$y$', fontsize=15)                               # y-axis label
bar = plt.colorbar(label=r"$\log_{\rm 10}(\rho/\rho_0)$")   # Colorbar
bar.set_label(r"$\log_{\rm 10}(\rho/\rho_0)$", fontsize=15) # Colorbar label
plt.title("Shearing Sheet - Density")                         # Plot title
plt.savefig("density-plot.pdf", dpi=300)                      # Save figure as a PDF file with set quality
plt.close()                                                   # Close figure after saving
