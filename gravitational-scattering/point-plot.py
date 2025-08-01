## Import critical libraries/modules
import pencil as pc             # Local module for handling Pencil Code outputs
import numpy as np              # Arrays
import matplotlib.pyplot as plt # Plotting

## Read in timeseries from Pencil Code output
ts       = pc.read.ts(datadir='data') # Pulling timeseries data
xp1, yp1 = ts.xpm, ts.ypm             # Setting x and y variables

## Check the grid size to apply to the plot
grid   = pc.read.grid()
x1, y1 = grid.x[3], grid.y[3]   # first corner coordinates
x2, y2 = grid.x[-4], grid.y[-4] # second corner coordinates

## Plot location of point masses at end of simulation
plt.plot(xp1, yp1, 'r+')                   # Plotting mean particle position
plt.xlabel(r'$x$', fontsize=15)            # x-axis label
plt.ylabel(r'$y$', fontsize=15)            # y-axis label
plt.axis('scaled')                         # Axes are scaled to match one another
plt.xlim(x1, x2)                           # Set x-axis limits
plt.ylim(y1, y2)                           # Set y-axis limits
plt.title("Shearing Sheet - Scattering")   # Plot title
plt.tight_layout()                         # Remove overlapping and clipping
plt.savefig("point-plot.pdf", dpi=300)     # Save figure as a PDF file with set quality
plt.close()                                # Close figure after saving

