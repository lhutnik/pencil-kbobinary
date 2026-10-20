## Import critical libraries/modules
import pencil as pc             # Local module for handling Pencil Code outputs
import numpy as np              # Arrays
import matplotlib.pyplot as plt # Plotting

## Read in VAR files from Pencil Code output
ff   = pc.read.var(trimall=True) # Pulling full VAR file output
rhop = ff.potself[0,:,:]         # Particle density at z=0 (z,y,x)

## Read in previous VAR files from Pencil Code output
ff2 = pc.read.var(trimall=True, datadir='~/pencil-kbobinary/gravitational-scattering/data')


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
bar.set_label(r"Self Potential", fontsize=15) # Colorbar label
#plt.title("")                      # Plot title
plt.savefig("potential.pdf", dpi=300)                      # Save figure as a PDF file with set quality
plt.close()                                                   # Close figure after saving

## Setup 1D plot
rangepot= np.linspace(-6+1/32,6-1/32, 192)
potself = ff.potself[7,95,:]              # Roughly middle of z, y ranges (midplane)
rangepot2= np.linspace(-2+1/32,2-1/32,64)
potself2 = ff2.potself[7,31,:]

## Plot potential at midplane
#plt.imshow(rhop,origin='lower',cmap='inferno',interpolation='none',
#		extent=[x1,x2,y1,2])                                  # Plotting as image
M = 1
G = 1
xrange = np.linspace(-2.,2.,500)
pointpot = -G*M/np.abs(xrange)

plt.plot(xrange, pointpot, linestyle='--', label='Point Source')
plt.plot(rangepot, potself, label='12x12 Grid')
plt.plot(rangepot2, potself2, label='4x4 Grid', color='magenta')
print('Minimum:',np.min(potself))

plt.xlabel(r'$x$', fontsize=15)                             # x-axis label
plt.ylabel(r'Self Potential', fontsize=15)                  # y-axis label
#bar = plt.colorbar()   # Colorbar
#bar.set_label(r"potself", fontsize=15) # Colorbar label
plt.ylim(-4,4)                                 # Potential limits (zoom in to view)
plt.legend()

plt.title("Self Potential Slice for $y=0$, $z=0$")          # Plot title
plt.savefig("potential2.pdf", dpi=300)                      # Save figure as a PDF file with set quality
plt.close()                                                 # Close figure after saving
