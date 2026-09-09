## Import critical libraries/modules
import pencil as pc             # Local module for handling Pencil Code outputs
import numpy as np              # Arrays
import matplotlib.pyplot as plt # Plotting

## Read in VAR files from Pencil Code output
#ff  = pc.read.var(trimall=True)  # Pulling full VAR file output
#pot = ff.potself[0,:,:]          # Self gravitational potential at z=0 (z,y,x)
var = pc.read.var(var_file='VAR0', trimall=True) # 
if hasattr(var, 'phi'):
    pot = var.phi
elif hasattr(var, 'potpointmass'):
    pot = var.potpointmass
else:
    raise ValueError("Potential field not found in this VAR file.")



## Check the grid size to apply to the plot
grid   = pc.read.grid()
x1, y1 = grid.x[3], grid.y[3]   # First corner coordinates
x2, y2 = grid.x[-4], grid.y[-4] # Second corner coordinates

## Plot 2D potential at midplane
plt.imshow(pot, origin='lower', cmap='inferno', interpolation='none', extent=[x1,x2,y1,y2]) # Plotting as image
plt.xlabel(r'$x$', fontsize=15)                 # x-axis label
plt.ylabel(r'$y$', fontsize=15)                 # y-axis label
bar = plt.colorbar()                            # Colorbar
bar.set_label(r"Self Potential", fontsize=15)   # Colorbar label
#plt.title("")                                  # Plot title
plt.savefig("potential-2dfigure.pdf", dpi=300)  # Save figure as a PDF file with set quality
plt.close()                                     # Close figure after saving

## Plot 1D potential line through point sources
rangepot = np.linspace(-4,4,256)
pot1d = ff.potself[0,0,:]
print(len(pot1d))
plt.plot(rangepot, pot1d, label='$256^2$ case')
print('Minimum:',np.min(potself))

plt.xlabel(r'$x$', fontsize=15)                             # x-axis label
plt.ylabel(r'Self Potential', fontsize=15)                  # y-axis label
#bar = plt.colorbar()   # Colorbar
#bar.set_label(r"potself", fontsize=15) # Colorbar label
#plt.ylim(-4,4)                                 # Potential limits (zoom in to view)
plt.legend()

plt.title("Self Potential Slice for $y=0$, $z=0$")   # Plot title
plt.savefig("potential-1dfigure.pdf", dpi=300)       # Save figure as a PDF file with set quality
plt.close()                                          # Close figure after saving
