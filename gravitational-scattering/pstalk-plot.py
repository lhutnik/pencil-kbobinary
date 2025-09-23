## Import critical libraries/modules
import os                       # File/directory checks
import pencil as pc             # Local module for handling Pencil Code outputs
import numpy as np              # Arrays
import matplotlib.pyplot as plt # Plotting

## Import time and pstalk information
timepath = os.path.join('.','data/tstalk.dat')
timelen = np.loadtxt(timepath)
#print(timelen)
nstalk  = int(timelen[1])      # Second element of timelen provides number of pstalk timesteps recorded
headerpath = os.path.join('.','output/file_0001.txt')
header  = np.loadtxt(headerpath)
header_len = len(header[0,:])  # Number of columns in each pstalk timestep
print(header_len)

## Create new array to input particle 1 data
parray1 = np.empty(shape=(nstalk-1, header_len))
for i in range(0, nstalk-1):   # For every pstalk save
    num = '{0:04}'.format(i+1) # Format each output file call e.g. 1 -> '0001'
    path = os.path.join(".",'output/file_'+num+'.txt')
    output = np.loadtxt(path)  # Load the output row from a file
    parray1[i,:] = output[0,:] # First row elements are substituted into the larger array
## Create new array to input particle 2 data
parray2 = np.empty(shape=(nstalk-1, header_len))
for i in range(0, nstalk-1):   # For every pstalk save
    num = '{0:04}'.format(i+1) # Format each output file call e.g. 1 -> '0001'
    path = os.path.join(".",'output/file_'+num+'.txt')
    output = np.loadtxt(path)  # Load the output row from a file
    parray2[i,:] = output[1,:] # First row elements are substituted into the larger array

## Check the grid size to apply to the plot
grid   = pc.read.grid()         # Read relevant information from grid data
x1, y1 = grid.x[3], grid.y[3]   # first corner coordinates
x2, y2 = grid.x[-4], grid.y[-4] # second corner coordinates

## Plot location of particles throughout simulation
fig, ax = plt.subplots(1, 1)                 # Figure and axes established
ax.plot(parray1[:,1], parray1[:,2], 'r.', label='aps='+str(parray1[0,4]))   # Plotting particle 1 position
ax.plot(parray2[:,1], parray2[:,2], 'b.', label='aps='+str(parray2[0,4]))   # Plotting particle 2 position
ax.set_xlabel(r'$x$', fontsize=15)           # x-axis label
ax.set_ylabel(r'$y$', fontsize=15)           # y-axis label
plt.axis('scaled')                           # Axes are scaled to match one another
ax.set_xlim(x1, x2)                          # Set x-axis limits
ax.set_ylim(y1, y2)                          # Set y-axis limits
#ax.set_title(r'$t={}$'.format(parray1[-1,0])) # Plot title
plt.legend()                                 # Add legend
plt.tight_layout()                           # Remove overlapping and clipping
plt.savefig('pstalk-figure.pdf', dpi=300)       # Save figure as a PDF file with set quality
plt.close()                                  # Close figure after saving

