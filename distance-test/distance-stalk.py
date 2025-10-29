## Import critical libraries/modules
import os                        # File/directory checks
import pencil as pc              # Local module for handling Pencil Code outputs
import numpy as np               # Arrays
import matplotlib.pyplot as plt  # Plotting
from helpers import colored_line # Helper function for plotting line with gradient
from helpers import printColor   # Add color to certain printed lines for readability


### (1) IMPORT SIMULATION DATA
## Import time information
timepath = os.path.join('.','data/tstalk.dat')
timelen = np.loadtxt(timepath)
nstalk = int(timelen[1]) # Second element of timelen provides number of pstalk timesteps recorded 
print('The number of pstalk files is:',nstalk)

## Import number of columns (variables) in output files
headerpath = os.path.join('.','output/file_0001.txt') # Most likely, the first file will exist
header  = np.loadtxt(headerpath)
header_len = len(header[0,:])  # Number of columns in each pstalk timestep
print('The number of columns is:',header_len)

## Create new array to input particle 1 data
parray1 = np.empty(shape=(627-1, header_len))
for i in range(0, 627-1):   # For every pstalk save
    num = '{0:04}'.format(i+1) # Format each output file call e.g. 1 -> '0001'
    path = os.path.join(".",'output1/file_'+num+'.txt')
    output = np.loadtxt(path)  # Load the output row from a file
    parray1[i,:] = output[0,:] # First row elements are substituted into the larger array
## Save massive particle data to .csv
np.savetxt('parray1.csv', parray1, delimiter=',', header='it, index, xp, yp, vx, vy, rhopswarm, aps')

## Create new array to input particle 2 data for r=1
parray21 = np.empty(shape=(627-1, header_len))
for i in range(0, 627-1):   # For every pstalk save
    num = '{0:04}'.format(i) # Format each output file call e.g. 1 -> '0001'
    path = os.path.join(".",'output1/file_'+num+'.txt')
    output = np.loadtxt(path)  # Load the output row from a file
    parray21[i,:] = output[1,:] # First row elements are substituted into the larger array
## Save near massless particle data to .csv
np.savetxt('parray2-1.csv', parray21, delimiter=',', header='it, index, xp, yp, vx, vy, rhopswarm, aps')

## Create new array to input particle 2 data for r=2
parray22 = np.empty(shape=(1777-1, header_len))
for i in range(0, 1777-1):   # For every pstalk save
    num = '{0:04}'.format(i) # Format each output file call e.g. 1 -> '0001'
    path = os.path.join(".",'output2/file_'+num+'.txt')
    output = np.loadtxt(path)  # Load the output row from a file
    parray22[i,:] = output[1,:] # First row elements are substituted into the larger array
## Save near massless particle data to .csv
np.savetxt('parray2-2.csv', parray22, delimiter=',', header='it, index, xp, yp, vx, vy, rhopswarm, aps')

## Create new array to input particle 2 data for r=2
parray23 = np.empty(shape=(2500-1, header_len))
for i in range(0, 2500-1):   # For every pstalk save
    num = '{0:04}'.format(i) # Format each output file call e.g. 1 -> '0001'
    path = os.path.join(".",'output3/file_'+num+'.txt')
    output = np.loadtxt(path)  # Load the output row from a file
    parray23[i,:] = output[1,:] # First row elements are substituted into the larger array
## Save near massless particle data to .csv
np.savetxt('parray2-3.csv', parray23, delimiter=',', header='it, index, xp, yp, vx, vy, rhopswarm, aps')


### (2) PLOT TRAJECTORIES OVER TIME
## Check the grid size to apply to the plot
grid = pc.read.grid() # Read relevant information from grid data
lx = grid.Lx          # Length of x-axis side
ly = grid.Ly          # Length of y-axis side

## Plot location of particles throughout simulation
fig, ax = plt.subplots(1, 1)                                                                     # Figure and axes established
ax.plot(parray1[:,2], parray1[:,3], 'r.', label=r'Sink $aps=$'+str(parray1[0,7]), zorder=3)      # Plotting particle 1 position
#min, max = np.max(parray2[:,0]), np.max(parray2[:,0])                                            # Minimum and maximum of time
#lines2 = colored_line(parray2[:,2], parray2[:,3], parray2[:,0], ax, linewidth=5, cmap='seismic') # Particle 2 position with time

ax.plot(parray21[:,2], parray21[:,3], label="Test Particle $r=1$",linestyle='solid',alpha=0.7,zorder=6)
ax.plot(parray22[:,2], parray22[:,3], label="Test Particle $r=2$",linestyle='solid',alpha=0.7,zorder=5)
ax.plot(parray23[:,2], parray23[:,3], label="Test Particle $r=3$",linestyle='solid',alpha=0.7,zorder=4)

## Point for starting positions of particles
ax.scatter([1,2,3],[0,0,0], c='k')

## Plot analytical solution provided input
#ax.plot(x_pos1, y_pos1, color='green', label='From Periapsis', zorder=5)                 # Analytical expectation from periapsis position and velocity in sim 
#ax.plot(x_pos2, y_pos2, '-.',color='magenta', label='From Initial Conditions', zorder=5) # Analytical expectation from initial conditons

## Add line indicating periapsis position vector
#ax.plot((0, x_coord), (0, y_coord), linewidth=2, linestyle='-', color='k', zorder=2)

## Add text indicating mass and initial velocity
#props = dict(boxstyle='round', facecolor='white', alpha=0.7)
#ax.text(0.44, .10, 'Expected Mass * '+"{:.2E}".format(mass_multiple),
#        transform=ax.transAxes, fontsize=8, verticalalignment='top', bbox=props, zorder=6)
#ax.text(0.44, .05, 'Initial $v$: $v_x$='+"{:.2E}".format(vx)+' $v_y$='+"{:.2E}".format(vy),
#        transform=ax.transAxes, fontsize=8, verticalalignment='top', bbox=props, zorder=6)

## Figure settings
#cbar = fig.colorbar(lines2, ax=ax, label='Time (code units)') # Time colorbar
#cbar.ax.invert_yaxis()                               # Flip color axis for better comparison with plot (time starts on top)
ax.set_xlabel(r'$x$', fontsize=15)                   # x-axis label
ax.set_ylabel(r'$y$', fontsize=15)                   # y-axis label
plt.axis('scaled')                                   # Axes are scaled to match one another
ax.set_xlim(-4, 4)                             # Set x-axis limits, assuming centered evenly on origin
ax.set_ylim(-4, 4)                             # Set y-axis limits, assuming centered evenly on origin
ax.set_title(r'Testing Circular Orbit at Different $r$') # Plot title
plt.legend(loc='best')                               # Add legend
ax.tick_params(axis='both', which='minor', length=0) # Set tick parameters (0 length)
ax.grid(which='major', alpha=0.5)                    # Applying grid based on major ticks
plt.tight_layout()                                   # Remove overlapping and clipping
plt.savefig('distance-figure.pdf', dpi=300)            # Save figure as a PDF file with set quality
plt.close()                                          # Close figure after saving




