## Import critical libraries/modules
import os                       # File/directory checks
#import pencil as pc             # Local module for handling Pencil Code outputs
import numpy as np              # Arrays
import matplotlib.pyplot as plt # Plotting

from matplotlib.collections import LineCollection
from matplotlib.ticker import AutoMinorLocator
from helpers import colored_line


## Import time information
timepath = os.path.join('.','data/tstalk.dat')
timelen = np.loadtxt(timepath)
nstalk  = int(timelen[1])      # Second element of timelen provides number of pstalk timesteps recorded
print('The number of pstalk files is:',nstalk)
## Import number of columns (variables) in output files
headerpath = os.path.join('.','output/file_0001.txt') # Most likely, the first file will exist
header  = np.loadtxt(headerpath)
header_len = len(header[0,:])  # Number of columns in each pstalk timestep
print('The number of columns is:',header_len)

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
x1, y1 = grid.x[3], grid.y[3]   # First corner coordinates
x2, y2 = grid.x[-4], grid.y[-4] # Second corner coordinates

## Line 1 to plot for analytic solution
range = np.arange(-np.pi,np.pi,0.01)    # Range of true anomalies to cover in radians
peri_angle = -0.704527                  # radians; periapsis angle from +x direction
a = -0.4455599810                       # Semi-major axis calculated separately
e = 2.16751617                          # Eccentricity calculated separately
r_array = (a*(1-e**2))/(1+e*np.cos(range-peri_angle))         # Array of radial distances
x_pos1, y_pos1 = r_array*np.cos(range), r_array*np.sin(range) # Polar coordinates to cartesian
#x_pos1 = np.cos(-peri_angle)*(x_pos1) + np.sin(-peri_angle)*(y_pos1)
#y_pos1 = -np.sin(-peri_angle)*(x_pos1) + np.cos(-peri_angle)*(y_pos1)

## Line 2 to plot for analytic solution
range = np.arange(-np.pi,np.pi,0.01)   # Range of true anomalies to cover in radians
peri_angle = -0.704527                 # radians; periapsis angle from +x direction
a = -0.89879412                        # Semi-major axis calculated separately
e = 1.7991716794                       # Eccentricity calculated separately
r_array = (a*(1-e**2))/(1+e*np.cos(range-peri_angle))         # Array of radial distances
x_pos2, y_pos2 = r_array*np.cos(range), r_array*np.sin(range) # Polar coordinates to cartesian
#x_pos2 = np.cos(-peri_angle)*(x_pos2) + np.sin(-peri_angle)*(y_pos2)
#y_pos2 = -np.sin(-peri_angle)*(x_pos2) + np.cos(-peri_angle)*(y_pos2)

# Provide a rotation transformation of our initial x and y coordinates 

## Plot location of particles throughout simulation
fig, ax = plt.subplots(1, 1)                                                   # Figure and axes established
ax.plot(parray1[:,2], parray1[:,3], 'r.', label=r'$aps=$'+str(parray1[0,5]))   # Plotting particle 1 position
min, max = np.max(parray2[:,0]), np.max(parray2[:,0])                          # Minimum and maximum of time
lines2 = colored_line(parray2[:,2], parray2[:,3], parray2[:,0], ax, linewidth=5, cmap='seismic') # Particle 2 position with time

## Plot analytical solution provided input
ax.plot(x_pos1, y_pos1, color='o', label='From Periapsis', zorder=10) # Solution directly from numerical results 
ax.plot(x_pos2, y_pos2, '-.',color='m', label='From Initial Conditions', zorder=10) # Analytical expectation based on periapsis from simulation

## Figure settings
cbar = fig.colorbar(lines2, ax=ax, label='Time (code units)') # Time colorbar
cbar.ax.invert_yaxis() # Flip color axis for better comparison with plot (time starts on top)
ax.set_xlabel(r'$x$', fontsize=15)           # x-axis label
ax.set_ylabel(r'$y$', fontsize=15)           # y-axis label
plt.axis('scaled')                           # Axes are scaled to match one another
ax.set_xlim(x1, x2)                          # Set x-axis limits
ax.set_ylim(y1, y2)                          # Set y-axis limits
ax.set_title(r'Gravitational Scattering Trajectory') # Plot title
plt.legend(loc='upper left')                 # Add legend
grid_scale = 4/32
ax.xaxis.set_minor_locator(AutoMinorLocator(grid_scale)) # Ticks for x axis
ax.yaxis.set_minor_locator(AutoMinorLocator(grid_scale)) # Ticks for y axis
ax.tick_params(axis='both', which='minor', length=0)   # Set tick parameters (0 length)
ax.grid(which='major', alpha=0.5)   # Applying grid based on minor ticks
plt.tight_layout()                           # Remove overlapping and clipping
plt.savefig('pstalk-figure.pdf', dpi=300)    # Save figure as a PDF file with set quality
plt.close()                                  # Close figure after saving

