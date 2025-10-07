## Import critical libraries/modules
import os                        # File/directory checks
import pencil as pc              # Local module for handling Pencil Code outputs
import numpy as np               # Arrays 
import matplotlib.pyplot as plt  # Plotting
from helpers import colored_line # Helper functions for plotting and printing

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
parray1 = np.empty(shape=(nstalk-1, header_len))
for i in range(0, nstalk-1):   # For every pstalk save
    num = '{0:04}'.format(i+1) # Format each output file call e.g. 1 -> '0001'
    path = os.path.join(".",'output/file_'+num+'.txt')
    output = np.loadtxt(path)  # Load the output row from a file
    parray1[i,:] = output[0,:] # First row elements are substituted into the larger array
## Save near massive particle data to .csv
np.savetxt('parray1.csv', parray1, delimiter=',', header='it, index, xp, yp, vx, vy, rhopswarm, aps')

## Create new array to input particle 2 data
parray2 = np.empty(shape=(nstalk-1, header_len))
for i in range(0, nstalk-1):   # For every pstalk save
    num = '{0:04}'.format(i+1) # Format each output file call e.g. 1 -> '0001'
    path = os.path.join(".",'output/file_'+num+'.txt')
    output = np.loadtxt(path)  # Load the output row from a file
    parray2[i,:] = output[1,:] # First row elements are substituted into the larger array
## Save near massless particle data to .csv
np.savetxt('parray2.csv', parray2, delimiter=',', header='it, index, xp, yp, vx, vy, rhopswarm, aps')


### (2) CALCULATE ORBITAL PARAMETERS OF HYPERBOLIC ORBIT
## Determine distance, velocity, and angle at periapsis
vp_array = np.loadtxt("parray2.csv", delimiter=',')       # Load array of positions and velocities over time
vp = 0                                                    # Minimum 'magnitude' of position vector
for i in range(len(vp_array[:,0])):                       # For all time indices...
    line_m = np.sqrt(vp_array[i,4]**2 + vp_array[i,5]**2) # Velocity vector magnitude 
    if line_m > vp:                                       # If higher, then this could be the periapsis 
        vp = line_m                                       # Record new highest velocity
        x_coord = vp_array[i,2]                           # x distance
        y_coord = vp_array[i,3]                           # y distance
        rp = np.sqrt(x_coord**2 + y_coord**2)             # Radial distance at periapsis
angle = np.arctan2(y_coord, x_coord)                      # tan(theta)=y/x=opp/adj
print('Maximum velocity:',vp)
print('Closest distance:',rp,'(Components: x=',x_coord,'y=',y_coord,')')
print('Origin to periapsis angle:',angle,'radians')

## From quantities in code units, determine characteristics of hyperbolic trajectory
rhopswarm   = 64.0 # Mass density set for sink particles
cell_side   = 256  # Grid points per dimension
side_length = 4.0  # Length of one side

G = 1/(4*np.pi)                              # Gravitational constant 
M = (rhopswarm) * (side_length/cell_side)**2 # Mass of particle at origin; set by column density and grid cell area
a1 = ((2/rp)-(vp**2)/(G*M))**(-1)            # Semi-major axis in code length units; solve from vis-viva equation
e1 = 1 - rp/a                                # Eccentricity
print('A representative hyperbolic trajectory has:')
print('Semi-major axis:',a1,'code length units')
print('Eccentricity:',e1)

## Alternatively, use initial conditions to determine orbital parameters of hyperbolic orbit
xp = 1.0   # Initial x position
yp = 1.99  # Initial y position 
vx = 0.0   # Initial x velocity
vy = -0.05 # Initial y velocity

## Calculate orbital parameters from expected mass and initial conditions 
M = (rhopswarm) * (side_length/cell_side)**2 # Column density is being used for 2D case
mu   = G * M                  # Gravitational parameter
r    = np.sqrt(xp**2 + yp**2) # Initial radial distance
v2   = vx**2 + vy**2          # Initial velocity^2 magnitude
a    = 1/(2/r - v2/mu)        # Semi-major axis from rewritten vis-viva equation
vinf = np.sqrt(-mu/a)         # Set r=infinity for velocity at infinity
h    = abs(vy*xp - vx*yp)     # Specific angular momentum at any point in trajectory
b    = h/vinf                 # Impact parameter or semi-minor axis
e    = np.sqrt(1+b**2/a**2)   # Eccentricity from semi-minor and semi-major axes
print('A representative hyperbolic trajectory for a 2D mass:')
print('Semi-major axis:',a,'code length units')
print('Eccentricity:',e)


### (3) SOLVE FOR POINTS ALONG TRAJECTORY
## Line 1 to plot for analytic solution
range = np.arange(-np.pi*0.8,np.pi*0.4,0.01) # Range of true anomalies to cover in radians
peri_angle = angle                           # radians; periapsis angle from +x direction
a = a1                                       # Semi-major axis calculated separately
e = e1                                       # Eccentricity calculated separately
r_array = (a*(1-e**2))/(1+e*np.cos(range-peri_angle))         # Array of radial distances
x_pos1, y_pos1 = r_array*np.cos(range), r_array*np.sin(range) # Polar coordinates to cartesian

## Line 2 to plot for analytic solution
range = np.arange(-np.pi*0.8,np.pi*0.4,0.01) # Range of true anomalies to cover in radians
peri_angle = angle                           # radians; periapsis angle from +x direction
a = a2                                       # Semi-major axis calculated separately
e = e2                                       # Eccentricity calculated separately
r_array = (a*(1-e**2))/(1+e*np.cos(range-peri_angle))         # Array of radial distances
x_pos2, y_pos2 = r_array*np.cos(range), r_array*np.sin(range) # Polar coordinates to cartesian


### (4) PLOT TRAJECTORIES OVER TIME
## Check the grid size to apply to the plot
grid   = pc.read.grid()         # Read relevant information from grid data
x1, y1 = grid.x[3], grid.y[3]   # First corner coordinates
x2, y2 = grid.x[-4], grid.y[-4] # Second corner coordinates

## Plot location of particles throughout simulation
fig, ax = plt.subplots(1, 1)                                                                     # Figure and axes established
ax.plot(parray1[:,2], parray1[:,3], 'r.', label=r'$aps=$'+str(parray1[0,5]))                     # Plotting particle 1 position
min, max = np.max(parray2[:,0]), np.max(parray2[:,0])                                            # Minimum and maximum of time
lines2 = colored_line(parray2[:,2], parray2[:,3], parray2[:,0], ax, linewidth=5, cmap='seismic') # Particle 2 position with time

## Plot analytical solution provided input
ax.plot(x_pos1, y_pos1, color='green', label='From Periapsis', zorder=10)                 # Analytical expectation from periapsis position and velocity in sim 
ax.plot(x_pos2, y_pos2, '-.',color='magenta', label='From Initial Conditions', zorder=10) # Analytical expectation from initial conditons

## Figure settings
cbar = fig.colorbar(lines2, ax=ax, label='Time (code units)') # Time colorbar
cbar.ax.invert_yaxis()                               # Flip color axis for better comparison with plot (time starts on top)
ax.set_xlabel(r'$x$', fontsize=15)                   # x-axis label
ax.set_ylabel(r'$y$', fontsize=15)                   # y-axis label
plt.axis('scaled')                                   # Axes are scaled to match one another
ax.set_xlim(x1, x2)                                  # Set x-axis limits
ax.set_ylim(y1, y2)                                  # Set y-axis limits
ax.set_title(r'Gravitational Scattering Trajectory') # Plot title
plt.legend(loc='upper left')                         # Add legend
ax.tick_params(axis='both', which='minor', length=0) # Set tick parameters (0 length)
ax.grid(which='major', alpha=0.5)                    # Applying grid based on major ticks
plt.tight_layout()                                   # Remove overlapping and clipping
plt.savefig('pstalk-figure.pdf', dpi=300)            # Save figure as a PDF file with set quality
plt.close()                                          # Close figure after saving




