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
parray1 = np.empty(shape=(nstalk-1, header_len))
for i in range(0, nstalk-1):   # For every pstalk save
    num = '{0:04}'.format(i+1) # Format each output file call e.g. 1 -> '0001'
    path = os.path.join(".",'output/file_'+num+'.txt')
    output = np.loadtxt(path)  # Load the output row from a file
    parray1[i,:] = output[0,:] # First row elements are substituted into the larger array
## Save massive particle data to .csv
np.savetxt('parray1.csv', parray1, delimiter=',', header='it, index, xp, yp, vx, vy, rhopswarm, aps')

## Create second array to input particle 1 data from smaller grid
parray1b = np.empty(shape=(nstalk-1, header_len))
for i in range(0, nstalk-1):   # For every pstalk save
    num = '{0:04}'.format(i+1) # Format each output file call e.g. 1 -> '0001'
    path = os.path.join(".",'output1/file_'+num+'.txt')
    output = np.loadtxt(path)  # Load the output row from a file
    parray1b[i,:] = output[0,:] # First row elements are substituted into the larger array


## Create new array to input particle 2 data
parray2 = np.empty(shape=(nstalk-1, header_len))
for i in range(0, nstalk-1):   # For every pstalk save
    num = '{0:04}'.format(i) # Format each output file call e.g. 1 -> '0001'
    path = os.path.join(".",'output/file_'+num+'.txt')
    output = np.loadtxt(path)  # Load the output row from a file
    parray2[i,:] = output[1,:] # First row elements are substituted into the larger array
## Save near massless particle data to .csv
np.savetxt('parray2.csv', parray2, delimiter=',', header='it, index, xp, yp, vx, vy, rhopswarm, aps')

## Create second array to input particle 2 data from smaller grid
parray2b = np.empty(shape=(nstalk-1, header_len))
for i in range(0, nstalk-1):   # For every pstalk save
    num = '{0:04}'.format(i) # Format each output file call e.g. 1 -> '0001'
    path = os.path.join(".",'output1/file_'+num+'.txt')
    output = np.loadtxt(path)  # Load the output row from a file
    parray2b[i,:] = output[1,:] # First row elements are substituted into the larger array


### (2) CALCULATE ORBITAL PARAMETERS OF HYPERBOLIC ORBIT
## Determine distance, velocity, and angle at periapsis
vp_array = np.loadtxt("parray2.csv", delimiter=',')       # Load array of positions and velocities over time
vp = 0                                                    # Minimum 'magnitude' of velocity vector
for i in range(len(vp_array[:,0])):                       # For all time indices...
    line_m = np.sqrt(vp_array[i,4]**2 + vp_array[i,5]**2) # Velocity vector magnitude 
    if line_m > vp:                                       # If higher, then this could be the periapsis 
        vp = line_m                                       # Record new highest velocity
        vx = vp_array[i,4]                                # x velocity
        vy = vp_array[i,5]                                # y velocity
        x_coord = vp_array[i,2]                           # x distance
        y_coord = vp_array[i,3]                           # y distance
        rp = np.sqrt(x_coord**2 + y_coord**2)             # Radial distance at periapsis
angle = np.arctan2(y_coord, x_coord)                      # tan(theta)=y/x=opp/adj
print('Maximum velocity:',vp,'(Components: vx=',vx,'vy=',vy,')')
print('Closest distance:',rp,'(Components: x=',x_coord,'y=',y_coord,')')
print('Origin to periapsis angle:',angle,'radians')

## From quantities in code units, determine characteristics of hyperbolic trajectory
grid   = pc.read.grid(trim=True) # Read from grid data
rhopswarm   = parray1[-1,6]      # Mass density set for central sink particle
dx = grid.dx                     # Grid width in x dimension
print('Cell width is:',dx)

G = 1/(4*np.pi)           # Gravitational constant
M = (rhopswarm) * (dx)**2 # Mass of particle at origin; set by column density and grid cell area

mass_multiple = 1.00      # As a test, multiply mass
M = (rhopswarm) * (dx)**2 * mass_multiple # TEST MASS
mu = G * M                # Gravitational parameter

a1 = (2/rp-(vp**2)/mu)**(-1) # Semi-major axis in code length units; solve from vis-viva equation
e1 = 1 - rp/a1               # Eccentricity
printColor('Hyperbolic trajectory from periapsis:','HEADER')
print('Semi-major axis:',a1,'code length units')
print('Eccentricity:',e1)

## Alternatively, use initial conditions to determine orbital parameters of hyperbolic orbit
xp = parray2[0,2] # Initial x position
yp = parray2[0,3] # Initial y position
vx = parray2[0,4] # Initial x velocity
vy = parray2[0,5] # Initial y velocity

## Calculate orbital parameters from expected mass and initial conditions
ri   = np.sqrt(xp**2 + yp**2) # Initial radial distance
vi2  = vx**2 + vy**2          # Initial velocity^2 magnitude
a2   = (2/ri - vi2/mu)**(-1)  # Semi-major axis from rewritten vis-viva equation
vinf = np.sqrt(-mu/a2)        # Set r=infinity for velocity at infinity
h    = abs(vy*xp - vx*yp)     # Specific angular momentum at any point in trajectory
b    = h/vinf                 # Impact parameter or semi-minor axis
e2   = np.sqrt(1+b**2/a2**2)  # Eccentricity from semi-minor and semi-major axes
e2   = np.sqrt(1 + (2*vi2*h**2)/mu**2 - ri*vi2/mu) # Eccentricity solved in a different fashion
printColor('Hyperbolic trajectory from initial conditions:','HEADER')
print('Semi-major axis:',a2,'code length units')
print('Eccentricity:',e2)


### (3) SOLVE FOR POINTS ALONG TRAJECTORY
f_inf = np.arccos(-1/e1)
range = np.arange(-1.4*np.pi,-0.2*np.pi,0.001) # Range of true anomalies to cover in radians
peri_angle = angle                            # radians; periapsis angle from +x direction

## Line 1 to plot for analytic solution (periapsis case)
a = a1  # Semi-major axis calculated separately
e = e1  # Eccentricity calculated separately
r_array = (a*(1-e**2))/(1+e*np.cos(range-peri_angle))         # Array of radial distances
x_pos1, y_pos1 = r_array*np.cos(range), r_array*np.sin(range) # Polar coordinates to cartesian
## Line 2 to plot for analytic solution (initial condition case)
a = a2  # Semi-major axis calculated separately
e = e2  # Eccentricity calculated separately
r_array = (a*(1-e**2))/(1+e*np.cos(range-peri_angle))         # Array of radial distances
x_pos2, y_pos2 = r_array*np.cos(range), r_array*np.sin(range) # Polar coordinates to cartesian


### (4) PLOT TRAJECTORIES OVER TIME
## Check the grid size to apply to the plot
grid = pc.read.grid() # Read relevant information from grid data
lx = grid.Lx          # Length of x-axis side
ly = grid.Ly          # Length of y-axis side

## Plot location of particles throughout simulation
fig, ax = plt.subplots(1, 1)                                                                     # Figure and axes established
ax.plot(parray1[:,2], parray1[:,3], 'r.', label=r'Sink $aps=$'+str(parray1[0,7]), zorder=3)           # Plotting particle 1 position
min, max = np.max(parray2[:,0]), np.max(parray2[:,0])                                            # Minimum and maximum of time
#lines2 = colored_line(parray2[:,2], parray2[:,3], parray2[:,0], ax, linewidth=5, cmap='seismic') # Particle 2 position with time

ax.plot(parray2[:,2], parray2[:,3], label="Test Particle (12x12)")
ax.plot(parray2b[:,2], parray2b[:,3], label="Test Particle (4x4)",color='magenta')

## Plot analytical solution provided input
#ax.plot(x_pos1, y_pos1, color='green', label='From Periapsis', zorder=5)                 # Analytical expectation from periapsis position and velocity in sim 
#ax.plot(x_pos2, y_pos2, '-.',color='magenta', label='From Initial Conditions', zorder=5) # Analytical expectation from initial conditons

## Add line indicating periapsis position vector
ax.plot((0, x_coord), (0, y_coord), linewidth=2, linestyle='-', color='k', zorder=2)

## Add text indicating mass and initial velocity
props = dict(boxstyle='round', facecolor='white', alpha=0.7)
ax.text(0.44, .10, 'Expected Mass * '+"{:.2E}".format(mass_multiple),
        transform=ax.transAxes, fontsize=8, verticalalignment='top', bbox=props, zorder=6)
ax.text(0.44, .05, 'Initial $v$: $v_x$='+"{:.2E}".format(vx)+' $v_y$='+"{:.2E}".format(vy),
        transform=ax.transAxes, fontsize=8, verticalalignment='top', bbox=props, zorder=6)

## Figure settings
#cbar = fig.colorbar(lines2, ax=ax, label='Time (code units)') # Time colorbar
#cbar.ax.invert_yaxis()                               # Flip color axis for better comparison with plot (time starts on top)
ax.set_xlabel(r'$x$', fontsize=15)                   # x-axis label
ax.set_ylabel(r'$y$', fontsize=15)                   # y-axis label
plt.axis('scaled')                                   # Axes are scaled to match one another
ax.set_xlim(-2, 2)                             # Set x-axis limits, assuming centered evenly on origin
ax.set_ylim(-2, 2)                             # Set y-axis limits, assuming centered evenly on origin
ax.set_title(r'Circular Test Orbit of Massless Particle') # Plot title
plt.legend(loc='best')                               # Add legend
ax.tick_params(axis='both', which='minor', length=0) # Set tick parameters (0 length)
ax.grid(which='major', alpha=0.5)                    # Applying grid based on major ticks
plt.tight_layout()                                   # Remove overlapping and clipping
plt.savefig('pstalk-figure.pdf', dpi=300)            # Save figure as a PDF file with set quality
plt.close()                                          # Close figure after saving




