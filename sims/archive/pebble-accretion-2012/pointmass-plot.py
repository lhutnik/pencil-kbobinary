## Import critical libraries 
import numpy as np                # Pi constant
import matplotlib.pyplot as plt  # Plotting
import astropy.constants as const # Physical constants
import pencil as pc

## Set constants
pi     = np.pi
AU     = const.au.cgs.value      # Astronomical unit in cm
#Mearth = const.M_earth.cgs.value # Mass of Earth in g
Msun   = const.M_sun.cgs.value   # Mass of Sun in g
G      = const.G.cgs.value       # Gravitational constant
amu    = const.u.cgs.value       # Atomic mass unit
kb     = const.k_B.cgs.value     # Boltzmann constant
Rgas   = const.R.cgs.value       # Gas constant
yr     = 3.1558149504e7          # Year in seconds


### IMPORT POINT MASS PARTICLE DATA
#sim = pc.get_sim()
ts = pc.read.ts()
xq1, xq2 = ts.xq1, ts.xq2
yq1, yq2 = ts.yq1, ts.yq2
vxq1, vxq2 = ts.vxq1, ts.vxq2


### PLOT TRAJECTORIES OVER TIME
## Check the grid size to apply to the plot
grid = pc.read.grid() # Read relevant information from grid data
lx = grid.Lx          # Length of x-axis side
ly = grid.Ly          # Length of y-axis side

## Plot location of particles throughout simulation
fig, ax = plt.subplots(1, 1)                                                                     # Figure and axes established
#ax.plot(parray1[:,2], parray1[:,3], 'r.', label=r'Sink $aps=$'+str(parray1[0,7]), zorder=3)           # Plotting particle 1 position
#min, max = np.max(parray2[:,0]), np.max(parray2[:,0])                                            # Minimum and maximum of time
#lines2 = colored_line(parray2[:,2], parray2[:,3], parray2[:,0], ax, linewidth=5, cmap='seismic') # Particle 2 position with time

ax.plot(xq1, yq1, label="$m=1$")
ax.plot(xq2, yq2, label="$m=1e-30$", color='magenta')

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
ax.set_xlim(-2, 2)                             # Set x-axis limits, assuming centered evenly on origin
ax.set_ylim(-2, 2)                             # Set y-axis limits, assuming centered evenly on origin
ax.set_title(r'Circular Test Orbit of Near Massless Point Mass') # Plot title
plt.legend(loc='best')                               # Add legend
ax.tick_params(axis='both', which='minor', length=0) # Set tick parameters (0 length)
ax.grid(which='major', alpha=0.5)                    # Applying grid based on major ticks
plt.tight_layout()                                   # Remove overlapping and clipping
plt.savefig('pointmass-figure.pdf', dpi=300)            # Save figure as a PDF file with set quality
plt.close()                                          # Close figure after saving
