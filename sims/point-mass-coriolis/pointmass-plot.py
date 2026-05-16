## Import critical libraries 
import numpy as np                # Pi constant
import matplotlib.pyplot as plt  # Plotting
import astropy.constants as const # Physical constants
import pencil as pc
import pencil_old as oldpc

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


### IMPORT SIMULATION DATA
#sim = pc.get_sim()
## Timeseries data
ts = pc.read.timeseries.ts()
t = ts.t
xq1, xq2 = np.array(ts.xq1), np.array(ts.xq2)     # x positions
yq1, yq2 = np.array(ts.yq1), np.array(ts.yq2)     # y positions
vxq1, vxq2 = np.array(ts.vxq1), np.array(ts.vxq2) # x velocities 
vyq1, vyq2 = np.array(ts.vyq1), np.array(ts.vyq2) # y velocities 

## Grid data 
grid = pc.read.grid() # Check the grid size to apply to the plot
lx = grid.Lx          # Length of x-axis side
ly = grid.Ly          # Length of y-axis side

## Point mass specific data
qvar = oldpc.read_qvar()
imass = qvar.mass
m1 = imass[0] 
m2 = imass[1]


### POSITION VS. TIME
fig, ax = plt.subplots(1, 1) # Initiate figure and axis

## Plot location of particles throughout simulation
ax.plot(xq1, yq1, label="Primary", alpha=0.5, color='blue')
ax.plot(xq2, yq2, label="Secondary", alpha=0.5, color='orange')

## Figure settings
ax.set_xlabel(r'$x$', fontsize=15)                   # x-axis label
ax.set_ylabel(r'$y$', fontsize=15)                   # y-axis label
ax.set_title(r'Position vs. Time of $M_1$={}, $M_2$={} Binary'.format(m1, m2)) # Plot title
plt.axis('scaled')                                   # Axes are scaled to match one another
ax.set_xlim(-2, 2)                                   # Set x-axis limits, assuming centered evenly on origin
ax.set_ylim(-2, 2)                                   # Set y-axis limits, assuming centered evenly on origin
plt.legend(loc='best')                               # Add legend
ax.tick_params(axis='both', which='minor', length=0) # Set tick parameters (0 length)
ax.grid(which='major', alpha=0.5)                    # Applying grid based on major ticks
plt.tight_layout()                                   # Remove overlapping and clipping
plt.savefig('pointmass-figure.pdf', dpi=300)         # Save figure as a PDF file with set quality
plt.close()                                          # Close figure after saving


### TOTAL ENERGY VS. TIME PLOT
fig, ax = plt.subplots(1, 1) # Figure and axes established

## Plot energy over time for the system of particles
G = 1 # Gravitational constant
r = np.sqrt((xq1 - xq2)**2 + (yq1 - yq2)**2)
U = - (G)*(m1*m2)/r # System gravitational potential energy
T = 0.5*m1*(vxq1**2 + vyq1**2) + 0.5*m2*(vxq2**2 + vyq2**2) # Sum of kinetic energies
E = U + T # Energy is the sum of kinetic and potential energy
ax.plot(t, E, color='purple')

## Figure settings
ax.set_ylabel(r'$E=T+U$', fontsize=15)                   # x-axis label
ax.set_xlabel(r'$t$', fontsize=15)                   # y-axis label
#ax.set_xlim(-2, 2)                             # Set x-axis limits, assuming centered evenly on origin
#ax.set_ylim(-2, 2)                             # Set y-axis limits, assuming centered evenly on origin
ax.set_title(r'Energy vs. Time of $M_1$={}, $M_2$={} Binary'.format(m1, m2)) # Plot title
ax.tick_params(axis='both', which='minor', length=0) # Set tick parameters (0 length)
ax.grid(which='major', alpha=0.5)                    # Applying grid based on major ticks
plt.tight_layout()                                   # Remove overlapping and clipping
plt.savefig('energy-figure.pdf', dpi=300)            # Save figure as a PDF file with set quality
plt.close()                                          # Close figure after saving


### TOTAL ANGULAR MOMENTUM VS. TIME PLOT
fig, ax = plt.subplots(1, 1) # Figure and axes established

## Plot energy over time for the system of particles
G = 1 # Gravitational constant
r = np.sqrt((xq1 - xq2)**2 + (yq1 - yq2)**2)
U = - (G)*(m1*m2)/r # System gravitational potential energy
T = 0.5*m1*(vxq1**2 + vyq1**2) + 0.5*m2*(vxq2**2 + vyq2**2) # Sum of kinetic energies
E = U + T # Energy is the sum of kinetic and potential energy
ax.plot(t, E, color='purple')

## Figure settings
ax.set_ylabel(r'$E=T+U$', fontsize=15)                   # x-axis label
ax.set_xlabel(r'$t$', fontsize=15)                   # y-axis label
#ax.set_xlim(-2, 2)                             # Set x-axis limits, assuming centered evenly on origin
#ax.set_ylim(-2, 2)                             # Set y-axis limits, assuming centered evenly on origin
ax.set_title(r'Energy vs. Time of $M_1$={}, $M_2$={} Binary'.format(m1, m2)) # Plot title
ax.tick_params(axis='both', which='minor', length=0) # Set tick parameters (0 length)
ax.grid(which='major', alpha=0.5)                    # Applying grid based on major ticks
plt.tight_layout()                                   # Remove overlapping and clipping
plt.savefig('energy-figure.pdf', dpi=300)            # Save figure as a PDF file with set quality
plt.close()                                          # Close figure after saving
