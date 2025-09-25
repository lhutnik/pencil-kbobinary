## Import critical libraries/modules
import os                       # File/directory checks
import pencil as pc             # Local module for handling Pencil Code outputs
import numpy as np              # Arrays
import matplotlib.pyplot as plt # Plotting

from matplotlib.collections import LineCollection
def colored_line(x, y, c, ax, **lc_kwargs):
    """
    Plot a line with a color specified along the line by a third value.

    It does this by creating a collection of line segments. Each line segment is
    made up of two straight lines each connecting the current (x, y) point to the
    midpoints of the lines connecting the current point with its two neighbors.
    This creates a smooth line with no gaps between the line segments.

    Parameters
    ----------
    x, y : array-like
        The horizontal and vertical coordinates of the data points.
    c : array-like
        The color values, which should be the same size as x and y.
    ax : Axes
        Axis object on which to plot the colored line.
    **lc_kwargs
        Any additional arguments to pass to matplotlib.collections.LineCollection
        constructor. This should not include the array keyword argument because
        that is set to the color argument. If provided, it will be overridden.

    Returns
    -------
    matplotlib.collections.LineCollection
        The generated line collection representing the colored line.
    """
    if "array" in lc_kwargs:
        warnings.warn('The provided "array" keyword argument will be overridden')

    # Default the capstyle to butt so that the line segments smoothly line up
    default_kwargs = {"capstyle": "butt"}
    default_kwargs.update(lc_kwargs)

    # Compute the midpoints of the line segments. Include the first and last points
    # twice so we don't need any special syntax later to handle them.
    x = np.asarray(x)
    y = np.asarray(y)
    x_midpts = np.hstack((x[0], 0.5 * (x[1:] + x[:-1]), x[-1]))
    y_midpts = np.hstack((y[0], 0.5 * (y[1:] + y[:-1]), y[-1]))

    # Determine the start, middle, and end coordinate pair of each line segment.
    # Use the reshape to add an extra dimension so each pair of points is in its
    # own list. Then concatenate them to create:
    # [
    #   [(x1_start, y1_start), (x1_mid, y1_mid), (x1_end, y1_end)],
    #   [(x2_start, y2_start), (x2_mid, y2_mid), (x2_end, y2_end)],
    #   ...
    # ]
    coord_start = np.column_stack((x_midpts[:-1], y_midpts[:-1]))[:, np.newaxis, :]
    coord_mid = np.column_stack((x, y))[:, np.newaxis, :]
    coord_end = np.column_stack((x_midpts[1:], y_midpts[1:]))[:, np.newaxis, :]
    segments = np.concatenate((coord_start, coord_mid, coord_end), axis=1)

    lc = LineCollection(segments, **default_kwargs)
    lc.set_array(c)  # set the colors of each segment

    return ax.add_collection(lc)


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

## Additional line to plot for analytic solution
range = np.arange(-np.pi,np.pi,0.02)   # Range of true anomalies to cover in radians
peri_angle = -0.6842318043085592       # radians; periapsis angle from +x direction
a = -0.065846062618                    # Semi-major axis calculated separately
e = 15.66145548                        # Eccentricity calculated separately
r_array = (a*(1-e**2))/(1+e*np.cos(range-peri_angle))       # Array of radial distances
x_pos, y_pos = r_array*np.cos(range), r_array*np.sin(range) # Polar coordinates to cartesian

## Plot location of particles throughout simulation
fig, ax = plt.subplots(1, 1)                                                   # Figure and axes established
ax.plot(parray1[:,2], parray1[:,3], 'r.', label=r'$aps=$'+str(parray1[0,5]))   # Plotting particle 1 position
min, max = np.max(parray2[:,0]), np.max(parray2[:,0])                          # Minimum and maximum of time
lines2 = colored_line(parray2[:,2], parray2[:,3], parray2[:,0], ax, linewidth=5, cmap='seismic') # Particle 2 position with time
## Plot analytical solution provided input
ax.plot(x_pos, y_pos, 'g', label='Analytical', zorder=10) # Analytical expectation based on periapsis from simulation

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
major_ticks = np.arange(-4, 4, 4/32)         # Ticks for grid
ax.set_xticks(major_ticks)                   # Applying ticks to x axis
ax.set_yticks(major_ticks)                   # Applying ticks to y axis
ax.grid(which='major',zorder=1, alpha=0.5)   # Applying grid based on major ticks
plt.tight_layout()                           # Remove overlapping and clipping
plt.savefig('pstalk-figure.pdf', dpi=300)    # Save figure as a PDF file with set quality
plt.close()                                  # Close figure after saving

