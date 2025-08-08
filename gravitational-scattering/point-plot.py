## Import critical libraries/modules
import os                       # File/directory checks
import pencil as pc             # Local module for handling Pencil Code outputs
import numpy as np              # Arrays
import matplotlib.pyplot as plt # Plotting

## Read in timeseries from Pencil Code output
#ts       = pc.read.ts(datadir='data') # Pulling timeseries data
#xp1, yp1 = ts.xpm, ts.ypm             # Setting x and y variables

## Input or generate array with particle information
if os.path.exists(os.path.join('.','outparray.csv')):       # If a .csv is already available, do not reproduce
    parray = np.loadtxt('outparray.csv', delimiter=',')     # Reuse previously made .csv from relevant data
else: 
    SIM = pc.get_sim()                                      # Pull simulation 
    varfiles = 'all' #'VAR###', 'last4', 'first3'           # Specifying which VAR/PVAR files to pull
    varlist = SIM.get_varlist(pos=varfiles, particle=False) # Pull time and cell information from VAR files
    pvarlist = SIM.get_varlist(pos=varfiles, particle=True) # Pull particle information from PVAR files

    parray = np.empty(shape=(len(pvarlist), 11))             # Empty array to insert data into
    for i, (f, p) in enumerate(zip(varlist, pvarlist)):
        ff = pc.read.var(datadir=SIM.datadir, var_file=f, quiet=True, trimall=False)
        pp = pc.read.pvar(datadir=SIM.datadir, pvarfile=p, quiet=True)
        #l_ipars = pp.ipars # Particle number
        l_px = pp.xp    # Particle x position                                 
        l_py = pp.yp    # Particle y position 
        #l_pz = pp.zp   # Particle z position 
        l_vx = pp.vpx   # Particle x velocity 
        l_vy = pp.vpy   # Particle y velocity 
        #l_vz = pp.vpz  # Particle z velocity
        l_aps = pp.aps  # Particle radii   
        t = ff.t
        parray[i] = [t, l_px[0], l_py[0], l_vx[0], l_vy[0], l_px[1], l_py[1], l_vx[1], l_vy[1], l_aps[0], l_aps[1]]
    print(parray)
    np.savetxt('outparray.csv', parray, delimiter=',', header='t, l_px[0], l_py[0], l_vx[0], l_vy[0], l_px[1], l_py[1], l_vx[1], l_vy[1], l_aps[0], l_aps[1]')

## Check the grid size to apply to the plot
grid   = pc.read.grid()         # Read relevant information from grid data
x1, y1 = grid.x[3], grid.y[3]   # first corner coordinates
x2, y2 = grid.x[-4], grid.y[-4] # second corner coordinates

## Plot location of particles throughout simulation
fig, ax = plt.subplots(1, 1)                 # Figure and axes established
ax.plot(parray[:,1], parray[:,2], 'r.', label='aps='+str(parray[-1,9]))    # Plotting particle 1 position
ax.plot(parray[:,5], parray[:,6], 'b.', label='aps='+str(parray[-1,10]))   # Plotting particle 2 position
ax.set_xlabel(r'$x$', fontsize=15)           # x-axis label
ax.set_ylabel(r'$y$', fontsize=15)           # y-axis label
plt.axis('scaled')                           # Axes are scaled to match one another
ax.set_xlim(x1, x2)                          # Set x-axis limits
ax.set_ylim(y1, y2)                          # Set y-axis limits
ax.set_title(r'$t={}$'.format(parray[-1,0])) # Plot title
plt.legend()                                 # Add legend
plt.tight_layout()                           # Remove overlapping and clipping
plt.savefig('point-plot.pdf', dpi=300)       # Save figure as a PDF file with set quality
plt.close()                                  # Close figure after saving

