## Import critical libraries/modules
import os
import pencil as pc             # Local module for handling Pencil Code outputs
import numpy as np              # Arrays
import matplotlib.pyplot as plt # Plotting

## Read in timeseries from Pencil Code output
ts       = pc.read.ts(datadir='data') # Pulling timeseries data
xp1, yp1 = ts.xpm, ts.ypm             # Setting x and y variables

SIM = pc.get_sim()


if os.path.exists(os.path.join('.','outparray.csv')):
    parray = np.loadtxt('outparray.csv', delimiter=',')
else: 
    varfiles = 'all' #'VAR###', 'last4', 'first3'
    varlist = SIM.get_varlist(pos=varfiles, particle=False)
    pvarlist = SIM.get_varlist(pos=varfiles, particle=True)

    parray = np.empty(shape=(len(pvarlist), 9))
    print(parray.shape)
    for i, (f, p) in enumerate(zip(varlist, pvarlist)):
        ff = pc.read.var(datadir=SIM.datadir, var_file=f, quiet=True, trimall=False)
        pp = pc.read.pvar(datadir=SIM.datadir, pvarfile=p, quiet=True) #flist=['px','pvx']
        #l_ipars = pp.ipars # Particle number
        l_px = pp.xp       # Particle x position                                 
        l_py = pp.yp       # Particle y position 
        #l_pz = pp.zp       # Particle z position 
        l_vx = pp.vpx     # Particle x velocity 
        l_vy = pp.vpy      # Particle y velocity 
        #l_vz = pp.vpz      # Particle z velocity
        t = ff.t
        parray[i] = [t, l_px[0], l_py[0], l_vx[0], l_vy[0], l_px[1], l_py[1], l_vx[1], l_vy[1]]
        #print(l_ipars)
        #print(l_px)
        #print(l_vx)
    print(parray)
    np.savetxt('outparray.csv', parray, delimiter=',', header='t, l_px[0], l_py[0], l_vx[0], l_vy[0], l_px[1], l_py[1], l_vx[1], l_vy[1]')

## Check the grid size to apply to the plot
grid   = pc.read.grid()
x1, y1 = grid.x[3], grid.y[3]   # first corner coordinates
x2, y2 = grid.x[-4], grid.y[-4] # second corner coordinates

## Plot location of point masses at end of simulation
fig, ax = plt.subplots(1, 1)               # Figure and axes established
ax.plot(parray[:,1], parray[:,2], 'ro')    # Plotting particle 1 position
ax.plot(parray[:,5], parray[:,6], 'b.')    # Plotting particle 1 position
ax.set_xlabel(r'$x$', fontsize=15)         # x-axis label
ax.set_ylabel(r'$y$', fontsize=15)         # y-axis label
plt.axis('scaled')                         # Axes are scaled to match one another
ax.set_xlim(x1, x2)                        # Set x-axis limits
ax.set_ylim(y1, y2)                        # Set y-axis limits
ax.set_title("Gravitational Scattering Test") # Plot title
plt.tight_layout()                         # Remove overlapping and clipping
plt.savefig("point-plot.pdf", dpi=300)     # Save figure as a PDF file with set quality
plt.close()                                # Close figure after saving

