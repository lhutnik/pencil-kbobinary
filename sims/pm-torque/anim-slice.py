## Import critical libraries/modules
import pencil as pc             # Local module for handling Pencil Code outputs
import pencil_old as oldpc      # Utilize old Pencil Code functions
import numpy as np              # Arrays
import matplotlib.pyplot as plt # Plotting
import matplotlib.animation as animation # Animation

### IMPORT SIMULATION DATA
## Timeseries data
ts = pc.read.timeseries.ts()
t = ts.t
xq1, xq2 = np.array(ts.xq1), np.array(ts.xq2)     # x positions
yq1, yq2 = np.array(ts.yq1), np.array(ts.yq2)     # y positions
#vxq1, vxq2 = np.array(ts.vxq1), np.array(ts.vxq2) # x velocities
#vyq1, vyq2 = np.array(ts.vyq1), np.array(ts.vyq2) # y velocities

## Read in VAR and QVAR files from Pencil Code output
#ff   = pc.read.var(trimall=True) # Pulling full VAR file output
#rho = ff.rho[0,:,:]              # Gas density at z=0, midplane; [z,y,x]
slices = pc.read.slices(field='rho',extension='xy') # Slice object with xy plane and times
#rho = rho_all[frame, :, :]

## Interpolate position data to slice timesteps
xq1_interp, xq2_interp = np.interp(slices.t, ts.t, xq1), np.interp(slices.t, ts.t, xq2)
yq1_interp, yq2_interp = np.interp(slices.t, ts.t, yq1), np.interp(slices.t, ts.t, yq2)

#qvar = oldpc.read_qvar()    # Pulling QVAR file output 
#imass = qvar.mass           # Particle masses 
#m1, m2 = imass[0], imass[1] 
m1, m2 = 0.5, 0.5
#ixp = qvar.xq               # Particle positions at recorded time
#iyp = qvar.yq
#p1, p2 = (ixp[0], iyp[0]), (ixp[1], iyp[1])
#t = qvar.t                  # Simulation time at recording

print(f"Slice data shape: {slices.xy.rho.shape}") 
print(f"Number of time steps: {len(slices.t)}")

## Simulation settings (set already in start.in, run.in) in code units
c_s = 5.6746e1
G = 1
a_hel = 1.1906e4
m_sol = 1.519e10

## Check the grid size to apply to the plot
grid   = pc.read.grid()
x1, y1 = grid.x[3], grid.y[3]   # First corner coordinates (avoiding ghost cells)
x2, y2 = grid.x[-4], grid.y[-4] # Second corner coordinates (avoiding ghost cells)
lx = grid.Lx          # Length of x-axis side
ly = grid.Ly          # Length of y-axis side
x1, x2 = -4, -4+lx
y1, y2 = -4, -4+ly

### PLOT DENSITY
fig, ax = plt.subplots(1, 1, layout="tight") # Initiate figure and axis
## Figure settings
ax.set_xlabel(r'$x$', fontsize=15)                   # x-axis label
ax.set_ylabel(r'$y$', fontsize=15)                   # y-axis label
#ax.set_title(r'$M_1$={} $M_2$={}, $t$={}'.format(m1, m2, f"{t:.5e}")) # Plot title
plt.axis('scaled')                                   # Axes are scaled to match one another
ax.set_xlim(x1, x2)                                  # Set x-axis limits, assuming centered evenly on origin
ax.set_ylim(y1, y2)                                  # Set y-axis limits, assuming centered evenly on origin
nx, ny = len(grid.x)-6 , len(grid.y)-6               # Determine number of cells in each dimension, subtracting 3 ghost cells from each side


## Initialize plot with empty data to be updated later
im = ax.imshow(np.zeros((nx, ny)), origin='lower', cmap='inferno', 
               interpolation='none', extent=[x1,x2,y1,y2], vmin=0.5, vmax=1.5)
p1 = (xq1_interp[0], yq1_interp[0])
p2 = (xq2_interp[0], yq2_interp[0])
scat1 = ax.scatter(xq1_interp[0], yq1_interp[0], label="Primary", color='blue', marker='.') # Position of primary plotted
scat2 = ax.scatter(xq2_interp[0], yq2_interp[0], label="Secondary", color='red', marker='.') # Position of secondary plotted
title = ax.set_title('')
com = ax.scatter(0, 0, color='k', marker='x') # Center of mass indicator

#plt.legend(loc='best')                               # Add legend
bar = plt.colorbar(im, extend='both')   # Colorbar
bar.set_label(r"$\rho_{code}-1$", fontsize=15) # Colorbar label
ax.tick_params(axis='both', which='minor', length=0) # Set tick parameters (0 length)
ax.grid(which='major', alpha=0.5)                    # Applying grid based on major ticks
plt.tight_layout()                                   # Remove overlapping and clipping

## Animation function
def slice_update(frame):
    global p1
    global p2
    p1, p2 = (xq1_interp[frame], yq1_interp[frame]), (xq2_interp[frame], yq2_interp[frame])
    scat1.set_offsets(p1) # Set new position of primary
    scat2.set_offsets(p2) # Set new position of secondary

    rho_frame = slices.xy.rho[frame, :, :] # Access xy attribute, then rho
    im.set_data(rho_frame) # Set density displayed at given frame

    current_t = slices.t[frame] # Recover simulation time when slice was taken
    title.set_text(r'$M_1$={}, $M_2$={}, $t$={:.5e}'.format(m1, m2, current_t)) # Update title
    return im, title, scat1, scat2

## Plot location of particles throughout simulation
#ax.scatter(p1[0], p1[1], label="Primary", color='blue')
#ax.scatter(p2[0], p2[1], label="Secondary", color='orange')

## Plot the Hill radius and its smoothing fraction
#r_Hill1 = a_hel*(m1/(3*(m1+m_sol)))**(1/3)
#r_Hill2 = a_hel*(m2/(3*(m2+m_sol)))**(1/3)
#r_Hillsys = a_hel*(1/(3*(1+m_sol)))**(1/3)
#hill1 = plt.Circle(p1, r_Hill1, color='blue', alpha=0.5, fill=False) 
#hill2 = plt.Circle(p2, r_Hill2, color='orange', alpha=0.5, fill=False) 
#hillsys = plt.Circle((0,0), r_Hillsys, color='green', alpha=0.5, fill=False, label='Mutual', linestyle=':') 
#ax.add_patch(hill1), ax.add_patch(hill2), ax.add_patch(hillsys)

## Plot smoothing fraction 
#frac1, frac2 = 0.12599, 0.12599
#r_smooth1 = r_Hill1 * frac1
#r_smooth2 = r_Hill2 * frac2
#smooth1 = plt.Circle(p1, r_smooth1, color='blue', alpha=0.5, fill=False, linestyle='--') 
#smooth2 = plt.Circle(p2, r_smooth2, color='orange', alpha=0.5, fill=False, linestyle='--') 
#ax.add_patch(smooth1), ax.add_patch(smooth2)

## Plot the Bondi radius
#r_Bondi1 = (2*G*m1)/(c_s**2)
#r_Bondi2 = (2*G*m2)/(c_s**2)
#bondi1 = plt.Circle(p1, r_Bondi1, color='purple', fill=False) 
#bondi2 = plt.Circle(p2, r_Bondi2, color='red', fill=False) 
#ax.add_patch(bondi1), ax.add_patch(bondi2)

## Save and close animation
ani = animation.FuncAnimation(fig, slice_update, frames=len(slices.t), blit=False)  # Create animation 
ani.save('binary_evolution.mp4', writer='ffmpeg', fps=10, dpi=300) # Save animation with set frame rate
plt.close()                                                        # Close figure after saving

