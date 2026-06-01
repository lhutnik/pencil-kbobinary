## Import critical libraries/modules
import pencil as pc             # Local module for handling Pencil Code outputs
import pencil_old as oldpc      # Utilize old Pencil Code functions
import numpy as np              # Arrays
import matplotlib.pyplot as plt # Plotting

## Read in VAR and QVAR files from Pencil Code output
#ff   = pc.read.var(trimall=True) # Pulling full VAR file output
#rho = ff.rho[0,:,:]              # Particle density at z=0; [z,y,x]

qvar = oldpc.read_qvar() # Pulling QVAR file output 
imass = qvar.mass # Particle masses 
m1, m2 = imass[0], imass[1] 
ixp = qvar.xq # Particle positions at recorded time
iyp = qvar.yq
p1, p2 = (ixp[0], iyp[0]), (ixp[1], iyp[1])

## Simulation settings (set already in start.in, run.in) in code units
c_s = 1#5.6476e1
G = 1
a_hel = 1.1906e4
m_sol = 1.519e10

## Check the grid size to apply to the plot
grid   = pc.read.grid()
x1, y1 = grid.x[3], grid.y[3]   # First corner coordinates (avoiding ghost cells)
x2, y2 = grid.x[-4], grid.y[-4] # Second corner coordinates (avoiding ghost cells)
lx = grid.Lx          # Length of x-axis side
ly = grid.Ly          # Length of y-axis side

### PLOT PARTICLE POSITION
fig, ax = plt.subplots(1, 1, layout="tight") # Initiate figure and axis

## Plot gas density across the simulation 
#plt.imshow(rho,origin='lower', cmap='inferno', interpolation='none', extent=[x1,x2,y1,y2]) # Plotting as image

## Plot location of particles throughout simulation
ax.scatter(p1[0], p1[1], label="Primary", color='blue')
ax.scatter(p2[0], p2[1], label="Secondary", color='orange')

## Plot the Hill radius and its smoothing fraction
r_Hill1 = a_hel*(m1/(3*(m1+m_sol)))**(1/3)
r_Hill2 = a_hel*(m2/(3*(m2+m_sol)))**(1/3)
r_Hillsys = a_hel*(1/(3*(1+m_sol)))**(1/3)
hill1 = plt.Circle(p1, r_Hill1, color='blue', alpha=0.5, fill=False) 
hill2 = plt.Circle(p2, r_Hill2, color='orange', alpha=0.5, fill=False) 
hillsys = plt.Circle((0,0), r_Hillsys, color='green', alpha=0.5, fill=False, label='Mutual') 
ax.add_patch(hill1), ax.add_patch(hill2), ax.add_patch(hillsys)

## Plot smoothing fraction 
frac = 0.25
r_smooth1 = a_hel*(m1/(3*(m1+m_sol)))**(1/3) * frac
r_smooth2 = a_hel*(m2/(3*(m2+m_sol)))**(1/3) * frac
smooth1 = plt.Circle(p1, r_smooth1, color='blue', alpha=0.5, fill=False, linestyle='--') 
smooth2 = plt.Circle(p2, r_smooth2, color='orange', alpha=0.5, fill=False, linestyle='--') 
ax.add_patch(smooth1), ax.add_patch(smooth2)

## Plot the Bondi radius
r_Bondi1 = (2*G*m1)/(c_s**2)
r_Bondi2 = (2*G*m2)/(c_s**2)
bondi1 = plt.Circle(p1, r_Bondi1, color='purple', fill=False) 
bondi2 = plt.Circle(p2, r_Bondi2, color='red', fill=False) 
#ax.add_patch(bondi1), ax.add_patch(bondi2)

## Figure settings
ax.set_xlabel(r'$x$', fontsize=15)                   # x-axis label
ax.set_ylabel(r'$y$', fontsize=15)                   # y-axis label
ax.set_title(r'Position vs. Time of $M_1$={}, $M_2$={} Binary'.format(m1, m2)) # Plot title
plt.axis('scaled')                                   # Axes are scaled to match one another
ax.set_xlim(-4, -4+lx)                               # Set x-axis limits, assuming centered evenly on origin
ax.set_ylim(-4, -4+ly)                               # Set y-axis limits, assuming centered evenly on origin
plt.legend(loc='best')                               # Add legend
ax.tick_params(axis='both', which='minor', length=0) # Set tick parameters (0 length)
ax.grid(which='major', alpha=0.5)                    # Applying grid based on major ticks
plt.tight_layout()                                   # Remove overlapping and clipping
plt.savefig('slice-plot.pdf', dpi=300)               # Save figure as a PDF file with set quality
plt.close()                                          # Close figure after saving

