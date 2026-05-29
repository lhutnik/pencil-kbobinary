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
c_s = 5.6476e1
G = 1
a_hel = 1.906e4
m_sol = 1.519e10

## Check the grid size to apply to the plot
grid   = pc.read.grid()
x1, y1 = grid.x[3], grid.y[3]   # First corner coordinates (avoiding ghost cells)
x2, y2 = grid.x[-4], grid.y[-4] # Second corner coordinates (avoiding ghost cells)

## Plot 2d density field at midplane
#plt.imshow(rho,origin='lower', cmap='inferno', interpolation='none', extent=[x1,x2,y1,y2]) # Plotting as image
fig, ax = plt.subplots()

## Plot position of particles at saved time 
ax.scatter(p1[0], p1[1], color='blue',label='Primary')
ax.scatter(p2[0], p2[1], color='orange',label='Secondary')

## Plot the Hill radius and its smoothing fraction
r_Hill1 = a_hel*(m1/(3*(m1+m_sol)))**(1/3)
r_Hill2 = a_hel*(m2/(3*(m2+m_sol)))**(1/3)
hill1 = plt.Circle(p1, r_Hill1, color='b', fill=False) 
hill2 = plt.Circle(p2, r_Hill2, color='b', fill=False) 
ax.add_patch(hill1), ax.add_patch(hill2)

## Plot the Bondi radius
r_Bondi1 = (2*G*m1)/(c_s**2)
r_Bondi2 = (2*G*m2)/(c_s**2)
bondi1 = plt.Circle(p1, r_Bondi1, color='r', fill=False) 
bondi2 = plt.Circle(p2, r_Bondi2, color='r', fill=False) 
ax.add_patch(bondi1), ax.add_patch(bondi2)

## Plot specific settings
plt.xlabel(r'$x$', fontsize=15)                             # x-axis label
plt.ylabel(r'$y$', fontsize=15)                             # y-axis label
plt.xlim(-4, 4)                                            # x limits
plt.ylim(-4, 4)                                            # y limits
#plt.axis('scaled')                                   # Axes are scaled to match one another
#bar = plt.colorbar(label=r"$\log_{\rm 10}(\rho/\rho_0)$")   # Colorbar
#bar.set_label(r"$\log_{\rm 10}(\rho/\rho_0)$", fontsize=15) # Colorbar label
plt.legend()

#plt.title("Shearing Sheet - Scattering")                      # Plot title
plt.savefig("slice-plot.pdf", dpi=300)                      # Save figure as a PDF file with set quality
plt.close()                                                   # Close figure after saving
