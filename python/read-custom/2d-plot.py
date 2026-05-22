## Import critical libraries/modules
import pencil as pc             # Local module for handling Pencil Code outputs
import numpy as np              # Arrays
import matplotlib.pyplot as plt # Plotting

## Read in VAR files from Pencil Code output
ts = 5000                        # Timestep 
ff   = pc.read.var(trimall=True) # Pulling full VAR file output
rho = ff.rho[0,:,:]              # Particle density at z=0; [z,y,x]

## Simulation settings (set already in start.in, run.in)
cs = 1
G = 1
a_helio = 20
m_helio = 10 
m1, m2 = 0.5, 0.5
p1, p2 = (0,1), (0,-1) # Particle positions

## Check the grid size to apply to the plot
grid   = pc.read.grid()
x1, y1 = grid.x[3], grid.y[3]   # first corner coordinates
x2, y2 = grid.x[-4], grid.y[-4] # second corner coordinates

## Plot 2d density field at midplane
plt.imshow(rho,origin='lower',cmap='inferno',interpolation='none',
		extent=[x1,x2,y1,2])                                  # Plotting as image
#ax.add_patch(circle2)

## Plot position of particles at timestep 
plt.plot(p1, label='Primary')
plt.plot(p1, label='Secondary')

## Plot the Hill radius and its smoothing fraction
r_Hill1 = a_helio*(m1/(3*(m1+m_helio)))**(1/3)
r_Hill2 = a_helio*(m2/(3*(m2+m_helio)))**(1/3)
hill1 = plt.Circle(p1, r_Hill1, color='b', fill=False) 
hill2 = plt.Circle(p1, r_Hill2, color='b', fill=False) 

## Plot the Bondi radius
r_Bondi1 = (2*G*m1)/(cs**2)
r_Bondi2 = (2*G*m2)/(cs**2)
bondi1 = plt.Circle(p1, r_Bondi1, color='r', fill=False) 
bondi2 = plt.Circle(p1, r_Bondi2, color='r', fill=False) 

## Plot specific settings
plt.xlabel(r'$x$', fontsize=15)                               # x-axis label
plt.ylabel(r'$y$', fontsize=15)                               # y-axis label
bar = plt.colorbar(label=r"$\log_{\rm 10}(\rho/\rho_0)$")   # Colorbar
bar.set_label(r"$\log_{\rm 10}(\rho/\rho_0)$", fontsize=15) # Colorbar label
plt.legend()

#plt.title("Shearing Sheet - Scattering")                      # Plot title
plt.savefig("density-plot.pdf", dpi=300)                      # Save figure as a PDF file with set quality
plt.close()                                                   # Close figure after saving
