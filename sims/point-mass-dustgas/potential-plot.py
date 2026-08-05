## Import critical libraries/modules
import pencil as pc             # Local module for handling Pencil Code outputs
import numpy as np              # Arrays
import matplotlib.pyplot as plt # Plotting

## Read in VAR files from Pencil Code output
#ff  = pc.read.var(trimall=True)  # Pulling full VAR file output
#pot = ff.potself[0,:,:]          # Self gravitational potential at z=0 (z,y,x)

## Read in point mass potential data saved as CSV files
q1 = np.loadtxt('omega2_output_mass1.csv')
q2 = np.loadtxt('omega2_output_mass2.csv')

rr = q1[:, 0]        # Radial distance from point mass
omega2_q1, omega2_q2 = q1[:, 1], q2[:, 1] # Omega^2 acceleration term

## Check the grid size to apply to the plot
#grid   = pc.read.grid()
#x1, y1 = grid.x[3], grid.y[3]   # First corner coordinates
#x2, y2 = grid.x[-4], grid.y[-4] # Second corner coordinates

## Plot 2D potential at midplane
#plt.imshow(pot, origin='lower', cmap='inferno', interpolation='none', extent=[x1,x2,y1,y2]) # Plotting as image
#plt.xlabel(r'$x$', fontsize=15)                 # x-axis label
#plt.ylabel(r'$y$', fontsize=15)                 # y-axis label
#bar = plt.colorbar()                            # Colorbar
#bar.set_label(r"Self Potential", fontsize=15)   # Colorbar label
#plt.title("")                                  # Plot title
#plt.savefig("potential-2dfigure.pdf", dpi=300)  # Save figure as a PDF file with set quality
#plt.close()                                     # Close figure after saving

## Plot 1D potential line through point sources
#rangepot = np.linspace(-4,4,256)
#pot1d = ff.potself[0,0,:]
print('Length of Omega^2 list:',len(omega2_q1))
#plt.plot(rangepot, pot1d, label='$256^2$ case')
G = 1.
m_p = 0.5
smoothing_length = 2.5e-1
omega2_ayl = []
for r in rr:
	if r <= smoothing_length:
		omega2 = G*(m_p)/(smoothing)**3 * (3*r/smoothing_length - 4)
		omega2_ayl.append(omega2)
	elif r > smoothing_length:
		omega2 = -G*(m_p)/(r)**3
		omega2_ayl.append(omega2)
plt.plot(rr, omega2_ayl, linestyle='--', label='Analytical', alpha=0.5)
plt.plot(rr, omega2_q1, marker='o', linestyle='-', label='Pencil Output', alpha=0.5)
plt.grid(True)
#print('Minimum:',np.min(potself))

plt.xlabel(r'$x$', fontsize=15)             # x-axis label
plt.ylabel(r'$\Omega^2$ Term', fontsize=15) # y-axis label
#bar = plt.colorbar()   # Colorbar
#bar.set_label(r"potself", fontsize=15) # Colorbar label
#plt.ylim(-4,4)                                 # Potential limits (zoom in to view)
plt.legend()

plt.title('Boley Smoothed Potential/Acceleration Factor') # Title of plot
plt.savefig("potential-1dfigure.pdf", dpi=300)       # Save figure as a PDF file with set quality
plt.close()                                          # Close figure after saving
