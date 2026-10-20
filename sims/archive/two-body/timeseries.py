## Import critical Python libraries/modules
import pencil as pc             # Local module for reading in Pencil Code outputs
import numpy as np              # Arrays
import matplotlib.pyplot as plt # Plotting 

## Import timeseries data from timeseries file
ts=pc.read.ts(datadir='data')

## Plot timeseries 
plt.plot(ts.t/2/np.pi,ts.rhopmax)             # Simple line plot between time and maximum particle density
plt.xlabel(r'$t/T_0$')                        # x-axis label
plt.ylabel(r'max($\rho_p$)/$\rho_0$')         # y-axis label
plt.yscale('log')                             # Plot y scale (logarithmic)
plt.title("Maximum Particle Density vs Time") # Plot title
plt.savefig("density.pdf", dpi=300)           # Save figure as a PDF file with set quality
plt.close()                                   # Close figure instance
