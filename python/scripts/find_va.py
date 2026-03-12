## Import critical libraries
import numpy as np

## Given outparray.csv from pvpr-plot.py, determine distance, velocity, and angle at periapsis
vp_array = np.loadtxt("outparray.csv", delimiter=',')     # Load array of positions and velocities over time
vp = 0                                                    # Minimum 'magnitude' of position vector
for i in range(len(vp_array[:,0])):                       # For all times indices
    line_m = np.sqrt(vp_array[i,3]**2 + vp_array[i,4]**2) # Velocity vector magnitude 
    if line_m > vp:                                       # If higher, then this could be the periapsis 
        vp = line_m
        x_coord = vp_array[i,1]
        y_coord = vp_array[i,2]
        rp = np.sqrt(x_coord**2 + y_coord**2)

angle = np.arctan(y_coord/x_coord)                        # tan(theta)=y/x=opp/adj

print('Maximum velocity:',vp)
print('Closest distance:',rp,'(Components: x=',x_coord,'y=',y_coord,')')
print('Origin to periapsis angle:',angle,'radians')

## From quantities in code units, determine characteristics of hyperbolic trajectory
rhopswarm   = 1.0
cell_side   = 32  # Grid cells per side
side_length = 4.0 # Length of one side


G = 1/(4*np.pi)       # Gravitational constant #0.01061032953945969
M = (side_length/cell_side)**3 * (rhopswarm) # Mass of particle at origin

a = ((2/rp)-(vp**2)/(G*M))**(-1) # Semi-major axis in code length units
e = 1 - rp/a                     # Eccentricity
print('A representative hyperbolic trajectory has:')
print('Semi-major axis:',a,'code length units')
print('Eccentricity:',e)
