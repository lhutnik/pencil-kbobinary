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

## Alternatively, use initial conditions and updated mass
xp = 1.0
yp = 1.99
xv = 0.0
yv = -0.05

M = (rhopswarm) * (side_length/cell_side)**2 # Column density is being used for 2D case?
mu = G*M                # Gravitational parameter

r = np.sqrt(xp**2 + yp**2) # Initial radial distance
v2 = xv**2 + yv**2      # Initial velocity^2 magnitude

a = 1/(2/r - v2/mu)  # Semi-major axis from rewritten vis-viva equation

vinf = np.sqrt(-mu/a)   # Set r=infinity for velocity at infinity

h = abs(xv*yp - yv*xp) # Specific angular momentum at any point in trajectory
b = h/vinf           # Impact paramter or semi-minor axis

e = np.sqrt(1+b**2/a**2) # Eccentricity from semi-minor and semi-major axes

print('A representative hyperbolic trajectory for a 2D mass:')
print('Semi-major axis:',a,'code length units')
print('Eccentricity:',e)

