## Import critical libraries 
import numpy as np                # Pi constant 
#import scipy as sp               # 
from math import sqrt             # Math functions
#import pylab as plt              #
#import matplotlib                # 
import astropy.constants as const # Constants

## Set constants 
pi     = np.pi
AU     = const.au.cgs.value
Mearth = const.M_earth.cgs.value
Msun   = const.M_sun.cgs.value
G      = const.G.cgs.value
amu    = const.u.cgs.value
kb     = const.k_B.cgs.value
Rgas   = const.R.cgs.value
yr     = 3.1558149504e7

sigma_mol = 2e-15
Gamma     = 1.4 # ~5/3
Gamma1    = 1./Gamma
mmol      = 2.3
Rgasmu    = Rgas/mmol
cp        = Gamma*Rgasmu/(Gamma-1)
cv        = cp/Gamma

## Input choices  
l3D = False # Setting 3D vs. 2D case
rr  = 20    # Distance of binary in AU
h   = 0.05  # Scale height ratio at distance rr
e   = 0.0   # Eccentricity of mutual binary orbit
Hill_frac = 0.05 # Fraction of Hill radius for furthest separation

## Solving for other variables given choices 
r = rr*AU                          # Radial distance in cm
H = h*r                            # Gas disk scale height 
Omega = sqrt(G*Msun)/r**1.5        # 
cs = Omega*H                       # Local sound speed
Temperature = cs**2/(cp*(Gamma-1)) # Gas temperature 
#print("Distance=",rr," AU")
#print("Aspect Ratio=",H/r)
#print("Scale Height=",H," cm=", H/AU, "AU")
#print("Sound speed=",cs/1e5," km/s")
#print("Temperature=",Temperature," K")
#print("Omega=",Omega," s^{-1}")

## Code units (set in start.in or other configuration files)
cs_code    = 1 
#H_code     = 1 # Set automatically if binary orbit is small
rho_code   = 1
Omega_code = 1

## Mass of our initial binary
Mpluto        = 1.309e25     # g
Mplanetesimal = 0.1 * Mpluto # g

## Orbit of our initial binary
r_Hill      = r*(Mplanetesimal/(3*Msun))**(1/3)
bin_sep     = r_Hill*Hill_frac      # Separation of binary at apoapsis

## Computing unit length, time
unit_length = 2 * bin_sep    # Twice binary separation to box edge
H_code = H/unit_length
print("Scale Height H (code): ",H_code)
unit_time     = 1/Omega * 1/Omega_code 
print("Unit length: ",unit_length)
print("Unit time: ",unit_time)

#unit_velocity = cs # Overspecified, but consistent
unit_velocity = unit_length/unit_time   
print("Unit velocity: ",unit_velocity)

## Choose the Toomre Q value (>1 for stability)
Q = 30 

## Choosing Q sets Sigma (surface density) and rho (density)
Sigma = cs*Omega/(G*Q*pi)
if (l3D): # If solving for 3D case
    rho = Sigma/(sqrt(2*pi)*H)
    unit_volume_density = rho/rho_code
    Sigma_code = rho_code * sqrt(2*pi) * H_code
    unit_column_density = Sigma/Sigma_code
    unit_mass = unit_volume_density * unit_length**3
else: # If solving for 2D case
    Sigma_code = rho_code # in 2D rho is actually column density
    unit_column_density = Sigma/Sigma_code
    rho = Sigma/(sqrt(2*pi)*H)
    # actual code units for volume density
    actual_rho_code = Sigma_code/(sqrt(2*pi)*H_code)
    unit_volume_density = rho/actual_rho_code
    unit_mass = unit_column_density * unit_length**2
Msun_code = Msun/unit_mass
print("Unit mass: ",unit_mass)
print("M_Sun (code): ",Msun_code)
#print("unit_volume_density=",unit_volume_density," g/cm^3")
#print("unit_column_density=",unit_column_density," g/cm^2")

## box volume 
#Lx_code = 2.
#Ly_code = 2.
#Lz_code = 2.

#Lx=Lx_code*unit_length
#Ly=Ly_code*unit_length
#Lz=Lz_code*unit_length

#box_area = Lx*Ly
#mbox = Sigma * box_area

#volume_code = Lx_code*Ly_code*Lz_code
#mbox_code = rho0_code * volume_code 
#unit_mass = mbox/mbox_code

## Determining code unit mass of planetesimals in simulation
Mplanetesimal_code = Mplanetesimal/unit_mass
print("Planetesimal mass (code): ",Mplanetesimal_code)

## Gravitational constant in code units
G_code = cs_code*Omega_code/(Q*pi*Sigma_code)
print("Gravitational constant (code): ",G_code)

## Hill radius and separation of planetesimals in code units
sep_code    = bin_sep/unit_length   # Binary separation in code units 
print("Planetesimal separation (code): ",sep_code)
## Particle positions along common axis at apoapsis
x1 = -(Mplanetesimal)/(Mplanetesimal+Mplanetesimal)*sep_code
x2 = (Mplanetesimal)/(Mplanetesimal+Mplanetesimal)*sep_code
print("Particle position at apoapsis:")
print("x1 =",x1)
print("x2 =",x2)
## [[We further want to calculate the separation from COM for differing masses]]

## Calculate relative velocity of planetesimals 
v_rel      = 0.5*sqrt(G*(Mplanetesimal+Mplanetesimal)*((2/bin_sep)-((1+e)/bin_sep)))
v_rel_code = v_rel/unit_velocity
v1         = -(Mplanetesimal)/(Mplanetesimal+Mplanetesimal)*v_rel_code # Velocity of particle 1 at apoapsis 
v2         = (Mplanetesimal)/(Mplanetesimal+Mplanetesimal)*v_rel_code  # Velocity of particle 2 at apoapsis 
print("Velocities of particles:")
print("v1 =",v1)
print("v2 =",v2)

## Calculate rhopswarm given planetesimal mass 
rhopswarm = Mplanetesimal_code * 1e2 # Black box conversion from mass to rhopswarm from scattering problems
print("rhopswarm of particles:")
print("rhopswarm1 =",rhopswarm)
print("rhopswarm2 =",rhopswarm)

