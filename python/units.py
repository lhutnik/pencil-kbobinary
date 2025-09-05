## Import critical libraries 
import numpy as np                # Pi constant 
#import scipy as sp               # 
from math import sqrt             # Math functions
#import pylab as plt              #
#import matplotlib                # 
import astropy.constants as const # Constants

## Set constants
pi     = np.pi
AU     = const.au.cgs.value      # Astronomical unit in cm
#Mearth = const.M_earth.cgs.value # Mass of Earth in g
Msun   = const.M_sun.cgs.value   # Mass of Sun in g
G      = const.G.cgs.value       # Gravitational constant
amu    = const.u.cgs.value       # Atomic mass unit
kb     = const.k_B.cgs.value     # Boltzmann constant
Rgas   = const.R.cgs.value       # Gas constant
yr     = 3.1558149504e7          # Year in seconds

sigma_mol = 2e-15
Gamma     = 1.4                    # For a diatomic gas; 5/3 for monatomic, 8/6 for non-linear triatomic
#Gamma1    = 1./Gamma
mmol      = 2.3                    #
Rgasmu    = Rgas/mmol              #
cp        = Gamma*Rgasmu/(Gamma-1) #
cv        = cp/Gamma               # 

## Input choices
l3D = False # Setting 3D vs. 2D case
rr  = 20    # Distance of binary in AU
h   = 0.05  # Scale height ratio at distance rr
e   = 0.0   # Eccentricity of mutual binary orbit
Q   = 30    # Toomre Q; >1 for gravitationally stable region
mass_ratio = 0.5    # Binary mass ratio
Hill_frac  = 0.05 # Fraction of Hill radius for furthest separation (apoapsis)

## Solving for other variables given choices 
r     = rr*AU                 # Radial distance in cm
H     = h*r                   # Gas disk scale height 
Omega = sqrt(G*Msun)/r**1.5   # Local angular velocity around Sun
cs    = Omega*H               # Local sound speed
T     = cs**2/(cp*(Gamma-1))  # Gas temperature 
#print("Distance=",rr," AU")
#print("Aspect Ratio=",H/r)
#print("Scale Height=",H," cm=", H/AU, "AU")
#print("Sound speed=",cs/1e5," km/s")
#print("Temperature=",T," K")
#print("Omega=",Omega," s^{-1}")

## Code units (set in start.in or other configuration files)
cs_code    = 1 
#H_code     = 1 # Set automatically if binary orbit is much smaller than 1 AU
rho_code   = 1
Omega_code = 1

## Mass of our initial binary
Mpluto   = 1.309e25      # g
Msystem  = 0.1 * Mpluto  # Sum of binary masses in g
Mplanet1 = Msystem/(1+mass_ratio) # Mass of primary in g
Mplanet2 = Msystem - Mplanet1     # Mass of secondary in g
Mplanetesimal = Msystem/2

## Orbit of our initial binary
r_Hill      = r*(Msystem/(3*Msun))**(1/3)
bin_sep     = r_Hill*Hill_frac      # Separation of binary at apoapsis

## Computing unit length, time
unit_length = bin_sep              # Twice binary separation to box edge
H_code = H/unit_length             # Code scale height set by unit length
unit_time = 1/Omega * 1/Omega_code # Unit time set by orbital angular velocity
print("Scale Height H (code): ",H_code)
print("Unit length: ",unit_length,"cm")
print("Unit time: ",unit_time,"s")

#unit_velocity = cs # Overspecified, but consistent
unit_velocity = unit_length/unit_time
print("Unit velocity: ",unit_velocity,"cm/s")

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
print("Unit mass: ",unit_mass,"g")
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
M1code = Mplanet1/unit_mass
M2code = Mplanet2/unit_mass
print("M1 (code): ",M1code)
print("M2 (code): ",M2code)

## Gravitational constant in code units
G_code = cs_code*Omega_code/(Q*pi*Sigma_code)
print("Gravitational constant (code): ",G_code)

## Hill radius and separation of planetesimals in code units
sep_code = bin_sep/unit_length   # Binary separation in code units 
print("Planetesimal separation (code): ",sep_code)
## Particle positions along common axis at apoapsis
x1 = -(Mplanet2)/(Msystem)*sep_code # Position of particle 1 at apoapsis
x2 = (Mplanet1)/(Msystem)*sep_code  # Position of particle 2 at apoapsis
print("Particle position at apoapsis:")
print("x1 =",x1)
print("x2 =",x2)

## Calculate relative velocity of planetesimals 
v_rel      = sqrt(G*(Msystem)*((2/bin_sep)-((1+e)/bin_sep)))
v_rel_code = v_rel/unit_velocity
v1         = -(Mplanet2)/(Msystem)*v_rel_code # Velocity of particle 1 at apoapsis 
v2         = (Mplanet1)/(Msystem)*v_rel_code  # Velocity of particle 2 at apoapsis 
print("Velocities of particles:")
print("v1 =",v1)
print("v2 =",v2)

## Calculate rhopswarm given planetesimal mass 
rhopswarm1 = M1code * 1e2 # Black box conversion from mass to rhopswarm
rhopswarm2 = M2code * 1e2
print("rhopswarm of particles:")
print("rhopswarm1 =",rhopswarm1)
print("rhopswarm2 =",rhopswarm2)

