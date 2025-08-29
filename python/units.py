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
Gamma=1.4
Gamma1=1./Gamma
mmol=2.3
Rgasmu=Rgas/mmol
cp=Gamma*Rgasmu/(Gamma-1)
cv=cp/Gamma

## Choices for units  
l3D=False # Setting 3D vs. 2D case
rr = 20   # Distance of binary in AU
h = 0.05  # Scale height ratio at distance rr

## Solving for other variables given choices 
r = rr*AU
H = h*r
Omega = sqrt(G*Msun)/r**1.5
cs = Omega*H
Temperature=cs**2/(cp*(Gamma-1)) 

print("Distance=",rr," AU")
print("Aspect Ratio=",H/r)
print("Scale Height=",H," cm=", H/AU, "AU")
print("Sound speed=",cs/1e5," km/s")
print("Temperature=",Temperature," K")
print("Omega=",Omega," s^{-1}")

## Code units (set in start.in or other configuration files)
cs_code    = 1 
H_code     = 1
rho_code   = 1
Omega_code = 1

## Computing unit length, time
unit_length   = H / H_code
unit_time     = 1/Omega * 1/Omega_code 

print("Unit length: ",unit_length," code units")
print("Unit time: ",unit_time," code units")


#unit_velocity = cs -- overspecified, but consistent
unit_velocity = unit_length/unit_time   
print("Unit velocity: ",unit_velocity," code units")

## Choose the Q value 
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

print("unit_mass=",unit_mass)
print("Msun in code units=",Msun_code)
print("unit_volume_density=",unit_volume_density," g/cm^3")
print("unit_column_density=",unit_column_density," g/cm^2")

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
Mpluto        = 1.309e25     # g
Mplanetesimal = 0.1 * Mpluto # g
Mplanetesimal_code = Mplanetesimal/unit_mass

print("planetesimal mass in physics units=",Mplanetesimal," g")
print("planetesimal mass in code units=",Mplanetesimal_code)

## Gravitational constant in code units
G_code = cs_code*Omega_code/(Q*pi*Sigma_code)
print("gravitational constant in code units:",G_code)

