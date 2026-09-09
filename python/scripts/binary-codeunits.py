## Import modules
import numpy as np

## Set constants
au = 1.49e13      # [cm]
mu = 2.3
gamma = 1.4
kb = 1.38e-16     # [CGS]
amu= 1.66e-24     # [g]
Msun = 2e33       # [g]
Mpluto = 1.3e25   # [g]
G = 6.68e-8       # [CGS]
yr = 3.1e7          # [s]

#######################################################################
# code units - set by the binary parameters

# vcirc_code      = 1     # circular velocity of binary
# Omegabin_code   = 1     # angular frequency of binary
# separation_code = 1     # separation of binary

#######################################################################

# Physical variables
r                 = 20 * au # [cm]; heliocentric distance
T                 = 20      # [K]; local temperature of gas disk

#######################################################################

## Choices that set the physical units

Mp                = 5e-3 * Mpluto  # mass of the binary
Rhill             = r*np.cbrt(Mp/Msun/3)
a                 = 0.01 * Rhill   # separation of the binary

# sanity check -- should give cs_code = 50 and Omegasun_code = 0.1
#a= 28317294633.91888  #cm
#Mp = 1.7160758597858634e+23  #g

#######################################################################

#
# Distance sets the binary period around the Sun. This compares to the binary period.
# Temperature sets the sound speed. Sound speed compares to the circular velocity.
#

print("Time quantities")
Omegasun=np.sqrt(G*Msun)/r**1.5 # Heliocentric orbital frequency
print("Omegasun=",Omegasun," 1/s")
print("Period sun=",2*np.pi/Omegasun/yr," yr \n")
Omegabin  = np.sqrt(G*Mp/a**3)
print("Omegabin=",Omegabin," 1/s")
print("Period bin=",2*np.pi/Omegabin/yr," yr \n")

print("Velocity quantities")
cs = np.sqrt(gamma*T*kb/mu/amu) # Sound speed
print("cs=",cs," cm/s")
vcirc     = Omegabin*a # Binary orbital velocity
print("vcirc=",vcirc," cm/s \n")

print("Code units, finally")
cs_code = cs/vcirc # Sound speed / velocity code unit
print("cs_code=",cs_code)
Omegasun_code = Omegasun/Omegabin # Heliocentric Omega / frequency code unit 
print("Omegasun_code=",Omegasun_code)


## Accretion variables & resolution
Nx=128 # Number of mesh points
Lx=8 # Domain side length
dx=Lx/Nx # Domain resolution

print("dx=",dx)

beta = 0.1 # Beta parameter
dv = beta * cs_code/2 # Random velocity of dust
rb = 1/dv**2 # Bondi radius
tb = rb/dv # Bondi time
ts = tb # Bondi time = stopping/friction time
r = 2*np.sqrt(rb*dv*ts) # Accretion radius

print("friction time, in code units=",ts)
print("Stokes number=",ts*Omegasun_code)
print("racc/dx",r/dx)
print("racc/Lx",r/Lx)

## Determine viscous and advective timesteps

alpha = 1e-4 # Sunyaez alpha parameter
H_code = cs_code/Omegasun_code # Scale height in code units for isothermal layers case
nu = alpha*cs_code*H_code # Viscosity

C = 0.5 #CFL number (somewhere between 0 and 1)
dt_visc = C * dx**2 / nu
dt_cs = C * dx / cs_code

print("visc, cs timesteps=",dt_visc,dt_cs)
