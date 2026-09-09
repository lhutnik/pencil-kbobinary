import numpy as np

# constants
au = 1.49e13
mu = 2.3
gamma=1.4
kb=1.38e-16
amu= 1.66e-24
Msun=2e33
Mpluto = 1.3e25 #g
G = 6.68e-8
yr=3.1e7

#######################################################################
# code units - set by the binary parameters

# vcirc_code      = 1     # circular velocity of binary
# Omegabin_code   = 1     # angular frequency of binary
# separation_code = 1     # separation of binary

#######################################################################

# physical variables
r                 = 20*au # distance from the sun
T                 = 20    # temperature in K

#######################################################################

# choices that set the physical units

Mp                = 5e-3*Mpluto  # mass of the binary
Rhill             = r*np.cbrt(Mp/Msun/3)
a                 = 0.01*Rhill  # separation of the binary

# sanity check -- should give cs_code = 50 and Omegasun_code = 0.1
#a= 28317294633.91888  #cm
#Mp = 1.7160758597858634e+23  #g

#######################################################################

#
# Distance sets the binary period around the Sun. This compares to the binary period.
# Temperature sets the sound speed. Sound speed compares to the circular velocity.
#

print("Time quantities")
Omegasun=np.sqrt(G*Msun)/r**1.5
print("Omegasun=",Omegasun," 1/s")
print("Period sun=",2*np.pi/Omegasun/yr," yr \n")
Omegabin  = np.sqrt(G*Mp/a**3)
print("Omegabin=",Omegabin," 1/s")
print("Period bin=",2*np.pi/Omegabin/yr," yr \n")

print("Velocity quantities")
cs = np.sqrt(gamma*T*kb/mu/amu)
print("cs=",cs," cm/s")
vcirc     = Omegabin*a
print("vcirc=",vcirc," cm/s \n")

print("Code units, finally")
cs_code = cs/vcirc
print("cs_code=",cs_code)
Omegasun_code = Omegasun/Omegabin
print("Omegasun_code=",Omegasun_code)


### accretion variables & resolution

Nx=128
Lx=8
dx=Lx/Nx

print("dx=",dx)

beta = 0.1

dv = beta * cs_code/2

rb = 1/dv**2

tb = rb/dv

ts = tb

r = 2*np.sqrt(rb*dv*ts)

print("friction time, in code units=",ts)
print("Stokes number=",ts*Omegasun_code)
print("racc/dx",r/dx)

print("racc/Lx",r/Lx)

#######

#viscosity -- timestep

alpha = 1e-4
H_code = cs_code/Omegasun_code
nu = alpha*cs_code*H_code

C=0.5
dt_visc = C*dx**2/nu
dt_cs = C*dx/cs_code

print("visc, cs timesteps=",dt_visc,dt_cs)
