# %%
## Import libraries 
import numpy as np                # Pi constant and math
from math import sqrt             # Math functions
import astropy.constants as const # Constants

# %%
## Constants
pi     = np.pi
AU     = const.au.cgs.value      # [cm]; Astronomical unit (AU)
Mearth = const.M_earth.cgs.value # [g]; Mass of Earth
Msun   = const.M_sun.cgs.value   # [g]; Mass of Sun
Mpluto = 1.3025e25               # [g]; Mass of Pluto
G      = const.G.cgs.value       # [cm^3 g^-1 s^-2]; Gravitational constant
amu    = const.u.cgs.value       # [g]; Atomic mass unit
kb     = const.k_B.cgs.value     # [ergs cm^-2 s^-1 K^-4]; Boltzmann constant
Rgas   = const.R.cgs.value       # [ergs K^-1 mol^-1]; Gas constant
yr     = 3.1558149504e7          # [s]; one year
gamma  = 1.4                     # Adiabatic index; diatomic=7/5, monatomic=5/3, non-linear triatomic=8/6
mmol   = 2.3                     # Mean molecular weight (proton masses); assumes bulk H2, He, and trace molecular gas
Rgasmu = Rgas/mmol               # Gas constant / mean molecular weight
cp     = gamma*Rgasmu/(gamma-1)  # Specific heat capacity at constant pressure
cv     = cp/gamma                # Specific heat capacity at constant volume

# %%
## INPUT CHOICES FOR PHYSICAL UNITS
rr = 20      # [AU]; heliocentric distance of binary 
r  = rr * AU # [cm]; heliocentric distance of binary
#h  = 0.05    # Scale height ratio at distance rr
#H  = h * r   # [cm]; gas disk scale height 
T  = 20      # [K]; local temperature of gas (isothermal case) 
Q  = 30      # Toomre Q; >1 for gravitationally stable region
alpha = 1e-4 # Shakura and Sunyaev alpha parameter
beta = 1e-1  # Beta parameter

## BINARY SETTINGS
mass_ratio = 1       # Binary mass ratio (0,1]; M2/M1 = f
Hill_frac  = 0.01    # Fraction of mutual Hill radius for furthest separation (apoapsis) [0.01,0.4]
e          = 0.0     # Eccentricity of mutual binary orbit [0,1)
Msystem  = 5e-3 * Mpluto              # [g]; Sum of binary component masses
Mplanet1 = Msystem / (1 + mass_ratio) # [g]; Mass of primary 
Mplanet2 = Msystem - Mplanet1         # [g]; Mass of secondary

## SIMULATION SETTINGS
Lx = 8     # Full x or y grid size in code units; assumes equal x,y scale
Nx = 128   # Number of grid points/cells per dimension
dx = Lx/Nx # Domain resolution
C = 0.5    # CFL number (0,1)

# %%
## Solve for binary and disk values
r_Hill  = r * (Msystem/(3*(Msun+Msystem)))**(1/3) # [cm]; Mutual hill radius of binary
bin_sep = r_Hill * Hill_frac                      # [cm]; Separation of binary at apoapsis
c_s = np.sqrt(gamma*T*kb/mmol/amu)                # [cm/s]; local sound speed

# %%
## Solve for code units of length, mass, and time based on the binary specifications 
unit_length = bin_sep # [cm]; initial binary separation defines the unit length
print("Unit length: ",f"{unit_length:.5e}","cm,",f"{unit_length/AU:.5e}",'AU')

Omega_bin = sqrt((G*Msystem) / (bin_sep)**3) # [Hz]; Keplerian frequency of binary
unit_time = 1 / Omega_bin                    # [s]; Unit time set by orbital period of binary; 2*pi = 1 orbital period
Omega_sun = sqrt((G*Msun) / (r)**3)          # [Hz]; Keplerian frequency around Sun
print("Unit time: ",f"{unit_time:.5e}","s,",f"{unit_time/(60*60*24):.5e}",'days,',f"{unit_time/yr:.5e}",'years')

unit_velocity = unit_length/unit_time
v_circ = Omega_bin * bin_sep # Binary orbital velocity
print("Unit velocity: ",f"{unit_velocity:.5e}","cm/s")
print("Circular velocity ",f"{v_circ:.5e}","cm/s")

unit_mass = Msystem # [g]; Mass of the binary system defines a unit mass; mass fractions follow naturally in these units
print("Unit mass: ",f"{unit_mass:.5e}","g")

# %%
## Solve for configuration parameters in code units
print("vvv Configuration Parameters vvv")

cs_code = c_s / unit_velocity # Sound speed in code units
print("Sound Speed (code): ",f"{cs_code:.5e}")

Omega_code = Omega_sun / Omega_bin # Heliocentric Keplerian frequency in code units
print("Omega (code): ",f"{Omega_code:.5e}")

H_code = cs_code / Omega_code # Scale height in code units
print("H (code): ",f"{H_code:.5e}")

sep_code = bin_sep/unit_length # Binary separation in code units
#print("Planetesimal separation (code): ",f"{sep_code:.5e}")
print("Semi-major axis (code): ", f"{sep_code/2:.5e}") # Distant from center of mass for equal mass system

nu = alpha*cs_code*H_code
print("Viscosity (code): ",f"{nu:.5e}")

# %%
## Solve for diagnostic parameters in code units
print("vvv Diagnostic Parameters vvv")

r_hill_code = r_Hill / unit_length # Hill radius in code units
print("R_Hill (code): ",f"{r_hill_code:.5e}")

msun_code = Msun / unit_mass # Heliocentric distance in code units
print("Solar mass (code): ",f"{msun_code:.5e}")

x1 = -(Mplanet2)/(Msystem)*sep_code # Position of particle 1 at apoapsis
x2 = (Mplanet1)/(Msystem)*sep_code  # Position of particle 2 at apoapsis
print("Initial positions along x-axis:")
print("x1 =",f"{x1:.4e}")
print("x2 =",f"{x2:.4e}")

v_rel      = sqrt(G*(Msystem)*((2/bin_sep)-((1+e)/bin_sep))) # Relative velocity of planetesimals 
v_rel_code = v_rel/unit_velocity
v1         = -(Mplanet2)/(Msystem)*v_rel_code # Velocity of particle 1 at apoapsis 
v2         = (Mplanet1)/(Msystem)*v_rel_code  # Velocity of particle 2 at apoapsis 
print("Initial tangential velocities:")
print("v1 =",f"{v1:.5e}")
print("v2 =",f"{v2:.5e}")

GM1 = G * Mplanet1 # [cm^3 s^-2]; GM of planetesimal 1 in physical units
GM1_code = GM1/((unit_length**3) * (unit_time**-2)) # Convert from physical to code units (cm^3 s^-2 in denominator)
G_code = GM1_code / (Mplanet1/unit_mass) # Divide by established mass of planetesimal 1 (somewhere between 0 and 1)
print("Gravitational constant (solving for M_sys=1): ",f"{G_code:.5e}")

dv = beta * cs_code/2 # Random velocity of dust
rb = 1/dv**2          # Bondi length 
tb = rb/dv            # Bondi time
ts = tb               # Setting the Bondi time equal to the friction time (ideally equal for accretion)
r_acc = 2*np.sqrt(rb*dv*ts) # Accretion radius
print("Friction time (code): ",f"{ts:.5e}")
print("Stokes number: ",f"{ts*Omega_code:.5e}")
print("r_acc/dx: ",r_acc/dx)
print("r_acc/Lx: ",r_acc/Lx)

dt_visc = C * dx**2 / nu   # Viscous timestep
dt_cs   = C * dx / cs_code # Advective timestep
print("Viscous timestep (code): ",f"{dt_visc:.5e}")
print("Advective timestep (code): ",f"{dt_cs:.5e}")
