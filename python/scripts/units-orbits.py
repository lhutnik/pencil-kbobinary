# %%
## Import critical libraries 
import numpy as np                # Pi constant and math
from math import sqrt             # Math functions
import astropy.constants as const # Constants

# %%
## Constants
pi     = np.pi
AU     = const.au.cgs.value      # [cm]; Astronomical unit (AU)
Mearth = const.M_earth.cgs.value # [g]; Mass of Earth
Msun   = const.M_sun.cgs.value   # [g]; Mass of Sun
Mpluto = 1.309e25                # [g]; Mass of Pluto
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
rr = 20      # [AU]; Heliocentric distance of binary 
r  = rr * AU # [cm]; Heliocentric distance of binary
h  = 0.05    # Scale height ratio at distance rr
H  = h * r   # [cm]; Gas disk scale height 
Q  = 30      # Toomre Q; >1 for gravitationally stable region
alpha = 1e-3 # Shakura and Sunyaev alpha parameter

## BINARY SETTINGS
mass_ratio = 1       # Binary mass ratio (0,1]; M2/M1 = f
Hill_frac  = 0.3     # Fraction of mutual Hill radius for furthest separation (apoapsis); [1,40]%
e          = 0.0     # Eccentricity of mutual binary orbit [0,1)
aps1       = 0.01    # Size of sink particle 1 in code units
aps2       = 0.01    # Size of sink particle 2 in code units

Msystem  = 1e-2 * Mpluto              # [g]; Sum of binary component masses
Mplanet1 = Msystem / (1 + mass_ratio) # [g]; Mass of primary 
Mplanet2 = Msystem - Mplanet1         # [g]; Mass of secondary

## SIMULATION SETTINGS
l3D         = False # 3D or not 3D run (2D)
grid_size   = 8     # Full x or y grid size in code units; assumes equal x,y scale
grid_points = 128   # Number of grid points/cells per dimension

## SET CODE UNITS (set in start.in or other configuration files)
cs_code    = 1
rho_code   = 1
#Omega_code = 1   # 2*pi time unit = 1 orbit of the binary
#H_code     = 1   # Set automatically by Omega and cs in code (cs_code/Omega_code)
Msystem_code = 1 # Pointmasses module assumes the system mass sums to 1 in code units

# %%
## Solve for binary values
r_Hill  = r * (Msystem/(3*(Msun+Msystem)))**(1/3) # [cm]; Mutual hill radius of binary
bin_sep = r_Hill * Hill_frac                      # [cm]; Separation of binary at apoapsis

# %%
## Solve for code units of length, mass, and time based on the binary specifications 
#unit_length = H # [cm]; Unit length is set by scale height
unit_length = bin_sep # [cm]; initial binary separation defines the unit length
print("Unit length: ",f"{unit_length:.4e}","cm,",f"{unit_length/AU:.4e}",'AU')

Omega_bin   = sqrt((G*Msystem) / (bin_sep)**3) # [Hz]; Keplerian frequency of binary
unit_time   = 1 / Omega_bin                      # [s]; Unit time set by orbital period of binary; 2*pi = 1 orbital period
print("Unit time: ",f"{unit_time:.4e}","s,",f"{unit_time/(60*60*24):.4e}",'days,',f"{unit_time/yr:.4e}",'years')

unit_velocity = unit_length/unit_time
print("Unit velocity: ",f"{unit_velocity:.4e}","cm/s")
#print("Sound speed: ",f"{cs:.4e}","cm/s") # This should be similar to the unit velocity

unit_mass = Msystem # [g]; Mass of the binary system defines a unit mass; mass fractions follow naturally in these units
print("Unit mass: ",f"{unit_mass:.4e}","g")

# %%
## Solve for other variables given inputs 
Omega_sol = sqrt(G*Msun) / r**1.5   # [s^-1]; Local angular frequency around Sun
cs    = Omega_sol * H               # [cm/s]; Local sound speed
T     = cs**2 / (cp*(gamma-1))      # [K]; Local gas temperature

## Solve for binary values in code units
M1_code = Msystem_code / (1 + mass_ratio) # [g]; Mass of primary 
M2_code = Msystem_code - M1_code  # [g]; Mass of secondary

## Solve for Omega around the central star in code units 
Omega_code = Omega_sol * unit_time # Omega / unit frequency
print("Omega (code): ",f"{Omega_code:.4e}")

## Solve for H in code length units
H_code = H / unit_length
print("H (code): ",f"{H_code:.4e}")

## Solve for Hill radius in code length units
r_hill_code = r_Hill / unit_length
print("R_Hill (code): ",f"{r_hill_code:.4e}")

## Solve for sound speed in code velocity units
cs_code = cs / unit_velocity
print("Sound Speed (code): ",f"{cs_code:.4e}")

## Solve for heliocentric distance in code length units
r_code = r / unit_length
print("Heliocentric r (code): ",f"{r_code:.4e}")

## Solve for heliocentric distance in code length units
msun_code = Msun / unit_mass
print("Solar mass (code): ",f"{msun_code:.4e}")

## Solve for fraction of Hill radius to smooth potential out to (covers 4 grid cells at least)
r_Hill1_code = r_code * (M1_code/(3*(M1_code+msun_code)))**(1/3)
r_Hill2_code = r_code * (M2_code/(3*(M2_code+msun_code)))**(1/3)
required_cells = 4 # Cells to cover over 1 smoothing radius
smooth_r = (grid_size/grid_points) * required_cells # Convert grid coverage to unit length
smooth_f1 = smooth_r / r_Hill1_code # Divide by r_hill to get fraction
smooth_f2 = smooth_r / r_Hill2_code # Divide by r_hill to get fraction
print("Smoothed potential Hill radius fractions: ",f"{smooth_f1:.4e}",f"{smooth_f2:.4e}")
print("Smoothed potential radius (code): ",f"{smooth_r:.4e}")

## Solve for rho_0 in code units 
Sigma = cs*Omega_sol/(G*Q*pi) # Determine Sigma (surface density) and rho (density) from Q
rho = Sigma/(sqrt(2*pi)*H)
rho_code = rho / (unit_mass * unit_length**(-3))
print("rho_0 (code): ",f"{rho_code:.4e}")

# %%
## Determine Sigma (surface density) and rho (density) from Q
Sigma = cs*Omega_sol/(G*Q*pi)

## Solve for unit mass
if (l3D): # If solving for 3D case,
    rho = Sigma/(sqrt(2*pi)*H)
    unit_volume_density = rho/rho_code
    Sigma_code = rho_code * sqrt(2*pi) * H_code
#    unit_column_density = Sigma/Sigma_code
#    unit_mass = unit_volume_density * unit_length**3
else: # If solving for 2D case,
    Sigma_code = rho_code # In 2D, rho is the column density
    unit_column_density = Sigma/Sigma_code
    rho = Sigma/(sqrt(2*pi)*H) # Solve for density in physical units
#    
#    ## Solve for unit mass given unit volume density
#    actual_rho_code = Sigma_code/(sqrt(2*pi)*H_code)
#    unit_volume_density = rho/actual_rho_code
#    unit_mass = unit_column_density * unit_length**2
#print("Unit mass: ",f"{unit_mass:.4e}","g")

# %%
## Solve for the gravitational constant in code units given the system mass must equal 1
GM1 = G * Mplanet1 # [cm^3 s^-2]; GM of planetesimal 1 in physical units
GM1_code = GM1/((unit_length**3) * (unit_time**-2)) # Convert from physical to code units (cm^3 s^-2 in denominator)
G_code = GM1_code / M1_code # Divide by established mass of planetesimal 1 (somewhere between 0 and 1)
print("Gravitational constant (solving for M_sys=1): ",f"{G_code:.6e}")

## Check whether another method of finding GM results in the same value
M1_code_u, M2_code_u = Mplanet1/unit_mass, Mplanet2/unit_mass # Code unit mass of planetesimals
#print("M1 (code): ",f"{M1_code_u:.4e}")
#print("M2 (code): ",f"{M2_code_u:.4e}")

G_code_u = cs_code * Omega_code / (Q*pi*Sigma_code) # Gravitational constant solved by taking other code units
#print("Gravitational constant (code alone): ",f"{G_code_u:.6e}")
#print("Gravitational constant (dividing GM_code/(M_sys=1)): ",f"{M1_code_u*G_code_u/M1_code:.6e}")


## Compare methods to find GM in code units
M1_code_u, M2_code_u = Mplanet1/unit_mass, Mplanet2/unit_mass # Code unit mass of planetesimals
#print("G*M1 (code only): ",f"{M1_code_u*G_code_u:.6e}")
#print("G*M1 (M_sys=1): ",f"{M1_code*G_code:.6e}")

#print("G*M2 (code only): ",f"{M2_code_u*G_code_u:.6e}")
#print("G*M2 (M_sys=1): ",f"{M2_code*G_code:.6e}")

# %%
## Solve for initial binary separation in code units
sep_code = bin_sep/unit_length # Binary separation in code units 
print("Planetesimal separation (code): ",f"{sep_code:.6e}")
print("Semi-major axis (code): ", f"{sep_code/2:.6e}") # Distant from center of mass for equal mass system

# %%
## Solve for initial positions and velocities assuming apoapsis 
## Particle positions along common axis at apoapsis
x1 = -(Mplanet2)/(Msystem)*sep_code # Position of particle 1 at apoapsis
x2 = (Mplanet1)/(Msystem)*sep_code  # Position of particle 2 at apoapsis
print("Initial positions along axis:")
print("x1 =",f"{x1:.4e}")
print("x2 =",f"{x2:.4e}")

## Calculate relative velocity of planetesimals 
v_rel      = sqrt(G*(Msystem)*((2/bin_sep)-((1+e)/bin_sep)))
v_rel_code = v_rel/unit_velocity
v1         = -(Mplanet2)/(Msystem)*v_rel_code # Velocity of particle 1 at apoapsis 
v2         = (Mplanet1)/(Msystem)*v_rel_code  # Velocity of particle 2 at apoapsis 
print("Initial tangential velocities:")
print("v1 =",f"{v1:.4e}")
print("v2 =",f"{v2:.4e}")

# %%
## Calculate rhopswarm given planetesimal mass (sink particle only case)
## Option 1: mass set by size of a single grid cell
#rhopswarm1 = M1_code * (grid_size/grid_points)**(-3)
#rhopswarm2 = M2_code * (grid_size/grid_points)**(-3)
#print("If set by grid scale, then rhopswarm of particles:")
#print("rhopswarm1 =",rhopswarm1)
#print("rhopswarm2 =",rhopswarm2)
## Option 2: mass set by size of sink particles (aps)
#rhopswarm1 = M1_code * (3/(4*pi)) * (aps1/2)**(-3)
#rhopswarm2 = M2_code * (3/(4*pi)) * (aps2/2)**(-3)
#print("If set by particle size, then rhopswarm of particles:")
#print("rhopswarm1 =",rhopswarm1)
#print("rhopswarm2 =",rhopswarm2)

# %%
cfl = (grid_size/grid_points)/cs_code # Minimum timestep to remain at a Courant number of 1
print("CFL Timestep: ",f"{cfl:.6e}")

nu = alpha*cs_code*H_code
print("Viscosity required (code): ",f"{nu:.6e}")


