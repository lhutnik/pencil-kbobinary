import numpy as np
import scipy as sp
import math
import pylab as plt

frac_hill=0.1

G=6.67384e-8
Msun=2e33 
AU=1.49597871e13
sigma_mol = 2e-15
amu = 1.6605402e-24

#
# UT 
#

x1=20.6/2 * 1e5
y1=19.9/2 * 1e5
z1=9.4/2  * 1e5

x2=15.4/2 * 1e5
y2=13.8/2 * 1e5
z2=9.8/2  * 1e5

vol1=4./3*math.pi*x1*y1*z1
vol2=4./3*math.pi*x2*y2*z2

a1 = np.cbrt(x1*y1*z1)
a2 = np.cbrt(x2*y2*z2)

print 'equivalent radius (km)=',a1/1e5,a2/1e5

a_bullet = np.array([a1,a2])

rho_bullet = 0.5  #cgs

m1 = vol1 * rho_bullet
m2 = vol2 * rho_bullet

print 'masses=',m1,m2

# m1+m2 = 1 

mu_body = 4e10 # rigidity 4e9 N/m2; 4e10 g/cm/s2 
Rbody = np.cbrt(a1**3+a2**3)
Mbody = m1+m2

klove =	3./2 / (1. + 19.*mu_body*Rbody/(2.*G*Mbody*rho_bullet))

print 'klove=',klove
Q = 100


r=45*AU

rhill = r*((m1+m2)/(3*Msun))**(1./3)
#r=20*AU
#r=8*AU

print 'r/AU=',r/AU
print 'rhill/AU, rhill/km=',rhill/AU,rhill/1e5

lnOmega2 = np.log(G*Msun)-3.*np.log(r) 
Omega = np.sqrt(np.exp(lnOmega2))

Rgas=8.314e7
gamma=1.4
gamma1=1./gamma
mmol=2.3
Rgasmu=Rgas/mmol
cp=gamma*Rgasmu/(gamma-1)
cv=cp/gamma

#
#  Minimum Mass Solar Nebula
#

Sigma0 = 1700.
T0=280.
p = 1.5
q = 0.5

Temperature = T0 * (r/AU)**(-q)
Sigma_gas = Sigma0 * (r/AU)**(-p)

cs2=Temperature*cp*(gamma-1.)

cs=np.sqrt(cs2)
H = cs/Omega

print 'Sigma_gas =',Sigma_gas,'  g/cm^2'
print 'T =',Temperature,'  K'
print 'cs=',cs,' cm/s'
print 'H=',H,' cm'
print 'H/AU=',H/AU
print 'H/r=',H/r

rho_g = Sigma_gas/(np.sqrt(2*math.pi)*H)

print 'rho_gas =',rho_g,' g/cm^3'

mfp = mmol*amu/(rho_g * sigma_mol)

print 'lambda (mean free path) =',mfp,' cm'

Kn = mfp/(2*a_bullet)

print 'Knudsen number=',Kn

knn=3*Kn

dlnpdlnr = p + q/2. + 1.5  
eta = 1./2 * dlnpdlnr * (H/r)**2
vk = Omega*r
print 'Center of mass velocity km/s',vk/1e5
ugas = eta * vk

a_bin = frac_hill*rhill
g0=G*(m1+m2)
period_bin = np.sqrt(4*math.pi**2/g0 * a_bin**3)
Omega_bin = 2*math.pi/period_bin
Omega_sun = np.sqrt(G*Msun/r**3)

vbin=Omega_bin*a_bin
print 'v_bin=',vbin,' cm/s'

uu = ugas+vbin

Ma_cm = ugas/cs

print 'headwind velocity =',ugas,' cm/s'
print 'Mach number=',Ma_cm

mu_visc = np.sqrt(8./math.pi) * rho_g * cs * mfp/3. 

print 'dynamical viscosity=',mu_visc,' g/cm/s'

Rey = 2*a_bullet*rho_g*uu/mu_visc 

print 'Reynolds number=',Rey

Cdstk = 24./Rey *(1+0.27*Rey)**0.43 + 0.47*(1-np.exp(-0.04*Rey**0.38)) 

rho_gas = Sigma_gas / np.sqrt(2*math.pi) / H 

tau = 4*mfp*rho_bullet/(3*rho_gas*Cdstk*cs) * 1./(Ma_cm*Kn)



#tau = 1./3 * a_bullet**2 * rho_bullet / cs * np.sqrt(math.pi/2.) * sigma_mol/ (mmol * amu) 

#St = np.sqrt(32*math.pi)/(knn*Ma_cm) * mfp * rho_bullet / Sigma_gas * (knn+1)**2/Cdstk

yr = 3.1e7
print 'period (yr) =',period_bin/yr
print 'tau/(1e7yr)=',tau/yr/1e7
print 'Stokes numbers rel to binary orbit (/1e7)=',tau*Omega_bin/1e7
print 'Stokes numbers rel to solar orbit (/1e6)=',tau*Omega_sun/1e6

#print 'St/1e5=',St/1e5

#
# Effective one body 
#
tau1=tau[0]
tau2=tau[1]
tau_eff=tau1*tau2*(m2+m1)/(tau2*m2 + tau1*m1)
tw = tau1*tau2/(tau1-tau2)
ugas_eff=ugas * tau_eff/tw
print 'tau_eff/(Myr):', tau_eff/yr/1e6
print 'tau_wind/(Myr):', tw/yr/1e6
print 'ugas (code units):', ugas/vbin
print 'ugas_eff (code units):', ugas_eff/vbin
print 'effective tau/1e7 (code units):',tau_eff*Omega_bin/1e7
print 'wind reduction factor',tau_eff/tw

print 'Omega_sun/Omega_bin=',Omega_sun/Omega_bin
print 'Omega_coriolis (code)=',Omega_sun/Omega_bin


#tau_quadratic = 16/15. * a_bullet * rho_bullet / rho_g / ugas
