
import math 
import numpy as np 

G = 6.67e-8 
msun = 1.99e33
AU = 1.49e13
a = 44.5813998 * AU

d1 = 19e5 
d2 = 14e5 

r1 = d1/2.
r2 = d2/2.

rho = 1. 
pi = math.pi 

m1= rho * 4.*pi/3 * r1**3
m2= rho * 4.*pi/3 * r2**3

mu = G*(m1+m2)

rhill = a * (1./3 * (m1+m2)/msun)**(1./3)

print 'm1=',m1/1e3,' kg, m2=',m2/1e3,' kg, ratio sum/sun',(m1+m2)/msun
print 'rhill=',rhill/AU,' AU',rhill/1e5,' km'

r = rhill/10

print 'distance=',r 
print 'distance/radii',r/(r1+r2)

T = np.sqrt(4*pi**2/mu * r**3)

#hour = 3600.
#day  = 24.
#year = 365.

year = 3.154e7

print 'period, years=',T/year

Omega = 2*pi/T 

a1 = m2/m1 * r 
a2 = (1-m2/m1) * r 

v2 = Omega * a2 
v1 = Omega * a1 

print 'orbital velocities',v2,v1,' cm/s' 

period_sun = np.sqrt(4*pi**2/(G*msun) * a**3)

Omega_sun = 2*pi/period_sun

vk = Omega_sun * a 

print 'Keplerian velocity=',vk/1e5,' km/s'
eta = 1e-2 

vgas = eta*vk 

print 'relaltive has velocity=',vgas,' cm/s'

vcode = vgas/(v2+v1)

print 'gas code velocity=',vcode
