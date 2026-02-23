
import math
import numpy as np

yr  = 3.1e7
#day = 3600*23.93419

##Earth-Moon

#period = 27.321661*day
#k2p=0.299
#Qp=12
#mp=5.9736e24    #earths' mass (kg)
#cp=6.378e6   #earth's radius (m)

#k2s=0.030 
#Qs=27
#ms=7.249e22
#cs=1.73753e6     #moon's radius (m)

#a=3.844e8  #orbital sma moon
#n=2*math.pi/period #mean motion moon

## Moon tide 
#adot_moon = 3*k2s/Qs * mp/ms * (cs/a)**5 * n * a
#taus = a/adot_moon
#print 'earth tide on Moon=',taus/yr/1e7,' 1e7 yr'

## Earth tide  
#adot_earth = 3*k2p/Qp * ms/mp * (cp/a)**5 * n * a
#taup = a/adot_earth
#print 'moon tide on earth=',taup/yr/1e10,' 1e10 yr'

## Retrograde (UT);
G=6.67e-11
mu=4e9
M=1.01e15
A=7.8e3
gc=G*M/A**2
rho=0.5e3
mu_tilde=19.*mu/(2*rho*gc*A)
k2=1.5/(1+mu_tilde)
Q=100
m1=M
m2=m1/2.
c1=A
a=6000e3
n = np.sqrt(G*(m1+m2))/a**1.5

adot = 3.*k2/Q * m1/m2 * (c1/a)**5 * a*n
tau = a/adot
print 'tide on Mu69=',tau/yr/1e6,' Myr'
