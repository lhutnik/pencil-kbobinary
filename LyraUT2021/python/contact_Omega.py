import numpy as np
import scipy as sp
import math
import pylab as plt
import pencil as pc
from mpl_toolkits.axes_grid1 import make_axes_locatable

G=6.67384e-8
Msun=1.99e33 
AU=1.49597871e13
sigma_mol = 2e-15
amu = 1.6605402e-24
pi=math.pi

nr=100

rsun=np.logspace(1.,2.,nr)*AU

x1=20.6/2 * 1e5
y1=19.9/2 * 1e5
z1=9.4/2  * 1e5

x2=15.4/2 * 1e5
y2=13.8/2 * 1e5
z2=9.8/2  * 1e5

vol1=4./3*pi*x1*y1*z1
vol2=4./3*pi*x2*y2*z2

a1 = np.cbrt(x1*y1*z1)
a2 = np.cbrt(x2*y2*z2)

abullet = np.array([a1,a2])

rho_bullet = 0.5

m1 = vol1 * rho_bullet
m2 = vol2 * rho_bullet
g0 = G*(m1+m2)

rhill = rsun * (1./3 * (m1+m2)/Msun)**(1./3)

#print 'm1=',m1/1e3,' kg, m2=',m2/1e3,' kg, ratio sum/sun',(m1+m2)/Msun
#print 'rhill=',rhill/AU,' AU',rhill/1e5,' km'

datadir='/Users/wlyra/pencil-code/f90/pencil-wlyra/UltimaThule/Apr19_NBody_Wind100_Mass21/data'
par=pc.read_param(datadir=datadir)
par2=pc.read_param(param2=True,datadir=datadir)
if (par.iprimary == 1):
    iprimary = 0
    isecondary = 1
else:
    iprimary = 1
    isecondary = 0

lnOmega2 = np.log(G*Msun)-3.*np.log(rsun) 
Omega_sun = np.sqrt(np.exp(lnOmega2))

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

p=1.5
q=0.5

Temperature = 280 * (rsun/AU)**(-q)
Sigma_gas = 1700 * (rsun/AU)**(-p)

cs2=Temperature*cp*(gamma-1.)

cs=np.sqrt(cs2)
H = cs/Omega_sun

#print 'Sigma_gas =',Sigma_gas,'  g/cm^2'
#print 'T =',Temperature,'  K'
#print 'cs=',cs,' cm/s'
#print 'H=',H,' cm'
#print 'H/AU=',H/AU
#print 'H/r=',H/rsun

rho_g = Sigma_gas/(np.sqrt(2*pi)*H)

#print 'rho_g =',rho_g,' g/cm^3'

mfp = mmol*amu/(rho_g * sigma_mol)

#print 'lambda (mean free path) =',mfp,' cm'

na = 2

Kn = np.zeros([nr,na])
for i in range(2): 
    Kn[:,i] = mfp/(2*abullet[i])

#print 'Knudsen number=',Kn

knn=3*Kn

dlnpdlnr = p + q/2. + 1.5  
eta = 1./2 * dlnpdlnr * (H/rsun)**2
vk = Omega_sun*rsun
ugas = eta * vk

Ma = ugas/cs

#print 'headwind velocity =',ugas,' cm/s'
#print 'Mach number=',Ma

mu_visc = np.sqrt(8./pi) * rho_g * cs * mfp/3. 

#print 'dynamical viscosity=',mu_visc,' g/cm/s'

Rey = np.zeros([nr,na])
for i in range(2): 
    Rey[:,i] = 2*abullet[i]*rho_g * ugas/mu_visc 

#

print 'Reynolds number=',Rey

Cdstk = 24/Rey + 3.6*Rey**(-0.313)

tau = np.zeros([nr,na])
for i in range(2):
    tau[:,i] = 4*mfp*rho_bullet/(3*rho_g*Cdstk[:,i]*cs) * 1./(Ma*Kn[:,i])
    #St[:,i] = np.sqrt(32.*pi)/(knn[:,i]*Ma) * mfp * rho_bullet / Sigma_gas * (knn[:,i]+1)**2/Cdstk[:,i]

#tau=St
tau_eff=tau[:,0]*tau[:,1]*(m2+m1)/(tau[:,1]*m2 + tau[:,0]*m1)
tw = tau[:,0]*tau[:,1]/(tau[:,0]-tau[:,1])
ugas_eff=ugas * tau_eff/tw

#plt.show()

nb=99
a_bin=np.zeros([nr,nb])
binary_separation = np.logspace(-4,0,nb)
for ir in range(nr):
    a_bin[ir,:] = binary_separation * rhill[ir]

#print 'distance=',a_bin
#print 'distance/radii',a_bin/(r1+r2)

period_bin = np.sqrt(4*pi**2/g0 * a_bin**3)
Omega_bin = 2*pi/period_bin

vbin = Omega_bin*a_bin

Omega_sun_Omega_bin = np.zeros([nr,nb])

for ib in range(nb):
    Omega_sun_Omega_bin[:,ib] = Omega_sun/Omega_bin[:,ib]

#period_sun = 2*pi/Omega_sun
yr=3.154e7
#print 'contact time=',contact
#plt.plot(rsun/AU,contact[:,50]/yr)
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(12,5))
img = ax1.contourf(rsun/AU,binary_separation,np.log10(np.transpose(Omega_sun_Omega_bin)),256,cmap="inferno")

#CS4=ax1.contour(rsun/AU,binary_separation,np.log10(np.transpose(contact)/yr),levels=[4.5,5,5.5,6,6.5],colors='white',linewidths=1)
#fmt = {}
#strs = [ r'$3\times 10^4$', r'$10^5$', r'$3\times 10^5$', r'$10^6$', r'$3\times 10^6$']
#for l,s in zip( CS4.levels, strs ):
#    fmt[l] = s
#ax1.clabel(CS4, fmt=fmt, colors='white', fontsize=14)


#robject = 2e6
#ax1.plot(rsun/AU,robject/rhill,color='black')
#ax1.annotate('20 km', [11,2e-3],rotation=-8,size=14,color='black')

#ax1.set_title("Binary hardening timescale")

ax1.set_xscale("log")
ax1.set_yscale("log")

ax1.set_xlabel(r'$r$ (AU)')
ax1.set_ylabel(r'$a_{\rm bin}/r_{\rm Hill}$')

divider = make_axes_locatable(ax1)
cim     = divider.append_axes("right", size = "2%", pad = 0.3)
#plt.colorbar(img,cax=cim,orientation="vertical",label=r'$\log_{10} [t_{\rm inspiral}$ (yr)]',ticks=[4.5,5.0,5.5,6.0,6.5])
plt.colorbar(img,cax=cim,orientation="vertical",label=r'$\log_{10} [\Omega_s/\Omega_b]$')#,ticks=[4.5,5.0,5.5,6.0,6.5])

awc = ugas*(1/tau[:,1]-1/tau[:,0])
rwish = np.sqrt(G*(m1+m2)/awc)
print rwish/rhill

f=binary_separation

ir=65

coriolis = np.zeros([nr,nb])
drag     = np.zeros([nr,nb])

for ir in range(100):
    coriolis[ir,:] = Omega_sun[ir] * vbin[ir,:]
    #gravity  = Omega_bin
    drag[ir,:]     = vbin[ir,:] / tau_eff[ir]

img = ax2.contourf(rsun/AU,binary_separation,np.log10(np.transpose(coriolis/drag)),256,cmap="inferno")

divider = make_axes_locatable(ax2)
cim2     = divider.append_axes("right", size = "2%", pad = 0.3)
plt.colorbar(img,cax=cim2,orientation="vertical")

#ax2.set_xscale("log")
#ax2.set_yscale("log")

plt.tight_layout()
plt.show()
#plt.savefig("hardening_timescale.png")
