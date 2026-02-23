import numpy as np
import scipy as sp
import math
import pylab as plt

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=MEDIUM_SIZE) # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE) # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE) # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE) # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE) # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE) # legend fontsize
plt.rc('figure', titlesize=SMALL_SIZE) # fontsize of the figure title


frac_hill=0.1

G=6.67384e-8
Msun=2e33 
AU=1.49597871e13
sigma_mol = 2e-15
amu = 1.6605402e-24
Rgas=8.314e7
gamma=1.4
gamma1=1./gamma
mmol=2.3
Rgasmu=Rgas/mmol
cp=gamma*Rgasmu/(gamma-1)
cv=cp/gamma
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

rho_bullet = 0.5

m1 = vol1 * rho_bullet
m2 = vol2 * rho_bullet

print 'masses=',m1,m2

nr=100
rr=np.logspace(-0.5,2.5,nr)

r=rr*AU

rhill = r*((m1+m2)/(3*Msun))**(1./3)

sma = rhill*frac_hill
Omega_bin = np.sqrt(G*(m1+m2))/sma**1.5
vbin = Omega_bin*sma

lnOmega2 = np.log(G*Msun)-3.*np.log(r) 
Omega = np.sqrt(np.exp(lnOmega2))
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

rho_g = Sigma_gas/(np.sqrt(2*math.pi)*H)
print 'rho_gas =',rho_g,' g/cm^3'

mfp = mmol*amu/(rho_g * sigma_mol)
print 'lambda (mean free path) =',mfp,' cm'

Kn=np.zeros([nr,2])
for i in range(2):
    Kn[:,i] = mfp/(2*a_bullet[i])

#dlnpdlnr = p + q/2. + 1.5  
#eta = 1./2 * dlnpdlnr * (H/r)**2
#vk = Omega*r

mu_visc = np.sqrt(8./math.pi) * rho_g * cs * mfp/3. 
rho_gas = Sigma_gas / np.sqrt(2*math.pi) / H 

ng=5
ugas = np.array([0,1e-3,1e-2,5e-2,1e-1])*1e5

tau=np.zeros([nr,ng,2])

yr = 3.15576e7
for ig in range(ng):
    uu = vbin + ugas[ig]        
    Ma=uu/cs
    for i in range(2): 
        Rey = 2*a_bullet[i]*rho_g*uu/mu_visc 
        Cdstk = 24./Rey *(1+0.27*Rey)**0.43 + 0.47*(1-np.exp(-0.04*Rey**0.38)) 
        tau[:,ig,i] = 4*mfp*rho_bullet/(3*rho_gas*Cdstk*cs) * 1./(Ma*Kn[:,i])        
#
# Effective one body 
#
    tau1=tau[:,ig,0]
    tau2=tau[:,ig,1]
    tau_eff=tau1*tau2*(m2+m1)/(tau2*m2 + tau1*m1)

    if (ig == 0 or ig== 3):
        plt.plot(r/AU,tau_eff/yr/1e6,label=np.str(np.int(ugas[ig]/1e2))+' m/s')
    else:
        plt.plot(r/AU,tau_eff/yr/1e6,label=np.str(np.int(ugas[ig]/1e2))+' m/s',linestyle='--')

plt.xscale('log')
plt.xlim([0.8,200])
plt.ylim([1e-3,1e2])
plt.yscale('log')
plt.legend(loc='lower right',shadow=True, fancybox=True)

plt.xlabel(r'$r$ (AU)')
plt.ylabel(r'$\tau_{\rm eff}$ (Myr)')
plt.title('Drag time vs distance')

#plt.axvline(45,linestyle=':',color='grey')
#plt.axhline(10,linestyle=':',color='grey')

plt.savefig("StoppingTimeDistance.png")
#plt.show()
