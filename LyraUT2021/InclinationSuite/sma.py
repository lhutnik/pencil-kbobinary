def quantities_individual(x,y,z,vx,vy,vz,OO):

    rad  =  np.sqrt(x**2 + y**2 + z**2)
    v2   =  vx**2 + vy**2 + vz**2
    sma  =  1./(2./rad - v2)
    hx   =  vz*y - vy*z #- OO*x*z
    hy   =  vx*z - vz*x #- OO*y*z
    hz   =  vy*x - vx*y #- OO*(x**2+y**2)
    hh   =  np.sqrt(hx**2+hy**2+hz**2)
    ecc  =  np.sqrt(1 - hh**2/sma)
    incl =  180/math.pi * np.arccos(hz/hh)

    return rad,hz,hh,sma,ecc,v2,incl

import numpy as np
import math
import pylab as plt
import pencil as pc

fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(13,4))

ts=pc.read_ts(datadir="./Incl00Ecc00/data")
par2=pc.read_param(param2=True,datadir="Incl00Ecc00/data")
OO=par2.omega_rotating
tau=par2.stokesnumber[1]
ugas=par2.ugas
rad,hz,hh,sma,ecc,v2,incl = quantities_individual(ts.xq2,ts.yq2,ts.zq2,ts.vxq2,ts.vyq2,ts.vzq2,OO) 

ax1.plot(ts.t,sma,color='black',label='Numerical')
sma_solution = sma[0]*np.exp(-2*ts.t/tau)
ax1.plot(ts.t,sma_solution,linestyle=':',color='red',label='Analytical')

mu=1. # G*(m1+m2)
ax2.plot(ts.t,ecc,color='black')#,label='Numerical')
ax2.set_ylim([0.,0.3])
ecc_solution = np.sin(3./2*ugas*np.sqrt(sma[0]/mu)*(1-np.exp(-ts.t[0:2000]/tau)))
ax2.plot(ts.t[0:2000],ecc_solution[0:2000],linestyle=':',color='green')#,label='Numerical')

#ecc_solution = np.sin(3./2*ugas*np.sqrt(sma[0]/mu)*(1-np.exp(-ts.t/tau)))
#ecc_solution = np.cos(3./2*ugas*np.sqrt(sma[0]/mu)*(1-np.exp(-ts.t/tau)) + np.arccos(ecc[0]))
#ax2.plot(ts.t,ecc_solution,linestyle=':',color='red',label='Analytical No Coriolis')

left, bottom, width, height = [0.515, 0.59, 0.135, 0.3]
sx2 = fig.add_axes([left, bottom, width, height])
sx2.plot(ts.t[0:1000],ecc[0:1000],color='black',label='Numerical')
ecc_solution = np.sin(3./2*ugas*np.sqrt(sma[0]/mu)*(1-np.exp(-ts.t[0:1000]/tau)))
sx2.plot(ts.t[0:1000],ecc_solution,linestyle=':',color='green',label='Analytical '+r' $\Omega_\odot=0$')
sx2.set_ylim([0.,0.4])
sx2.legend()

ax3.plot(ts.t,hh,color='black',label='Numerical')
ax3.plot(ts.t,hh[0]*np.exp(-ts.t/tau),linestyle=':',color='red',label='Analytical')
hh_solution = np.exp(-ts.t[0:2000]/tau) * (hh[0]-1 + np.cos(3./2*sma[0]*ugas*(1-np.exp(-ts.t[0:2000]/tau))))  
ax3.plot(ts.t[0:2000],hh_solution,linestyle=':',color='green',label='Analytical '+r' $\Omega_\odot=0$')
left, bottom, width, height = [0.755, 0.225, 0.135, 0.3]
sx3 = fig.add_axes([left, bottom, width, height])
sx3.plot(ts.t[0:1000],hh[0:1000],color='black')
sx3.plot(ts.t[0:1000],hh_solution[0:1000],linestyle=':',color='green') #,label='Analytical '+r' $\Omega_\odot=0$')
sx3.set_ylim([0.965,1.0])
#sx3.legend()


#ax3.plot(ts.t,hh_solution,linestyle=':',color='red',label='Analytical No Coriolis')

ax1.set_xlabel(r'$t$ (yr)')
ax1.set_ylabel(r'$a$')
ax1.set_title("Semimajor Axis")
ax1.legend()
#ax1.set_ylim([0,1])

ax2.set_xlabel(r'$t$ (yr)')
ax2.set_ylabel(r'$e$')
ax2.set_title("Eccentricity")
#ax2.legend()
#ax2.set_ylim([0,1])

ax3.set_xlabel(r'$t$ (yr)')
ax3.set_ylabel(r'$h$')
ax3.set_title("Angular Momentum")
ax3.legend()
ax3.set_ylim([0.15,1])

plt.tight_layout()
plt.show()
#plt.savefig("Binary.png")
