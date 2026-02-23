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


def quantities_pair(x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2,OO):

    x=x1-x2
    y=y1-y2
    z=z1-z2
    vx=vx1-vx2
    vy=vy1-vy2
    vz=vz1-vz2

    rad  =  np.sqrt(x**2 + y**2 + z**2)
    v2   =  vx**2 + vy**2 + vz**2
    sma  =  1./(2./rad - v2)
    hx   =  vz*y - vy*z
    hy   =  vx*z - vz*x
    hz   =  vy*x - vx*y
    hh   =  np.sqrt(hx**2+hy**2+hz**2)
    ecc  =  np.sqrt(1 - hh**2/sma)
    incl =  180/math.pi * np.arccos(hz/hh)
    return rad,hz,hh,sma,ecc,v2,incl



import numpy as np
import math
import pylab as plt
import pencil as pc

fig, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(1,5,figsize=(12,3))

ts1=pc.read_ts(datadir="./Incl90Ecc03/data")
ts2=pc.read_ts(datadir="./Incl90Ecc03_nogas/data")
ts3=pc.read_ts(datadir="./Incl90Ecc03_nogas_masses/data")
par2=pc.read_param(param2=True,datadir="Incl90Ecc03/data")

OO=par2.omega_rotating
tau=par2.stokesnumber[1]
ugas=par2.ugas
rad1,hz1,hh1,sma1,ecc1,v21,incl1 = quantities_individual(ts1.xq2,ts1.yq2,ts1.zq2,ts1.vxq2,ts1.vyq2,ts1.vzq2,OO) 
rad2,hz2,hh2,sma2,ecc2,v22,incl2 = quantities_individual(ts2.xq2,ts2.yq2,ts2.zq2,ts2.vxq2,ts2.vyq2,ts2.vzq2,OO) 
rad3,hz3,hh3,sma3,ecc3,v23,incl3 = quantities_pair(ts3.xq1,ts3.yq1,ts3.zq1,ts3.vxq1,ts3.vyq1,ts3.vzq1,
                                                   ts3.xq2,ts3.yq2,ts3.zq2,ts3.vxq2,ts3.vyq2,ts3.vzq2,OO)

time1 = ts1.t/tau
time2 = ts2.t/tau
time3 = ts3.t/tau

ax1.plot(time1,sma1,color='red',label='Gas')
ax1.plot(time2,sma2,color='black',label='No Gas')

#sma_solution = sma1[0]*np.exp(-2*ts1.t/tau)
#ax1.plot(time1,sma_solution,linestyle=':',color='red',label='Analytical with Gas')

mu=1. # G*(m1+m2)
ax2.plot(time1,ecc1,color='red',label='Gas')
ax2.plot(time2,ecc2,color='black',label='No Gas')
mu=1.0

ax3.plot(time1,hh1,color='red',label='Gas')
ax3.plot(time2,hh2,color='black',label='No Gas')

ax4.plot(time1,incl1,color='red',label='Gas')
ax4.plot(time2,incl2,color='black',label='No Gas')

ax5.plot(time1,np.sqrt(v21),color='red',label='Gas')
ax5.plot(time2,np.sqrt(v22),color='black',label='No Gas')

ax1.set_xlabel(r'$t/\tau$')
ax1.set_ylabel(r'$a/a_0$')
ax1.set_title("Semimajor Axis")
ax1.legend()

ax2.set_xlabel(r'$t/\tau$')
ax2.set_ylabel(r'$e$')
ax2.set_title("Eccentricity")
#ax2.legend()

ax3.set_xlabel(r'$t/\tau$')
ax3.set_ylabel(r'$h/h_0$')
ax3.set_title("Angular Momentum")
#ax3.legend()

ax4.set_xlabel(r'$t/\tau$')
ax4.set_ylabel(r'$I$')
ax4.set_title("Inclination")
#ax4.legend()

ax5.set_xlabel(r'$t/\tau$')
ax5.set_ylabel(r'$u/u_0$')
ax5.set_title("Velocity")
#ax5.legend()

plt.tight_layout()
#plt.show()
plt.savefig("GasNoGas.png")
#plt.savefig("GasNoGas_Masses.png")
