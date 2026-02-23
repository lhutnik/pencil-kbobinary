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

fig, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(10,6))

print "reading ecc 09"
ts=pc.read_ts(datadir="./Incl00Ecc09/data")
par2=pc.read_param(param2=True,datadir="Incl00Ecc00/data")
OO=par2.omega_rotating

rad,hz,hh,sma,ecc,v2,incl = quantities_individual(ts.xq2,ts.yq2,ts.zq2,ts.vxq2,ts.vyq2,ts.vzq2,OO) 
ax1.plot(ts.t,incl,color='red',label="e=0.9")
ax1.set_xlabel(r'$t$ (yr)')
ax1.set_ylabel(r'$i$')
ax1.set_title("Inclination")

ax2.plot(ts.t,ecc,color='red')
ax2.set_xlabel(r'$t$ (yr)')
ax2.set_ylabel(r'$e$')
ax2.set_title("Eccentricity")

ax3.plot(ts.t,hh,color='red')
ax3.set_xlabel(r'$t$ (yr)')
ax3.set_ylabel(r'$h$')
ax3.set_title("Total Ang Mom")

ax4.plot(ts.t,sma,color='red')
ax4.set_xlabel(r'$t$ (yr)')
ax4.set_ylabel(r'$a$')
ax4.set_title("Semimajor axis")

ax5.plot(ts.t,rad,color='red')
ax5.set_xlabel(r'$t$ (yr)')
ax5.set_ylabel(r'$r$')
ax5.set_title("Separation")

ax6.plot(ts.t,np.sqrt(v2),color='red')
ax6.set_xlabel(r'$t$ (yr)')
ax6.set_ylabel(r'$v$')
ax6.set_title("Velocity")

print "reading ecc 07"
ts=pc.read_ts(datadir="./Incl00Ecc07/data")
rad,hz,hh,sma,ecc,v2,incl = quantities_individual(ts.xq2,ts.yq2,ts.zq2,ts.vxq2,ts.vyq2,ts.vzq2,OO)
ax1.plot(ts.t,incl,color='blue',label="e=0.7")
ax2.plot(ts.t,ecc,color='blue')
ax3.plot(ts.t,hh,color='blue')
ax4.plot(ts.t,sma,color='blue')
ax5.plot(ts.t,rad,color='blue')
ax6.plot(ts.t,np.sqrt(v2),color='blue')

print "reading ecc 05"
ts=pc.read_ts(datadir="./Incl00Ecc05/data")
rad,hz,hh,sma,ecc,v2,incl = quantities_individual(ts.xq2,ts.yq2,ts.zq2,ts.vxq2,ts.vyq2,ts.vzq2,OO)
ax1.plot(ts.t,incl,color='green',label="e=0.5")
ax2.plot(ts.t,ecc,color='green')
ax3.plot(ts.t,hh,color='green')
ax4.plot(ts.t,sma,color='green')
ax5.plot(ts.t,rad,color='green')
ax6.plot(ts.t,np.sqrt(v2),color='green')

print "reading ecc 03"
ts=pc.read_ts(datadir="./Incl00Ecc03/data")
rad,hz,hh,sma,ecc,v2,incl = quantities_individual(ts.xq2,ts.yq2,ts.zq2,ts.vxq2,ts.vyq2,ts.vzq2,OO)
ax1.plot(ts.t,incl,color='orange',label="e=0.3")
ax2.plot(ts.t,ecc,color='orange')
ax3.plot(ts.t,hh,color='orange')
ax4.plot(ts.t,sma,color='orange')
ax5.plot(ts.t,rad,color='orange')
ax6.plot(ts.t,np.sqrt(v2),color='orange')

print "reading ecc 00"
ts=pc.read_ts(datadir="./Incl00Ecc00/data")
rad,hz,hh,sma,ecc,v2,incl = quantities_individual(ts.xq2,ts.yq2,ts.zq2,ts.vxq2,ts.vyq2,ts.vzq2,OO)
ax1.plot(ts.t,incl,color='black',label="e=0.0")
ax2.plot(ts.t,ecc,color='black')
ax3.plot(ts.t,hh,color='black')
ax4.plot(ts.t,sma,color='black')
ax5.plot(ts.t,rad,color='black')
ax6.plot(ts.t,np.sqrt(v2),color='black')

ax5.set_yscale("log")
ax4.set_yscale("log")
ax1.legend()

plt.tight_layout()
#plt.savefig("incl150_suite.png")
plt.show()
