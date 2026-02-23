def quantities_individual(x,y,z,vx,vy,vz,OO):

    rad  =  np.sqrt(x**2 + y**2 + z**2)
    v2   =  vx**2 + vy**2 + vz**2
    sma  =  1./(2./rad - v2)
    #hh   =  np.sqrt((vz*y - vy*z)**2 + (vx*z - vz*x)**2 + (vy*x-vx*y)**2)
    hx   =  vz*y - vy*z - OO*x*z
    hy   =  vx*z - vz*x - OO*y*z
    hz   =  vy*x - vx*y - OO*(x**2+y**2)
    hh   =  np.sqrt(hx**2+hy**2+hz**2)
    ecc  =  np.sqrt(1 - hh**2/sma)
    incl =  180/math.pi * np.arccos((vy*x-vx*y)/hh)
    
    return rad,hz,hh,sma,ecc,v2,incl

import pylab as plt
import numpy as np
import pencil as pc
import math

fig, ((ax1,ax2,axB,ax3,ax4,ax5),(ax6,ax7,axC,ax8,ax9,axA)) = plt.subplots(2,6,figsize=(13,6))

ts=pc.read_ts()
t = ts.t /2/math.pi * 5.5 #/ 1e6

par2=pc.read_param(param2=True)
if (par2.lcoriolis_force):
    OO=par2.omega_rotating
else:
    OO=0.0
    
rad,hz,hh,sma,ecc,v2,incl = quantities_individual(ts.xq2,#-ts.xq1,
                                               ts.yq2,#-ts.yq1,
                                               ts.zq2,#-ts.zq1,
                                               ts.vxq2,#-ts.vxq1,
                                               ts.vyq2,#-ts.vyq1,
                                               ts.vzq2,OO)#-ts.vzq1)

par=pc.read_param()
if (par.iprimary == 1):
    iprimary = 0
    isecondary = 1
else:
    iprimary = 1
    isecondary = 0

ax1.plot(t,incl,color='red',label='Numerical')
ax1.set_xlabel(r'$t$ (yr)')
ax1.set_ylabel(r'$i$')
ax1.set_title("Inclination")

ax2.plot(t,ecc,color='red',label='Numerical')
ax2.set_xlabel(r'$t$ (yr)')
ax2.set_ylabel(r'$e$')
ax2.set_title("Eccentricity")

axB.plot(t,hz,color='red',label='Numerical')
axB.set_xlabel(r'$t$ (yr)')
axB.set_ylabel(r'$h_z$')
axB.set_title("Vertical Ang Mom")

axC.plot(t,hh,color='red',label='Numerical')
axC.set_xlabel(r'$t$ (yr)')
axC.set_ylabel(r'$h$')
axC.set_title("Total Ang Mom")

#ax3.plot(ts.xq2-ts.xq1,ts.yq2-ts.yq1,color="black")

ax6.plot(t,sma,color='red',label='Numerical')
ax6.set_xlabel(r'$t$ (yr)')
ax6.set_ylabel(r'$a$')
ax6.set_title("Semimajor axis")

x=ts.xq2
y=ts.yq2
z=ts.zq2
vx=ts.vxq2
vy=ts.vyq2
vz=ts.vzq2

ax3.plot(t,x)
ax3.set_xlabel("t")
ax3.set_ylabel("x")
#ax3.set_title("Orbit")

ax4.plot(t,y)
ax4.set_xlabel("t")
ax4.set_ylabel("y")
#ax4.set_title("Orbit")

ax5.plot(t,z)
ax5.set_xlabel("t")
ax5.set_ylabel("z")
#ax5.set_title("Orbit")

ax8.plot(t,vx)
ax8.set_xlabel("t")
ax8.set_ylabel("vx")
#ax8.set_title("Orbit")

ax9.plot(t,vy)
ax9.set_xlabel("t")
ax9.set_ylabel("vy")
#ax9.set_title("Orbit")

axA.plot(t,vz)
axA.set_xlabel("t")
axA.set_ylabel("vz")
#axA.set_title("Orbit")

ax7.plot(t,rad,color='red',label='Numerical')
ax7.set_xlabel(r'$t$ (yr)')
ax7.set_ylabel(r'$r$')
ax7.set_title("Separation")

#ax6.plot(t,np.sqrt(v2),color='red',label='Numerical')
#ax6.set_xlabel(r'$t$ (Myr)')
#ax6.set_ylabel(r'$v$')
#ax6.set_title("Velocity")

#plt.savefig("draworbit.png")
plt.tight_layout()
plt.show()

