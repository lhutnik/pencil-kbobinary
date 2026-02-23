
def quantities_individual(x,y,vx,vy):

    rad  =  np.sqrt(x**2 + y**2)
    phi  =  np.arctan2(y,x)
    vrad =  np.cos(phi)*vx + np.sin(phi)*vy
    vphi = -np.sin(phi)*vx + np.cos(phi)*vy
    

    return rad,vrad,vphi

def quantities_pair(x1,y1,vx1,vy1,x2,y2,vx2,vy2):

    x=x1-x2
    y=y1-y2
    vx=vx1-vx2
    vy=vy1-vy2
    
    rad  =  np.sqrt(x**2 + y**2)
    v2   =  vx**2 + vy**2
    sma  =  1./(2./rad - v2)
    phi  =  np.arctan2(y,x)
    vr   =  np.cos(phi)*vx + np.sin(phi)*vy
    vphi = -np.sin(phi)*vx + np.cos(phi)*vy
    Lphi =  vphi * rad
    ecc  =  np.sqrt(1 - Lphi**2/sma)
    
    return rad,vr,vphi,Lphi,sma,ecc

import math
import pencil as pc 
import pylab as plt 
import numpy as np

plt.rcParams['axes.formatter.useoffset'] = False

fig, ((ax1,ax2,ax3,axA),(ax4,ax5,ax6,axB)) = plt.subplots(2,4,figsize=(12,6))

#fig, ((ax1,ax2),(ax3,axA)) = plt.subplots(2,2,figsize=(8,6))

ts1=pc.read_ts()
par=pc.read_param(param2=True)
time = ts1.t/2/math.pi

r1,vrad1,vphi1   = quantities_individual(ts1.xq1,ts1.yq1,ts1.vxq1,ts1.vyq1)
r2,vrad2,vphi2   = quantities_individual(ts1.xq2,ts1.yq2,ts1.vxq2,ts1.vyq2)
rad,vr,vphi,lphi,sma,ecc = quantities_pair      (ts1.xq1,ts1.yq1,ts1.vxq1,ts1.vyq1,
                                         ts1.xq2,ts1.yq2,ts1.vxq2,ts1.vyq2)

# From Dermott & Murray
m2 = 0.5
m1 = 0.5
a1 = m2/(m1+m2) * sma
a2 = m1/(m1+m2) * sma
h1 = lphi*(m2/(m1+m2))**2
h2 = lphi*(m1/(m1+m2))**2

tau=par.stokesnumber[0]

solution_a = sma[0]*np.exp(-2*ts1.t/tau)
ax1.plot(time,sma,color='red',label='Numerical')
#ax1.plot(time,rad  ,color='red' ,linestyle='-.',label='Separation')
ax1.plot(time,solution_a,color='blue',linestyle='-.',label='Analytical (circ)')
ax1.set_yscale("log")

solution_l = lphi[0]*np.exp(-ts1.t/tau)
ax2.plot(time,lphi,color='red',label='Numerical')
ax2.plot(time,solution_l,color='blue',linestyle='-.',label='Analytical')
ax2.set_yscale("log")

ax3.plot(time,ecc  ,color='red')

solution_v = vphi[0]*np.exp(ts1.t/tau)
axA.plot(time,vphi,color='red',label='Numerical')
axA.plot(time,solution_v,color='blue',linestyle='-.',label='Analytical (circ)')

ax4.plot(time,r1   ,color='red' ,label='sep primary')
#ax4.plot(time,r2   ,color='green',label='sep secondary')
ax4.plot(time,a1   ,color='blue',linestyle='-.',label='sma primary')
#ax4.plot(time,a2   ,color='orange',linestyle='-.',label='sma secondary')
ax4.set_yscale("log")
solution_r = r1[0]*np.exp(-2*ts1.t/tau)
ax4.plot(time,solution_r,color='green',linestyle='--',label='Numerical')

hh1 = r1 * vphi1
hh2 = r2 * vphi2

#ax5.plot(time,vphi1,color='red',label='primary lphi frac bin')
#ax5.plot(time,vphi2,color='green',label='secondary')
u = par.ugas
ax5.plot(time,vphi1,color='red')
solution_vphi = vphi1[0] + u/tau  * ts1.t 
ax5.plot(time,solution_vphi,color='green',linestyle='--')
solution_vphi = vphi1[0] * np.exp(ts1.t/tau) 

ax5.plot(time,solution_vphi,color='green',linestyle='--')
ax5.plot(time,vphi1[0] + u*(1 - np.exp(-ts1.t/tau)),color='cyan',linestyle='--')

#ax5.plot(time,vphi,color='blue',label='primary direct')
#ax5.plot(time,,color='orange',label='secondary')

axB.plot(time,vrad1,color='red',label='primary')
axB.plot(time,vrad2,color='green',label='secondqry')

ax6.plot(time,hh1,color='blue')
ax6.plot(time,h1,color='red')


omega = vphi/rad
T = 2*math.pi/omega

solution1 = hh1[0] * np.exp(-ts1.t/tau) + r1[0]*T[0]*u/(2*math.pi*tau)*ts1.t * np.exp(-2.5*ts1.t/tau)
sol_early = hh1[0] + r1[0]*u/(tau)  * ts1.t * np.exp(-2.0*ts1.t/tau)
#solution2 = hh1[0] * np.exp(-ts1.t/tau) - 
sol_late  = hh1[0] * np.exp(-ts1.t/tau) 
ax6.plot(time,solution1,linestyle='--',color='magenta')

ax6.plot(time,sol_early,linestyle='--',color='cyan')
#ax6.plot(time,solution2,linestyle='-.',color='green')
ax6.plot(time,sol_late,linestyle='--',color='green')

#ax6.set_xlim([0,300])
#ax6.set_ylim([1e-3,10])
#ax6.set_yscale('log')

#ax6.plot(time,vphi1,color='red',label='primary')
#ax6.plot(time,vphi2,color='green',label='secondqry')

ax1.set_ylabel(r'$r$')
ax4.set_ylabel(r'$r$')
ax3.set_ylabel(r'$e$')
ax1.set_xlabel(r'$t/T_0$')
ax4.set_xlabel(r'$t/T_0$')
ax3.set_xlabel(r'$t/T_0$')
ax2.set_ylabel(r'$a$')
axB.set_ylabel(r'$v_r$')
#ax6.set_ylabel(r'$v_\phi$')
#ax6.set_ylabel(r'$v_\phi$')
ax2.set_xlabel(r'$t/T_0$')
axB.set_xlabel(r'$t/T_0$')
ax6.set_xlabel(r'$t/T_0$')

ax1.set_title("Semimajor axis - (equal mass, St=10$^5$)")
ax2.set_title("Angular momentum")
ax3.set_title("Eccentricity")
ax4.set_title("Semimajor axis (rel to CM)")
ax5.set_title("Angular momentum (rel to CM)")
axB.set_title("Radial Velocity (rel to CM)")
axA.set_title("Azimuthal Velocity")
#ax6.set_title("Azimuthal Velocity (rel to CM)")
ax6.set_title("Angular momentum (rel to CM)")

#solution = r1_rel[0]*np.exp(-ts1.t/(2*1e5))

#ax1.plot(time1,r_rel,color='red')
#ax1.plot(time1,solution,color='black',linestyle='dashed')
ax1.legend(loc='best',shadow=True, fancybox=True)
ax1.legend(loc='best',shadow=True, fancybox=True)
ax2.legend(loc='best',shadow=True, fancybox=True)
ax3.legend(loc='best',shadow=True, fancybox=True)
axA.legend(loc='best',shadow=True, fancybox=True)
ax4.legend(loc='best',shadow=True, fancybox=True)
ax5.legend(loc='best',shadow=True, fancybox=True)
ax6.legend(loc='best',shadow=True, fancybox=True)
axB.legend(loc='best',shadow=True, fancybox=True)

ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()
ax5.grid()
ax6.grid()
axA.grid()

#ax1.set_xlim([0,1000])
#ax2.set_xlim([0,1000])
#ax3.set_xlim([0,1000])

#ax1.set_ylim([0.85,1.01])
#ax2.set_ylim([0.43,0.51])
#ax3.set_ylim([0.43,0.51])

plt.tight_layout()
plt.show()
#plt.savefig('EqualMass_St1e5_nogas_analytical.png')
