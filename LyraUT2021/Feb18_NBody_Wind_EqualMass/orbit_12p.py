
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

plt.rcParams['axes.formatter.useoffset'] = False

fig, ((ax1,ax2,ax3,ax4),(ax5,ax6,ax7,ax8),(ax9,axA,axB,axC)) = plt.subplots(3,4,figsize=(17,8),sharex=True)
#fig, ((ax1,ax2,ax3,ax4),(ax5,ax6,ax7,ax8),(ax9,axA,axB,axC)) = plt.subplots(3,4,figsize=(14,6),sharex=True)

ts1=pc.read_ts()
par=pc.read_param()
par2=pc.read_param(param2=True)
time = ts1.t/2/math.pi

r1,vrad1,vphi1   = quantities_individual(ts1.xq1,ts1.yq1,ts1.vxq1,ts1.vyq1)
r2,vrad2,vphi2   = quantities_individual(ts1.xq2,ts1.yq2,ts1.vxq2,ts1.vyq2)
rad,vr,vphi,lphi,sma,ecc = quantities_pair      (ts1.xq1,ts1.yq1,ts1.vxq1,ts1.vyq1,
                                         ts1.xq2,ts1.yq2,ts1.vxq2,ts1.vyq2)

# From Dermott & Murray
m2 = par.pmass[1]
m1 = 1-m2
a1 = m2/(m1+m2) * sma
a2 = m1/(m1+m2) * sma
h1 = lphi*(m2/(m1+m2))**2
h2 = lphi*(m1/(m1+m2))**2

tau=par2.stokesnumber[0]

solution_a = sma[0]*np.exp(-2*ts1.t/tau)
ax1.plot(time,sma,color='red',label=r'$a$')
ax1.plot(time,rad,color='green',linestyle='-.',label=r'$r$')
ax1.plot(time,solution_a,color='blue',linestyle=':',label=r'$a_s=a_0 e^{-2t/\tau}$')
ax1.set_yscale("log")

solution_l = lphi[0]*np.exp(-ts1.t/tau)
ax2.plot(time,lphi,color='red',label=r'$r v_\phi$')
ax2.plot(time,solution_l,color='blue',linestyle=':',label=r'$L_s=L_0 e^{-t/\tau}$')
ax2.set_yscale("log")

ax3.plot(time,ecc  ,color='red')

solution_v = vphi[0]*np.exp(ts1.t/tau)
ax4.plot(time,vphi,color='red',label=r'$v_\phi$')
ax4.plot(time,solution_v,color='blue',linestyle=':',label=r'$v_{\phi,s}=v_{\phi,0} e^{t/\tau}$')

nn=1

ax5.plot(time,a1   ,color='red',label=r'$a_1$')
ax5.plot(time,r1   ,color='green',linestyle='-.',label=r'$r_1$')
ax5.plot(time,solution_a * m2/(m1+m2),color='blue',linestyle=':',label=r'$a_s (m_2/M)$')
ax5.set_yscale("log")
hh1 = r1 * vphi1
ax6.plot(time[::nn],hh1[::nn],color='red',label=r'$r_1 v_{\phi 1}$')
ax6.plot(time,h1,color='blue',linestyle=':',label=r'$L (m_2/M)^2$')
#ax6.set_yscale("log")
ax7.plot(time[::nn],vrad1[::nn],color='red')
ax8.plot(time[::nn],vphi1[::nn],color='red',label=r'$v_{\phi 1}$')
ax8.plot(time,solution_v*(m2/(m1+m2)),color='blue',linestyle=':',label=r'$v_{\phi,s} (m_2/M)$')

ax9.plot(time,a2   ,color='red',label=r'$a_2$')
ax9.plot(time,r2   ,color='green',linestyle='-.',label=r'$r_2$')
ax9.plot(time,solution_a*(m1/(m1+m2)),color='blue',linestyle=':',label=r'$a_s (m_1/M)$')
ax9.set_yscale("log")
hh2 = r2 * vphi2
axA.plot(time[::nn],hh2[::nn],color='red',label=r'$r_2 v_{\phi 2}$')
axA.plot(time,h2,color='blue',linestyle=':',label=r'$L (m_1/M)^2$')
#axA.set_yscale("log")
axB.plot(time[::nn],vrad2[::nn],color='red')
axC.plot(time[::nn],vphi2[::nn],color='red',label=r'$v_{\phi 2}$')
axC.plot(time,solution_v*(m2/(m1+m2)),color='blue',linestyle=':',label=r'$v_{\phi,s} (m_1/M)$')

ax1.set_ylabel(r'$a$')
ax2.set_ylabel(r'$L$')
ax3.set_ylabel(r'$e$')
ax4.set_ylabel(r'$v_\phi$')

ax5.set_ylabel(r'$a_1$')
ax6.set_ylabel(r'$L_1$')
ax7.set_ylabel(r'$v_{r 1}$')
ax8.set_ylabel(r'$v_{\phi 1}$')

ax9.set_ylabel(r'$a_2$')
axA.set_ylabel(r'$L_2$')
axB.set_ylabel(r'$v_{r 2}$')
axC.set_ylabel(r'$v_{\phi 2}$')

#ax1.set_xlabel(r'$t/T_0$')
#ax2.set_xlabel(r'$t/T_0$')
#ax3.set_xlabel(r'$t/T_0$')
#ax4.set_xlabel(r'$t/T_0$')
#ax5.set_xlabel(r'$t/T_0$')
#ax6.set_xlabel(r'$t/T_0$')
#ax7.set_xlabel(r'$t/T_0$')
#ax8.set_xlabel(r'$t/T_0$')
ax9.set_xlabel(r'$t/T_0$')
axA.set_xlabel(r'$t/T_0$')
axB.set_xlabel(r'$t/T_0$')
axC.set_xlabel(r'$t/T_0$')


ax1.set_title("Semimajor axis - "+r'$(m_1=m_2, \tau\Omega$='+str(int(tau))+")")
ax2.set_title("Angular momentum")
ax3.set_title("Eccentricity")
ax4.set_title("Azimuthal Velocity")

ax5.set_title("Semimajor axis (rel to CM)")
ax6.set_title("Angular momentum (rel to CM)")
ax7.set_title("Radial Velocity (rel to CM)")
ax8.set_title("Azimuthal Velocity (rel to CM)")

#ax9.set_title("Semimajor axis (rel to CM)")
#axA.set_title("Angular momentum (rel to CM)")
#axB.set_title("Radial Velocity (rel to CM)")
#axC.set_title("Azimuthal Velocity (rel to CM)")

#solution = r1_rel[0]*np.exp(-ts1.t/(2*1e5))

a0 = 6000.
sep = 20. 

contact = sep/a0

ax1.axhline(contact,linestyle=':',color='grey',label='contact')
ax5.axhline(contact/2,linestyle=':',color='grey',label='contact')
ax9.axhline(contact/2,linestyle=':',color='grey',label='contact')


ax1.legend(loc='best',shadow=True, fancybox=True)
ax2.legend(loc='best',shadow=True, fancybox=True)
#ax3.legend(loc='best',shadow=True, fancybox=True)
ax4.legend(loc='best',shadow=True, fancybox=True)
ax5.legend(loc='best',shadow=True, fancybox=True)
ax6.legend(loc='best',shadow=True, fancybox=True)
#ax7.legend(loc='best',shadow=True, fancybox=True)
ax8.legend(loc='best',shadow=True, fancybox=True)
ax9.legend(loc='best',shadow=True, fancybox=True)
axA.legend(loc='best',shadow=True, fancybox=True)
#axB.legend(loc='best',shadow=True, fancybox=True)
axC.legend(loc='best',shadow=True, fancybox=True)

ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()
ax5.grid()
ax6.grid()
ax7.grid()
ax8.grid()
ax9.grid()
axA.grid()
axB.grid()
axC.grid()

#ax7.set_ylim([-0.008,0.008])
#axB.set_ylim([-0.008,0.008])

plt.tight_layout()
#plt.show()
plt.savefig('EqualMass_St1e3_ugas100_analytical_12p.png')
