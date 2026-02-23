
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

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rcParams['axes.formatter.useoffset'] = False

fig, (ax6,axB,ax5) = plt.subplots(1,3,figsize=(10,4))

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
u = par.ugas

hh1 = r1 * vphi1
#hh2 = r2 * vphi2

ax6.plot(time,hh1,color='red')
#ax6.plot(time,h1,color='red')
omega = vphi/rad
T = 2*math.pi/omega

#solution1 = hh1[0] * np.exp(-ts1.t/tau) + r1[0]*T[0]*u/(2*math.pi*tau)*ts1.t * np.exp(-2.5*ts1.t/tau)
#ax6.plot(time,solution1,linestyle='--',color='magenta')

sol_early_1 = hh1[0] + r1[0]*u/(tau)  * ts1.t * np.exp(-2.0*ts1.t/tau)
sol_early_2 = hh1[0] - r1[0]*u/(tau)  * ts1.t * np.exp(-2.0*ts1.t/tau)
sol_late  = hh1[0] * np.exp(-ts1.t/tau) 
ax6.plot(time,sol_early_1,linestyle='--',color='magenta')
ax6.plot(time,sol_early_2,linestyle='--',color='magenta')
ax6.plot(time,sol_late,linestyle='--',color='green')

axins = inset_axes(ax6,width='50%',height='50%',loc=7)#, width=XX, height=XX, loc=XX)
axins.plot(time[0:2000],hh1[0:2000],color='red')
axins.plot(time[0:2000],sol_early_1[0:2000],linestyle='--',color='magenta')
axins.plot(time[0:2000],sol_early_2[0:2000],linestyle='--',color='magenta')
axins.plot(time[0:2000],sol_late[0:2000],linestyle='--',color='green')
x0, x1, y0, y1 = 0, 10, -4, 4
axins.axis([x0, x1, y0, y1])
mark_inset(ax6, axins, loc1=2, loc2=3, fc="black", ec="0.0")

axB.plot(time,vrad1,color='red')
sol_up=vrad1[0] + u*(1 - np.exp(-ts1.t/tau))
sol_do=vrad1[0] - u*(1 - np.exp(-ts1.t/tau))
axB.plot(time,sol_up,color='cyan',linestyle='--')
axB.plot(time,sol_do,color='cyan',linestyle='--')
axins = inset_axes(axB,width='50%',height='50%',loc=7)#, width=XX, height=XX, loc=XX)
axins.plot(time[0:2000],vrad1[0:2000],color='red')
axins.plot(time[0:2000],sol_up[0:2000],linestyle='--',color='cyan')
axins.plot(time[0:2000],sol_do[0:2000],linestyle='--',color='cyan')
x0, x1, y0, y1 = 0.0, 10, -7, 7
axins.axis([x0, x1, y0, y1])
mark_inset(axB, axins, loc1=2, loc2=3, fc="black", ec="0.0")

ax5.plot(time,vphi1,color='red')
solution_vphi = vphi1[0] * np.exp(ts1.t/tau) 
ax5.plot(time,solution_vphi,color='green',linestyle='--')
sol_up=vphi1[0] + u*(1 - np.exp(-ts1.t/tau))
sol_do=vphi1[0] - u*(1 - np.exp(-ts1.t/tau))
ax5.plot(time,sol_up,color='cyan',linestyle='--')
ax5.plot(time,sol_do,color='cyan',linestyle='--')
axins = inset_axes(ax5,width='50%',height='50%',loc=7)#, width=XX, height=XX, loc=XX)
axins.plot(time[0:2000],vphi1[0:2000],color='red')
axins.plot(time[0:2000],sol_up[0:2000],linestyle='--',color='cyan')
axins.plot(time[0:2000],sol_do[0:2000],linestyle='--',color='cyan')
axins.plot(time[0:2000],solution_vphi[0:2000],color='green',linestyle='--')
x0, x1, y0, y1 = 0.0, 10, -7, 7
axins.axis([x0, x1, y0, y1])
mark_inset(ax5, axins, loc1=2, loc2=3, fc="black", ec="0.0")

axB.set_title("Radial Velocity")
ax6.set_title("Angular momentum")
ax5.set_title("Azimuthal Velocity")

axB.set_ylabel(r'$v_{r 1}$')
ax6.set_ylabel(r'$L_1$')
ax5.set_ylabel(r'$v_{\phi 1}$')

#ax5.legend(loc='best',shadow=True, fancybox=True)
#ax6.legend(loc='best',shadow=True, fancybox=True)
#axB.legend(loc='best',shadow=True, fancybox=True)

ax5.grid()
ax6.grid()
axB.grid()

ax5.set_xlabel(r'$t/T_0$')
ax6.set_xlabel(r'$t/T_0$')
axB.set_xlabel(r'$t/T_0$')



plt.tight_layout()
#plt.show()
plt.savefig('EqualMass_St1e5_wind_inlets.png')
