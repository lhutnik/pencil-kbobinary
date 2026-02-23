import math
import pencil as pc 
import pylab as plt 
import numpy as np

plt.rcParams['axes.formatter.useoffset'] = False
#plt.rcParams['agg.path.chunksize'] = 10000000

#fig,((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9),(axA,axB,axC))=plt.subplots(4,3,figsize=(10,10))

fig,((ax1,ax2,axA),(ax3,ax4,axB)) = plt.subplots(2,3,figsize=(12,8))
#fig,(ax2) = plt.subplots(1,1,figsize=(6,4))

ts1=pc.read_ts(datadir='./NBody_Wind/data')
ts2=pc.read_ts(datadir='./NBody_NoWind/data')

par=pc.read_param(datadir='./NBody_Wind/data',param2=True)

xprimary_nowind       =  ts2.xq1
yprimary_nowind       =  ts2.yq1
xsecondary_nowind     =  ts2.xq2
ysecondary_nowind     =  ts2.yq2

vxprimary_nowind       =  ts2.vxq1
vyprimary_nowind       =  ts2.vyq1
vxsecondary_nowind     =  ts2.vxq2
vysecondary_nowind     =  ts2.vyq2

xrel_nowind  = xsecondary_nowind - xprimary_nowind
yrel_nowind  = ysecondary_nowind - yprimary_nowind
vxrel_nowind  = vxsecondary_nowind - vxprimary_nowind
vyrel_nowind  = vysecondary_nowind - vyprimary_nowind

r_rel_nowind = np.sqrt(xrel_nowind**2 + yrel_nowind**2)

v2_rel_nowind     =  vxrel_nowind**2 + vyrel_nowind**2
a_rel_nowind       =  1./(2./r_rel_nowind - v2_rel_nowind)
phi_rel_nowind    =  np.arctan2(yrel_nowind,xrel_nowind)
vphi_rel_nowind   = -np.sin(phi_rel_nowind)*vxrel_nowind + np.cos(phi_rel_nowind)*vyrel_nowind
angular_momentum_nowind = vphi_rel_nowind * r_rel_nowind
e_rel_nowind       =  np.sqrt(1 - angular_momentum_nowind**2/a_rel_nowind)

where_are_NaNs = np.isnan(e_rel_nowind)
e_rel_nowind[where_are_NaNs] = 0

#--
xprimary_wind       =  ts1.xq1
yprimary_wind       =  ts1.yq1
xsecondary_wind     =  ts1.xq2
ysecondary_wind     =  ts1.yq2

vxprimary_wind       =  ts1.vxq1
vyprimary_wind       =  ts1.vyq1
vxsecondary_wind     =  ts1.vxq2
vysecondary_wind     =  ts1.vyq2

xrel_wind  = xsecondary_wind - xprimary_wind
yrel_wind  = ysecondary_wind - yprimary_wind
vxrel_wind  = vxsecondary_wind - vxprimary_wind
vyrel_wind  = vysecondary_wind - vyprimary_wind

r_rel_wind = np.sqrt(xrel_wind**2 + yrel_wind**2)

v2_rel_wind     =  vxrel_wind**2 + vyrel_wind**2
a_rel_wind       =  1./(2./r_rel_wind - v2_rel_wind)
phi_rel_wind    =  np.arctan2(yrel_wind,xrel_wind)
vphi_rel_wind   = -np.sin(phi_rel_wind)*vxrel_wind + np.cos(phi_rel_wind)*vyrel_wind
angular_momentum_wind = vphi_rel_wind * r_rel_wind
e_rel_wind       =  np.sqrt(1 - (angular_momentum_wind)**2/a_rel_wind)
#--



time_wind = ts1.t/2/math.pi
time_nowind = ts2.t/2/math.pi

n=1

ax1.plot(time_nowind[::n],a_rel_nowind[::n],color='red',label='no wind')
ax1.plot(time_wind[::n],a_rel_wind[::n],color='blue',label='wind')
ax1.set_ylabel("SMA")

ax3.plot(time_nowind[::n],e_rel_nowind[::n],color='red',label='no wind')
ax3.plot(time_wind[::n],e_rel_wind[::n],color='blue',label='wind')
ax3.set_ylabel(r'$e=\sqrt{1-v_\phi^2r^2/a}$')
ax3.set_title("Eccentricity")
ax3.set_xlabel(r'$t/T_0$')
ax3.set_ylim([0,1])

axA.plot(time_wind,vphi_rel_wind,color='red',label='no wind')
axA.plot(time_nowind,vphi_rel_nowind,color='blue',label='no wind')

#axA.plot(time_wind[::n],e_rel_wind[::n],color='blue',label='wind')
#axA.set_ylabel(r'$e=\sqrt{1-v_\phi^2r^2/a}$')
#axA.set_title("Eccentricity")
#axA.set_xlabel(r'$t/T_0$')
#axA.set_ylim([0,1])


#ax2.plot(time_wind[::n]  ,r_rel_wind[::n],color='black',linestyle='solid',label=r'$u_{\rm gas}$='+str(int(par.ugas)))
#dt = par.dt
#istride=par.it1
#orbit = 2*math.pi
#dt_per_orbit = int(orbit/dt)
#nt=int(ts1.t[len(ts1.t)-1]/2/math.pi)
#tmean=np.zeros(nt)
#sepmean=np.zeros(nt)
#apoastron=np.zeros(nt)
#periastron=np.zeros(nt)
#for i in range(nt):
#    it1 =  i*dt_per_orbit/istride
#    it2 =  (i+1)*dt_per_orbit/istride
#    print i,it1,it2,ts1.t[it1]/2/math.pi,ts1.t[it2]/2/math.pi
#    tmean[i]=np.mean(time_wind[it1:it2])
#    sepmean[i]=np.mean(r_rel_wind[it1:it2])
#    apoastron[i]=np.max(r_rel_wind[it1:it2])
#    periastron[i]=np.min(r_rel_wind[it1:it2])
#ax2.plot(tmean,sepmean,color='blue',linestyle='solid',label=r'$u_{\rm gas}=$'+str(int(par.ugas))+', orb average')
#ax2.plot(tmean,apoastron,color='black',linestyle='solid',label=r'$u_{\rm gas}=$'+str(int(par.ugas))+', apoastron')
#ax2.plot(tmean,periastron,color='grey',linestyle='solid',label=r'$u_{\rm gas}=$'+str(int(par.ugas))+', periastron')

ax2.plot(time_nowind[::n],r_rel_nowind[::n],color='red',linestyle='solid',label=r'$u_{\rm gas}$=0')

rbox = np.convolve(r_rel_wind, np.ones((n,))/n, mode='valid')
tbox = np.convolve(time_wind, np.ones((n,))/n, mode='valid')

apocenter_wind = a_rel_wind * (1+e_rel_wind)
pericenter_wind = a_rel_wind * (1-e_rel_wind)

ax2.plot(tbox,rbox,color='blue',linestyle='solid',label=r'$u_{\rm gas}=$'+str(int(par.ugas))+', separation')
ax2.plot(time_wind[::n],apocenter_wind[::n],color='black',linestyle='solid',label=r'$u_{\rm gas}=$'+str(int(par.ugas))+', apoastron')
ax2.plot(time_wind[::n],pericenter_wind[::n],color='grey',linestyle='solid',label=r'$u_{\rm gas}=$'+str(int(par.ugas))+', periastron')


ax2.set_title(r'St=1e5 linear, $m_1/m_2$=2.5, $r_0/r_{\rm hill}$=0.1')
ax2.set_ylabel('Orbital Separation')
ax2.set_yscale('log')
ax2.axhline(1./300,linestyle='dotted',color='gray',label='25 km')
#ax2.axhline(10,linestyle='dotted',color='gray',label=r'$r_{\rm Hill}$')
ax2.set_xlabel(r'$t/T_0$')
ax2.set_ylim([1./400,2.])

ax1.set_xlim([0,25421])
ax2.set_xlim([0,25421])
ax3.set_xlim([0,25421])
ax4.set_xlim([0,25421])
#ax2.set_ylim([0.39,0.47])

ax4.plot(time_wind[::n]  ,  angular_momentum_wind[::n],color='blue',linestyle='solid',label='wind')
ax4.plot(time_nowind[::n],angular_momentum_nowind[::n],color='red',linestyle='solid',label='no wind')
ax4.set_title('Angular momentum')
ax4.set_ylabel(r'$v_\phi r$')
ax4.set_xlabel(r'$t/T_0$')
ax4.set_ylim([0,1.1])

ax1.legend(loc='best',shadow=True, fancybox=True)
ax2.legend(loc='best',shadow=True, fancybox=True)
ax3.legend(loc='best',shadow=True, fancybox=True)
ax4.legend(loc='best',shadow=True, fancybox=True)

plt.grid()
plt.show()
#plt.savefig('Wind_vs_NoWind_St1e5_ugas100_eccentricity.png')
