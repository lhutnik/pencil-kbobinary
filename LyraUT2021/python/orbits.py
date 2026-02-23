import math
import pencil as pc 
import pylab as plt 
import numpy as np

plt.rcParams['axes.formatter.useoffset'] = False
#plt.rcParams['agg.path.chunksize'] = 10000000

#fig,((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9),(axA,axB,axC))=plt.subplots(4,3,figsize=(10,10))

fig,((ax1,ax2,axA),(ax3,ax4,axB)) = plt.subplots(2,3,figsize=(10,6))
#fig,(ax2) = plt.subplots(1,1,figsize=(6,4))

#ts1=pc.read_ts(datadir='./NBody_Wind_EqualMass/data')
#par=pc.read_param(datadir='./NBody_Wind_EqualMass/data',param2=True)

ts1=pc.read_ts()
par=pc.read_param(param2=True)

xprimary       =  ts1.xq1
yprimary       =  ts1.yq1
xsecondary     =  ts1.xq2
ysecondary     =  ts1.yq2

vxprimary       =  ts1.vxq1
vyprimary       =  ts1.vyq1
vxsecondary     =  ts1.vxq2
vysecondary     =  ts1.vyq2

xrel  = xsecondary - xprimary
yrel  = ysecondary - yprimary
vxrel  = vxsecondary - vxprimary
vyrel  = vysecondary - vyprimary

r_rel = np.sqrt(xrel**2 + yrel**2)

v2_rel     =  vxrel**2 + vyrel**2
a_rel       =  1./(2./r_rel - v2_rel)
phi_rel    =  np.arctan2(yrel,xrel)
vphi_rel   = -np.sin(phi_rel)*vxrel + np.cos(phi_rel)*vyrel
angular_momentum = vphi_rel * r_rel
e_rel       =  np.sqrt(1 - angular_momentum**2/a_rel)

where_are_NaNs = np.isnan(e_rel)
e_rel[where_are_NaNs] = 0

time = ts1.t/2/math.pi

n=1

ax1.plot(time[::n],a_rel[::n],color='red')
ax1.set_ylabel("SMA")

solution = a_rel[0]*np.exp(-time/(2*1e4))

ax1.plot(time[::n],a_rel[::n],color='red')
ax1.plot(time,solution,color='black',linestyle='dashed')


ax3.plot(time[::n],e_rel[::n],color='red')
ax3.set_ylabel(r'$e=\sqrt{1-v_\phi^2r^2/a}$')
ax3.set_title("Eccentricity")
ax3.set_xlabel(r'$t/T_0$')
ax3.set_ylim([0,1])

ax2.plot(time[::n],r_rel[::n],color='red',linestyle='solid',label=r'$u_{\rm gas}$=0')

#rbox = np.convolve(r_rel_wind, np.ones((n,))/n, mode='valid')
#tbox = np.convolve(time_wind, np.ones((n,))/n, mode='valid')
#apocenter_wind = a_rel_wind * (1+e_rel_wind)
#pericenter_wind = a_rel_wind * (1-e_rel_wind)
#ax2.plot(tbox,rbox,color='blue',linestyle='solid',label=r'$u_{\rm gas}=$'+str(int(par.ugas))+', separation')
#ax2.plot(time_wind[::n],apocenter_wind[::n],color='black',linestyle='solid',label=r'$u_{\rm gas}=$'+str(int(par.ugas))+', apoastron')
#ax2.plot(time_wind[::n],pericenter_wind[::n],color='grey',linestyle='solid',label=r'$u_{\rm gas}=$'+str(int(par.ugas))+', periastron')

ax2.set_title(r'St=1e5 linear, $m_1/m_2$=2.5, $r_0/r_{\rm hill}$=0.1')
ax2.set_ylabel('Orbital Separation')
ax2.set_yscale('log')
ax2.axhline(1./300,linestyle='dotted',color='gray',label='25 km')
#ax2.axhline(10,linestyle='dotted',color='gray',label=r'$r_{\rm Hill}$')
ax2.set_xlabel(r'$t/T_0$')
ax2.set_ylim([1./400,2.])

#ax1.set_xlim([0,25421])
#ax2.set_xlim([0,25421])
#ax3.set_xlim([0,25421])
#ax4.set_xlim([0,25421])
#ax2.set_ylim([0.39,0.47])

ax4.plot(time[::n],angular_momentum[::n],color='red',linestyle='solid')
ax4.set_title('Angular momentum')
ax4.set_ylabel(r'$v_\phi r$')
ax4.set_xlabel(r'$t/T_0$')
ax4.set_ylim([0,1.1])

ax1.legend(loc='best',shadow=True, fancybox=True)
ax2.legend(loc='best',shadow=True, fancybox=True)
ax3.legend(loc='best',shadow=True, fancybox=True)
ax4.legend(loc='best',shadow=True, fancybox=True)

phi1    =  np.arctan2(yprimary,xprimary)
vphi1   = -np.sin(phi1)        *vxprimary   + np.cos(phi1)        *vyprimary
phi2    =  np.arctan2(ysecondary,xsecondary)
vphi2   = -np.sin(phi2)*vxsecondary + np.cos(phi2)*vysecondary

axA.plot(ts1.t/2/math.pi,vphi1,label='primary')
axA.plot(ts1.t/2/math.pi,vphi2,linestyle='dashed',label='secondary')

plt.grid()
plt.tight_layout()
plt.show()
#plt.savefig('Wind_St1e5_ugas100_EqualMass.png')
