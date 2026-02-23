import math
import pencil as pc 
import pylab as plt 
import numpy as np

plt.rcParams['axes.formatter.useoffset'] = False
#plt.rcParams['agg.path.chunksize'] = 10000000

#fig,((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9),(axA,axB,axC))=plt.subplots(4,3,figsize=(10,10))

#fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,8))
fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,4))

fig,(ax2) = plt.subplots(1,1,figsize=(6,4))

#ts=pc.read_ts(filename='time_series.190105.1111')
#ts=pc.read_ts(filename='time_series_1.dat',datadir='.')

#ts1=pc.read_ts(filename='time_series_dt1e-2_St5e6_1.dat',datadir='.')
#ts2=pc.read_ts(filename='time_series_dt1e-2_St5e6_2.dat',datadir='.')

ts1=pc.read_ts(filename='time_series_dt1e-2_St5e6_quad_1.dat',datadir='.')
ts2=pc.read_ts(filename='time_series_dt1e-2_St5e6_quad_2.dat',datadir='.')
ts3=pc.read_ts(filename='time_series_dt1e-2_St5e6_quad_3.dat',datadir='.')
ts4=pc.read_ts(filename='time_series_dt1e-2_St5e6_quad_4.dat',datadir='.')

#ts=pc.read_ts()
#par=pc.read_param(param2=True)

xstar = np.concatenate((ts1.xq1,ts2.xq1,ts3.xq1,ts4.xq1))
ystar = np.concatenate((ts1.yq1,ts2.yq1,ts3.yq1,ts4.yq1))

xplanet = np.concatenate((ts1.xq2,ts2.xq2,ts3.xq2,ts4.xq2))
yplanet = np.concatenate((ts1.yq2,ts2.yq2,ts3.yq2,ts4.yq2))

time = np.concatenate((ts1.t,ts2.t,ts3.t,ts4.t))

#print 'plotting 1'
#ax1.plot(xstar,ystar,label='primary',color='orange')
#ax1.plot(xplanet,yplanet,label='secondary',color='red')
#ax1.set_title(r'St=5e6, $u_{\rm gas}$=0, $m_1/m_2$=2.5')
#ax1.set_ylabel('Y')
#ax1.set_xlabel('X')
#print 'plotting 1 done'

#phi=np.linspace(0,2*math.pi,100)

rstar   = np.sqrt(xstar**2 + ystar**2)
rplanet = np.sqrt(xplanet**2 + yplanet**2)

print 'plotting 2'
n=1000

ax2.plot(time[::n]/2/math.pi,rplanet[::n],color='red',label='secondary')
ax2.plot(time[::n]/2/math.pi,rstar[::n],color='orange',label='primary')
ax2.plot(time/2/math.pi,rplanet-rstar,color='black',linestyle='dashed',label='separation')
ax2.set_title(r'St=5e6 linear, $u_{\rm gas}$=0, $m_1/m_2$=2.5, $r_0/r_{\rm hill}$=0.1')
ax2.set_ylabel(r'$r$')
ax2.set_yscale('log')
#ax2.axhline(1./300,linestyle='dotted',color='gray',label='25 km')
ax2.set_xlabel(r'$t/T_0$')
print 'plotting 2 done'

#ax3.plot(ts.t/2/math.pi,rplanet-rstar,color='black')
#ax3.set_title('Orbital Decay')
#ax3.set_ylabel(r'$r_2-r_1$')
#ax3.set_xlabel(r'$t/T_0$')

#mplanet = 0.285714285714
#mstar = 1- mplanet

#xcm = (mstar * xstar + mplanet * xplanet)/(mstar + mplanet)
#ycm = (mstar * ystar + mplanet * yplanet)/(mstar + mplanet)

#print 'plotting 3'
#ax3.plot(xcm,ycm)
#ax3.set_title('Center of Mass')
#ax3.set_ylabel('X')
#ax3.set_xlabel('Y')

#plt.axes().set_aspect('equal')
#ax1.set_aspect('equal')
#ax3.set_aspect('equal')
#print 'plotting 3 - done'

#ax1.legend(loc='best',shadow=True, fancybox=True)
print 'prepare legend'
ax2.legend(loc='best',shadow=True, fancybox=True)
print 'legend done'

#print 'prepare layout'
#plt.tight_layout()
#print 'layout done'
print 'prepare show'
plt.grid()
#plt.show()
print 'done showing'
plt.savefig('BinaryDecay_St5e6_ugas0_q2_quadratic.png')
