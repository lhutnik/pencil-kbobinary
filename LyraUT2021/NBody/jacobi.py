import math
import pencil as pc 
import pylab as plt 
import numpy as np

plt.rcParams['axes.formatter.useoffset'] = False

#fig,((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9),(axA,axB,axC))=plt.subplots(4,3,figsize=(10,10))
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(20,10))

ts=pc.read_ts()

#ax1.plot(ts.t,ts.xq1)
#ax2.plot(ts.t,ts.xq2)
#ax3.plot(ts.t,ts.xq3)

#ax4.plot(ts.t,ts.yq1)
#ax5.plot(ts.t,ts.yq2)
#ax6.plot(ts.t,ts.yq3)

#ax7.plot(ts.t,ts.vxq1)
#ax8.plot(ts.t,ts.vxq2)
#ax9.plot(ts.t,ts.vxq3)

#axA.plot(ts.t,ts.vyq1)
#axB.plot(ts.t,ts.vyq2)
#axC.plot(ts.t,ts.vyq3)

xstar = ts.xq1*np.cos(ts.yq1)
ystar = ts.xq1*np.sin(ts.yq1)

xplanet = ts.xq2*np.cos(ts.yq2)
yplanet = ts.xq2*np.sin(ts.yq2)

xparticle = ts.xq3*np.cos(ts.yq3)
yparticle = ts.xq3*np.sin(ts.yq3)

ax1.plot(xstar,ystar,label='star')
ax1.plot(xplanet,yplanet,label='planet')
ax1.plot(xparticle,yparticle,label='test particle')
ax1.set_title('Orbits')
ax1.set_ylabel('Y')
ax1.set_xlabel('X')

phi=np.linspace(0,2*math.pi,100)
ax1.plot(ts.xq3[0]*np.cos(phi),ts.xq3[0]*np.sin(phi),linestyle='dotted',color='black')

rpar2 = xparticle**2 + yparticle**2

v2 = ts.vxq3**2 + ts.vyq3**2

mu2 = 1e-2
mu1 = 1-mu2

r1 = np.sqrt((xparticle-xstar)**2 + (yparticle-ystar)**2)
r2 = np.sqrt((xparticle-xplanet)**2 + (yparticle-yplanet)**2)

#jacobi_synodic = rpar2 + 2*(mu1/r1 + mu2/r2) - v2

ksi = xparticle
eta = yparticle

vxparticle =  ts.vxq3 * np.cos(ts.yq3) - ts.vyq3 * np.sin(ts.yq3)
vyparticle =  ts.vxq3 * np.sin(ts.yq3) + ts.vyq3 * np.cos(ts.yq3)

ksidot = vxparticle
etadot = vyparticle

jacobi_inertial = 2*(mu1/r1 + mu2/r2) - v2 + 2*(ksi*etadot - eta*ksidot)

ax2.plot(ts.t/2/math.pi,(jacobi_inertial-jacobi_inertial[0])/jacobi_inertial[0])

ax2.set_title('Jacobi constant')
ax2.set_ylabel(r'$(J-J_0)/J_0$')
ax2.set_xlabel(r'$t/T_0$')

#plt.axes().set_aspect('equal')
ax1.set_aspect('equal')
ax1.legend()

plt.show()
#plt.savefig('jacobi1e-4.png')
