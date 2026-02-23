
import numpy as np
import pylab as plt

figure, (ax1,ax2) = plt.subplots(1, 2,figsize=(8,4),sharex=True)

f0           = np.loadtxt("./InclinationSuite/AllIn/Incl90/timeseries.dat",skiprows=1)
f1           = np.loadtxt("./InclinationSuite/AllIn/Incl91/timeseries.dat",skiprows=1)
f2           = np.loadtxt("./InclinationSuite/AllIn/Incl92/timeseries.dat",skiprows=1)
f3           = np.loadtxt("./InclinationSuite/AllIn/Incl93/timeseries.dat",skiprows=1)
f4           = np.loadtxt("./InclinationSuite/AllIn/Incl94/timeseries.dat",skiprows=1)
f5           = np.loadtxt("./InclinationSuite/AllIn/Incl95/timeseries.dat",skiprows=1)
f6           = np.loadtxt("./InclinationSuite/AllIn/Incl96/timeseries.dat",skiprows=1)
#f7           = np.loadtxt("./InclinationSuite/AllIn/Incl97/timeseries.dat",skiprows=1)


it          = f0[:,0]
t           = f0[:,1]
dt          = f0[:,2]
a           = f0[:,3]
e           = f0[:,4]
h           = f0[:,5]
cosI        = f0[:,6]
hhcte       = f0[:,7]
period_yr   = f0[:,8]
Omega1_hr   = f0[:,9]
Omega2_hr   = f0[:,10]

yr=3.15576008e7
Myr=1e6*yr
hour=3600.

ax1.set_title("Inclination")
ax1.set_ylabel(r'$I$')
ax1.plot(t,np.degrees(np.arccos(cosI)),color='black')
ax1.plot(f1[:,1],np.degrees(np.arccos(f1[:,6])),color='red')
ax1.plot(f2[:,1],np.degrees(np.arccos(f2[:,6])),color='blue')
ax1.plot(f3[:,1],np.degrees(np.arccos(f3[:,6])),color='orange')
ax1.plot(f4[:,1],np.degrees(np.arccos(f4[:,6])),color='cyan',linestyle='--')
ax1.plot(f5[:,1],np.degrees(np.arccos(f5[:,6])),color='magenta',linestyle='--')
ax1.plot(f6[:,1],np.degrees(np.arccos(f6[:,6])),color='green',linestyle='--')
ax1.axhline(98,linestyle=':',color='black')
#ax3.plot(f7[:,1],np.degrees(np.arccos(f7[:,6])),color='yellow')

ax2.set_title("Pericenter")
ax2.set_ylabel(r'$a(1-e)$')
ax2.set_yscale("log")
ax2.plot(t,a*(1-e),color='black',label=r'$90^\circ$')
ax2.plot(f1[:,1],f1[:,3]*(1-f1[:,4]),color='red',label=r'$91^\circ$')
ax2.plot(f2[:,1],f2[:,3]*(1-f2[:,4]),color='blue',label=r'$92^\circ$')
ax2.plot(f3[:,1],f3[:,3]*(1-f3[:,4]),color='orange',label=r'$93^\circ$')
ax2.plot(f4[:,1],f4[:,3]*(1-f4[:,4]),color='cyan',label=r'$94^\circ$',linestyle='--')
ax2.plot(f5[:,1],f5[:,3]*(1-f5[:,4]),color='magenta',label=r'$95^\circ$',linestyle='--')
ax2.plot(f6[:,1],f6[:,3]*(1-f6[:,4]),color='green',label=r'$96^\circ$',linestyle='--')
ax2.axhline(30,linestyle=':',color='black')
ax2.legend(loc='best',fancybox=True,shadow=True) 
#ax6.plot(f7[:,1],f7[:,3]*(1-f7[:,4]),color='yellow')

ax1.set_xlabel("t (Myr)")
ax2.set_xlabel("t (Myr)")

plt.tight_layout()
plt.savefig("inclination_kozai_MU69.png")
#plt.show()
