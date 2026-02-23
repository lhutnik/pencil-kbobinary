
import numpy as np
import pylab as plt

figure, ((ax1,ax2,ax3)) = plt.subplots(1, 3,figsize=(12,3.5),sharex=True)

f0           = np.loadtxt("./InclinationSuiteEcc0.4/frachill01/KozaiOnly/Incl99/timeseries.dat",skiprows=1)
f1           = np.loadtxt("./InclinationSuiteEcc0.4/frachill01/KozaiTidal/Incl99/timeseries.dat",skiprows=1)
f2           = np.loadtxt("./InclinationSuiteEcc0.4/frachill01/KozaiTidalJ2/Incl99/timeseries.dat",skiprows=1)
f3           = np.loadtxt("./InclinationSuiteEcc0.4/frachill01/AllIn/Incl99/timeseries.dat",skiprows=1)

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

ax1.set_title("Semimajor Axis")                             
ax1.set_ylabel(r'$a$(km)')
ax1.plot(t,a,color='green',label='No dissipation')
ax1.plot(f1[:,1],f1[:,3],color='blue',label=r'$Q$ tidal')
ax1.plot(f2[:,1],f2[:,3],color='orange',label=r'$Q+J_2$')
ax1.plot(f3[:,1],f3[:,3],color='red',label=r'$Q+J_2+$ drag')
ax1.legend(loc='best',fancybox=True,shadow=True) 

ax2.set_title("Inclination")
ax2.set_ylabel(r'$I$')
ax2.plot(t,np.degrees(np.arccos(cosI)),color='green')
ax2.plot(f1[:,1],np.degrees(np.arccos(f1[:,6])),color='blue')
ax2.plot(f2[:,1],np.degrees(np.arccos(f2[:,6])),color='orange')
ax2.plot(f3[:,1],np.degrees(np.arccos(f3[:,6])),color='red')

ax3.set_title("Pericenter")
ax3.set_ylabel(r'$a(1-e)$(km)')
ax3.set_yscale("log")
ax3.plot(t,a*(1-e),color='green')
ax3.plot(f1[:,1],f1[:,3]*(1-f1[:,4]),color='blue')
ax3.plot(f2[:,1],f2[:,3]*(1-f2[:,4]),color='orange')
ax3.plot(f3[:,1],f3[:,3]*(1-f3[:,4]),color='red')
ax3.axhline(30,linestyle=':',color='black')

ax1.set_xlabel("t (Myr)")
ax2.set_xlabel("t (Myr)")
ax3.set_xlabel("t (Myr)")

#ax1.set_xlim([0,0.05])
#ax2.set_xlim([0,0.05])
#ax3.set_xlim([0,0.05])
#ax4.set_xlim([0,0.05])
#ax5.set_xlim([0,0.05])
#ax6.set_xlim([0,0.05])
#ax7.set_xlim([0,0.05])
#ax8.set_xlim([0,0.05])
#ax9.set_xlim([0,0.05])

plt.tight_layout()
plt.savefig("incl99_qJ2drag.png")
#plt.show()
