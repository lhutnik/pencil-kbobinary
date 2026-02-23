
import numpy as np
import pylab as plt

figure, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3, 3,figsize=(12,6),sharex=True)

f0           = np.loadtxt("./InclinationSuite/frachill004/AllIn/Incl90/timeseries.dat",skiprows=1)
f1           = np.loadtxt("./InclinationSuite/frachill004/KozaiOnly/Incl90/timeseries.dat",skiprows=1)
#f2           = np.loadtxt("./InclinationSuite/KozaiTidalJ2/Incl99/timeseries.dat",skiprows=1)
#f3           = np.loadtxt("./InclinationSuite/AllIn/Incl99/timeseries.dat",skiprows=1)
#f4           = np.loadtxt("./InclinationSuite/KozaiOnly/Incl94/timeseries.dat",skiprows=1)
#f5           = np.loadtxt("./InclinationSuite/KozaiOnly/Incl95/timeseries.dat",skiprows=1)
#f6           = np.loadtxt("./InclinationSuite/KozaiOnly/Incl96/timeseries.dat",skiprows=1)
#f7           = np.loadtxt("./InclinationSuite/KozaiOnly/Incl97/timeseries.dat",skiprows=1)


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

ax1.set_title("Angular Momentum")
ax1.set_ylabel(r'$h/h_0$')
ax1.plot(t,h/h[0],color='black')
ax1.plot(f1[:,1],f1[:,5]/f1[:,5][0],color='red',linestyle='--')
#ax1.plot(f2[:,1],f2[:,5]/f2[:,5][0],color='blue')
#ax1.plot(f3[:,1],f3[:,5]/f3[:,5][0],color='orange')
#ax1.plot(f4[:,1],f4[:,5]/f4[:,5][0],color='cyan')
#ax1.plot(f5[:,1],f5[:,5]/f5[:,5][0],color='magenta')
#ax1.plot(f6[:,1],f6[:,5]/f6[:,5][0],color='green')
#ax1.plot(f7[:,1],f7[:,5]/f7[:,5][0],color='yellow')


#ax2.set_title("Semimajor Axis")                             
ax2.set_ylabel(r'$a (km)$')

ax2.plot(t,a,color='black',label='All In')
ax2.plot(f1[:,1],f1[:,3],color='red',label='No Dissipation',linestyle='--')

#ax2.plot(t,a,color='black',label='No dissipation')
#ax2.plot(f1[:,1],f1[:,3],color='red',label=r'$Q$ tidal',linestyle='--')
#ax2.plot(f2[:,1],f2[:,3],color='blue',label=r'$Q+J_2$')
#ax2.plot(f3[:,1],f3[:,3],color='orange',label=r'$Q+J_2+$ drag')
#ax2.plot(f4[:,1],f4[:,3],color='cyan')
#ax2.plot(f5[:,1],f5[:,3],color='magenta')
#ax2.plot(f6[:,1],f6[:,3],color='green')
#ax2.plot(f7[:,1],f7[:,3],color='yellow')
ax2.legend(loc='best',fancybox=True,shadow=True) 


ax3.set_title("Inclination")
ax3.set_ylabel("I")
ax3.plot(t,np.degrees(np.arccos(cosI)),color='black')
ax3.plot(f1[:,1],np.degrees(np.arccos(f1[:,6])),color='red',linestyle='--')
#ax3.plot(f2[:,1],np.degrees(np.arccos(f2[:,6])),color='blue')
#ax3.plot(f3[:,1],np.degrees(np.arccos(f3[:,6])),color='orange')
#ax3.plot(f4[:,1],np.degrees(np.arccos(f4[:,6])),color='cyan')
#ax3.plot(f5[:,1],np.degrees(np.arccos(f5[:,6])),color='magenta')
#ax3.plot(f6[:,1],np.degrees(np.arccos(f6[:,6])),color='green')
#ax3.plot(f7[:,1],np.degrees(np.arccos(f7[:,6])),color='yellow')

ax4.set_title("Eccentricity")
ax4.set_ylabel(r'$e$')
ax4.plot(t,e ,color='black')
ax4.plot(f1[:,1],f1[:,4],color='red',linestyle='--')
#ax4.plot(f2[:,1],f2[:,4],color='blue')
#ax4.plot(f3[:,1],f3[:,4],color='orange')
#ax4.plot(f4[:,1],f4[:,4],color='cyan')
#ax4.plot(f5[:,1],f5[:,4],color='magenta')
#ax4.plot(f6[:,1],f6[:,4],color='green')
#ax4.plot(f7[:,1],f7[:,4],color='yellow')

#ax5.set_title("Kozai const H")
ax5.set_ylabel(r'$\cos{I} \sqrt{1-e^2}$')
ax5.plot(t,hhcte,color='black')
ax5.plot(f1[:,1],f1[:,7],color='red',linestyle='--')
#ax5.plot(f2[:,1],f2[:,7],color='blue')
#ax5.plot(f3[:,1],f3[:,7],color='orange')
#ax5.plot(f4[:,1],f4[:,7],color='cyan')
#ax5.plot(f5[:,1],f5[:,7],color='magenta')
#ax5.plot(f6[:,1],f6[:,7],color='green')
#ax5.plot(f7[:,1],f7[:,7],color='yellow')

ax6.set_title(r'$a(1-e)$')
ax6.set_ylabel("Pericenter")
ax6.set_yscale("log")
ax6.plot(t,a*(1-e),color='black')
ax6.plot(f1[:,1],f1[:,3]*(1-f1[:,4]),color='red',linestyle='--')
#ax6.plot(f2[:,1],f2[:,3]*(1-f2[:,4]),color='blue')
#ax6.plot(f3[:,1],f3[:,3]*(1-f3[:,4]),color='orange')
#ax6.plot(f4[:,1],f4[:,3]*(1-f4[:,4]),color='cyan')
#ax6.plot(f5[:,1],f5[:,3]*(1-f5[:,4]),color='magenta')
#ax6.plot(f6[:,1],f6[:,3]*(1-f6[:,4]),color='green')
#ax6.plot(f7[:,1],f7[:,3]*(1-f7[:,4]),color='yellow')

ax7.set_title("Periods")
ax7.set_ylabel(r'$T$ (hr)')
ax7.set_yscale("log")
ax7.legend(loc='best',fancybox=True,shadow=True)
ax7.plot(t,period_yr*yr/hour,color='black',label='orbit')
ax7.plot(t,Omega1_hr,color='black',linestyle='--',label='spin 1')           
ax7.plot(t,Omega2_hr,color='black',linestyle=':',label='spin 2') 

ax8.set_title("Timestep")
ax8.set_ylabel(r'$\Delta t/T_{\rm binary}$')
ax8.plot(t,dt*Myr/(period_yr*yr),color='black')
ax8.plot(f1[:,1],f1[:,2]*Myr/(f1[:,8]*yr),color='red',linestyle='--')
#ax8.plot(f2[:,1],f2[:,2]*Myr/(f2[:,8]*yr),color='blue')
#ax8.plot(f3[:,1],f3[:,2]*Myr/(f3[:,8]*yr),color='orange')
#ax8.plot(f4[:,1],f4[:,2]*Myr/(f4[:,8]*yr),color='cyan')
#ax8.plot(f5[:,1],f5[:,2]*Myr/(f5[:,8]*yr),color='magenta')
#ax8.plot(f6[:,1],f6[:,2]*Myr/(f6[:,8]*yr),color='green')
#ax8.plot(f7[:,1],f7[:,2]*Myr/(f7[:,8]*yr),color='yellow')

ax7.set_xlabel("t (Myr)")
ax8.set_xlabel("t (Myr)")
ax9.set_xlabel("t (Myr)")

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
#plt.savefig("kozai_timesteps_RK3.png")
plt.show()
