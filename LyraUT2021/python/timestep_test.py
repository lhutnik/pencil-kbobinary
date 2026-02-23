
import numpy as np
import pylab as plt

figure, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2, 3,figsize=(12,6),sharex=True)

#f0           = np.loadtxt('workfile_KozaiOnly_I100_E04.dat',skiprows=1)
#f1           = np.loadtxt('workfile_KozaiDiss_I100_E04.dat',skiprows=1)
#f2           = np.loadtxt('workfile_KozaiDissJ2_I100_E04.dat',skiprows=1)
#f3           = np.loadtxt('workfile_allin_I100_E04.dat',skiprows=1)

f0           = np.loadtxt('Kozai_I95_E04_dt1e-2_all/workfile.dat',skiprows=1)
f1           = np.loadtxt('Kozai_I95_E04_dt1e-3_all/workfile.dat',skiprows=1)
f2           = np.loadtxt('Kozai_I95_E04_dt1e-4_all/workfile.dat',skiprows=1)
f3           = np.loadtxt('Kozai_I95_E04_dt1e-5_all/workfile.dat',skiprows=1)
#f3           = np.loadtxt(open('./workfile.dat','rt').readlines()[:-1],skiprows=1)
#f0           = np.loadtxt('workfile_KozaiOnly_I80_E04.dat',skiprows=1)
#f1           = np.loadtxt('workfile_KozaiOnly_I80_E04_dt1e-2.dat',skiprows=1)
#f2           = np.loadtxt('workfile_1e-4.dat',skiprows=1)



it          = f0[:,0]
t           = f0[:,1]
dt          = f0[:,2]
a           = f0[:,3]
e           = f0[:,4]
h           = f0[:,5]
inclination = f0[:,6]
omega       = f0[:,7]
hhcte       = f0[:,8]
ffcte       = f0[:,9]
tf          = f0[:,10]
n           = f0[:,11]
Omega1      = f0[:,12]
Omega2      = f0[:,13]
#hzcte       = f0[:,14]

#it,t/Myr, dt/Myr, a/km, e, h,np.degrees(inclination), np.degrees(omega), hhcte, ffcte, tf/yr, twopi/n/yr, twopi/Omega1/yr, twopi/Omega2/yr)

ax1.set_ylabel(r'$h/h_0$')
#ax1.set_xlabel("t (Myr)")
ax1.set_title("Angular Momentum")
                             
ax2.set_ylabel(r'$a (km)$')
#ax2.set_xlabel("t (Myr)")
ax2.set_title("Semimajor Axis")

ax3.set_ylabel("I")
#ax3.set_xlabel("t (Myr)")
ax3.set_title("Inclination")

ax4.set_ylabel(r'$e$')
ax4.set_xlabel("t (Myr)")
ax4.set_title("Eccentricity")

ax5.set_ylabel(r'$\cos{I} \sqrt{1-e^2}$')
ax5.set_xlabel("t (Myr)")
ax5.set_title("Kozai const H")

##ax6.set_title(r'$-2 - 3e^2 + (3 + 12e^2 - 15e^2\cos^2(\omega))\sin^2I$')
#ax6.set_title(r'$F_{\rm const}$')
##ax6.set_xlabel("t (Myr)")
#ax6.set_ylabel("Kozai const F")

ax6.set_title(r'$a(1-e)$')
ax6.set_xlabel("t (Myr)")
ax6.set_ylabel("Pericenter")

ax1.plot(t,h/h[0],color='black',linestyle=':')
ax1.plot(f1[:,1],f1[:,5]/f1[:,5][0],color='red',linewidth=2.5)
ax1.plot(f2[:,1],f2[:,5]/f2[:,5][0],color='blue')
ax1.plot(f3[:,1],f3[:,5]/f3[:,5][0],color='green',linestyle='--')

#ax2.plot(t,a,color='black',label='No dissipation')
#ax2.plot(f1[:,1],f1[:,3],color='red',label=r'$Q$ tidal')
#ax2.plot(f2[:,1],f2[:,3],color='blue',label=r'$Q+J_2$')
##ax2.plot(f3[:,1],f3[:,3],color='green',label=r'$Q+J_2+$ drag')
#ax2.legend(loc='best',fancybox=True,shadow=True)

ax2.plot(t,a,color='black',label=r'$dt=10^{-2}T_{\rm orb}$',linestyle=':')
ax2.plot(f1[:,1],f1[:,3],color='red',label=r'$dt=10^{-3}T_{\rm orb}$',linewidth=2.5)
ax2.plot(f2[:,1],f2[:,3],color='blue',label=r'$dt=10^{-4}T_{\rm orb}$')
ax2.plot(f3[:,1],f3[:,3],color='green',label=r'$dt=10^{-5}T_{\rm orb}$',linestyle='--')
#ax2.plot(f3[:,1],f3[:,3],color='green',label=r'$Q+J_2+$ drag')
ax2.legend(loc='best',fancybox=True,shadow=True)

ax6.plot(t,a*(1-e),color='black',linestyle=':')
ax6.plot(f1[:,1],f1[:,3]*(1-f1[:,4]),color='red',linewidth=2.5)
ax6.plot(f2[:,1],f2[:,3]*(1-f2[:,4]),color='blue')
ax6.plot(f3[:,1],f3[:,3]*(1-f3[:,4]),color='green',linestyle='--')

#ax2.plot(f2[:,1],f2[:,3],color='red')
#ax2.plot(f2[:,1],f2[:,3]*(1-f2[:,4]))#,color='green')
#ax2.plot(f2[:,1],f2[:,3]*(1+f2[:,4]),color='blue')
#ax2.plot(f3[:,1],f3[:,3],color='red')
#ax2.plot(f3[:,1],f3[:,3]*(1-f3[:,4]))#,color='green')
#ax2.plot(f3[:,1],f3[:,3]*(1+f3[:,4]),color='blue')
#ax2.set_yscale("log")
ax6.set_yscale("log")

ax3.plot(t,inclination,color='black',linestyle=':')
ax3.plot(f1[:,1],f1[:,6],color='red',linewidth=2.5)
ax3.plot(f2[:,1],f2[:,6],color='blue')
ax3.plot(f3[:,1],f3[:,6],color='green',linestyle='--')
        
ax4.plot(t,e ,color='black',linestyle=':')
ax4.plot(f1[:,1],f1[:,4],color='red',linewidth=2.5)
ax4.plot(f2[:,1],f2[:,4],color='blue')
ax4.plot(f3[:,1],f3[:,4],color='green',linestyle='--')
    
ax5.plot(t,hhcte,color='black',linestyle=':')
ax5.plot(f1[:,1],f1[:,8],color='red',linewidth=2.5)
ax5.plot(f2[:,1],f2[:,8],color='blue')
ax5.plot(f3[:,1],f3[:,8],color='green',linestyle='--')

#ax6.plot(t,ffcte)
plt.tight_layout()
plt.savefig("kozai_timesteps_RK3_all_in.png")
#plt.show()
