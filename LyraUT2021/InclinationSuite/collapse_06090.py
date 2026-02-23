def quantities_individual(x,y,z,vx,vy,vz,OO):

    rad  =  np.sqrt(x**2 + y**2 + z**2)
    v2   =  vx**2 + vy**2 + vz**2
    sma  =  1./(2./rad - v2)
    hx   =  vz*y - vy*z #- OO*x*z
    hy   =  vx*z - vz*x #- OO*y*z
    hz   =  vy*x - vx*y #- OO*(x**2+y**2)
    hh   =  np.sqrt(hx**2+hy**2+hz**2)
    ecc  =  np.sqrt(1 - hh**2/sma)
    incl =  180/math.pi * np.arccos(hz/hh)

    return rad,hz,hh,sma,ecc,v2,incl

import numpy as np
import math
import pylab as plt
import pencil as pc

nelements=5
ninclination=3

fig, axs = plt.subplots(ninclination,nelements,figsize=(12,6),sharex=True)
#plt.annotate("I=0",xy=(0.1,0.1))
#axs=axs.ravel()

#    (ax1,ax2,ax3,ax4),(ax5,ax6,ax7,ax8),(ax9,ax10,ax11,ax12)) = plt.subplots(3,4,figsize=(12,6),sharex=True)

#plt.subplots_adjust(left=0.3, right=0.97, top=0.9, bottom=0.1, hspace=0.25)

#for i in range(3):
#    for j in range(4):
#        axs[i,j].minorticks_on()
#        axs[i,j].tick_params(which='both',direction='in',top=True,right=True, labelsize=16)
#        axs[i,j].tick_params(which='major',length=10)
#        axs[i,j].tick_params(which='minor',length=5)
#plt.tick_params(axis='both', which='major',labelsize=18)

tau=1e5

par2=pc.read_param(param2=True,datadir="Incl00Ecc00/data")
OO=par2.omega_rotating

for i in range(nelements):
    for j in range(ninclination): 
        axs[j,i].set_xlim([0,1])
        axs[j,i].minorticks_on()
        axs[j,i].tick_params(which='both',direction='in',top=True,right=True)
        
for i in range(nelements): 
    axs[2,i].set_xlabel(r'$t/\tau$')
#
axs[0,0].set_ylabel(r'$I \quad (I_0=0^\circ)$')
axs[1,0].set_ylabel(r'$I \quad (I_0=60^\circ)$')
axs[2,0].set_ylabel(r'$I \quad (I_0=90^\circ)$')
for j in range(ninclination): 
    axs[j,1].set_ylabel(r'$a/a_0$')
    axs[j,2].set_ylabel(r'$h/h_0$')
    axs[j,3].set_ylabel(r'$e$')
    axs[j,4].set_ylabel(r'$u/u_0$')
#
axs[0,0].set_title("Inclination")
axs[0,1].set_title("Semimajor axis")
axs[0,2].set_title("Angular Momentum")
axs[0,3].set_title("Eccentricity")
axs[0,4].set_title("Velocity")
#
window = 10000
weights = np.repeat(1.0, window)/window
#
angle_str=['00','60','90']
ecc_array=[0.9,0.7,0.5,0.3,0]
ecc_str=['09','07','05','03','00']

color=['red','blue','green','orange','black']

for iangle in range(ninclination):
    for iecc in range(nelements):
        datadir="./Incl"+angle_str[iangle]+"Ecc"+ecc_str[iecc]+"/data"
        print datadir
        ts=pc.read_ts(datadir=datadir)
        rad,hz,hh,sma,ecc,v2,incl = quantities_individual(ts.xq2,ts.yq2,ts.zq2,ts.vxq2,ts.vyq2,ts.vzq2,OO)
        #print iangle,iecc,color[iecc],str(eccarray[iecc])
        time=np.convolve(ts.t/tau,weights,'valid')
        axs[iangle,0].plot(time,np.convolve(incl        ,weights,'valid'),color=color[iecc])
        axs[iangle,1].plot(time,np.convolve(sma         ,weights,'valid'),color=color[iecc],label="e="+str(ecc_array[iecc]))
        axs[iangle,2].plot(time,np.convolve(hh          ,weights,'valid'),color=color[iecc])
        axs[iangle,3].plot(time,np.convolve(ecc         ,weights,'valid'),color=color[iecc])
        axs[iangle,4].plot(time,np.convolve(np.sqrt(v2) ,weights,'valid'),color=color[iecc])        
        
#axs5.set_yscale("log")
axs[0,1].legend(loc='best',shadow='True',fancybox='True')

plt.tight_layout()
plt.savefig("inclsuite_0_60_90deg.png")
#print 'showing'
#plt.show()

