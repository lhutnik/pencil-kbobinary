
import numpy as np
import pylab as plt

figure, (ax1,ax2) = plt.subplots(1, 2,figsize=(6,2.5),sharex=True)

f0           = np.loadtxt("./TestRuns/ConservationKozaiOnly_d1e-3/timeseries.dat",skiprows=1)

t           = f0[:,1]
a           = f0[:,3]
hhcte       = f0[:,7]

ax1.set_title("Semimajor Axis")                             
ax1.set_ylabel(r'$1-a/a_0$')
ax1.set_xlim([0,1])
ax1.plot(t,1-a/a[0],color='black')

ax2.set_title("Kozai Constant")
#ax2.set_ylabel(r'$\cos{I} \sqrt{1-e^2}$')
ax2.set_ylabel(r'$1-L_z/L_{z0}$')
ax2.set_xlim([0,1])
ax2.plot(t,1-hhcte/hhcte[0],color='black')

ax1.set_xlabel("t (Myr)")
ax2.set_xlabel("t (Myr)")

plt.tight_layout()
#plt.savefig("kozai_conservation.png")
plt.show()
