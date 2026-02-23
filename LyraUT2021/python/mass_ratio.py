import numpy as np
import scipy as sp
import math
import pylab as plt

m1 = np.linspace(0.5,1.0,10000)

tau1 = m1**(2./3)
tau2 = (1-m1)**(2./3)

tw = tau1 * tau2/(tau1 - tau2)

plt.plot(1/m1-1,tw)
plt.yscale('log')
plt.xlabel(r'$m_2/m_1$')
plt.ylabel(r'$\tau_1\tau_2/(\tau_1+\tau_2)$')
plt.title("Wind timescale vs mass ratio")
plt.ylim([0,5e3])
plt.axvline(0.540477641253,linestyle='dashed',color='black')
plt.annotate(r'$MU69$',[0.55,3e-3],size=14)

plt.savefig("mu69_massratio.png")
#plt.show()
