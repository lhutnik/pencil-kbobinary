
#import numpy as np
#import pylab as plt

import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.dates as mdates
from matplotlib.transforms import Transform
from matplotlib.ticker import (
    AutoLocator, AutoMinorLocator)

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

#plt.rc('font', size=MEDIUM_SIZE) # controls default text sizes
#plt.rc('axes', titlesize=MEDIUM_SIZE) # fontsize of the axes title
#plt.rc('axes', labelsize=MEDIUM_SIZE) # fontsize of the x and y labels
#plt.rc('xtick', labelsize=SMALL_SIZE) # fontsize of the tick labels
#plt.rc('ytick', labelsize=SMALL_SIZE) # fontsize of the tick labels
#plt.rc('legend', fontsize=MEDIUM_SIZE) # legend fontsize
#plt.rc('figure', titlesize=BIGGER_SIZE) # fontsize of the figure title

x=np.logspace(1,4,100)
y=np.degrees(np.arccos(np.sqrt(6./5 * 1/x)) )

rhill = 4.3e4
rb=30.
x_hill = rhill/rb

fig,ax=plt.subplots(constrained_layout=True)

ax.set_xscale("log")
ax.plot(x,180-y,color='black')
ax.axvline(x_hill,linestyle='--')
ax.annotate(r'$r_{\rm Hill}$ for MU69',xy=(x_hill+200,105))

ax.set_ylabel(r'$I_{\rm crit}$')
ax.set_xlabel(r'$a/R_b$')
ax.set_xlim([1e1,1e4])
ax.set_title("Critical Inclination vs Semimajor Axis")

def xin(x):
    return x*30.

def xout(x):
    return x/30.

secax = ax.secondary_xaxis('top', functions=(xin, xout))
secax.set_xlabel(r'$a_0$ (km)')

plt.show()
#plt.savefig("icrit.png")
