
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

ecc=0.1

#plt.rc('font', size=MEDIUM_SIZE) # controls default text sizes
#plt.rc('axes', titlesize=MEDIUM_SIZE) # fontsize of the axes title
#plt.rc('axes', labelsize=MEDIUM_SIZE) # fontsize of the x and y labels
#plt.rc('xtick', labelsize=SMALL_SIZE) # fontsize of the tick labels
#plt.rc('ytick', labelsize=SMALL_SIZE) # fontsize of the tick labels
#plt.rc('legend', fontsize=MEDIUM_SIZE) # legend fontsize
#plt.rc('figure', titlesize=BIGGER_SIZE) # fontsize of the figure title

x=np.logspace(1,4,100)
rba=1./x
y1=np.degrees(np.arccos(np.sqrt(6./5 * 1/x)) )
y=np.degrees(np.arccos(np.sqrt(3./5 * (1-(1-rba)**2))))


rhill = 42994.72935601677
rb=30.
x_hill = rhill/rb

#fig,ax=plt.subplots(constrained_layout=True)

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,6))

ax1.set_xscale("log")
ax1.plot(x,180-y,color='black')
ax1.plot(x,180-y1,color='black',linestyle=':')
ax1.axvline(x_hill,linestyle='--')
ax1.annotate(r'$r_{\rm Hill}$ for MU69',xy=(x_hill+200,105))

ax1.set_ylabel(r'$I_{\rm crit}$')
ax1.set_xlabel(r'$a/R_b$')
ax1.set_xlim([1e1,1e4])
ax1.set_title("Critical Inclination vs Semimajor Axis - Kozai only")


ax2.set_xscale("log")
ax2.plot(x,180-y,color='black')
ax1.plot(x,180-y1,color='black',linestyle=':')
ax2.axvline(x_hill,linestyle='--')
ax2.annotate(r'$r_{\rm Hill}$ for MU69',xy=(x_hill+200,105))

ax2.set_ylabel(r'$I_{\rm crit}$')
ax2.set_xlabel(r'$a/R_b$')
ax2.set_xlim([1e1,1e4])
ax2.set_title("Critical Inclination vs Semimajor Axis - Full Model")


xfrachill04=.4*rhill/rb
if (ecc==0.4): 
    ax1.plot(np.repeat(xfrachill04,4),list(range(90,94)),'o',color='blue',markersize=4)
    ax1.plot(np.repeat(xfrachill04,1),[94],'o',mfc='none',color='blue',markersize=4)
elif (ecc==0.1):
    ax1.plot(np.repeat(xfrachill04,3),list(range(90,93)),'o',color='green',markersize=4)
    ax1.plot(np.repeat(xfrachill04,1),[93],'o',mfc='none',color='green',markersize=4)
else:
    print("fried chicken")
#
#
#
xfrachill02=.2*rhill/rb
if (ecc==0.4): 
    ax1.plot(np.repeat(xfrachill02,5),list(range(90,95)),'o',color='blue',markersize=4)
    ax1.plot(np.repeat(xfrachill02,1),[95],'o',mfc='none',color='blue',markersize=4)
elif (ecc==0.1):
    ax1.plot(np.repeat(xfrachill02,4),list(range(90,94)),'o',color='green',markersize=4)
    ax1.plot(np.repeat(xfrachill02,1),[94],'o',mfc='none',color='green',markersize=4)
else:
    print("fried chicken")
#
ax2.plot(np.repeat(xfrachill02,9),list(range(90,99)),'o',color='red',markersize=4)
ax2.plot(np.repeat(xfrachill02,1),[99],'o',mfc='none',color='red',markersize=4)
#

xfrachill01=.1*rhill/rb
if (ecc==0.4): 
    ax1.plot(np.repeat(xfrachill01,7),list(range(90,97)),'o',color='blue',markersize=4)
    ax1.plot(np.repeat(xfrachill01,1),[97],'o',mfc='none',color='blue',markersize=4)
elif (ecc==0.1):
    ax1.plot(np.repeat(xfrachill01,6),list(range(90,96)),'o',color='green',markersize=4)
    ax1.plot(np.repeat(xfrachill01,1),[96],'o',mfc='none',color='green',markersize=4)
else:
    print("fried chicken") 
#
ax2.plot(np.repeat(xfrachill01,10),list(range(90,100)),'o',color='red',markersize=4)
ax2.plot(np.repeat(xfrachill01,1),[100],'o',mfc='none',color='red',markersize=4)
#
#
#
xfrachill004=.04*rhill/rb
if (ecc==0.4): 
    ax1.plot(np.repeat(xfrachill004,10),list(range(90,100)),'o',color='blue',markersize=4)
    ax1.plot(np.repeat(xfrachill004,1),[100],'o',mfc='none',color='blue',markersize=4)
elif (ecc==0.1):
    ax1.plot(np.repeat(xfrachill004,9),list(range(90,99)),'o',color='green',markersize=4)
    ax1.plot(np.repeat(xfrachill004,1),[99],'o',mfc='none',color='green',markersize=4)
else:
    print("fried chicken") 
#
#
#
xfrachill002=.02*rhill/rb
ax1.plot(np.repeat(xfrachill002,14),list(range(90,104)),'o',color='blue',markersize=4)
ax1.plot(np.repeat(xfrachill002,1),[104],'o',mfc='none',color='blue',markersize=4)
#
#
#
xfrachill001=.01*rhill/rb
ax1.plot(np.repeat(xfrachill001,20),list(range(90,110)),'o',color='blue',markersize=4)
ax1.plot(np.repeat(xfrachill001,1),[110],'o',mfc='none',color='blue',markersize=4)



def xin(x):
    return x*30.

def xout(x):
    return x/30.

secax1 = ax1.secondary_xaxis('top', functions=(xin, xout))
secax1.set_xlabel(r'$a_0$ (km)')

secax2 = ax2.secondary_xaxis('top', functions=(xin, xout))
secax2.set_xlabel(r'$a_0$ (km)')

plt.show()
#plt.savefig("CriticalInclination_vs_SMA.png")
