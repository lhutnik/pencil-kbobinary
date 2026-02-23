
#inclination 95: collapse; ecc=0.4; 0.1 or 1 period

# dt (orb); tmax=0.01
#       10;  af=1803.1739723597138
#        1;  af=1803.454989860288
#      0.1; af=1803.1403182272027
#     0.01; af=1803.091746952749

def XYZJ(J2,r,a,n,e,Oe,Oq,Oh,OO):
    
    fac = 1.5 * J2 * (r/a)**2 * n/(1-e**2)**2/OO**2
    X = fac * Oe*Oh
    Y = fac * Oq*Oh
    Z = fac * (2*Oh**2 - Oe**2 - Oq**2)*.5

    return X,Y,Z

def XYZVW(Q,klove,m1,m2,G,r,a,mu,n,e,Oe,Oq,Oh):
    
    #tv = 3./2. * (1+klove)**2/klove * Q * n * r**3/(G*m1)
    #tf = tv/9. * (a/r)**8  * m1**2/((m1+m2)*m2) * (1+klove)**(-2.)
    tf = 1./6 * Q/klove * m1/m2 * 1./np.sqrt(G*(m1+m2)) * 1./r**5 * a**(13./2)
    V = 9./tf * ((1+15./4*e**2 + 15./8*e**4 + 5./64*e**6)/(1-e**2)**(13./2) - 11./18*Oh/n*(1+3./2*e**2+1./8*e**4)/(1-e**2)**5)
    W = 1./tf * ((1+15./2*e**2 + 45./8*e**4 + 5./16*e**6)/(1-e**2)**(13./2) - Oh/n*(1+3.*e**2+3./8*e**4)/(1-e**2)**5)
        
    X = -m2*klove*(r/a)**5/(2*mu*n) * Oh*Oe/(1-e**2)**2 - Oq*(1+9./2*e**2+5./8*e**4)/(2*n*tf*(1-e**2)**5)
    Y = -m2*klove*(r/a)**5/(2*mu*n) * Oh*Oq/(1-e**2)**2 + Oe*(1+3./2*e**2+1./8*e**4)/(2*n*tf*(1-e**2)**5)
    Z =  m2*klove*(r/a)**5/(2*mu*n) * ((2*Oh**2-Oe**2-Oq**2)/(2*(1-e**2)**2) + 15*G*m2/a**3 * (1+3./2*e**2+1./8*e**4)/(1-e**2)**5)

    return X,Y,Z,V,W

def get_orbit(e,h,G,Mbody):
    #e = np.sqrt(ev[0]**2+ev[1]**2+ev[2]**2)
    #h = np.sqrt(hv[0]**2+hv[1]**2+hv[2]**2)
    a = h**2/(G*Mbody*(1-e**2))
    n = np.sqrt(G*Mbody)/a**1.5 
    q = a*(1-e)

    return a,n,q

#def get_unitvector(av,a):
#    ahat = av/a
#    return ahat

def cross(a,b):
    c=np.zeros(3) # h cross e 
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = a[2]*b[0] - a[0]*b[2]
    c[2] = a[0]*b[1] - a[1]*b[0]
    return c

def dot(a,b):
    c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    return c

def get_eqh_components(A,x,y,z):

        Asize = np.sqrt(np.sum(A**2))
        
        A_dot_x = dot(A/Asize,x)
        A_dot_y = dot(A/Asize,y)
        A_dot_z = dot(A/Asize,z)
#        
        alpha = np.arccos(A_dot_z)
        cosa  = np.cos(alpha)
        sina  = np.sin(alpha)
#        
        beta = math.atan2(A_dot_y,A_dot_x)
        cosb = np.cos(beta)
        sinb = np.sin(beta)
#        
        Ax = Asize*sina*cosb
        Ay = Asize*sina*sinb
        Az = Asize*cosa

        return Ax,Ay,Az

import numpy as np
import pylab as plt
import math
import sys

one  = np.double(1.)
zero = np.double(0.)

twopi      = 2*math.pi
alpha_ts   = np.double([0.   , -5./9.  ,-153./128.])
beta_ts    = np.double([1./3., 15./16. ,   8./15. ])
#

figure, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3, 3,figsize=(12,6))

#figure, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2, 3,figsize=(12,6))

yr = 3.15576008e7*one
Myr = yr*1e6
hour = 3600.
            
itmax     = 100000000
#itdiagnos = 100
itdiagnos = 1000
tmax      = Myr
isave     = 10000

frac_hill=0.1*one
G=6.67384e-8*one
Msun=2e33*one 
AU=1.49597871e13*one

km = 1e5*one

#
# UT 
#

x1=20.6/2 * 1e5
y1=19.9/2 * 1e5
z1=9.4/2  * 1e5

x2=15.4/2 * 1e5
y2=13.8/2 * 1e5
z2=9.8/2  * 1e5

#x1=100. * km
#y1=100. * km
#z1=100. * km

#x2=100. * km
#y2=100. * km
#z2=100. * km

vol1=4./3*math.pi*x1*y1*z1
vol2=4./3*math.pi*x2*y2*z2

r1 = np.cbrt(x1*y1*z1)
r2 = np.cbrt(x2*y2*z2)

print 'equivalent radius (km)=',r1/km,r2/km

rho_bullet = 0.5  #cgs
#rho_bullet = 1.0*one

m1 = vol1 * rho_bullet
m2 = vol2 * rho_bullet

I1=0.2*m1*(x1**2+y1**2)
I2=0.2*m2*(x2**2+y2**2)

print 'masses=',m1,m2
print 'inertia moments=',I1,I2

mu_body = 4e10*one # rigidity 4e9 N/m2; 4e10 g/cm/s2 
Rbody = np.cbrt(r1**3+r2**3)
Mbody = m1+m2
GM = G*Mbody
GM1 = 1./GM
sqrtGM = np.sqrt(GM)

klove1 =	3./2 / (1. + 19.*mu_body*r1/(2.*G*m1*rho_bullet))
klove2 =	3./2 / (1. + 19.*mu_body*r2/(2.*G*m2*rho_bullet))
#klove1 =	r1/1e10
#klove2 =	r1/1e10

print 'klove=',klove1,klove2
Q1 = 100*one
Q2 = 100*one

asun=45*AU
rhill = asun*((m1+m2)/(3*Msun))**(1./3)
a_bin = frac_hill*rhill
#a_bin = 100*r1
#a_bin = 1000*km 
n_bin = np.sqrt(G*(m1+m2))/a_bin**1.5

klove = np.mean([klove1,klove2])
Q= np.mean([Q1,Q2])

J21 = 1./10 * (x1**2 + y1**2 - 2*z1**2) / r1**2
J22 = 1./10 * (x2**2 + y2**2 - 2*z2**2) / r2**2

print "J2=",J21,J22

tt = list()
qq = list()
QQ = list()
aa = list()
ee = list()
hh = list()
ee0 = list()
ee1 = list()
ee2 = list()
hh0 = list()
hh1 = list()
hh2 = list()
qq0 = list()
qq1 = list()
qq2 = list()
OO1e= list()
OO1q = list()
OO1h = list()
OO2e = list()
OO2q = list()
OO2h = list()
ii = list()
hconst = list()
fconst = list()
ttf = list()
OO1 = list()
OO2 = list()
oomm = list()
nn = list()
dtt = list()

#figure, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3, 3,figsize=(12,6))

ax1.plot([],[], 'o')
ax1.set_ylabel("h")
ax1.set_xlabel("t")
ax1.set_title("Angular Momentum")
ax1.set_autoscaley_on(True)
ax1.ticklabel_format(useOffset=False)
                             
ax2.plot([],[], 'o')
ax2.set_ylabel("a")
ax2.set_xlabel("t")
ax2.set_title("Semimajor Axis")
ax2.set_autoscaley_on(True)
ax2.ticklabel_format(useOffset=False)

ax3.plot([],[], 'o')
ax3.set_ylabel("I")
ax3.set_xlabel("t")
ax3.set_title("Inclination")
ax3.set_autoscaley_on(True)
ax3.ticklabel_format(useOffset=False)

ax4.plot([],[], 'o')
ax4.set_ylabel("e")
ax4.set_xlabel("t")
ax4.set_title("Eccentricity")
ax4.set_autoscaley_on(True)
ax4.ticklabel_format(useOffset=False)

ax5.plot([],[], 'o')
ax5.set_ylabel("hconst")
ax5.set_xlabel("t")
ax5.set_title("Hconst")
ax5.set_autoscaley_on(True)
ax5.ticklabel_format(useOffset=False)

ax6.plot([],[], 'o')
ax6.set_ylabel("fconst")
ax6.set_xlabel("t")
ax6.set_title("Fconst")
ax6.set_autoscaley_on(True)
ax6.ticklabel_format(useOffset=False)

ax7.plot([],[], 'o')
ax7.set_ylabel(r'$P(yr)$')
ax7.set_xlabel("t")
ax7.set_title("Spin/Orbit")
ax7.set_autoscaley_on(True)
ax7.ticklabel_format(useOffset=False)
ax7.set_yscale("log")

ax8.plot([],[], 'o')
ax8.set_ylabel(r'$t_f$')
ax8.set_xlabel("t")
ax8.set_title("Tidal Dissipation")
ax8.set_autoscaley_on(True)
ax8.ticklabel_format(useOffset=False)

ax9.plot([],[], 'o')
ax9.set_ylabel(r'$\omega$')
ax9.set_xlabel("t")
ax9.set_title("omega")
ax9.set_autoscaley_on(True)
ax9.ticklabel_format(useOffset=False)

#
#  Eccentricity
#

a=a_bin
e=0.4*one
h=np.sqrt(a*(1-e**2)*GM)
t=zero
n=sqrtGM/a**1.5 

ev = np.zeros(3)
ev[0]=one
ev[1]=zero
ev[2]=zero

hv = np.zeros(3)
hv[0]=zero
hv[1]=zero
hv[2]=one

qv = np.zeros(3)
#qv[0]=0
#qv[1]=1.
#qv[2]=0.
qv=cross(hv,ev)

# 15 hr period
O1 = np.zeros(3)
O1[0]=0.
O1[1]=0.
O1[2]=2*math.pi/(15*hour) 

O2=O1

mu=m1*m2/(m1+m2)
nsun = np.sqrt(G*Msun)/asun**1.5
esun=0.0417249
Hsun = np.sqrt(G*Msun*asun*(1-esun**2))
Pout=2*math.pi/nsun

#
# inclination 
#

I = np.radians(95.*one)
sinI=np.sin(I)
cosI=np.cos(I)

omega = np.radians(zero)
beta = np.radians(270.*one)-omega
cosb = np.cos(beta)
sinb = np.sin(beta)

H=np.zeros(3)
H[0]=sinI*cosb
H[1]=sinI*sinb
H[2]=cosI

Hcte = cosI * np.sqrt(1-e**2)
Fcte = -2 - 3*e**2 + (3 + 12*e**2 - 15*e**2*np.cos(omega)**2)*sinI**2

h0=h
tau_drag = 10*Myr

f = open('workfile.dat','w+')
print 'it --- t (Myr) --- dt(Myr) --- a(km) --- e --- h --- i --- omega --- hconst --- fconst --- tf (yr) --- Pbin (yr) --- Spin1(yr) --- Spin2(yr) --- Lz '
f.write('it --- t (Myr) --- dt(Myr) --- a(km) --- e --- h --- i --- omega --- hconst --- fconst --- tf (yr) --- Pbin (yr) --- Spin1(yr) --- Spin2(yr) --- Lz \n')

lJ2_tidal=True
ldissipation=True
lspinorbit=True
lorbit_drag=True

for it in np.arange(itmax+1):
#
    n          = sqrtGM/a**1.5 
    period     = 2*math.pi/n
#
    dt         = 1e-2*period  * np.sqrt((1-e)/(1+e))
    dt_beta_ts = [i * dt for i in beta_ts]
#
    for itsub in np.arange(0,2):
        if (itsub == 0):
            dhdt   = zero
            dedt   = zero
            dhvdt  = zero*np.zeros(3)
            devdt  = zero*np.zeros(3)
            dO1dt  = zero*np.zeros(3)
            dO2dt  = zero*np.zeros(3)
            ds     = 0.
        else:
            dhdt   = alpha_ts[itsub]*dhdt
            dedt   = alpha_ts[itsub]*dedt
            dhvdt  = alpha_ts[itsub]*dhvdt
            devdt  = alpha_ts[itsub]*devdt
            dO1dt  = alpha_ts[itsub]*dO1dt
            dO2dt  = alpha_ts[itsub]*dO2dt
            ds     = alpha_ts[itsub]*ds
#
        ds=ds+one
#
        a,n,q = get_orbit(e,h,G,Mbody)
        qv = cross(hv,ev)
#        
        O1e,O1q,O1h=get_eqh_components(O1,ev,qv,hv)
        O2e,O2q,O2h=get_eqh_components(O2,ev,qv,hv)
        Omega1=np.sqrt(np.sum(O1**2))
        Omega2=np.sqrt(np.sum(O2**2))
#        
        XX1 = 0.
        XX2 = 0.
        YY1 = 0.
        YY2 = 0.
        ZZ1 = 0.
        ZZ2 = 0.    
        VV1 = 0.
        VV2 = 0.
        WW1 = 0.
        WW2 = 0. 
        if (lJ2_tidal==True): 
            XJ1,YJ1,ZJ1 = XYZJ(J21,r1,a,n,e,O1e,O1q,O1h,Omega1)
            XJ2,YJ2,ZJ2 = XYZJ(J22,r2,a,n,e,O2e,O2q,O2h,Omega2)
            XX1 = XX1 + XJ1
            XX2 = XX2 + XJ2
            YY1 = YY1 + YJ1
            YY2 = YY2 + YJ2
            ZZ1 = ZZ1 + ZJ1
            ZZ2 = ZZ2 + ZJ2    
#
        if (ldissipation==True):
            X1,Y1,Z1,V1,W1 = XYZVW(Q1,klove1,m1,m2,G,r1,a,mu,n,e,O1e,O1q,O1h)
            X2,Y2,Z2,V2,W2 = XYZVW(Q2,klove2,m2,m1,G,r2,a,mu,n,e,O2e,O2q,O2h)

            VV1 = VV1 + V1
            VV2 = VV2 + V2            
            WW1 = WW1 + W1
            WW2 = WW2 + W2
            XX1 = XX1 + X1
            XX2 = XX2 + X2
            YY1 = YY1 + Y1
            YY2 = YY2 + Y2
            ZZ1 = ZZ1 + Z1
            ZZ2 = ZZ2 + Z2

        See=0.
        Sqq=0.
        Shh=0.
        Sqh=0.
        Seh=0.
        Seq=0.
        if (lspinorbit==True):
#            
            #C = Msun*nsun**2/(4*(Mbody+Msun)*n*(1-e**2)**0.5*(1-esun**2)**1.5)
        
            Pin=2*math.pi/n
            tau = 2.*Pout**2/(3*math.pi*Pin) * (Mbody+Msun)/Msun * (1-esun**2)**1.5
            C = 1./(3*tau) * (1-e**2)**(-0.5)
#
            He,Hq,Hh = get_eqh_components(H,ev,qv,hv)
#        
            See=C*(1-3*He**2)
            Sqq=C*(1-3*Hq**2)
            Shh=C*(1-3*Hh**2)
            Sqh=-3*C*Hq*Hh
            Seh=-3*C*He*Hh
            Seq=-3*C*He*Hq
#
        Wd=0.
        if (lorbit_drag==True):
            Wd=1./tau_drag
#
        dedt = dedt - e*(VV1 + VV2      + 5*(1-e**2)*Seq)
        dhdt = dhdt - h*(WW1 + WW2 + Wd - 5*   e**2 *Seq)
#        
        devdt = devdt +                                 (ZZ1 + ZZ2 + (1-e**2)*(4*See-Sqq)   )*qv - (YY1 + YY2 + (1-e**2)*Sqh)*hv
        dhvdt = dhvdt + (YY1 + YY2 + (1-e**2)*Sqh)*ev - (XX1 + XX2 +          (4*e**2+1)*Seh)*qv
#
        muhI1=mu*h/I1
        muhI2=mu*h/I2
#        
        dO1dt = dO1dt + muhI1 * (-YY1*ev + XX1*qv + WW1*hv) 
        dO2dt = dO2dt + muhI2 * (-YY2*ev + XX2*qv + WW2*hv) 
#
###########################################################################################
#
        #tf_tmp = 1./6 * Q/klove * m1/m2 * 1./np.sqrt(G*(m1+m2)) * 1./r1**5 * a_tmp**(13./2)
        #dadt_tmp = dadt_tmp - a_tmp/(tf_tmp*(1-e_tmp)**(15./2))
        #dedt_tmp = dedt_tmp - 1/(tf_tmp*(1-e_tmp)**(13./2))
        #a_tmp = a_tmp + dt_beta_ts[itsub]*dadt_tmp
        #e_tmp = e_tmp + dt_beta_ts[itsub]*dedt_tmp
#
###########################################################################################
#
        e = e + dt_beta_ts[itsub]*dedt
        h = h + dt_beta_ts[itsub]*dhdt
        ev = ev + dt_beta_ts[itsub]*devdt
        hv = hv + dt_beta_ts[itsub]*dhvdt
        
        O1 = O1 + dt_beta_ts[itsub]*dO1dt
        O2 = O2 + dt_beta_ts[itsub]*dO2dt
                
        t = t + dt_beta_ts[itsub]*ds
#
        if (it==isave):
            np.savez("var",t=t,e=e,h=h,
                     ev0=ev[0],ev1=ev[1],ev2=ev[2],
                     hv0=hv[0],hv1=hv[1],hv2=hv[2],
                     O10=O1[0],O11=O1[1],O12=O1[2],
                     O20=O2[0],O21=O2[1],O22=O2[2])
        
    if ((it % itdiagnos) == 0):
        a,n,q = get_orbit(e,h,G,Mbody)
        qv=cross(hv,ev)
        cosI = dot(H,hv)
        inclination = np.arccos(cosI)
        sinI  = np.sin(inclination)
        beta = math.atan2(dot(H,qv),dot(H,ev))
        cosb = np.cos(beta)
        sinb = np.sin(beta)        
        omega = np.radians(270.) - beta
        while omega > twopi:
            omega = omega - twopi
        while omega < 0:
            omega = omega + twopi            
        cosom=np.cos(omega)
        tf = 1./6 * Q1/klove1 * m1/m2 * 1./np.sqrt(G*(m1+m2)) * 1./r1**5 * a**(13./2)
        hhcte = np.sqrt(1-e**2)*cosI
        ffcte = (-2 - 3*e**2 + (3 + 12*e**2 - 15*e**2*cosom**2)*sinI**2 )
        totzangmom=dot(mu*h*hv + I1*O1 + I2*O2,H) 
        
        print it,t/Myr, dt/Myr, a/km, e, h, np.degrees(inclination), np.degrees(omega), hhcte, ffcte, tf/yr, twopi/n/yr, twopi/Omega1/hour, twopi/Omega2/hour, totzangmom
        f.write("%d %E %E %E %E %E %E %E %E %E  %E %E %E %E %E\n"
                %(it,t/Myr, dt/Myr, a/km, e, h,np.degrees(inclination), np.degrees(omega), hhcte, ffcte, tf/yr, twopi/n/yr, twopi/Omega1/hour, twopi/Omega2/hour,totzangmom))

        tt.append(t/Myr)
        qq.append(q/km)
        QQ.append(a*(1+e)/km)
        hh.append(h/h0)
        aa.append(a/km)        
        ee.append(e)
        ii.append(np.degrees(inclination))
        #oomm.append(np.degrees(omega))
        oomm.append(totzangmom)
#        
        ee0.append(ev[0])
        ee1.append(ev[1])
        ee2.append(ev[2])
#        
        hh0.append(hv[0])
        hh1.append(hv[1])
        hh2.append(hv[2])
#
        qq0.append(qv[0])
        qq1.append(qv[1])
        qq2.append(qv[2])
#        
        OO1e.append(O1[0])
        OO1q.append(O1[1])
        OO1h.append(O1[2])
        OO1.append(twopi/np.sqrt(np.sum(O1**2))/yr)
#
        OO2e.append(O2[0])
        OO2q.append(O2[1])
        OO2h.append(O2[2])
        OO2.append(twopi/np.sqrt(np.sum(O2**2))/yr)
#
        hconst.append(hhcte)
        fconst.append(ffcte)
        ttf.append(tf/yr)
        nn.append(twopi/n/yr)
        dtt.append(dt)

    if ((it == itmax) or t > tmax or (a*(1-e) < 30*km)):

        print 'it,itmax=',it,itmax
        print 't,tmax=',t,tmax
        print 'q (km) =',a*(1-e)/km
    
        #ax1.plot(tt,hh0,color='blue' ,label="h0",linestyle=":")    
        #ax1.plot(tt,hh1,color='green',label="h1",linestyle=":")   
        #ax1.plot(tt,hh2,color='red'  ,label="h2",linestyle=":")
        ax1.plot(tt,hh)
               
        ax2.plot(tt,aa,color='red')
        ax2.plot(tt,qq,color='green')
        ax2.plot(tt,QQ,color='blue')
        ax2.set_yscale("log")
        
        ax3.plot(tt,ii,color='red')
        
        #ax4.plot(tt,ee0,color='blue' ,label="e0",linestyle=":")
        #ax4.plot(tt,ee1,color='green',label="e1",linestyle=":")
        #ax4.plot(tt,ee2,color='red'  ,label="e2",linestyle=":")
        ax4.plot(tt,ee ,color='black',label="total")
        #ax4.set_yscale("log")
    
        ax5.plot(tt,hconst)

        ax6.plot(tt,fconst)

        #ax7.plot(tt,OO1e,color='blue' ,linestyle=":")
        #ax7.plot(tt,OO1q,color='green',linestyle=":")
        #ax7.plot(tt,OO1h,color='red'  ,linestyle=":")
        ax7.plot(tt,OO1 ,color='black')           
        
        #ax7.plot(tt,OO2e,color='blue' ,linestyle="-.")
        #ax7.plot(tt,OO2q,color='green',linestyle="-.")
        #ax7.plot(tt,OO2h,color='red'  ,linestyle="-.")
        ax7.plot(tt,OO2 ,color='black',linestyle='--')           
        ax7.plot(tt,nn,color='cyan')

        ax8.plot(tt,ttf)

        ax9.plot(tt,dtt)
                        
        #ax1.legend(shadow=True,loc='best',fancybox=True)
        #ax4.legend(shadow=True,loc='best',fancybox=True)
        #ax8.legend(shadow=True,loc='best',fancybox=True)
        #ax7.legend(shadow=True,loc='best',fancybox=True)
        
        plt.tight_layout()
        plt.show()
        f.close()
        sys.exit()
        
        #plt.draw()
        #plt.pause(0.0001)
        #plt.clf()
