subroutine XYZJ(J2,r,a,n,e,Oe,Oq,Oh,OO,X,Y,Z)
!
! This subroutine calculates the tidal factors associated with 
! the permanent J2 quadrupole moment, as derived by Ragozzine &
! Brown (2009). 
!
  real :: fac
  real, intent(in):: J2,r,a,n,e,Oe,Oq,Oh,OO
  real, intent(out):: X,Y,Z
  
  fac = 1.5 * J2 * (r/a)**2 * n/(1-e**2)**2/OO**2
  X = fac * Oe*Oh
  Y = fac * Oq*Oh
  Z = fac * (2*Oh**2 - Oe**2 - Oq**2)*.5
!
endsubroutine XYZJ
!***********************************************************************  
subroutine XYZVW(Q,klove,m1,m2,G,r,a,mu,n,e,Oe,Oq,Oh,X,Y,Z,V,W)
!
! This subroutine calculates the tidal factors associated with 
! the induced quadrupole moment, as in the model of Eggleton et al. (2001),
! adapted by Fabrycky 2007 for the planetary case of tidal dissipation. 
!
  real :: tf
  real, intent(in):: Q,klove,m1,m2,G,r,a,mu,n,Oe,Oq,Oh,e
  real, intent(out):: X,Y,Z,V,W
!
  tf = 1./6 * Q/klove * m1/m2 * 1./sqrt(G*(m1+m2)) * 1./r**5 * a**(13./2)
  V = 9./tf * ((1+15./4*e**2 + 15./8*e**4 + 5./64*e**6)/(1-e**2)**(13./2) - 11./18*Oh/n*(1+3./2*e**2+1./8*e**4)/(1-e**2)**5)
  W = 1./tf * ((1+15./2*e**2 + 45./8*e**4 + 5./16*e**6)/(1-e**2)**(13./2) - Oh/n*(1+3.*e**2+3./8*e**4)/(1-e**2)**5)
!       
  X = -m2*klove*(r/a)**5/(2*mu*n) * Oh*Oe/(1-e**2)**2 - Oq*(1+9./2*e**2+5./8*e**4)/(2*n*tf*(1-e**2)**5)
  Y = -m2*klove*(r/a)**5/(2*mu*n) * Oh*Oq/(1-e**2)**2 + Oe*(1+3./2*e**2+1./8*e**4)/(2*n*tf*(1-e**2)**5)
  Z =  m2*klove*(r/a)**5/(2*mu*n) * ((2*Oh**2-Oe**2-Oq**2)/(2*(1-e**2)**2) + 15*G*m2/a**3 * (1+3./2*e**2+1./8*e**4)/(1-e**2)**5)
!
endsubroutine XYZVW
!***********************************************************************  
subroutine get_orbit(e,h,G,Mbody,a,n,q)
!
! This subroutine takes eccentricity and angular momentum, outputting 
! semimajor axis, period, and pericenter distance. 
!
  real, intent(in):: e,h,G,Mbody
  real, intent(out):: a,n,q
!  
  a = h**2/(G*Mbody*(1-e**2))
  n = sqrt(G*Mbody)/a**1.5 
  q = a*(1-e)
!
endsubroutine get_orbit
!***********************************************************************    
subroutine cross(a,b,c)
!  
! Cross product.
!
  real, dimension(3), intent(in):: a,b
  real, dimension(3), intent(out):: c
!  
  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)
!  
endsubroutine cross
!***********************************************************************   
subroutine dot(a,b,c)
!  
! Dot product.
!
  real, dimension(3), intent(in):: a,b
  real, intent(out) :: c
!  
  c = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
!  
endsubroutine dot
!***********************************************************************      
subroutine get_eqh_components(A,x,y,z,Ax,Ay,Az)
!
! This subroutine calculates the components of a given vector A in the 
! coordinate system defined by orthogonal vectors x, y, and z. 
!
  real :: Asize, A_dot_x, A_dot_y, A_dot_z
  real :: alpha, sina, cosa
  real :: beta,  sinb, cosb
!  
  real, dimension(3), intent(in) :: A
  real, dimension(3), intent(in) :: x,y,z
  real, intent(out):: Ax,Ay,Az
! 
  Asize = sqrt(sum(A**2))
!        
  call dot(A/Asize,x,A_dot_x)
  call dot(A/Asize,y,A_dot_y)
  call dot(A/Asize,z,A_dot_z)
!        
  alpha = acos(A_dot_z)
  cosa  = cos(alpha)
  sina  = sin(alpha)
!        
  beta = atan2(A_dot_y,A_dot_x)
  cosb = cos(beta)
  sinb = sin(beta)
!        
  Ax = Asize*sina*cosb
  Ay = Asize*sina*sinb
  Az = Asize*cosa
!
endsubroutine get_eqh_components
!***********************************************************************
subroutine tidalJ(a,n,e,&
           r1,O1e,O1q,O1h,Omega1,J21,&
           r2,O2e,O2q,O2h,Omega2,J22,&
           XX1,XX2,YY1,YY2,ZZ1,ZZ2)
!
! This subroutine adds the contribution of the permanent quadrupole to 
! the X, Y, and Z tidal parameters. 
!
  real, intent(in) :: a,n,e
  real, intent(in) :: r1,O1e,O1q,O1h,Omega1,J21
  real, intent(in) :: r2,O2e,O2q,O2h,Omega2,J22
  real, intent(inout) :: XX1,XX2,YY1,YY2,ZZ1,ZZ2
  real :: XJ1,YJ1,ZJ1,XJ2,YJ2,ZJ2
!
  call XYZJ(J21,r1,a,n,e,O1e,O1q,O1h,Omega1,XJ1,YJ1,ZJ1)
  call XYZJ(J22,r2,a,n,e,O2e,O2q,O2h,Omega2,XJ2,YJ2,ZJ2)
!  
  XX1 = XX1 + XJ1; XX2 = XX2 + XJ2
  YY1 = YY1 + YJ1; YY2 = YY2 + YJ2
  ZZ1 = ZZ1 + ZJ1; ZZ2 = ZZ2 + ZJ2
!
endsubroutine tidalJ
!***********************************************************************
subroutine dissipation(G,a,mu,n,e,&
           Q1,klove1,m1,r1,O1e,O1q,O1h,&
           Q2,klove2,m2,r2,O2e,O2q,O2h,&
           XX1,XX2,YY1,YY2,ZZ1,ZZ2,VV1,VV2,WW1,WW2)
!
! This subroutine adds the contribution of the induced quadrupole to 
! the V and W dissipation functions, as well as to the 
! X, Y, and Z tidal parameters. 
!
  real, intent(in) :: G,a,mu,n,e
  real, intent(in) :: Q1,klove1,m1,r1,O1e,O1q,O1h
  real, intent(in) :: Q2,klove2,m2,r2,O2e,O2q,O2h
  real, intent(inout) :: XX1,XX2,YY1,YY2,ZZ1,ZZ2,VV1,VV2,WW1,WW2
  real :: X1,Y1,Z1,V1,W1,X2,Y2,Z2,V2,W2
!
  call XYZVW(Q1,klove1,m1,m2,G,r1,a,mu,n,e,O1e,O1q,O1h,X1,Y1,Z1,V1,W1)
  call XYZVW(Q2,klove2,m2,m1,G,r2,a,mu,n,e,O2e,O2q,O2h,X2,Y2,Z2,V2,W2)
!
  XX1 = XX1 + X1; XX2 = XX2 + X2
  YY1 = YY1 + Y1; YY2 = YY2 + Y2
  ZZ1 = ZZ1 + Z1; ZZ2 = ZZ2 + Z2
  VV1 = VV1 + V1; VV2 = VV2 + V2            
  WW1 = WW1 + W1; WW2 = WW2 + W2
!  
endsubroutine dissipation
!***********************************************************************
subroutine thirdbody(Mbody,Msun,Pout,pi,esun,HHH,&
           n,e,ev,qv,hv,See,Sqq,Shh,Sqh,Seh,Seq)
!            
! This subroutine calculates the coupling of the solar orbit to the
! binary orbit, using the model of Eggleton et al. (2001), accurate to
! the quadrupole term. 
!
  real, intent(in) :: Mbody,Msun,Pout,pi,esun,n,e
  real, dimension(3), intent(in) :: HHH,ev,qv,hv
  real, intent(out) :: See,Sqq,Shh,Sqh,Seh,Seq
!  
  real :: Pin,tau,C
  real :: He,Hq,Hh
!
  Pin=2*pi/n
  tau = 2.*Pout**2/(3*pi*Pin) * (Mbody+Msun)/Msun * (1-esun**2)**1.5
  C = 1./(3*tau) * (1-e**2)**(-0.5)
!
  call get_eqh_components(HHH,ev,qv,hv,He,Hq,Hh)
!        
! Calculate the components of the Sij orbit coupling tensor. 
!
  See=C*(1-3*He**2)
  Sqq=C*(1-3*Hq**2)
  Shh=C*(1-3*Hh**2)
  Sqh=-3*C*Hq*Hh
  Seh=-3*C*He*Hh
  Seq=-3*C*He*Hq
!
endsubroutine thirdbody
!***********************************************************************
subroutine evolve_quantities(e,h,ev,hv,O1,O2,t,&
     dedt,dhdt,devdt,dhvdt,dO1dt,dO2dt,ds,dt_,levolve_spin)
!
! This subroutine adds the derivatives to the evolution quantities, 
! and advances the timestep. 
!
  real, intent(in) :: dt_,ds,dedt,dhdt
  real, dimension(3), intent(in) :: devdt,dhvdt,dO1dt,dO2dt
  real, intent(inout) :: e,h,t
  real, dimension(3), intent(inout) :: ev,hv,O1,O2
  logical, intent(in) :: levolve_spin
!
  e = e + dt_*dedt
  h = h + dt_*dhdt
  ev = ev + dt_*devdt
  hv = hv + dt_*dhvdt
!        
  if (levolve_spin) then 
    O1 = O1 + dt_*dO1dt    
    O2 = O2 + dt_*dO2dt
  endif
!                
  t = t + dt_*ds
!
endsubroutine evolve_quantities  
!***********************************************************************
!***********************************************************************          
program kozai_tidal
!
! Main program. This program solves the fourteen equations of the
! Kozai cycles plus tidal friction, with the addition of a planetary 
! permanent J2, and our modest addition of gas drag. 
!  
!
! Constants
!
  real :: pi=3.141592653589793238
  real :: yr=3.15576008d7
  real :: G=6.67430d-8
  real :: Msun=1.98847d33
  real :: AU=1.49597871d13
  real :: km=1d5
  real :: hour=3600.
!
! UT Parameters  
! 
  real :: x1km=20.6
  real :: y1km=19.9
  real :: z1km= 9.4
  real :: x2km=15.4
  real :: y2km=13.8
  real :: z2km= 9.8
  real :: rho_bullet=0.5 ! g/cm3
  real :: mu_body = 4d10 ! rigidity 4e9 N/m2; 4e10 g/cm/s2
  real :: Q1=100.
  real :: Q2=100.
  real :: asun_AU=45.
  real :: frac_hill=0.1
  real :: e=0.4
  real :: t=0.0

  real, dimension(3) :: ev=(/1.,0.,0./)
  real, dimension(3) :: qv=(/0.,1.,0./)  
  real, dimension(3) :: hv=(/0.,0.,1./)
  real, dimension(3) :: O1=(/0.,0.,0./)
  real, dimension(3) :: O2=(/0.,0.,0./)
  real, dimension(3) :: HHH =(/0.,0.,0./)
  real :: O1_hour = 15.
  real :: esun=0.0417249
  real :: I0_degree=95.
  real :: omega=0.0
  real :: tau_drag_Myr=10.
!
  real :: twopi
  real :: Myr,tmaxMyr=10,tmax
  real :: Rbody,Mbody,GM,GM1,sqrtGM
  real :: asun,rhill,a_bin,n_bin
  real :: J21,J22
  real :: a,h,n,mu,nsun,Hsun,Pout
  real :: II,sinI,cosI
  real :: beta,cosb,sinb
  real :: Hcte,h0,tau_drag
  real :: klove1,klove2,Q,klove
  real :: x1_body,y1_body,z1_body,vol1,r1,I1,m1
  real :: x2_body,y2_body,z2_body,vol2,r2,I2,m2
!
  integer :: itmax=100000000
  integer :: itdiagnos=1000
  logical :: lJ2_tidal=.true.
  logical :: ldissipation=.true.
  logical :: lthirdbody=.true.
  logical :: lorbit_drag=.true.
  logical :: levolve_spin
  real    :: dt_period=1d-2
!  
! Runtime parameters
!
  integer :: it,itsub,j
  real :: ds
  real :: period,dt
  real :: dhdt,dedt
  real, dimension(3) :: devdt,dhvdt,dO1dt,dO2dt
  real :: O1e,O1q,O1h,O2e,O2q,O2h
  real :: Omega1,Omega2
  real :: XX1,XX2,YY1,YY2,ZZ1,ZZ2   
  real :: VV1,VV2,WW1,WW2
  real :: See,Sqq,Shh,Sqh,Seh,Seq
  real :: Wd
  real :: muhI1,muhI2
!
  real, dimension(3) :: alpha_ts
  real, dimension(3) :: beta_ts
  real, dimension(3) :: dt_beta_ts
!
  namelist /input/ I0_degree, frac_hill, e, itmax, itdiagnos, lJ2_tidal, ldissipation, lthirdbody,&
       lorbit_drag, dt_period,x1km,y1km,z1km,x2km,y2km,z2km,rho_bullet,mu_body,Q1,Q2,asun_AU,esun,&
       O1_hour,O2_hour,tmaxMyr
!  
  open(2,file='input.in')
  read(2,nml=input)
  close(2)
!
  Myr  = yr*1d6
  tmax = tmaxMyr*Myr
  twopi = 2*pi
  alpha_ts   = (/0.   , -5./9.  ,-153./128./)
  beta_ts    = (/ 1./3., 15./16. ,   8./15. /)
!
! UT 
!
  x1_body=x1km/2 * km
  y1_body=y1km/2 * km
  z1_body=z1km/2 * km
  x2_body=x2km/2 * km
  y2_body=y2km/2 * km
  z2_body=z2km/2 * km
!
  vol1=4./3 *pi*x1_body*y1_body*z1_body
  vol2=4./3 *pi*x2_body*y2_body*z2_body
!
  r1 = (x1_body*y1_body*z1_body)**(1./3)
  r2 = (x2_body*y2_body*z2_body)**(1./3)
!
  print*,'equivalent radius (km)=',r1/km,r2/km
!
  m1 = vol1 * rho_bullet
  m2 = vol2 * rho_bullet
!
  I1=0.2*m1*(x1_body**2+y1_body**2)
  I2=0.2*m2*(x2_body**2+y2_body**2)
!
  print*,'masses=',m1,m2
  print*,'inertia moments=',I1,I2
!
!  One-body equivalent
!
  Rbody = (r1**3+r2**3)**(1./3)
  Mbody = m1+m2
  GM = G*Mbody
  GM1 = 1./GM
  sqrtGM = sqrt(GM)
!
!  Love number
!  
  klove1 = 3./2 / (1. + 19.*mu_body*r1/(2.*G*m1*rho_bullet))
  klove2 = 3./2 / (1. + 19.*mu_body*r2/(2.*G*m2*rho_bullet))
  print*,'klove=',klove1,klove2
  klove = .5*(klove1+klove2)
  Q = .5*(Q1+Q2)
!
  asun=asun_AU*AU
  rhill = asun*((m1+m2)/(3*Msun))**(1./3)
  a_bin = frac_hill*rhill
  n_bin = sqrt(G*(m1+m2))/a_bin**1.5
!
  J21 = 1./10 * (x1_body**2 + y1_body**2 - 2*z1_body**2) / r1**2
  J22 = 1./10 * (x2_body**2 + y2_body**2 - 2*z2_body**2) / r2**2

  print*,"J2=",J21,J22
!
!  Eccentricity
!
  a=a_bin
  h=sqrt(a*(1-e**2)*GM)
  n=sqrtGM/a**1.5 
!
  call cross(hv,ev,qv)
!
! Spin period
!
  O1(3)=2*pi/(O1_hour*hour) 
  O2=O1
  Omega1=sqrt(sum(O1**2))
  Omega2=sqrt(sum(O2**2))
  levolve_spin = lJ2_tidal.or.ldissipation
!
  mu=m1*m2/(m1+m2)
  nsun = sqrt(G*Msun)/asun**1.5
  Hsun = sqrt(G*Msun*asun*(1-esun**2))
  Pout=2*pi/nsun
!
! inclination 
!
  II = I0_degree * pi/180.
  sinI=sin(II)
  cosI=cos(II)
!
  beta = 3*pi/2 - omega  !beta=270-beta
  cosb = cos(beta)
  sinb = sin(beta)
  HHH(1)=sinI*cosb
  HHH(2)=sinI*sinb
  HHH(3)=cosI
!
  Hcte = cosI * sqrt(1-e**2)
!
  h0=h
  tau_drag = tau_drag_Myr * Myr
!
  open(unit=1,file="timeseries.dat",status="replace")
  print*,             'it --- t (Myr) --- dt(Myr) --- a(km) --- e --- h --- i --- hconst --- Pbin (yr) --- Spin1(hr) --- Spin2(hr) '
  write(unit=1,FMT=*) 'it --- t (Myr) --- dt(Myr) --- a(km) --- e --- h --- i --- hconst --- Pbin (yr) --- Spin1(hr) --- Spin2(hr) '
!
  timestepping: do it=1,itmax
!
    n          = sqrtGM/a**1.5 
    period     = 2*pi/n
    ! scale the timestep by the velocity at pericenter
    dt         = dt_period*period*sqrt((1-e)/(1+e)) 
    dt_beta_ts = dt*beta_ts
!
    stageing: do itsub=1,3
      if (itsub == 0) then
        dhdt   = 0.0 ; dedt   = 0.0
        dhvdt  = 0.0 ; devdt  = 0.0
        dO1dt  = 0.0 ; dO2dt  = 0.0
        ds     = 0
      else
        dhdt   = alpha_ts(itsub)*dhdt ;  dedt  = alpha_ts(itsub)*dedt
        dhvdt  = alpha_ts(itsub)*dhvdt; devdt  = alpha_ts(itsub)*devdt
        dO1dt  = alpha_ts(itsub)*dO1dt; dO2dt  = alpha_ts(itsub)*dO2dt
        ds     = alpha_ts(itsub)*ds
      endif
!      
      ds=ds+1
!
      call get_orbit(e,h,G,Mbody,a,n,q)
      call cross(hv,ev,qv)
!
      if (levolve_spin) then 
        call get_eqh_components(O1,ev,qv,hv,O1e,O1q,O1h)
        call get_eqh_components(O2,ev,qv,hv,O2e,O2q,O2h)      
        Omega1=sqrt(sum(O1**2))
        Omega2=sqrt(sum(O2**2))
        muhI1=mu*h/I1
        muhI2=mu*h/I2
      endif
!
      XX1 = 0.; XX2 = 0.; YY1 = 0.; YY2 = 0.; ZZ1 = 0.; ZZ2 = 0.
      VV1 = 0.; VV2 = 0.; WW1 = 0.; WW2 = 0.; Wd  = 0.
      See = 0.; Sqq = 0.; Shh = 0.; Sqh = 0.; Seh = 0.; Seq = 0.
!
      if (lJ2_tidal) call tidalJ(a,n,e,&
           r1,O1e,O1q,O1h,Omega1,J21,&
           r2,O2e,O2q,O2h,Omega2,J22,&
           XX1,XX2,YY1,YY2,ZZ1,ZZ2)
!
      if (ldissipation) call dissipation(G,a,mu,n,e,&
           Q1,klove1,m1,r1,O1e,O1q,O1h,&
           Q2,klove2,m2,r2,O2e,O2q,O2h,&
           XX1,XX2,YY1,YY2,ZZ1,ZZ2,VV1,VV2,WW1,WW2)
!
      if (lthirdbody) call thirdbody(Mbody,Msun,Pout,pi,esun,HHH,&
           n,e,ev,qv,hv,See,Sqq,Shh,Sqh,Seh,Seq)
!
      if (lorbit_drag) Wd=1./tau_drag
!
!  Evolution equations      
!
      dedt = dedt - e*(VV1 + VV2      + 5*(1-e**2)*Seq)
      dhdt = dhdt - h*(WW1 + WW2 + Wd - 5*   e**2 *Seq)
!        
      do j=1,3
        devdt(j) = devdt(j) + (ZZ1 + ZZ2 + (1-e**2)*(4*See-Sqq))*qv(j) - (YY1 + YY2 + (1-e**2)*Sqh)*hv(j)
        dhvdt(j) = dhvdt(j) + (YY1 + YY2 + (1-e**2)*Sqh)*ev(j) - (XX1 + XX2 + (4*e**2+1)*Seh)*qv(j)
!
        if (levolve_spin) then 
          dO1dt(j) = dO1dt(j) + muhI1 * (-YY1*ev(j) + XX1*qv(j) + WW1*hv(j)) 
          dO2dt(j) = dO2dt(j) + muhI2 * (-YY2*ev(j) + XX2*qv(j) + WW2*hv(j))
        endif
!        
      enddo
!
      call evolve_quantities(e,h,ev,hv,O1,O2,t,dedt,dhdt,devdt,dhvdt,dO1dt,dO2dt,ds,dt_beta_ts(itsub),levolve_spin)
!
    enddo stageing
!    
    if (mod(it,itdiagnos) == 0) then 
      call dot(HHH,hv,cosI)
      Hcte = sqrt(1-e**2)*cosI
      print*,             it,t/Myr, dt/Myr, a/km, e, h, cosI, Hcte,twopi/n/yr,twopi/Omega1/hour,twopi/Omega2/hour
      write(unit=1,FMT=*) it,t/Myr, dt/Myr, a/km, e, h, cosI, Hcte,twopi/n/yr,twopi/Omega1/hour,twopi/Omega2/hour
    endif
!
    if ((it == itmax) .or. (t > tmax) .or. (a*(1-e) < 30*km)) then
      print*,'it,itmax=',it,itmax
      print*,'t (Myr),tmax (Myr)=',t/Myr,tmax/Myr
      print*,'q (km) =',a*(1-e)/km
      close(unit=1)
      stop
    endif
!
  enddo timestepping
!   
endprogram kozai_tidal
