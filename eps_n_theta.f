      subroutine eps_n_theta(npol,Eb,w,q2,tm,um,u_min,thetacm,eps)

c     To calculate model theta_pq in CM and epsilon. This subroutine is largely
c     based on theta_cm.f function, which in turn is based Jochen's script.

      implicit none

      integer npol
      real Eb,w,q2,tm,um,thetacm,eps

      real thetacm1, thetacm2, thetacm3
      real thetacm_t, thetacm_u


      REAL s,omega,q,tmin,umin,u_min
      REAL p1cm,p3cm,e1cm,e3cm,p1lab
     
      REAL p4cm,e4cm

      REAL m2,m3,m4
      REAL m12,m22,m32,m42

      real mp,mp2,mpi,mpi2,mn,mn2

      real mw, mw2, pi

      parameter (mp=.93827231)   !mp
      parameter (mp2=.88035493)  !mp^2
      parameter (mpi=.13956995)  !mpi
      parameter (mpi2=.01947977) !mpi^2
      parameter (mn=.93956563)   !mn
      parameter (mn2=.88278357)  !mn^2
      parameter (pi=3.14159265)  !mn^2

      
ccc/*--------------------------------------------------*/
ccc Modification for omega 
ccc WL
ccc Date: 31/March/2016

      parameter (mw=.78265)      !mw
      parameter (mw2=.61254)     !mw2


ccc
ccc/*--------------------------------------------------*/

      parameter (m3=mw)
      parameter (m32=mw2)

c      if(npol.gt.0) then
c         m2=mp
c         m22=mp2
c         m4=mn
c         m42=mn2
c      else
c         m2=mn
c         m22=mn2
c         m4=mp
c         m42=mp2
c      end if

c       W = 2.48

c      q2 = 1.75

c       W=2.47
c      q2=1.6 
c       q2=2.35 


      m2=mp
      m22=mp2
      m4=mp
      m42=mp2
       


      s=w*w
      omega=(s+q2-m22)/(2*m2)
      q=sqrt(q2+omega**2)

      print*, q2

*     m12=q2    !error?
      m12=-q2   !mass squared of virtual photon.

      e1cm=(s+m12-m22)/(2*w)
      e3cm=(s+m32-m42)/(2*w)

      e4cm=(s+m42-m32)/(2*w)

      p1lab=q
      p1cm=p1lab*m2/w

      p3cm=sqrt(e3cm*e3cm-m32)
      p4cm=sqrt(e4cm*e4cm-m42)


      tmin=-((e1cm-e3cm)**2-(p1cm-p3cm)**2)

      umin=-((e1cm-e4cm)**2-(p1cm-p4cm)**2)




c      print*, 't mininum check in eps n theta check (eps_n_theta.f): ', tm, tmin
c
c      if (tm.ge.tmin .AND. (tm-tmin)/(4*p1cm*p3cm) .lt. 1) then
c         thetacm_t=2*asin(sqrt((tm-tmin)/(4*p1cm*p3cm)))
c      else
c         thetacm_t= 3.1415926
c         print*, 'eps_n_theta: *** tm=',tm,' <  tmin=',tmin,' ! ***'
c      endif



c      thetacm2 = acos((m32 - m12 - 2*m22  + 2*(e1cm + m2) * e4cm  -  2*e1cm * m2)/(2 * p1cm * p4cm)) 
c      thetacm2 = acos((m32 - m12 - 2*m22  + 2*(e1cm + m2) * e4cm  -  2*e1cm * m2)/(2 * p1cm * p4cm)) 

      thetacm2 = acos((2*e1cm*e4cm - m12 - m42 - um)/(2 * p1cm * p4cm)) 
      thetacm3 = acos((2*e1cm*e3cm - m12 - m32 - tm)/(2 * p1cm * p3cm)) 

      thetacm1 = 2*asin(sqrt((um-umin)/(4* p1cm * p4cm)))

      print*, "Momentum Check: ", p1cm, p3cm, p4cm
      print*, "t and u Check ", tm, um
      print*, "Min Check and x: ", tmin, umin, q2/(w**2 + q2 -mp2)

      print*, "Angle Check: ", thetacm, thetacm1, thetacm + thetacm1
      print*, "Angle Check other: ", thetacm3, thetacm2, thetacm2+thetacm3

c      stop





      if (um.ge.umin .AND. (um-umin)/(4*p1cm*p4cm) .lt. 1) then
         thetacm_u=2*asin(sqrt((um-umin)/(4*p1cm*p4cm)))
      else
         thetacm_u= 0.0
         print*, 'eps_n_theta: *** tm=',tm,' <  tmin=',tmin,' ! ***'
         print*, 'eps_n_theta: *** um=',um,' <  umin=',umin,' ! ***'
      endif


      u_min = umin

      thetacm = pi - thetacm_u

     
      
      print*, Eb, omega, q2, 2*(q2+omega**2), 4.*Eb*(Eb-omega)-q2
c      print*,1.+2.*(q2+omega**2)/(4.*Eb*(Eb-omega)-q2)

      eps=1.+2.*(q2+omega**2)/(4.*Eb*(Eb-omega)-q2)
      eps=1./eps





c       print*,"asdasdasda ",  q2, w, Eb, eps, thetacm, (tm-tmin)/(4*p1cm*p3cm), 
c      *   um/(4*p1cm*p3cm)
c    
c 

c      print*, w, q2, umin
c
c      stop

c      print*, thetacm, eps




c      write(*,'(a13,7(F8.5,1x))')
c     *     'eps_n_theta: ',w,q2,t,tmin,thetacm,eps,omega

      end
