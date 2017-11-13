      program calc_xsect


      implicit none

c This script computes the experimental cross section using the ratio
c DATA/MC(Yields * CS(MC)

c      character*2 prv_it
c      common prv_it

c      integer q2_bin
c      integer u_bin, phi_bin
c      common u_bin, phi_bin
      
c     Get number of the previous iteration.
      
c      if(iargc().ne.1) then
c         print*,'*** usage: calc_xsect prv_it ***'
c         stop
c      end if
c      call getarg(1,prv_it)

c     Calculate unseparated cross-sections. Now settings are for the piplus data (+)




c      stop

      call xsect(+1,1.60,0.32)
      call xsect(+1,1.60,0.59)

      call xsect(+1,2.45,0.27)
      call xsect(+1,2.45,0.55)

      stop
      end

*=======================================================================

      subroutine xsect(npol_set,q2_set,eps_set)

      implicit none

      integer npol_set
      real q2_set,eps_set

      integer kset,nbin

      character*80 r_fn, kin_fn, xunsep_fn
      character*2 pol

      integer it,ip
      real Eb,eps

      real one

      integer nt,nphi
c      parameter (nt=3,nphi=8)
c      parameter (nt=2,nphi=10)

      real r,dr,w,dw,q2,dq2,tt,dtt,th_pos,th_cm
      real tm,tmn,tmx
      real um
      real eps_mod,th_mod,x_mod,um_min
      real x_real,dx_real

      integer ipol
      real th_pq

      real phi

      real, Dimension(10) :: u_bin_boundary

      real q2_bin

      integer u_bin, phi_bin
c      common  u_bin, phi_bin

      character*80:: line

      integer j


      include './kin_xyz.inc'








c   /*--------------------------------------------------*/
c   Read the u and phi bins 

      open (unit = 22, file = "../u_bin_interval", action='read')
      read (22,*) q2_bin, u_bin, phi_bin

      nt = u_bin
      nphi = phi_bin 




c      read (22, '(A)') line  
c      read (22, '(A)') line  
c      print*,  line


c      read (line, end=20) 

c      print*,  trim(line)

c      read(line, *) u_bin_boundary(0), u_bin_boundary(1), u_bin_boundary(2)

c      read(line, *) (u_bin_boundary(j), j = 1, 3)

c      print*, u_bin_boundary(1)
c       read (22,*) 
c       read (22,*) q2_bin, u_bin, phi_bin
c       read (22,*) 
c 
c       print*,  u_bin, phi_bin

c      do j = 1, 3
c         print*, u_bin_boundary(j)
c      end do


      if(q2_set.eq.1.6) then

         read (22, '(A)') line  
         read(line, *) (u_bin_boundary(j), j = 1, u_bin+1)

c        u_bin_boundary = (/ 0.0, 0.12,  0.20, 0.40/)
c        u_bin_boundary = (/0.0, 0.10, 0.17, 0.32/)
 
      elseif(q2_set.eq.2.45) then
 
c        u_bin_boundary = (/ 0.0, 0.212, 0.33, 0.60/)
c        u_bin_boundary = (/0.0, 0.19, 0.30, 0.50/)
 
         read (22,*) 
         read (22,*) 
         read (22, '(A)') line  
         read(line, *) (u_bin_boundary(j), j = 1,  u_bin+1)

      endif

      
      


c      stop

      close(22)

      ipol=0
      q2=0.
      eps=0.
      tmn=0.
      tmx=0.
      kset=0
      open(55,file='./list.settings.omega')
      do while(ipol.ne.npol_set.or.q2.ne.q2_set.or.eps.ne.eps_set)
         read(55,*) ipol,q2,eps,th_pq,tmn,tmx,nbin,kset
c         write(6,2)ipol,q2,eps,th_pq,tmn,tmx,nbin,kset
c 2       format(i5,5f10.5,2i5)
      end do
      close(55)
      write(6,3)tmn,tmx,kset
 3    format(' tmn, tmx: ',2f10.5,'  for kset=',i5)
      if(tmn.eq.0..or.tmx.eq.0.) 
     *     stop '*** setting is not found in list.settings.fpi2'

      Eb=Ebeam(kset)

      if(npol_set.lt.0) then
         pol='mn'
      else
         pol='pl'
      end if

      write(6,4)Eb,q2,eps,pol
 4    format(' xsect: Eb=',f8.5,'   at Q2=',f7.4,
     *     '  eps=',f6.4,'  pol=',a2)

c     construct ratio data file name.

      write(r_fn,10) pol,nint(q2*100),nint(eps*100)
 10   format('averages/aver.',a2,'_',i3.3,'_',i2,'.dat')
      print*,'xsect: r_fn=',r_fn

      open(51,file=r_fn)

c     construct kinematics data file name.

      write(kin_fn,20) nint(q2*100)
 20   format('averages/avek.',i3.3,'.dat')
      print*,'xsect: kin_fn=',kin_fn

      open(52,file=kin_fn)

*     construct output file name.
      write(xunsep_fn,30) pol,nint(q2_set*100),nint(eps_set*100)
 30   format('xsects/x_unsep.',a2,'_',i3.3,'_',i2)
      print*,'xsect: xunsep_fn=',xunsep_fn
c      pause

      open(61,file=xunsep_fn,status='replace')

      
      print*, ""
      print*, ""
      print*, ""
      print*, "Check for the reading file "
      print*, ""

      nbin = u_bin


      do it=1,nbin

c         tm=tmn+(it-0.5)*(tmx-tmn)/nbin


         um = (u_bin_boundary(it) + u_bin_boundary(it+1)) / 2

c         print *, "11112222  " , nbin, u_bin_boundary(1), u_bin_boundary(2), u_bin_boundary(3), um 

c         print *, nbin, u_bin_boundary(it), u_bin_boundary(it+1)  

c         stop


         read(52,*) w,dw,q2,dq2,tt,dtt,th_pos
         write(6,32) w,dw,q2,dq2,tt,dtt,th_pos
 32      format('xsect: ',7f10.4)


         th_cm=th_pos


         tm = tt


c         print *,  w,dw,q2,dq2,tt,dtt,th_pos

c        stop

         th_cm=th_cm*3.14159D0/180.D0
         

         print*, " "
         print*, " ", th_cm
         print*, " "


         do ip=1,nphi

c            phi=(ip-0.5)*2.*3.14159/nphi + 10 * 3.14159 / 180
c            phi=(ip-0.5)*2.*3.14159/nphi + 25.7 * 3.14159 / 180
            phi=(ip-0.5)*2.*3.14159/nphi


            if (phi .le. 0.0) then
                
             phi= phi + 2 * 3.14159   
                
            else if (phi .gt. 2 * 3.14159 ) then

             phi= phi - 2 * 3.14159   

c            else if (phi .eq. 360) then
c             phi = 0

            end if

            

  
c            phi=(ip-1)*2.*3.14159/nphi
            read(51,*) r,dr

c            print *, "ratio check: ", r, dr


c            stop

            
c             print *, ip, r, dr, it, nphi, nbin
c             print *, "~~ ", w, dw, q2, dq2, th_pos
c             print *, "== ", tm, phi
c             print *, "-- ", npol_set, Eb, q2_set, w, q2, tm, phi



c            print*, q2_set, um

c            stop


c            call xmodel_ghm(npol_set,Eb,q2_set,w,q2,um,tm,phi,
c     *           eps_mod,th_mod,x_mod)

            call xmodel_ghm_two_model(npol_set,Eb,q2_set,w,q2,um,tm,phi,
     *           eps_mod,th_mod,x_mod,um_min)


c           stop


c            if (q2_set.eq.1.6) then
c                
c               print*, "Q2=1.6 parameterization is used! "
c
c                call xmodel_ghm_160(npol_set,Eb,q2_set,w,q2,um,tm,phi,
c     *           eps_mod,th_mod,x_mod)
c
cc                stop
c
c            else if (q2_set.eq.2.45) then
c
c               print*, "Q2=2.45 parameterization is used! "
c
c                call xmodel_ghm_245(npol_set,Eb,q2_set,w,q2,um,tm,phi,
c     *           eps_mod,th_mod,x_mod)
c
c            else
c                print*, "No parameterization is aviliable for Q2=", q2_set
c                stop
c            endif
            
c            print*,"testtest" , q2_set 
c            stop

c            call xmodel_ghm(npol_set,Eb,q2_set,w,q2,um,tm,phi,
c     *           eps_mod,th_mod,x_mod)


c            call xmodel_ghm(npol_set,Eb,q2_set,w,q2,um,tm,phi,
c     *           eps_mod,th_mod,x_mod)



c            stop


cc /*--------------------------------------------------*/
cc angle check
cc angle check is taken off temporarily 
cc WL 31/March/2016


c           if (abs(th_mod-th_cm).gt.1.e-4) then
c              write(6,*)' Angle error ',th_mod,th_cm
c              stop
c           endif
c
c

c /*--------------------------------------------------*/

c ratio is data/simc - see GH logbook, p.55

               x_real=x_mod*r
               dx_real=x_mod*dr/r




             if (x_real.eq.0.0) then
                dx_real = -1000
             endif
            


c             dx_real=x_mod*dr

c            x_real=x_mod
c            dx_real=0.01




            print*, "JJJJJJJJ  ratio check ", r, x_mod, x_real 



cc ratio is simc/data - see GH logbook, p.55
c            if (r.ne.0) then
c               x_real=x_mod/r
c               dx_real=x_mod*dr/r**2
c            else
c               x_real=0.
c               dx_real=0.
c            endif


c           print*, "ratio check ", r, x_mod, x_real 


c   /*--------------------------------------------------
c   This is where data output to file happens


            
            print*, "kkkkkk  ", x_mod, eps_mod, th_mod 

            write(61,40) x_real,dx_real,x_mod,eps_mod,
     *           th_mod*180./3.14159,phi*180./3.14159,tm,um,um_min,w,q2
 40         format(3G15.5,f8.5,2f7.2,5f8.5)

         end do                 !phi


c         stop 
         
c         print*, x_real,dx_real,x_mod,eps_mod, 80./3.14159,phi*180./3.14159,tm,um,w,q2

c         stop

c        Write out kinematics for Henk.
         if(npol_set.gt.0) write(99,'(5f8.3,2x,2f6.2)')
     *   w,q2,eps_mod,th_mod*180./3.14159,tm,eps_set,q2_set

      end do                    !t

      close(51)
      close(52)
      close(61)
      print*,' '

      end

*=======================================================================

      subroutine xmodel(npol_set,Eb,q2_set,w,q2,tm,phi,
     *     eps_mod,th_mod,x_mod)


      integer npol_set
      real Eb,q2_set,w,q2,tm,phi,eps_mod,th_mod, thetacm, x_mod

      real*8 sig
      real*8 t_gev,tprime_gev,q2_gev,w_gev,eps,tp
      real*8 lambda0_sq,lambdapi_sq,alphapi_t,alphapi_0
      real*8 a,b,c,d,e,f,fpi,fpi235
      real*8 m_pi0sq
      real*8 phicm, pi

      pi = 3.1415926

      call eps_n_theta(npol_set, Eb, w, q2, tm, thetacm, eps_mod)
      
c       print*, "check   ", thetacm
c 
c       stop

c      print*,"*****", npol_set, Eb, w, q2, tm, thetacm, eps_mod


      phicm = phi

      tprime_gev = tm

c      print*, "aaaaaaaaaaaaaaaaaaaaaaaa    " , tprime_gev

      tp = abs(tprime_gev)      ! just to make sure it's positive

c      if (abs(t_gev)<tp) then
c         write(6,*)' invalid -t>-tprime error',abs(t_gev),tp
c         stop
c      endif

      m_pi0sq= Mpi02/1.e6

* Normally, this would equal m_rho**2, but Laget has adjusted the value
* to reproduce the real photon data.
      lambda0_sq= 0.462         !GeV^2

* Pion saturating Regge trajectory.
* For t>(m_pi**2), the trajectory takes the usual form, 0.7(t-m_pi**2).
* As t --> -infty, the trajectory goes asymptotically to -1.
* Laget does not mention which function he uses for the asymptotic
* behavior, but I find that taking the hyperbolic tangent of the usual
* trajectory expression gives a curve which resembles the one in
* Laget's paper.
      if (t_gev>m_pi0sq) then
         alphapi_t=0.7*(t_gev-m_pi0sq)
      else
         alphapi_t=tanh(0.7*(t_gev-m_pi0sq))
      endif

* Use these instead of the usual Q^2 scaling of the response functions.
      alphapi_0=0.7*(0.-m_pi0sq) 
      lambdapi_sq=lambda0_sq*((1.+alphapi_t)/(1.+alphapi_0))**2
      fpi    = 1./(1.+q2_gev/lambdapi_sq)
      fpi235 = 1./(1.+2.35/lambdapi_sq)
      
* Fit parameters to t-dependence of Laget's p(e,e'p)omega response
* functions at Q^2=2.35, W=2.47 [ub/GeV^2].  Before fitting, I first
* divided the response functions by (fpi/fpi235)^2 and any sin(thetacm)
* factors.
      a = 0.16675 
      b = 0.89524 
      c = 3.5991
      d = 6.4562
      e = -9.6199
      f = 5.8319
      sigL = a*exp(-b*tp)+c*exp(-d*(tp**0.5))+e*exp(-f*(tp**0.25))

c      print*, fpi, fpi235

      sigL = sigL *(fpi/fpi235)**2

c      print*,"asdasd", sigL


      a = -0.12286
      b = 0.56383
      c = 1.4673
      d = 2.1988
      e = 0.65170
      f = 18.501
      sigT = a*exp(-b*tp)+c*exp(-d*(tp**0.5))+e*exp(-f*(tp**2))
      sigT = sigT *(fpi/fpi235)**2

      a = 0.46397
      b = 4.9179 
      c = 3.3023
      d = 3.1741
      sigLT = a*exp(-b*tp)+c*exp(-d*(tp**0.5))
      sigLT = sigLT*(fpi/fpi235)**2*sin(thetacm)
* Laget uses -sqrt(e(1+e)) instead of +sqrt(2e(1+e))
      sigLT = -sigLT/sqrt(2.)

      a = -0.26497
      b = 2.7655 
      c = 2.8034
      d = 3.8586
      sigTT = a*exp(-b*tp)+c*exp(-d*(tp**0.5))
      sigTT= sigTT*(fpi/fpi235)**2*(sin(thetacm))**2

* Since I have nothing better to go on, for now I assume W scales as
* 1/(W^2-mp^2)^2.
      wfactor=(2.47**2-mp**2)**2/(w**2-mp**2)**2

       
c      print*, ":::::::::: ", wfactor  


      sigL = sigL*wfactor
      sigT = sigT*wfactor
      sigLT= sigLT*wfactor
      sigTT= sigTT*wfactor

      sig = sigt +eps*sigl +eps*cos(2.*phicm)*sigtt
     *     +sqrt(2.*eps*(1.+eps))*cos(phicm)*siglt
 

c      print*, 88888888, pi, sig

     

      sig = sig/2./pi/1.d+06      !dsig/dtdphicm in microbarns/MeV^2/rad

c      sig = sig/2./pi      !dsig/dtdphicm in microbarns/GeV^2/rad


      x_mod = sig     



      th_mod=thetacm

c      if (phi.lt.0.3) then
         write(6,102) eps_mod,tm,sigL,sigT,sigTT,sigLT, sig
 102     format('xmodel: eps=',f5.3,' t=',f5.3,' sigL=',f7.2,' sigT=',f6.2,
     1        ' sigTT=',f5.2,' sigLT=',f5.2,' x_mod=',f10.6)
c     endif

      end

c /*--------------------------------------------------*/


      subroutine xmodel_ghm(npol_set,Eb,q2_set,w_gev,q2_gev,u_gev,tm,phicm,
     *     eps,th_mod,x_mod)

      integer npol_set
      real Eb, q2_set, w, q2, phi, th_mod
      real thetacm, x_mod, u_gev
      real q2_gev,w_gev,eps, tm, phicm

      real*8 sig
c      real*8 q2_gev,w_gev,eps,tp
      real*8 a,b,c,d,e,f,fpi,fpi235
      real*8 m_pi0sq
      real*8 pi

      pi = 3.1415926

      call eps_n_theta(npol_set, Eb, w_gev, q2_gev, tm, thetacm, eps)
      
c       print*, "check   ", thetacm
c 
c       stop

c      print*,"*****", npol_set, Eb, w, q2, tm, thetacm, eps_mod


c      eps = eps_mod
c      phicm = phi
c      q2_gev = q2
c      w_gev = w



c     /*--------------------------------------------------*/

c      up = abs(u_gev)      ! just to make sure it's positive

c       a = -0.00025
c       b =  0.051536
c       c =  1.0
c *      d =  0.138
c       d =  14.4
c 
c       sigl = a + b*log(q2_gev) * exp( c + d*log(q2_gev)* up)
c 
c       a = 0.0027333
c       b = 0.051536
c       c = 1.0
c *      d = 0.1863
c       d = 11.54
c 
c       sigt = a + b*log(q2_gev) * exp( c + d*log(q2_gev)* up)
c     
c       print*, q2_set, q2_gev, up
c       print*, c + d*log(q2_gev)* up
c       print*, log(q2_gev), sigt
c *      stop
c 
c 
c 
c       a = 0.00158
c       b = 0.00424 
c       c = 1.04
c 
c       siglt = a + b * exp( c * up)
c 
c       a = 0.00024 
c       b = 0.001 
c       c = 1.1789
c 
c       sigtt = a + b * exp( c * up)

      up = abs(u_gev)      ! just to make sure it's positive


c      a = 0.0495
c      b = 7
c      c = 0.0
c
c      sigl = a * exp( -b * up) + c
c
c      a = 0.088
c      b = 5 
c      c = 0.00022
c
c      sigt = a * exp( -b * up) + c
c
c      a = 0.0015;
c      b = 1.04;
c      c = 0.000158;
c
c      siglt = a * exp( -b * up) + c
c
c      a = 0.001;
c      b = 1.1789;
c      c = 0.00024;
c
c      sigtt = a * exp( -b * up) + c
c


c /*--------------------------------------------------
c // Sigma T
      a = 0.088
      b = 5
      c = 0.00022


      sigt = a * exp( -b * up) + c


c /*--------------------------------------------------
c // Sigma L  
      a = 0.0495
      b = 7
      c = 0.0

      sigl = a * exp( -b * up) + c


c /*--------------------------------------------------
c // Sigma TT  

      a = 0.0;

      sigtt = a * sin(thetacm) *sin(thetacm)

c /*--------------------------------------------------
c // Sigma LT  

      a = 0.0;

      siglt = a * sin(thetacm)

* Since I have nothing better to go on, for now I assume W scales as
* 1/(W^2-mp^2)^2.
      wfactor=(2.47**2-mp**2)**2/(w_gev**2-mp**2)**2

      sigl = sigl*wfactor
      sigt = sigt*wfactor
      siglt= siglt*wfactor
      sigtt= sigtt*wfactor


c      print*, eps, phicm
c     stop


      sig = sigt +eps*sigl +eps*cos(2.*phicm)*sigtt
     *     +sqrt(2.*eps*(1.+eps))*cos(phicm)*siglt

c      sig = sig/2./pi/1.d+06      !dsig/dtdphicm in microbarns/MeV^2/rad

      x_mod = sig     

      th_mod=thetacm

c      if (phi.lt.0.3) then
         write(6,102) eps,up,sigL,sigT,sigTT,sigLT, sig
 102     format('xmodel: eps=',f5.3,' u=',f5.3,' sigL=',f7.2,' sigT=',f6.2,
     1        ' sigTT=',f5.2,' sigLT=',f5.2,' x_mod=',f10.6)
c     endif

      end



c /*--------------------------------------------------*/
c /*--------------------------------------------------*/

      subroutine xmodel_ghm_two_model(npol_set,Eb,q2_set,w_gev,q2_gev,u_gev,tm,phicm,
     *     eps,th_mod,x_mod,um_min)

      implicit none

      integer npol_set
      real Eb, q2_set, w, q2, phi, th_mod
      real thetacm, x_mod, u_gev, Mp, m_p, up, u_min, um_min
      real q2_gev,w_gev,eps, tm, phicm

      real*8 sig
c      real*8 q2_gev,w_gev,eps,tp
      real*8 a,b,c,d,e,f,g,h,fpi,fpi235
      real*8 m_pi0sq
      real*8 pi
      real*8 t0, t1, t2, t3
      real*8 l0, l1, l2, l3
      real*8 lt0,lt1,lt2,lt3
      real*8 tt0,tt1,tt2,tt3
      real*8 sigt,sigl,siglt,sigtt,wfactor

      parameter (Mp=938.27231) ! Proton Mass in MeV
      parameter (m_p=Mp/1000)  ! Proton Mass in GeV
      parameter (pi = 3.1415926)
      
      up = abs(u_gev)      ! just to make sure it's positive

      call eps_n_theta(npol_set, Eb, w_gev, q2_gev, tm, up, u_min, thetacm, eps)
      
      um_min = u_min

      if (up.lt. u_min) then
         u_gev = u_min
         up = abs(u_gev)
      endif

      print*,"w_gev   =   ",  w_gev 
      print*,"q2_gev  =   ",  q2_gev
      print*,"up      =   ",  up
      print*,"phicm   =   ",  phicm
      print*,"thetacm =   ",  thetacm
      print*,"eps     =   ",  eps

cc/*--------------------------------------------------*/
cc/*--------------------------------------------------*/

      if (q2_set.eq.1.6) then

c        print*, "Q2=1.60 parameterization is used" 
        t0  =            7.73587                   
        t1  =            -7.9672                    
        t2  =             0.0000                    
        t3  =             0.0000                    
        l0  =            13.2553                             
        l1  =           -47.2633                               
        l2  =             0.0000                                 
        l3  =             0.0000                                     
        lt0 =            -0.3439                     
        lt1 =             5.9217                     
        lt2 =             0.0000                    
        lt3 =             0.0000                    
        tt0 =             8.1221                    
        tt1 =          -139.8422                    
        tt2 =             0.0000                    
        tt3 =             0.0000                   
              
              
              

      else if (q2_set.eq.2.45) then

c        print*, "Q2=2.45 parameterization is used" 

        t0  =          6.16527  
        t1  =          -4.2124  
        t2  =           0.0000  
        t3  =           0.0000  
        l0  =          12.2546       
        l1  =         -29.8629       
        l2  =           0.0000       
        l3  =           0.0000       
        lt0 =          -0.3620  
        lt1 =           3.1028  
        lt2 =           0.0000  
        lt3 =           0.0000  
        tt0 =          -7.4032  
        tt1 =          63.4705  
        tt2 =           0.0000  
        tt3 =           0.0000  

      else
          print*, "No parameterization is aviliable for Q2=", q2_set
          stop
      endif


c        t0  =  0.05 ;
c        t1  =  0.2  ;
c        t2  =  0.25 ;
c        t3  =  -0.9 ;
c        l0  =  0.7  ;
c        l1  =  -2.7 ;
c        l2  =  -0.5 ;
c        l3  =  2.5  ;
c        lt0 =  0.05 ;
c        lt1 =  -0.3 ;
c        lt2 =  -0.15;
c        lt3 =  0.7  ;
c        tt0 =  0.1  ;
c        tt1 =  -1.8 ;
c        tt2 =  -0.3 ;
c        tt3 =   3   ;

c     /*--------------------------------------------------*/
c     // Sigma T
c      sigt = t1 * exp( -t2 * up) + t3;

c      sigt = t0 + t1 * up + t2 * log(q2_gev) + t3 * up * log(q2_gev)


       

       sigt = t0 / sqrt(q2_gev) + t1 * up / sqrt(q2_gev) + t2 / sqrt(q2_gev) + t3 * up / sqrt(q2_gev)

c      sigt = t0 
       
c     /*--------------------------------------------------*/
c     // Sigma L
c      sigl = l1* (exp(-l2*(up-l3)) + l2*(up-l3));
c      sigl = l1 + l2*up;
c      sigl= l1 * exp( l2 * up ) + l3 / up
c      sigl = l1 * (up-l2)**2 + l3

c      sigl = l0 + l1 * up + l2 * log(q2_gev) + l3 * up * log(q2_gev)
      sigl = l0/(q2_gev*q2_gev) + l1 * up/(q2_gev*q2_gev) + l2 / q2_gev + l3 * up / q2_gev


c     /*--------------------------------------------------*/
c     // Sigma LT  
      siglt = (lt0/q2_gev + lt1 * up/q2_gev + lt2 / q2_gev + lt3 * up / q2_gev) * sin(thetacm)

c     /*--------------------------------------------------*/
c     // Sigma TT  
      sigtt = (tt0/q2_gev + tt1 * up/q2_gev + tt2 / q2_gev + tt3 * up / q2_gev) * sin(thetacm) * sin(thetacm)
c      sigtt = (tt1 * exp( tt2 * up ) + tt3 / up) * sin(thetacm)  * sin(thetacm)

c      print*, "theta check ", thetacm, sin(thetacm)
c      stop

* Since I have nothing better to go on, for now I assume W scales as
* 1/(W^2-mp^2)^2.
      wfactor= 1 / ((w_gev**2-m_p**2)**2)

      sigl = sigl*wfactor
      sigt = sigt*wfactor
      siglt= siglt*wfactor
      sigtt= sigtt*wfactor

c      print*, eps, phicm
c     stop

      sig = sigt + eps*sigl + sqrt(2.*eps*(1.+eps))*cos(phicm)*siglt 
     * + eps*cos(2.*phicm)*sigtt

c      sig = sig/2./pi/1.d+06      !dsig/dtdphicm in microbarns/MeV^2/rad



      x_mod = sig     !2pi * dsig/dtdphicm in microbarns/GeV^2

      th_mod = thetacm



c      if (phi.lt.0.3) then
        

c      write(6,102) eps,up,sigT,sigL,sigLT,sigTT,sig,wfactor,m_p,w_gev
c 102     format('xmodel: eps=',f5.3,' u=',f5.3,' sigT=',f7.5,' sigL=',f7.5,
c     1        ' sigLT=',f10.8,' sigTT=',f10.8,' x_mod=',f15.14,
c     1   ' wfactor=', f10.6, ' mp=',f5.4, ' w_gev=',f6.4)
c

      print*, t0, t1 , q2_gev, up, w_gev, wfactor, eps, sigT, sigL, sig 
c      stop


      print*, "asdasdasdas ", u_min 


      write(6,102) eps,up,sigT,sigL,sigLT,sigTT,sig,u_min
 102     format('xmodel: eps=',f5.3,' u=',f5.3,' sigT=',f7.5,' sigL=',f7.5,
     1        ' sigLT=',f10.8,' sigTT=',f10.8,' x_mod=',f15.14, 
     1        ' umin=', f10.8)





      



      end










c /*--------------------------------------------------*/
c /*--------------------------------------------------*/
c /*--------------------------------------------------*/


      subroutine xmodel_ghm_160(npol_set,Eb,q2_set,w_gev,q2_gev,u_gev,tm,phicm,
     *     eps,th_mod,x_mod)

      integer npol_set
      real Eb, q2_set, w, q2, phi, th_mod
      real thetacm, x_mod, u_gev, mp
      real q2_gev,w_gev,eps, tm, phicm

      real*8 sig
c      real*8 q2_gev,w_gev,eps,tp
      real*8 a,b,c,d,e,f,g,h,fpi,fpi235
      real*8 m_pi0sq
      real*8 pi
      real*8 l1,l2,l3
      real*8 t1,t2,t3
      real*8 lt1,lt2,lt3
      real*8 tt1,tt2,tt3

      pi = 3.1415926
      mp = .93827231
      
      up = abs(u_gev)      ! just to make sure it's positive

      call eps_n_theta(npol_set, Eb, w_gev, q2_gev, tm, up, thetacm, eps)


c      a =  0.18881E+00
c      b =  0.14443E+02
c      c =  0.90201E-02
c      d = -0.75302E-01
c      e = -0.10392E+01 
c      f =  0.11590E+00
c      g = 0.0;
c      h = 0.0;    

c      a =  0.20217E+00;
c      b =  0.17986E+02;
c      c =  0.60942E-02;
c      d =  0.0;
c      e =  0.0;
c      f =  0.11590E+00;
c      g =  0.21197E+00;
c      h = -0.24716E-01;
c
cc /*--------------------------------------------------
cc // Sigma T
c      sigt = a * exp( -b * up) + c
c
cc /*--------------------------------------------------
cc // Sigma L
c      sigl = d * exp( -e * up) + f
c
cc /*--------------------------------------------------
cc // Sigma TT  
c      sigtt = g * sin(thetacm) *sin(thetacm)
c
cc /*--------------------------------------------------
cc // Sigma LT  
c      siglt = h * sin(thetacm)
c


c       /*--------------------------------------------------*/
       t1  =    0.1     
       t2  =   -5.0       
       t3  =    0.0002  
       l1  =    0.05    
       l2  =    0.6     
       l3  =    0.02    
       lt1 =    0.000000      
       lt2 =    0.000000      
       lt3 =    0.000000      
       tt1 =    0.000000      
       tt2 =    0.000000      
       tt3 =    0.000000      


cc /*--------------------------------------------------
cc // Sigma T
c      sigt = t1 * exp( -t2 * up) + t3
c
cc /*--------------------------------------------------
cc // Sigma L
c      sigl = l1 * exp( l2 * up) + l3 
c
cc /*--------------------------------------------------
cc // Sigma TT  
c      sigtt = tt1 * exp( -tt2 * up) * sin(thetacm) * sin(thetacm)
c
cc /*--------------------------------------------------
cc // Sigma LT  
c      siglt = lt1 * exp( -lt2 * up) * sin(thetacm)



c     /*--------------------------------------------------*/
c     // Sigma T
c      sigt = t1 * exp( -t2 * up) + t3;
      sigt = t1 * exp(t2 * up) + t3
       
c     /*--------------------------------------------------*/
c     // Sigma L
c      sigl = l1* (exp(-l2*(up-l3)) + l2*(up-l3));
 
c      sigl = l1 + l2*up;

c      sigl= l1 * exp( l2 * up ) + l3 / up

      sigl = l1 * (up-l2)**2 + l3




c     /*--------------------------------------------------*/
c     // Sigma LT  
      siglt = (lt1 * exp( lt2 * up ) + lt3 / up) * sin(thetacm);
c      siglt = 0


c     /*--------------------------------------------------*/
c     // Sigma TT  
      sigtt = tt1 * sin(thetacm) * sin(thetacm);
      
c      sigtt = 0 

c      print*, "theta check ", thetacm, sin(thetacm)

c      stop



* Since I have nothing better to go on, for now I assume W scales as
* 1/(W^2-mp^2)^2.
      wfactor=(2.47**2-mp**2)**2/(w_gev**2-mp**2)**2

      sigl = sigl*wfactor
      sigt = sigt*wfactor
      siglt= siglt*wfactor
      sigtt= sigtt*wfactor


c      print*, eps, phicm
c     stop


      sig = sigt +eps*sigl +eps*cos(2.*phicm)*sigtt
     *     +sqrt(2.*eps*(1.+eps))*cos(phicm)*siglt

c      sig = sig/2./pi/1.d+06      !dsig/dtdphicm in microbarns/MeV^2/rad

      x_mod = sig     

      th_mod=thetacm

c      if (phi.lt.0.3) then
         write(6,102) eps,up,sigT,sigL,sigLT,sigTT, sig
 102     format('xmodel: eps=',f5.3,' u=',f5.3,' sigT=',f7.2,' sigL=',f6.2,
     1        ' sigLT=',f5.2,' sigTT=',f5.2,' x_mod=',f10.6)
c     endif

      end






c /*--------------------------------------------------*/

      subroutine xmodel_ghm_245(npol_set,Eb,q2_set,w_gev,q2_gev,u_gev,tm,phicm,
     *     eps,th_mod,x_mod)

      integer npol_set
      real Eb, q2_set, w, q2, phi, th_mod
      real thetacm, x_mod, u_gev
      real q2_gev,w_gev,eps, tm, phicm

      real*8 sig
c      real*8 q2_gev,w_gev,eps,tp
      real*8 a,b,c,d,e,f,g,h,fpi,fpi235
      real*8 l1,l2,l3
      real*8 t1,t2,t3
      real*8 lt1,lt2,lt3
      real*8 tt1,tt2,tt3
      real*8 m_pi0sq
      real*8 pi

      pi = 3.1415926

      
      up = abs(u_gev)      ! just to make sure it's positive

      call eps_n_theta(npol_set, Eb, w_gev, q2_gev, tm, up, thetacm, eps)

c      a =  0.82472E-01
c      b =  0.71359E+01
c      c =  0.13529E-02
c      d =  0.15738E+00
c      e =  0.69794E+01 
c      f =  0.23247E-02
c      g =  0.0
c      h =  0.0

c       /*--------------------------------------------------*/
c       a =  0.93161E-01;
c       b =  0.78498E+01;
c       c =  0.15340E-02;
c       d =  0.0;
c       e =  0.0; 
c       f =  0.0;
c       g =  0.0;
c       h =  0.0;
c
c
c
cc /*--------------------------------------------------
cc // Sigma T
c      sigt = a * exp( -b * up) + c
c
cc /*--------------------------------------------------
cc // Sigma L
c      sigl = d * exp( -e * up) + f
c
cc /*--------------------------------------------------
cc // Sigma TT  
c      sigtt = g * sin(thetacm) *sin(thetacm)
c
cc /*--------------------------------------------------
cc // Sigma LT  
c      siglt = h * sin(thetacm)



c       /*--------------------------------------------------*/
c       t1  =   0.27099E+00
c       t2  =   0.86639E+01
c       t3  =   0.00000E+00
c       l1  =   0.0           
c       l2  =   0.0           
c       l3  =   0.0           
c       tt1 =   0.0           
c       tt2 =   0.0           
c       lt1 =   0.0           
c       lt2 =   0.0           



       t1  =    0.1     
       t2  =   -5.0     
       t3  =    0.0002  
       l1  =    0.05    
       l2  =    0.6     
       l3  =    0.02    
       lt1 =    0.000000
       lt2 =    0.000000
       lt3 =    0.000000
       tt1 =    0.000000
       tt2 =    0.000000
       tt3 =    0.000000
                                  
                                  
                           
cc /*-----------------------------------------------------------
cc // Sigma T            gma T
c      sigt = t1 * exp( -= t1 * ext2 * up) + t3
c                                 
cc /*-----------------------------------------------------------
cc // Sigma L            gma L
c      sigl = l1 * exp( l= l1 * ex2 * up) + l3 
c                                 
cc /*-----------------------------------------------------------
cc // Sigma TT           gma TT  
c      sigtt = tt1 * exp( -tt2 * up) * sin(thetacm) * sin(thetacm)
c
cc /*--------------------------------------------------
cc // Sigma LT  
c      siglt = lt1 * exp( -lt2 * up) * sin(thetacm)
c

c
cc     /*--------------------------------------------------*/
cc     // Sigma T
c      sigt = t1 * exp( -t2 * up) + t3;
c       
cc     /*--------------------------------------------------*/
cc     // Sigma L
c      sigl = l1* (exp(-l2*(up-l3)) + l2*(up-l3));
c 
c
cc     /*--------------------------------------------------*/
cc     // Sigma LT  
c      siglt = (lt1 * up + lt2) * sin(thetacm);
c      
cc     /*--------------------------------------------------*/
cc     // Sigma TT  
c      sigtt = tt1 * sin(thetacm) * sin(thetacm);



c     /*--------------------------------------------------*/
c     // Sigma T
c      sigt = t1 * exp( -t2 * up) + t3;

      sigt = t1 * exp(t2 * up) + t3
       
c     /*--------------------------------------------------*/
c     // Sigma L
c      sigl = l1* (exp(-l2*(up-l3)) + l2*(up-l3));
c      sigl = l1 + l2*up;
c      sigl= l1 * exp( l2 * up ) + l3 / up

      sigl = l1 * (up-l2)**2 + l3

c     /*--------------------------------------------------*/
c     // Sigma LT  
      siglt = (lt1 * exp( lt2 * up ) + lt3 / up) * sin(thetacm);
c      siglt = 0


c     /*--------------------------------------------------*/
c     // Sigma TT  
      sigtt = tt1 * sin(thetacm) * sin(thetacm);




      print*, "aaa " , siglt 

c      stop

* Since I have nothing better to go on, for now I assume W scales as
* 1/(W^2-mp^2)^2.
      wfactor=(2.47**2-mp**2)**2/(w_gev**2-mp**2)**2

      sigl = sigl*wfactor
      sigt = sigt*wfactor
      siglt= siglt*wfactor
      sigtt= sigtt*wfactor


c      print*, eps, phicm
c     stop


      sig = sigt +eps*sigl +eps*cos(2.*phicm)*sigtt
     *     +sqrt(2.*eps*(1.+eps))*cos(phicm)*siglt

c      sig = sig/2./pi/1.d+06      !dsig/dtdphicm in microbarns/MeV^2/rad

      x_mod = sig     

      th_mod=thetacm

c      if (phi.lt.0.3) then
         write(6,102) eps,up,sigT,sigL,sigLT,sigTT, sig
 102     format('xmodel: eps=',f5.4,' u=',f5.4,' sigT=',f9.6,' sigL=',f9.6,
     1        ' sigLT=',f9.6,' sigTT=',f9.6,' x_mod=',f9.6)
c     endif

        
c      stop

      end






*=======================================================================

      include 'eps_n_theta.f'

