C---------------------------------------------------
C                           Spectrum of perturbations: P(k), sigma_8,V(50h-1Mpc)
C     Anatoly Klypin December 1997
C                                                       Holtzman           appoximation  
C                                             or       BBKS-style + Hu&Sugiyama approx
C                            Needs cdm.fit file -- parameters of models and fitting 
C-----------------------------------------------------
C                            Om   = Omega_0                Omc   = Omega_cdm_present
C                            Omb = Omega_baryon     Omnu =Omega_neutrino_present
C                            Ocurv = Omega_curvature_present
C                            hsmall = H_0/100km/s/Mpc = the Hubble constant
C                            ns    = slope of the power spectum (ns=1 for Harr-Zeld)
C                            AEXPN = expansion parameter = 1/(1+z), z =redshift
C                            All calculations are done in real Mpc,
C                            Only final outputs are scaled to  h^{-1} 
C-----------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NtabM = 100000)
      PARAMETER (Nm    = 100)
      PARAMETER (Na    = 10)
      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Ocurv,Par(6),AEXPN,Sn,ns,hsmall
      COMMON / Coeff/   deltac22,deltac,rf,Uklow,Ukup,GrowthDen_0
      COMMON /PTABLE/xkt(0:NtabM),Pkt(0:NtabM),StepK,alog0,iFlag,Ntab 
      CHARACTER         Name*16,Answer*1,Header*45,StartFile*16
      REAL*8              INTG,ns,INPUT,P2,P,Ptophat,Pgauss,Vtophat
      DIMENSION       aNdens(Nm,Na),aMasses(Nm),aRedshifts(Na)
      EXTERNAL           INTG,P2,P,Ptophat,Pgauss,Vtophat
      DATA                StartFile/'InStart.dat'/
      WRITE (*,*) '================< Enter File Name for output file: '
      READ  (*,'(a)') Name
      OPEN(1,file=Name)
      OPEN(10,file='number'//Name)
      OPEN(2,file='cdm.fit') !  parameters of models
      
C                                                     Constants
      PI     =3.1415926535
      AEXPN =1.
      WRITE (*,'(A,$)') ' Enter sigma_8 and slope n: '
      READ  (*,*) Sigma8,ns
      CALL MODEL(Line)
      H      =100.*hsmall     ! Hubble constant
      DO a=0.1,1.,0.1
         z =1./a-1.
         CALL AGE(t0,GrowthDen,GrowthVel,a)
         write (*,20) t0,z,a,GrowthDen,GrowthVel
         write (1,20) t0,z,a,GrowthDen,GrowthVel
         write (10,20) t0,z,a,GrowthDen,GrowthVel
 20      format(' Age=',f8.3,' z=',f6.2,' a=',f6.3,
     .          ' GrowthRateDen=',f7.4,
     .          ' GrowthRateVelocity=',f7.4)
      ENDDO
      AEXPN =1.
      a =AEXPN
      z =1./a-1.
      write (*,*) ' a=',a,' Aexpn=',AEXPN
      CALL AGE(t0,GrowthDen_0,GrowthVel_0,a)
         write (*,20) t0,z,a,GrowthDen,GrowthVel
         write (1,20) t0,z,a,GrowthDen,GrowthVel
         write (10,20) t0,z,a,GrowthDen,GrowthVel
 
      rf    =8./hsmall        ! top-hat for sigma_8
      Sn    = (Sigma8)**2 / 
     &          INTG(Ptophat,1.d-5,5.d+1) ! normalization of P(k)
c                                             bias parameter      
      Sigma8 = sqrt(Sn*(
     &          INTG(Ptophat,1.d-5,5.d+1)))
      write (*,10) hsmall,Sigma8,ns  ! check if all is self-consistent
      write (1,10) hsmall,Sigma8,ns
      write (10,10) hsmall,Sigma8,ns
 10   Format(' Hubble=',f7.3,' Sigma_8 =',G11.3,
     &       ' Slope n=',f7.2)

c                       Evolution of Sigma(mass,a),Fraction_collapsed(a)
      write(1,15)'    Z       a'
 15   format(a,T17,3(3x,'Mass',7x,'Sigma',
     &       '    FracMass'))
      Do i=0,30
         z =i/3.
         AEXPN    = 1./(1.+z)
         aMass1   = 1.e10
         ss1      = SigMass(aMass1)
         fr_mass1 = 0.5*derfc(1.68d0/sqrt(2.)/ss1)
         aMass2   = 1.e11
         ss2      = SigMass(aMass2)
         fr_mass2 = 0.5*derfc(1.68d0/sqrt(2.)/ss2)
         aMass3   = 1.e12
         ss3      = SigMass(aMass3)
         fr_mass3 = 0.5*derfc(1.68d0/sqrt(2.)/ss3)

         write (1,'(2f8.3,3(g11.3,f9.3,g11.3))')z,AEXPN,
     &          aMass1,ss1,fr_mass1,
     &          aMass2,ss2,fr_mass2,aMass3,ss3,fr_mass3
      EndDo 

C                      Mass function at different redshifts
      Do i=1,Nm
         aMasses(i) = 1.d10*10.**(i/float(Nm)*4.)
      EndDo 
      Do i=1,Na
         aRedshifts(i) = (i/float(Na)*10.)
      EndDo 
      aRedshifts(1) =0.

      Do j=1,Na
         z     = aRedshifts(j)
         AEXPN = 1./(1.+z) 
      Do i=1,Nm
         aMass = aMasses(i)
         aNdens(i,j) = aNumber(aMass)
      EndDo 
      EndDo 
      write (10,'(2a)')' M(sun/h)','  Redshifts'
      write (10,'(11x,16g11.3)') aRedshifts 
      Do i=1,Nm
         write (10,'(g11.3,16g11.3)')
     &    aMasses(i)*hsmall,(aNdens(i,j)/hsmall**3,j=1,Na)
      EndDo 


      rf =50./hsmall           ! top-hat for bulk velocity
      Veloc = GrowthVel*H*
     &           sqrt(Sn*(INTG(Vtophat,1.d-5,1.d+0)))      
      write (*,30) Veloc,rf*hsmall
      write (1,30) Veloc,rf*hsmall
 30   Format(' Bulk Velocity =',F8.2,
     .       'km/s for radius top-hat=',f8.1,'h-1Mpc')
      write (1,*) ' k/h        Pk*h^3   Power Spectrum at z=',1/AEXPN-1.
      DO i=1,110
        w =1.e-4*10.**(i/20.)
        Pk0  =P(w)*Sn*hsmall**3*(2.*pi**2) 
        write (1,14) w/hsmall,Pk0
      ENDDO
 14   Format(2G11.4,5G12.4)
C..........................................................Power spectrum and  
C                                                         Amplitude of fluctuations in a Box at redshift z
      WRITE (*,'(A,$)')
     &         'Enter Box(h^-1Mpc), Nparticles(1D) and redshift:'
      READ  (*,*) Box, NROW, z
      Box   =Box/hsmall           ! scale it to real Mpc
      Uklow =2.*pi/Box            ! frequency range for integrals
      Ukup  =2.*pi/Box*(NROW/2.)
      AEXPN =1./(1.+z)            ! expansion parameter
      a     =AEXPN
      CALL AGE(t0,GrowthDen,GrowthVel,a)
      sigma =(GrowthDen/GrowthDen_0)*
     &             sqrt(Sn*INTG(P2,Uklow/sqrt(2.),Ukup))
      WRITE (*,40) z,sigma,NROW,Box*hsmall
 40   format('  z=',f8.3,' delta\rho/rho in box=',f9.5,/
     .       5X,'Particles=',i4,' Box=',f8.2)
      write (1,*) ' k/h        Pk*h^3   Power Spectrum at z=',1/AEXPN-1.
      DO i=1,130
        w =1.e-4*10.**(i/20.)
        Pk0  =(GrowthDen/GrowthDen_0)**2*P(w)*Sn*hsmall**3*(2.*pi**2) 
        write (1,14) w/hsmall,Pk0
      ENDDO
C............................................................................................      
      WRITE (*,'(A,$)') ' Do you need a file to run PMstart (Yes/No) '
      READ  (*,'(a)') Answer
      IF(Answer.eq.'Y'.or.Answer.eq.'y')Then ! prepare input file
         OPEN(30,file=StartFile)        !  parameters of a run
         write(30,'(A)') 'Yes           Read parameters from cdm.fit'
         HEADER='N=128x256 L=20h-1CDMsig=0.239 z=30-----------'
         write (*,*) '------ Enter Header up to 45 characters'
         write (*,*) '       Example:'
         write (*,'(A)') HEADER 
         read  (*,'(A)') HEADER
         write (30,'(A)') HEADER
         write (30,*) AEXPN,'    Expansion Parameter'
         ASTEP =INPUT(' Step of the expansion parameter  =')
         write (30,*) ASTEP,'    Step in dAEXPN    '
         write (30,*) sigma,'    DRho/rho in box   '
         write (30,*) ns,   '    Slope of P(k)     '
         write (30,*) Box*hsmall,'    Box in  Mpc/h   '
         write (30,*) Ocurv,'    Omega_curvature   '
         Nseed =INPUT(' Enter random seed integer 1-2^31 =')
         write (30,*) Nseed,'    Random seed       '
         write (30,*) Line,'     Line number in cdm.fit'
         Write (*,'(A,$)') ' Enter Min and Max number of mass species ='
         read (*,*) Nmin_lev,Nmax_lev       
         write (30,*) Nmin_lev,Nmax_lev,' Min and Max number of species'
         write (*,*) ' Results were written to file: ',StartFile
         CLOSE (30)
      ENDIF
      STOP
	   END
c-----------------------------------  Pcold(k)*k^2
      REAL*8 FUNCTION P2(WK)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Ocurv,Par(6),AEXPN,Sn,ns,hsmall
      COMMON / Coeff/   deltac22,deltac,rf,Uklow,Ukup,GrowthDen_0
      Real*8 ns
        P2=WK**2*P(WK)
      RETURN
      END

C---------------------------------------------
C                                                                 Power spectrum
C                                                                wk = k= in real Mpc
      REAL*8 FUNCTION P(wk)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NtabM = 100000)
      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Ocurv,Par(6),AEXPN,Sn,ns,hsmall
      COMMON / FACT/    qqscaleb
      COMMON /PTABLE/xkt(0:NtabM),Pkt(0:NtabM),StepK,alog0,iFlag,Ntab  

      Real*8 ns
      If(iFlag.gt.0)Then       ! use approximations
      If(Par(6).ne.0.)Then   ! Holtzman approx
        sk= sqrt(wk)
        P = wk**ns /
     .         (1.+sk*(Par(2)
     .            +sk*(Par(3)
     .            +sk*(Par(4)
     .            +sk*(Par(5) )))) )**(2.*Par(6))
      Else                   ! BBKS + Sugiyama approx
c        Gammaeff =Om*hsmall/exp(Omb*(1.+sqrt(hsmall/0.5)/Om))
c        Q = wk/hsmall/Gammaeff
        Q = wk*qqscaleb
        
        P = wk**ns*( LOG(1.+Par(1)*Q)/(Par(1)*Q) )**2/
     .          sqrt(1.+Q*(Par(2)+
     .                  Q*(Par(3)**2+
     .                  Q*(Par(4)**3+
     .                  Q*(Par(5)**4) ))))
      EndIf
      Else                ! use table for P(k)
         P = Ppk(wk)
      EndIf 
      RETURN
      END
C-------------------------------------------------
C                                  Age of the Universe: t0 (z=0)
C                                  Expansion parameter: a =1/(1+z)
C                                  Growth_Rate_Density at a: GrowthDen
C                                     GrowthDen =\delta\rho/\rho
C                                     normalized to GrowthDen =a for Omega=1
C                                     GrowthDen < 1 at a=1 for Lcdm or OCDM
C                                  Growth_Rate_velocity at a: GrowthVel
C                                     GrowthVel = V_peculiar= (a/delta)(d delta/d a)
      SUBROUTINE AGE(t0,GrowthDen,GrowthVel,a)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (zero =1.d-12)
      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Ocurv,Par(6),AEXPN,Sn,ns,hsmall
      Common /OmegRat/ OmLOm0,OmcOm0
      Real*8  INTG,ns
      External Hnorm,Hage,Hgrow

      Oml    =Max(Abs(1.-Om-Ocurv),zero)
      OmLOm0 =Oml/Om
      OmcOm0 =Ocurv/Om
      t0     =9.766/hsmall/sqrt(Om)*INTG(Hage,zero,a)
      z      =1./a-1.
      Hubble = hsmall*sqrt(Om/a**3)*Hnorm(a)
         ww  =INTG(Hgrow,zero,a)
      GrowthVel = -(1.5+OmcOm0*a)/Hnorm(a)**2+
     &                 sqrt(a**5)/Hnorm(a)**3/ww
      GrowthDen =2.5*Hnorm(a)/sqrt(a**3)*ww
c         write (*,20) t0,z,a,Hubble,GrowthDen,GrowthVel
c         write (1,20) t0,z,a,Hubble,GrowthDen,GrowthVel
c 20      format(' Age=',f8.3,' z=',f6.2,' a=',f6.3,
c     .          ' Hubble/100=',f8.3,' GrowthRateDen=',f7.4,
c     .          ' GrowthRateVelocity=',f7.4)
      Return
      End
c-----------------------------------  P(k)*k^2*gauss(rf)
      REAL*8 FUNCTION Pgauss(WK)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Ocurv,Par(6),AEXPN,Sn,ns,hsmall
      COMMON / Coeff/   deltac22,deltac,rf,Uklow,Ukup,GrowthDen_0
      Real*8 ns
        Pgauss=WK**2*P(WK)*EXP(-(rf*WK)**2)
      RETURN
      END
c-----------------------------------  P(k)*k^4*gauss(rf)
      REAL*8 FUNCTION P4gauss(WK)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Ocurv,Par(6),AEXPN,Sn,ns,hsmall
      COMMON / Coeff/   deltac22,deltac,rf,Uklow,Ukup,GrowthDen_0
      Real*8 ns
        P4gauss=WK**4*P(WK)*EXP(-(rf*WK)**2)
      RETURN
      END
C-------------------------------------- P*k^2*Top-Hat Filter
      REAL*8 FUNCTION Ptophat(wk)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Ocurv,Par(6),AEXPN,Sn,ns,hsmall
      COMMON / Coeff/   deltac22,deltac,rf,Uklow,Ukup,GrowthDen_0
      Real*8 ns
	   IF (wk.lt.1.d-4) THEN
            Ptophat =wk**ns*wk**2
        ELSE
            X      =wk*rf
	         TOPHAT =( (SIN(X)-x*COS(X))*3./X**3 )**2
            Ptophat=P(wk)*wk**2*TOPHAT
        ENDIF
      RETURN
      END
C-------------------------------------- Number-density of halos
C                              Gives M*dN/dM in units Mpc**-3
C                              aMass is in real Msun units
c                              Calls function Fn(sigma)
      REAL*8 FUNCTION aNumber(aMass)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (pi = 3.1415926535d0)
      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Ocurv,Par(6),AEXPN,Sn,ns,hsmall
      COMMON / Coeff/   deltac22,deltac,rf,Uklow,Ukup,GrowthDen_0
      REAL*8            INTG,ns,INPUT,Ptophat
      EXTERNAL          INTG,Ptophat

      rf = (aMass/(1.154e12*Om*hsmall**2))**0.33333333
      a  = AEXPN

      CALL AGE(t0,GrowthDen,GrowthVel,a)
      Sig =(GrowthDen/GrowthDen_0)*
     &           sqrt(Sn*INTG(Ptophat,1.d-4,2.d+2))
      aM1 = aMass*1.03
      rf  = (aM1/(1.154e12*Om*hsmall**2))**0.33333333
      Sig1=(GrowthDen/GrowthDen_0)*
     &           sqrt(Sn*INTG(Ptophat,1.d-4,2.d+2))

      aM0 = aMass/1.03
      rf  = (aM0/(1.154e12*Om*hsmall**2))**0.33333333
      Sig0=(GrowthDen/GrowthDen_0)*
     &           sqrt(Sn*INTG(Ptophat,1.d-4,2.d+2))

      derivative = (Sig0-Sig1)/Sig/(aM1-aM0)
      factor     = Fn(Sig)
      aNumber    = 2.75e11*(Om*hsmall**2)*
     &             derivative * factor

c      write (*,'(9g12.4)')aMass,Sig,Sig0,Sig1,
c     &          derivative,factor,aNumber 
      if(aNumber.lt.0..or.derivative.lt.0.)
     &    write (*,'(a,12g12.4)')' Error PS:',aMass,Sig,Sig1,Sig0,
     &           derivative,factor,aNumber,aM0,aM1
      RETURN
      END
C-------------------------------------- Auxiliary for aNumber()
C                              Gives dependance on sigma(M)
C
      REAL*8 FUNCTION Fn(sigma)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (pi = 3.1415926535d0)
      PARAMETER (Deltac = 1.68)

      x  = Deltac/sigma
      Fn =sqrt(2./pi)*x*exp(-0.5*x**2)

      RETURN
      END
 
C-------------------------------------- rms fluctuations
C                                       for mass aMass
C                              aMass is in real Msun units
      REAL*8 FUNCTION SigMass(aMass)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (pi = 3.1415926535d0)
      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Ocurv,Par(6),AEXPN,Sn,ns,hsmall
      COMMON / Coeff/   deltac22,deltac,rf,Uklow,Ukup,GrowthDen_0
      REAL*8            INTG,ns,INPUT,Ptophat
      EXTERNAL          INTG,Ptophat

            rf= (aMass/(1.154e12*Om*hsmall**2))**0.33333333
            a     =AEXPN
      CALL AGE(t0,GrowthDen,GrowthVel,a)
      SigMass =(GrowthDen/GrowthDen_0)*
     &           sqrt(Sn*INTG(Ptophat,1.d-5,2.d+2))

      RETURN
      END

C-------------------------------------- P*Top-Hat Filter
      REAL*8 FUNCTION Vtophat(wk)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Ocurv,Par(6),AEXPN,Sn,ns,hsmall
      COMMON / Coeff/   deltac22,deltac,rf,Uklow,Ukup,GrowthDen_0
      Real*8 ns
   	IF (wk.lt.1.d-4) THEN
            Vtophat =wk**ns
        ELSE
            X=wk*rf
	    TOPHAT=( (SIN(X)-x*COS(X))*3./X**3 )**2
            Vtophat=P(wk)*TOPHAT
        ENDIF
      RETURN
      END
C-------------------------------------- Linear Pair-wise V
      REAL*8 FUNCTION Vpairwise(wk)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Ocurv,Par(6),AEXPN,Sn,ns,hsmall
      COMMON / Coeff/   deltac22,deltac,rf,Uklow,Ukup,GrowthDen_0
      Real*8 ns
            X        =wk*rf
            Vpairwise=P(wk)*(1.-sin(X)/X)
      RETURN
      END
C________________________________________Read parameters of the model
C                                 Om       = Omega_0
C                                 Omb     = Omega_baryon
C                                 Omnu  = Omega_neutrino
C                                 hsmall  = hubble constant
C                                 Par      = fitting parameters
C                                 Par(6) = 0  --- bbks-style+Hu&Sugiyama
C                                 Par(6) ne 0 --- holtzman-style
      SUBROUTINE MODEL(Line)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NtabM = 100000)
      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Ocurv,Par(6),AEXPN,Sn,ns,hsmall
      COMMON / FACT/    qqscaleb
      COMMON /PTABLE/xkt(0:NtabM),Pkt(0:NtabM),StepK,alog0,iFlag,Ntab  

      Real*8 ns,INPUT
      Character             Header*79      

      Line =INPUT(' Enter Line Number in  cdm.fit=   ')
      If(Line.gt.0)Then
         iFlag = Line
      Read(2,'(A)') Header
      If(Line.gt.2)Then
         Do i=1,Line -2
            Read (2,*) a
         EndDo
      EndIf
      Read (2,*) Om,Omb,Omc,Omnu,hsmall,Par
       
      Ocurv =INPUT(' Enter Omega curvature at z=0  =>')
      theta = 2.726/2.7  ! = T_cmb/2.7
      Ob0   =Omb/Om      ! Hu&Sugiyama fitting
      Omh2  =Om*hsmall**2
      a1    =(46.9*Omh2)**0.670*(1.+(32.1*Omh2)**(-0.532))
      a2    =(12.0*Omh2)**0.424*(1.+(45.*Omh2)**(-0.582))
      alpha =1./a1**Ob0/a2**(Ob0**3)
c      qqscaleb = theta**2/(1.-Ob0)**0.60/sqrt(alpha)/Omh2
      qqscaleb = 1./Omh2*exp(Omb*(1.+sqrt(2*hsmall)/Om))

      Write (*,20)Om,Omb,Omc,Omnu,hsmall,Par
      Write (1,20)Om,Omb,Omc,Omnu,hsmall,Par
      Write (10,20)Om,Omb,Omc,Omnu,hsmall,Par
 20   Format(' Model: Om_0=',F5.3,' Om_baryon=',F5.3,
     .       ' Om_cold=',F5.3,' O_nu=',F5.3,
     .       ' hsmall=',F4.2,/8x,'Parameters=',(6G9.3))
      Else
         iFlag = Line
         If(iFlag.lt.-3)Stop ' Error: model Line # less than -3'
         If(iFlag.eq.-1)Then  ! W.Hu table
            hsmall = 0.7
            Omb = 0.024/hsmall**2
            Om   = 0.3
            Ocurv = 0.
            Omc  = Om
         EndIf
         If(iFlag.eq.-2)Then  ! Eis.Hu approximation
            hsmall = 0.7
            Omb = 0.04
            Om   = 0.3
            Ocurv = 0.
            Omc  = Om
         EndIf
         If(iFlag.eq.-3)Then  ! Eis.Hu approximation: WMAP
            hsmall = 0.73
            Omb = 0.04
            Om   = 0.24
            Ocurv = 0.
            Omc  = Om
         EndIf
 
      CALL ReadTable(hsmall)

      Write (*,22)Om,Omb,hsmall
      Write (1,22)Om,Omb,hsmall
 22   Format(' Model: Om_0=',F5.3,' Om_baryon=',F5.3,
     .       ' hsmall=',F4.2, ' <== Using table for P(k)')
      EndIf
      Return
      End
C-------------------------------------------- Simpson integration
      REAL*8 FUNCTION INTG(FUNC,A,B)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (EPS=5.0d-4, JMAX=22)
      EXTERNAL FUNC
      OST=-1.D30
      OS= -1.D30
      ST =0.
      DO 11 J=1,JMAX
        CALL TRAPZD(FUNC,A,B,ST,J)
        INTG=(4.0d0*ST-OST)/3.0d0
        IF (ABS(INTG-OS).Le.EPS*ABS(OS)) RETURN
        OS=INTG
        OST=ST
11    CONTINUE
      WRITE (*,*)'Integration did not converge'
      RETURN
      END
C----------------------------------------------
      SUBROUTINE TRAPZD(FUNCC,A,B,S,N)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
        SAVE IT
        EXTERNAL FUNCC
      IF (N.EQ.1) THEN
        S=0.5d0*(B-A)*(FUNCC(A)+FUNCC(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5D0*DEL
        SUM=0.0D0
        DO 11 J=1,IT
          SUM=SUM+FUNCC(X)
          X=X+DEL
11      CONTINUE
        S=0.5D0*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END

C-------------------------------------------- Simpson integration
      REAL*8 FUNCTION INTG1(FUNC,A,B)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (EPS=1.0d-3, JMAX=22)
      EXTERNAL FUNC
      OST=-1.d30
      OS= -1.d30
      ST =0.
      DO 11 J=1,JMAX
        CALL TRAPZD1(FUNC,A,B,ST,J)
        INTG1=(4.0D0*ST-OST)/3.0D0
        IF (ABS(INTG1-OS).Le.EPS*ABS(OS)) RETURN
        OS=INTG1
        OST=ST
11    CONTINUE
c      WRITE (1,*)'Integration did not converge'
      RETURN
      END
C----------------------------------------------
      SUBROUTINE TRAPZD1(FUNCC,A,B,S,N)
C---------------------------------------
	IMPLICIT Real*8 (A-H,O-Z)
        SAVE IT
        EXTERNAL FUNCC
      IF (N.EQ.1) THEN
        S=0.5D0*(B-A)*(FUNCC(A)+FUNCC(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5D0*DEL
        SUM=0.0D0
        DO 11 J=1,IT
          SUM=SUM+FUNCC(X)
          X=X+DEL
11      CONTINUE
        S=0.5D0*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END
C---------------------------------- Read in variables      
      REAL*8 FUNCTION INPUT(text)
C-----------------------------------------
      Character text*(*)
          write (*,'(A,$)')text
          read (*,*) x
          INPUT =x
      Return
      End
C-----------------------------------------
      REAL*8 FUNCTION sinhh(x)
C-----------------------------------------
      IMPLICIT REAL*8 (A-H,P-Z)
         sinhh =0.5*(exp(MIN(x,33.d+0))-exp(-MIN(x,33.d+0)))
      RETURN
      END
C---------------------------------------
      REAL*8 FUNCTION Hh(x)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
         Hh =(sqrt(x/(1.d+0 +x**3)))**3 
      Return
      End
C---------------------------------------
      REAL*8 FUNCTION Hnorm(x)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      Common /OmegRat/ OmLOm0,OmcOm0
         Hnorm =sqrt(1.d+0 
     &                  +OmLOm0*x**3
     &                  +OmcOm0*x )
      Return
      End
C---------------------------------------
      REAL*8 FUNCTION Hage(x)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
         Hage =sqrt(x)/Hnorm(x) 
      Return
      End
C---------------------------------------
      REAL*8 FUNCTION Hgrow(x)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
         Hgrow =(sqrt(x)/Hnorm(x))**3 
      Return
      End
C---------------------------------------
      REAL*8 FUNCTION Ppk(x)
c                     interpolate table with p(k)
c                       x is in real 1/Mpc
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NtabM = 100000)
      COMMON /PTABLE/xkt(0:NtabM),Pkt(0:NtabM),StepK,alog0,iFlag,Ntab  

      If(x.ge.xkt(Ntab))Then  ! slope is ns =-3
         Ppk =Pkt(Ntab)/(x/xkt(Ntab))**3
         Return
      EndIf
      If(x.lt.xkt(1))Then                ! slope is ns=1
         Ppk =Pkt(1)*(x/xkt(1))
         Return
      EndIf
      ind = INT((log10(x)-alog0)/StepK) +1
      dk  = xkt(ind+1)-xkt(ind)
        Ppk   = (Pkt(ind)*(xkt(ind+1)-x)+Pkt(ind+1)*(x-xkt(ind)))/dk
c        write (*,'(5x," k=",g12.4," table=",2g12.4,"  P=",3g12.4)') 
c     &    x,xkt(ind),xkt(ind+1),
c     &    Ppk,Pkt(ind),Pkt(ind+1)
      Return
      End
C---------------------------------------
      SUBROUTINE  ReadTable
c                     interpolate table with p(k)
C---------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NtabM = 100000)
      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Ocurv,Par(6),AEXPN,Sn,ns,hsmall
      COMMON /PTABLE/xkt(0:NtabM),Pkt(0:NtabM),StepK,alog0,iFlag,Ntab  
      REAL*8              ns
      If(iFlag.eq.-1)Then  ! W.Hu table      
        Open(50,file='WHUpk.dat')
        twop = 2.*3.14145926**2

        Ntab =0
 10     read(50,*,end=30,err=30)xx,pp
        Ntab =Ntab +1
        xkt(Ntab) =xx   *hsmall             ! scale to real Mpc
        Pkt(Ntab) = pp    /xx**3*twop
        goto 10
 30     write (*,*)  ' Read ',Ntab,' lines from file WHUpk.dat'
        Goto 50
      EndIf 
      If(iFlag.eq.-2)Then  ! Eis.Hu table      
           Open(50,file='pk_EHu_Om0.3_Omb0.043_s8=0.9.dat')
           write (*,*)'  Read data from file',
     &                ' pk_EHu_Om0.3_Omb0.043_s8=0.9.dat'
      Else
           Open(50,file='pk_EHu_WMAP06.dat')
           write (*,*)'  Read data from file',
     &                ' pk_EHu_WMAP06.dat'

      EndIf 
           Do i=1,9              ! read header
              read(50,'(a)',end=100,err=100) Line
           EndDo 
           Ntab =0
 12        read(50,*,end=32,err=32)xx,pp
           Ntab =Ntab +1
           xkt(Ntab) =xx   *hsmall             ! scale to real Mpc
           Pkt(Ntab) = pp                      !  Pk
           goto 12
 32        write (*,*)  
           write (*,*)  ' Read ',Ntab,' lines from file ...'
           write (*,*)  ' hsmall =',hsmall


 50     If(Ntab.le.1)stop 'wrong table for p(k)'
        StepK =log10(xkt(2)/xkt(1))
        alog0   = log10(xkt(1))
        
        Do k =2,Ntab-1      ! test that spacing is the same
           ss =log10(xkt(k+1)/xkt(k))
c           write (80,'(5g12.4)') xkt(k),xkt(k+1)-xkt(k),ss,
c     &            abs(ss/StepK-1.) 
           If(abs(ss/StepK-1.).gt.2.e-2)Then
              write (*,*)  ' error in K spacing. k=',k,xkt(k+1),xkt(k)
              STOP
           EndIf
        EndDo 
        close(50)
      Return
 100    write (*,*)  ' Error reading table of P(k)'
        stop
      End
