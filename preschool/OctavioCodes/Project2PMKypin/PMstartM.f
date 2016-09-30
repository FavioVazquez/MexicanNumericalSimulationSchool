C _______________________ START 3-D SIMULATIONS
C
C                         Klypin, August 1993
C
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMinitial.h'	  
C		   NROW = number of particles in 1 dimension
C		   NGRID= number of cells        in 1 dimension
C		   Npage =  number of particles per page = NROW**2
C


      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Par(6),ns,qqscaleb,QSCALE
      COMMON / KINENRG/ SKINE,SX,SY,SZ,SX2,SY2,SZ2
      Character  Hd*5,Tail*4
      DATA PI		/3.1415926535/

      Hd  ='PMcrs'
      Tail='.DAT'
        Open(30,file='TestDataSerial.dat')       
      CALL InitValues(NBYTE,SCALEL)
      write (*,*) ' InitValues is done. N bytes=',NBYTE
      Max_lev_abs = INT(ALOG(REAL(Lblock))/ALOG(2.)+1.1)
      write (*,*) ' Absolute Limits on mass resolution levels: ',
     &                  1,Max_lev_abs
c      NspecM   = 1  ! Nspecies -1
c      Nmin_lev = 1  ! Nspecies
      write (*,*) ' Actual Max resolution Level=',Nmax_lev,
     &               ' Min Level=',Nmin_lev
      write (16,*) ' Actual Max resolution Level=',Nmax_lev,
     &               ' Min Level=',Nmin_lev
      If(Nmax_lev.gt.Max_lev_abs)Then
         write (*,*) ' Wrong number of Levels for mass refinement (',
     &                 Nmax_lev ,') Block size =',Lblock,
     &                ' Use Absolute value=',Max_lev_abs
         Nmax_lev =Max_lev_abs
      Endif 

      OPEN(UNIT=9,FILE='PMcrd.DAT',form='unformatted')
      write (9)             ! this clears old header file
      CLOSE (9)
      CALL RDTAPE                  ! this opens files on disk

      OPEN(UNIT=16,FILE='RESNEW.DAT')
      write (16,*) ' Seed=',Nseed
       AEXP0 = AEXPN
       AEXPV = AEXPN - ASTEP/2.
C                                    set mask for multiple mass resolution
      CALL ZeroMASK
        If(Nmax_lev.ne.Nmin_lev)Then
            Open(30,file='SelCoords.DAT',status='unknown')
            CALL SetMASK
        EndIf
C			                   get a realization of spectrum in /GRID/
       NRAND =Nseed
       CALL SPECTR(ALPHA,NRAND)
C			                        get the displacement vector by FFT
       IFOUR = INT(ALOG(FLOAT(NROW))/ALOG(2.)+0.5)
       Write (*,*) ' Ifour =',IFOUR,' Nrow=',NROW
       CALL VECTOR(IFOUR)
	    write (*,*) ' SPECTR is done. IFOUR = ',IFOUR,' Alpha=',ALPHA
        Fact =  sqrt(Om0+Oml0*AEXPV**3+Ocurv*AEXPV)
       VCONS = -ALPHA/(2.*PI/NGRID)*(AEXPV/AEXP0)*SQRT(AEXPV)*Fact
       XCONS =	 ALPHA/(2.*PI/NGRID)*(AEXPN/AEXP0)
           write (*,*) ' Scaling Constants: (x)=', XCONS ,' (v)=',VCONS 

       QFACT =  FLOAT(NGRID)/FLOAT(NROW)
                      SKINE = 0.
                      SX    = 0.
                      SY    = 0.
                      SZ    = 0.
                      SX2   = 0.
                      SY2   = 0.
                      SZ2   = 0.

C			           find x,v and write to file, page by page
       Do i=1,Max_lev_abs
          lspecies(i) =0
       Enddo 
       Nspecies =0
         Icurrent    =0                         ! current particle
            Do kSpec =Max_lev_abs,1,-1
                CALL BLOCKS(XCONS,VCONS,kSpec,Icurrent)   
            Enddo 
           Vscale =ScaleL*100.*hubble/NGRID
      CALL WriteData(Vscale,AexpV,Wtotal)
      CALL Species_Change(Wtotal_o)      ! rearrange species                     
            EKIN = 0.5*SKINE/AEXPV**2       
           WRITE (*,'('' Ekin='',E12.4,'' Weght per cell='',g12.5,
     .                   '' (must be 1)'')')   EKIN,Wtotal/NGRID**3
           WRITE (16,'('' Ekin='',E12.4,'' Weght per cell='',g12.5,
     .                   '' (must be 1)'')')   EKIN,Wtotal/NGRID**3
C			        
      CALL WRTAPE   ! write header and control data
C                      write pt.dat file: time-step for particles
         Nparticles =lspecies(Nspecies)
      Do ic1 =1,Nparticles
         XPt(ic1) =astep
      Enddo 
      open ( 60 , file = 'pt.dat' , form = 'unformatted' )
      write(60) (XPt(ic1),ic1=1,Nparticles)
      write (*,*) ' pt=', (XPt(ic1),ic1=1,10)
      STOP
      END
C--------------------------------------------------
C                                                     sqrt(Power spectrum)
C                                                           k = (2pi/L) wk
C--------------------------------------------------
      FUNCTION TRUNF(WK)      
      INCLUDE 'PMparameters.h'
      PARAMETER (NSPEC =NROW/2)
      PARAMETER (NtabM = 100000)
      COMMON /PTABLE/xkt(0:NtabM),Pkt(0:NtabM),Ntab,StepK,alog0,iFlag  
      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Par(6),ns,qqscaleb,QSCALE
      Real                 ns,k
      IF (WK.GE.FLOAT(NSPEC)) THEN
	       TRUNF =0.
	       RETURN
      ENDIF
        k = QSCALE*wk

      If(iFlag.gt.0)Then       ! use approximations
      If(Par(6).ne.0.)Then   ! Holtzman approx
        sk= sqrt(k)
        TRUNF = k**(ns/2.) /
     .         (1.+sk*(Par(2)
     .            +sk*(Par(3)
     .            +sk*(Par(4)
     .            +sk*(Par(5) )))) )**(Par(6))
      Else                   ! BBKS + Sugiyama approx
c        Gammaeff =Om*hsmall/exp(Omb*(1.+sqrt(hsmall/0.5)/Om))
c        Q = wk/hsmall/Gammaeff
        Q = k*qqscaleb
        TRUNF = k**(ns/2.)* LOG(1.+Par(1)*Q)/(Par(1)*Q)/
     .          sqrt(sqrt(1.+Q*(Par(2)+
     .                  Q*(Par(3)**2+
     .                  Q*(Par(4)**3+
     .                  Q*(Par(5)**4) )))))
      EndIf
      Else               ! use table for P(k)
         TRUNF = sqrt(Ppk(k))
      EndIf  
      RETURN
      END
C---------------------------------------
      FUNCTION Ppk(x)
c                     interpolate table with p(k)
c                       x is in real 1/Mpc
C---------------------------------------
      PARAMETER (NtabM = 100000)
      COMMON /PTABLE/xkt(0:NtabM),Pkt(0:NtabM),Ntab,StepK,alog0,iFlag  

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
      SUBROUTINE  ReadTable(hsmall)
c                     interpolate table with p(k)
C---------------------------------------
      PARAMETER (NtabM = 100000)
      COMMON /PTABLE/xkt(0:NtabM),Pkt(0:NtabM),Ntab,StepK,alog0,iFlag
      character*80 Line
  
      If(iFlag.eq.-1)Then  ! W.Hu table      
           Open(50,file='WHUpk.dat')
           twop = 2.*3.14145926**2

           Ntab =0
 10        read(50,*,end=30,err=30)xx,pp
           Ntab =Ntab +1
           xkt(Ntab) =xx   *hsmall             ! scale to real Mpc
           Pkt(Ntab) = pp    /xx**3*twop   ! k^3Pk -> Pk
           goto 10
 30        write (*,*)  
           write (*,*)  ' Read ',Ntab,' lines from file WHUpk.dat'
           write (*,*)  ' hsmall =',hsmall
           Goto 50
      EndIf 
      If(iFlag.eq.-2)Then  ! Eis.Hu table      
           Open(50,file='pk_EHu_Om0.3_Omb0.043_s8=0.9.dat')
      Else
           Open(50,file='pk_EHu_WMAP06.dat')
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
           write (*,*)  ' Read ',Ntab,' lines from file pk_EHu_Om0.3...'
           write (*,*)  ' hsmall =',hsmall


 50   If(Ntab.le.1)stop 'wrong table for p(k)'
        StepK =log10(xkt(2)/xkt(1))
        alog0   = log10(xkt(1))
        
        Do k =2,Ntab-1      ! test that spacing is the same
           ss =log10(xkt(k+1)/xkt(k))
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


C*********************************************************************
C	   Power spectrum for CDM. Wk-  wavenumber(Mpc^-1, h=0.5)
C                         P = wk *T^2(k)
C                         QSCALE = 2pi/Box_size,  QS =hubble**(-2)
c      FUNCTION TRUNF(WK)
c      COMMON / TRUNCOM/ QSCALE,g,gy
c      COMMON / TRUNCDM/ QS,A1,A2,A3,A4,A5,A32,A43,A54
c       Q = QSCALE*QS*WK
c       TRUNF = sqrt( WK*
c     +		 (LOG(1.+A1*Q)/(A1*Q))**2/
c     +		 SQRT(1.+Q*(A2 +Q*(A32   +Q*(A43+Q*A54)))) )
c
c      RETURN
c      END

C*********************************************************************
C			  INITIALIZE CONSTANTS:
C			      Scalel    =length of square in MPC
C			      AEXPN =expansion factor (current)
C			      ASTEP	  =increment of expansion factor: AEXPN =AEXP_0+N*ASTEP
C                     PARTW	=weight of a particle: total 'mass' per cell must be unity
c                                          even for LCDM or OCDM 
C			      AMPLT	=amplitude of perturbations inside the simulation box
C			      NBYTE	=record length in bytes
C
      SUBROUTINE InitValues(NBYTE,SCALEL)
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMinitial.h'
      PARAMETER (NtabM = 100000)
      COMMON /PTABLE/xkt(0:NtabM),Pkt(0:NtabM),Ntab,StepK,alog0,iFlag  
      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Par(6),ns,qqscaleb,QSCALE
      COMMON / FERMI  / Vscale
      Real                 ns      
      Character            Answer*1
      Real                 INPUT
      DATA PI		         /3.1415926535/

      Do i=1,100
         extras(i) =0.
      enddo 
       Write (*,*) 'Would you like to use a model provided in cdm.fit '
       Write (*,*) 'or your own model?' 
       Write (*,'(A,$)') ' Enter Yes for cdm.fit;  No for your model ='
       Read  (*,'(A)')  Answer
       HEADER='N=128x256 L=20h-1CDMsig=0.239 z=30-----------'
       write (*,*) '------ Enter Header for the run up to 45 characters'
       write (*,*) '       Example:'
       write (*,'(A)') HEADER                                   
       read  (*,'(A)') HEADER
       write (*,'(A)') HEADER
       AEXPN =INPUT(' Initial expansion parameter (0-1)=')
       ASTEP =INPUT(' Step of the expansion parameter  =')
       AMPLT =INPUT(' Amplitude of density fluctuations=')
       ns    =INPUT(' Slope of the Power spectrum     n=')
       SCALEL=INPUT(' Box size (Mpc/h)                 =')
       Ocurv =INPUT(' Omega_curvature   at present     =')

       Nseed =INPUT(' Random seed (Integer  1-2^31)    =') 
       If(Nseed.le.0)Then
          lux = 2               ! level for luxury
          extras(80) = 1.
          extras(81) = lux
          If(Nseed.eq.0)Nseed = 121071
          Nseed = ABS(Nseed)
          CALL rluxgo(lux,Nseed,0,0) ! initialize luxury
       Else
          extras(80) = 0.
       EndIf
c                       extras(80) is used as switch for
c                                        random numbers generators
c                         = 0 - old Randd+Gauss
c                         = 1 - luxury+Gauss3
c                        extras(81) = lux - level of luxury
          
       extras(100) = SCALEL                   ! store the box size
       If(Answer.eq.'Y' .or. Answer.eq.'y')Then
          CALL MODEL(hsmall)
          hubble=hsmall
          Om0   =Om
          Oml0  =1. -Om0 -Ocurv
       Else
       write (*,*) ' You use your cosmological model.'
       write (*,*) ' Be sure you provide routine TRUNF'
       hubble=INPUT(' The Hubble constant (h=H/100)    =')
       Om0   =INPUT(' Omega_0 matter    at present     =')
       Oml0  =INPUT(' Omega_lambda      at present     =')
       EndIf
       SCALEL=SCALEL/hubble   ! scale it to real megaparsecs      
       NBYTE = NPAGE*6*4

       Write (*,'(A,$)') ' Enter Min and Max number of mass species ='
       read (*,*) Nmin_lev,Nmax_lev
             W     = (FLOAT(NGRID)/FLOAT(NROW))**3
            PARTW = W 
       If(Nmin_lev.lt.Nmax_lev)Then    ! multiple mass resolution
              write(*,*) 
              write(*,*) ' You use multiple mass resolution'
       endif 
       If(Nmin_lev.gt.Nmax_lev)Then ! change the order of levels
          i = Nmax_lev
          Nmax_lev = Nmin_lev
          Nmin_lev = i
       endif

       ISTEP = 0
       TINTG = 0.
       AU0   = 0.
       AEU0  = 0.
       EKIN  = 0.
       EKIN1 = 0.
       EKIN2 = 0.
       NROWC = NROW
       NGRIDC= NGRID
	     QSCALE = 2.*PI/SCALEL
	     QS     = hubble**(-2)
        write (*,'(/3x,3(a,g12.3))')  ' Qscale=',QSCALE,
     &                                   ' Box=',SCALEL,' slope=',ns   

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
      SUBROUTINE MODEL(hsmall)
C---------------------------------------
      INCLUDE 'PMparameters.h'
      COMMON / TRUNCOM/ Om,Omb,Omc,Omnu,Par(6),ns,qqscaleb,QSCALE
      PARAMETER (NtabM = 100000)
      COMMON /PTABLE/xkt(0:NtabM),Pkt(0:NtabM),Ntab,StepK,alog0,iFlag  
      Real                 ns,INPUT
      Character             Header1*79

      Line1 =INPUT(' Enter Line Number in cdm.fit     =')
      If(Line1.gt.0)Then
         iFlag = Line1
      OPEN(2,file='cdm.fit')
      Read(2,'(A)') Header1
      If(Line1.gt.2)Then
         Do i=1,Line1-2
            Read (2,*) a
         EndDo
      EndIf
      Read (2,*) Om,Omb,Omc,Omnu,hsmall,Par
      CLOSE (2) 
      
      theta = 2.726/2.7  ! = T_cmb/2.7
      Ob0   =Omb/Om      ! Hu&Sugiyama fitting
      Omh2  =Om*hsmall**2
      a1    =(46.9*Omh2)**0.670*(1.+(32.1*Omh2)**(-0.532))
      a2    =(12.0*Omh2)**0.424*(1.+(45.*Omh2)**(-0.582))
      alpha =1./a1**Ob0/a2**(Ob0**3)
      qqscaleb = theta**2/(1.-Ob0)**0.60/sqrt(alpha)/Omh2

      Write (*,20)Om,Omb,Omc,Omnu,hsmall,Par
c      Write (16,20)Om,Omb,Omc,Omnu,hsmall,Par
 20   Format(' Model: Om_0=',F5.3,' Om_baryon=',F5.3,
     .       ' Om_cold=',F5.3,' O_nu=',F5.3,
     .       ' hsmall=',F4.2,/8x,'Parameters=',(6G9.3))
      Else                     ! read spectrum from table
         iFlag = Line1
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
C------------------------------------------------
C                              Clear MASK for multiple mass resolution                             
C                               Initiat to level Nmin_lev
      SUBROUTINE ZeroMask
C------------------------------------------------
      INCLUDE 'PMparameters.h'	    
      INCLUDE 'PMinitial.h'	  

         DO  MK3 = 1,Nblocks
         DO  MK2 = 1,Nblocks
         DO  MK1 = 1,Nblocks
	        iMask(MK1,MK2,MK3) =Nmin_lev
         ENDDO
         ENDDO
         ENDDO
      RETURN
      END
C------------------------------------------------
C                       Set MASK for multiple mass resolution  
C                               Read particles with lagrangian coordinates
C                                                to set the mask                           
C                               Block   = 1 = low resolution
C                                       = 2 = higher resolution
C                                       = 3 = even higher resolution
C               After the mask is set, make maximum resolution Nmax_lev
C
      SUBROUTINE SetMask
C------------------------------------------------
      INCLUDE 'PMparameters.h'	  
      INCLUDE 'PMinitial.h'	  
      Dimension nLevel(10)
      Character*70 Line

      Do i=1,5
         Read(30,'(A)') Line
         Write(*,'(A)') Line
      EnDDo
      Nlines =0
 10   Read(30,*,end=50,err=50) i,xw,yw,zw,
     +            vwx,vwy,vwz,im,jm,km
           Nlines =Nlines +1
c           write (*,*) ' Block=',im,jm,km
           If(im.lt.0.or.im.gt.Nblocks)Then
              write (*,*) ' wrong x_block: ',im,jm,km,' Line=',Nlines
              Stop
           EndIf
            If(jm.lt.0.or.jm.gt.Nblocks)Then
              write (*,*) ' wrong y_block: ',im,jm,km,' Line=',Nlines
              Stop
           EndIf
            If(km.lt.0.or.km.gt.Nblocks)Then
              write (*,*) ' wrong z_block: ',im,jm,km,' Line=',Nlines
              Stop
           EndIf
c             If(iMask(im,jm,km).ne.0)Then
c              write (*,*) ' Occupied block: ',im,jm,km,' Line=',Nlines
c              write (*,*) '                 Block=',iMask(im,jm,km)
c              Stop
c           EndIf
           
           iMask(im,jm,km) =Nmax_lev    ! mark high resolution block
               Do ik =-1,1
               k =km +ik
               If(k.le.0)k=k+Nblocks
               If(k.gt.Nblocks)k=k-Nblocks
               Do ij =-1,1
                  j =jm +ij
                  If(j.le.0)j=j+Nblocks
                  If(j.gt.Nblocks)j=j-Nblocks
                  Do ii =-1,1
                      i =im +ii
                      If(i.le.0)i=i+Nblocks
                      If(i.gt.Nblocks)i=i-Nblocks
                      ind =abs(ik)+abs(ij)+abs(ii) ! only seven blocks in crest
                                ! are marked out
c                      If(iMask(i,j,k).le.0.or.iMask(i,j,k).gt.Nspecies)
c     +                     Then
c                         Stop
c                      Endif 
                      If(ind.le.1)Then
                          iMask(i,j,k)=Nmax_lev ! change mask
                      Endif 
                  ENDDO
                ENDDO
                ENDDO
            GoTo 10
C                      Change the mask: iteratevely add layers of lower resolution
 50        write(*,*) ' Total Lines read=',Nlines
       DO kSpec =Nmax_lev-1,2,-1     !
         DO  MK3 = 1,Nblocks     ! find blocks adjasent to filled blocks
         DO  MK2 = 1,Nblocks     ! mark them with kSpec =Nspecies-1, -2, ..1
         DO  MK1 = 1,Nblocks
            iM =      iMask(MK1,MK2,MK3)
            If(iM.eq.kSpec+1)Then          ! change to kSpec if not high res
               Do ik =-1,1
               k =MK3 +ik
               If(k.le.0)k=k+Nblocks
               If(k.gt.Nblocks)k=k-Nblocks
               Do ij =-1,1
                  j =MK2 +ij
                  If(j.le.0)j=j+Nblocks
                  If(j.gt.Nblocks)j=j-Nblocks
                  Do ii =-1,1
                     i =MK1 +ii
                     If(i.le.0)i=i+Nblocks
                     If(i.gt.Nblocks)i=i-Nblocks
                     If(iMask(i,j,k).lt.kSpec)iMask(i,j,k)=kSpec ! change mask
                  ENDDO
               ENDDO
               ENDDO
            EndIf
         ENDDO
         ENDDO
         ENDDO
       Enddo 

       Do i=1,Max_lev_abs            ! count marked blocks
          nLevel(i)=0
       Enddo 
         DO  MK3 = 1,Nblocks
         DO  MK2 = 1,Nblocks
         DO  MK1 = 1,Nblocks
            iM =      iMask(MK1,MK2,MK3)
            If(iM.le.Nmax_lev.and.iM.ge.Nmin_lev)Then
               nLevel(iM) = nLevel(iM) +1 
            Else
               write (*,*) 'wrong mask: ',iM, ' block=',MK1,MK2,MK3
               write (*,*)'      it should be within limits:', 
     &                                Nmin_lev,Nmax_lev
               Stop
            EndIf
         ENDDO
         ENDDO
         ENDDO

         Ntot =0
       Do i=Max_lev_abs,1,-1
          icurrent =nLevel(i)*8**(i-1)
         Write (*,*) ' Level=',i,' Number of Blocks=',nLevel(i),
     &                   ' N_Particles=', icurrent
         Ntot = Ntot + icurrent
       Enddo   
         Write (*,*) ' Total expected particles =',Ntot
      RETURN
      END
C------------------------------------------------
C                             Make a realization of spectrum of perturbations:
C                             ALPHA  = normalization factor for displacements
C                             NRAND = seed for random numbers 
      SUBROUTINE SPECTR(ALPHA,NRAND)
C------------------------------------------------
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMinitial.h'
      INCLUDE 'luxury.h'
      REAL*8     SUMM, WI3,WJ3,WK3,WD,TS,TRX
C						    Set spectrum
      SUMM = 0.
      iRandom =INT(extras(80)+1.e-5) ! set type of Random Numbers
      iFlag =0                                                       ! 0 - old GAUSS 1-GAUSS3(ranlux)
C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (MK3,MK2,MK1 )
      DO  MK3 = 1,NROW
      DO  MK2 = 1,NROW
      DO  MK1 = 1,NROW
	     GRX(MK1,MK2,MK3) =0.
	     GRY(MK1,MK2,MK3) =0.
	     GRZ(MK1,MK2,MK3) =0.
      ENDDO
      ENDDO
      ENDDO
c                    Do NOT parallelize: random numbers
      DO  k = 1,NROW   
         ksign = -1
         kz    = k +NSPEC
         k3    = k - 1
         If(k3.GT.NSPEC)k3=k3-NSPEC
           IF(kz.gt.NROW)THEN
              kz    =kz -NROW
              ksign =1
           ENDIF
           WK3= K3**2
c	        write (*,*) '  K=',k,' K3=',K3,' Nspec=',NSPEC
           DO  j = 1,NROW
              jsign = -1
              jz    = j +NSPEC
              j3    = j - 1
              If(j3.GT.NSPEC)j3=j3-NSPEC
                IF(jz.gt.NROW)THEN
                   jz    =jz -NROW
                   jsign =1
                ENDIF
                WJ3= j3**2
                DO  i = 1,NROW
                   isign = -1
                   iz    = i +NSPEC
                   i3    = i - 1
                   If(i3.GT.NSPEC)i3=i3-NSPEC
                     IF(iz.gt.NROW)THEN
                        iz    =iz -NROW
                        isign =1
                     ENDIF
                     WI3= i3**2

	       IF(i3+j3+k3.EQ.0)THEN
		         GRX(1,1,1) = 0.
		         GRY(1,1,1) = 0.
		         GRZ(1,1,1) = 0.
	       ELSE
		  WD = WI3 + WJ3 + WK3
		  WK = SQRT(WD)
                      If(iRandom.eq.0)Then
                           TS = TRUNF(WK) * GAUSS(NRAND)
                      Else
                           TS = TRUNF(WK) * GAUSS3()
                      EndIf
		  TRX = TS / WD
		  GRX(iz, j, k) = i3 * TRX *isign
		  GRY( i,jz, k) = j3 * TRX *jsign
		  GRZ( i, j,kz) = k3 * TRX *ksign
		  SUMM = SUMM + TS**2
	       ENDIF
	    ENDDO
	   ENDDO
      ENDDO

      IF(SUMM.LE.0.)WRITE (*,*) ' Error!!! Summ over spectrum = 0'
       ALPHA = AMPLT / SQRT(SUMM) *sqrt(8.)
      Write (*,'(10x,a20,3g12.4)') 'SPECTR ==>  SUMM=', SUMM, ALPHA 

      RETURN
      END
C------------------------------------------------
C		                        Define coordinates and velosities for
C                                  particles with resolution = Indx
C                                  Indx =1 = high resolution
C                                  Indx =2 = medium resolution
C                                  Icurrent = number of particles
      SUBROUTINE BLOCKS(XCONS,VCONS,Indx,Icurrent)
C------------------------------------------------
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMinitial.h'	  
      Dimension nLevel(10)
      Real*8   sDispl
       Do i=1,  Max_lev_abs          ! number of particles at
          nLevel(i)=0                  ! different levels
       Enddo 
       Boff = Lblock/2 +0.5
       sDispl = 0.

      QFACT = FLOAT(NGRID)/FLOAT(NROW)
        DO  MK3 = 1,Nblocks
            M3   = (MK3-1)*Lblock   ! part of offset
        DO  MK2 = 1,Nblocks
            M2   = (MK2-1)*Lblock   ! part of offset
        DO  MK1 = 1,Nblocks
            M1   = (MK1-1)*Lblock   ! part of offset
	      iM =iMask(MK1,MK2,MK3)
C----------------------------------------
c                     M1 + Boff  = lagrangian coordinate of 'mean' particle
c                     Lblock**3 = total number of averaged cells
          If(iM.eq.Indx)Then 
          If(iM.eq.1)Then  ! Low resolution: one particle
c              DX =0.
c              DY =0.
c              DZ =0.
c              Do kp=M3+1,M3+Lblock   !  sum of displacements
c              Do jp=M2+1,M2+Lblock
c              Do ip=M1+1,M1+Lblock
c                 DX = DX +GRX(ip,jp,kp)
c                 DY = DY +GRY(ip,jp,kp)
c                 DZ = DZ +GRZ(ip,jp,kp)
c              EndDo
c              EndDo
c              EndDo
c              DX = DX/Lblock**3               ! mean displacement
c              DY = DY/Lblock**3
c              DZ = DZ/Lblock**3
                 
                 kp = M3+(Lblock+1)/2
                 jp = M2+(Lblock+1)/2
                 ip = M1+(Lblock+1)/2
                 DX = GRX(ip,jp,kp)
                 DY = GRY(ip,jp,kp)
                 DZ = GRZ(ip,jp,kp)

                Q3    = QFACT*(M3+Boff-1.) +1. ! scale them to box size
	        Q2    = QFACT*(M2+Boff-1.) +1.
	        Q1    = QFACT*(M1+Boff-1.) +1.
               Icurrent = Icurrent +1
               If(Icurrent.gt.Nmaxpart)Then
                    write (*,*)'Too many particles:',Icurrent,' block=',
     &                                                     MK1,MK2,MK3
                    STOP
               EndIf
               nLevel(Indx) =nLevel(Indx) +1
	         XPt(Icurrent) = Q1 - XCONS*DX  + 1.E-4  +0.5 
	         YPt(Icurrent) = Q2 - XCONS*DY  + 1.E-4  +0.5
	         ZPt(Icurrent) = Q3 - XCONS*DZ  + 1.E-4  +0.5
	         VXt(Icurrent) =          VCONS*DX
	         VYt(Icurrent) =          VCONS*DY
	         VZt(Icurrent) =          VCONS*DZ

c            write (*,400) Q1,Q2,Q3,DX,DY,DZ,Icurrent
c 400        format(6g11.4,i8)
           Else
C------------------------------------------
           If(iM.eq.Max_lev_abs)Then  ! Highest resolution: 
              Do kp=M3+1,M3+Lblock 
              Do jp=M2+1,M2+Lblock 
              Do ip=M1+1,M1+Lblock 
                 DX = GRX(ip,jp,kp)     !  find displacements
                 DY = GRY(ip,jp,kp)
                 DZ = GRZ(ip,jp,kp)
                 Q3    = QFACT*(kp-1.) +1. ! scale them to box size
	           Q2    = QFACT*(jp-1.) +1.
	           Q1    = QFACT*(ip-1.) +1.
                 Icurrent = Icurrent +1
                 If(Icurrent.gt.Nmaxpart)Then
                    write (*,*)'Too many particles:',Icurrent,' block=',
     &                                                     MK1,MK2,MK3
                    STOP
                 EndIf
               nLevel(Indx) =nLevel(Indx) +1
	         XPt(Icurrent) = Q1 - XCONS*DX  + 1.E-4  +0.5 
	         YPt(Icurrent) = Q2 - XCONS*DY  + 1.E-4  +0.5
	         ZPt(Icurrent) = Q3 - XCONS*DZ  + 1.E-4  +0.5
	         VXt(Icurrent) =          VCONS*DX
	         VYt(Icurrent) =          VCONS*DY
	         VZt(Icurrent) =          VCONS*DZ
                 sDispl        = sDispl +DX**2+DY**2+Dz**2
               EndDo
              EndDo
              EndDo

           Else
C-------------------------------------------
c            Medium resolution: 
c                   KBlock = 2**(iM-1)        = number of particles in 1D per block
c                   Nbcells= Lblock/KBlock = number of cells per particle
c                   (MK-1)*Lblock +(i-1)*Nbcells = offset of cells in 1D
c                   (MK-1)*Lblock +(i-1)*Nbcells +Nbcells/2 +1/2 = lagrangian
c                                                                 coordinate of particle
c                    M1   = (MK1-1)*Lbloc
              KBlock = 2**(iM-1)           ! number of particles in 1D per block
              Nbcells = Lblock/KBlock ! number of cells per particle in 1D
              Xlag     = Nbcells/2 +0.5 ! offset for lagrangian coordinates
              Nbcells3 =Nbcells**3      ! total number of cells per particle
c       write (*,*)' Lblock=',Lblock,' Nblocks=',Nblocks,
c     +                ' Npart/block/1D=',KBlock

              Do i =1,KBlock                          ! loop over particles
                 i1      =M1+(i-1)*Nbcells      ! offset for cells in X
              	  Q1  = QFACT*(i1+Xlag-1.) +1. ! scale it to box size
              Do j =1,KBlock
                 j1      =M2+(j-1)*Nbcells       ! offset for cells in X
 	              Q2  = QFACT*(j1+Xlag-1.) +1.
              Do k=1,KBlock
                 k1      =M3+(k-1)*Nbcells       ! offset for cells in X 
                 Q3    = QFACT*(k1+Xlag-1.) +1. 
c                 DX =0.
c                 DY =0.
c                 DZ =0.
c                 Nncells =0
c                 Do  kp=k1+1,k1+Nbcells !  sum of displacements
c                 Do  jp=j1+1,j1+Nbcells
c                 Do  ip=i1+1,i1+Nbcells
c                    Nncells =Nncells +1
c                    If(ip.gt.NROW.or.ip.lt.1)Then
c                       write (*,*) ' error'
c                    endif 
c                    If(jp.gt.NROW.or.jp.lt.1)Then
c                       write (*,*) ' error'
c                    endif 
c                    If(kp.gt.NROW.or.kp.lt.1)Then
c                       write (*,*) ' error'
c                    endif 
c                    DX = DX +GRX(ip,jp,kp)
c                    DY = DY +GRY(ip,jp,kp)
c                    DZ = DZ +GRZ(ip,jp,kp)
c                 EndDo
c                 EndDo
c                 EndDo
c                 DX = DX/ Nbcells3              ! mean displacement
c                 DY = DY/ Nbcells3
c                 DZ = DZ/ Nbcells3
c                    If(Nncells.ne.Nbcells3)Then
c                       write (*,*) ' error: wrong number of cells=',
c     +                      Nncells,Nbcells3
c                       stop
c                    endif 
c                    If(Nncells.ne.Nbcells3)Then
c                       write (*,*) ' error: wrong number of cells=',
c     +                      Nncells,Nbcells3
c                       stop
c                    endif 

                 
                 kp = k1+(Nbcells+1)/2
                 jp = j1+(Nbcells+1)/2
                 ip = i1+(Nbcells+1)/2
                 DX = GRX(ip,jp,kp)
                 DY = GRY(ip,jp,kp)
                 DZ = GRZ(ip,jp,kp)

                 Icurrent = Icurrent +1
                 If(Icurrent.gt.Nmaxpart)Then
                    write (*,*)'Too many particles:',Icurrent,' block=',
     &                                                     MK1,MK2,MK3
                    STOP
                 EndIf
                nLevel(Indx) =nLevel(Indx) +1
	          XPt(Icurrent) = Q1 - XCONS*DX  + 1.E-4  +0.5 
	          YPt(Icurrent) = Q2 - XCONS*DY  + 1.E-4  +0.5
	          ZPt(Icurrent) = Q3 - XCONS*DZ  + 1.E-4  +0.5
	          VXt(Icurrent) =          VCONS*DX
	          VYt(Icurrent) =          VCONS*DY
	          VZt(Icurrent) =          VCONS*DZ
             EndDo
             EndDo
             EndDo
           EndIf
           EndIf

      EndIf         ! end iM = Indx 
      EndDo
      EndDo
      EndDo

      If(nLevel(indx).gt.0)Then
               Nspecies = Nspecies +1
                lspecies(Nspecies)   =Icurrent
                wspecies(Nspecies) =PARTW*
     &                                        Lblock**3/8**(Indx-1)
           write (*, '(10x,a20,g12.4)')' RMS 3d diplacement=',
     &             XCONS*sqrt(sDispl/nLevel(indx))
           write (*,'(10x,a20,i4,a,i10,a,6i10)') ' Species=',Nspecies,
     &                    ' Particles=',Icurrent,' on  levels:',
     &                     (nLevel(i),i=1,Max_lev_abs   )
      endif 

      RETURN
      END
C------------------------------------------------
C                                              				   FFT of the spectrum
      SUBROUTINE VECTOR(IFOUR)
C------------------------------------------------
      INCLUDE 'PMparameters.h'	  
      INCLUDE 'PMinitial.h'	  
      dimension Uf(marr) , Vf(marr)
      dimension Qi(mf67) , jndx(mf67)

      jb1 =4
      jq   =IFOUR
      CALL SETF67(jb1,jq,jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

C				       x-component
C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( k,j,i,Uf,Vf,jp,jsl,ll1,k2,k4 )
      DO K=1,NROW
	   DO J=1,NROW
	    DO I=1,NROW
	       Uf(I) =GRX(I,J,K)
	    ENDDO
	    CALL FOUR67(jb1,jq,jp,jsl,ll1,k2,k7,Qi,jndx,Uf,Vf)
	    DO I=1,NROW
	       GRX(I,J,K) = Vf(I)
	    ENDDO
	   ENDDO
      ENDDO
      write (*,*) ' 1st is done'
C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( k,j,i,Uf,Vf,jp,jsl,ll1,k2,k4 )
      DO K=1,NROW
	   DO I=1,NROW
	    DO J=1,NROW
	       Uf(J) =GRX(I,J,K)
	    ENDDO
	    CALL FOUR67(jb1,jq,jp,jsl,ll1,k2,k7,Qi,jndx,Uf,Vf)
	    DO J=1,NROW
	       GRX(I,J,K) = Vf(J)
	    ENDDO
	   ENDDO
      ENDDO
C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( k,j,i,Uf,Vf,jp,jsl,ll1,k2,k4 )
      DO J=1,NROW
	   DO I=1,NROW
	    DO K=1,NROW
	       Uf(K) =GRX(I,J,K)
	    ENDDO
	    CALL FOUR67(jb1,jq,jp,jsl,ll1,k2,k7,Qi,jndx,Uf,Vf)
	    DO K=1,NROW
	       GRX(I,J,K) = Vf(K)/8.
	    ENDDO
	   ENDDO
      ENDDO
				       write (*,*) ' x is done'
C				       Y-component
C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( k,j,i,Uf,Vf,jp,jsl,ll1,k2,k4 )
      DO K=1,NROW
	   DO J=1,NROW
	    DO I=1,NROW
	       Uf(I) =GRY(I,J,K)
	    ENDDO
	    CALL FOUR67(jb1,jq,jp,jsl,ll1,k2,k7,Qi,jndx,Uf,Vf)
	    DO I=1,NROW
	       GRY(I,J,K) = Vf(I)
	    ENDDO
	   ENDDO
      ENDDO
C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( k,j,i,Uf,Vf,jp,jsl,ll1,k2,k4 )
      DO K=1,NROW
	   DO I=1,NROW
	    DO J=1,NROW
	       Uf(J) =GRY(I,J,K)
	    ENDDO
	    CALL FOUR67(jb1,jq,jp,jsl,ll1,k2,k7,Qi,jndx,Uf,Vf)
	    DO J=1,NROW
	       GRY(I,J,K) = Vf(J)
	    ENDDO
	   ENDDO
      ENDDO
C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (  k,j,i,Uf,Vf,jp,jsl,ll1,k2,k4 )
      DO J=1,NROW
	   DO I=1,NROW
	    DO K=1,NROW
	       Uf(K) =GRY(I,J,K)
	    ENDDO
	    CALL FOUR67(jb1,jq,jp,jsl,ll1,k2,k7,Qi,jndx,Uf,Vf)
	    DO K=1,NROW
	       GRY(I,J,K) = Vf(K)/8.
	    ENDDO
	   ENDDO
      ENDDO
				       write (*,*) ' y is done'
C				       Z-component
C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( k,j,i,Uf,Vf,jp,jsl,ll1,k2,k4 )
      DO K=1,NROW
	   DO J=1,NROW
	    DO I=1,NROW
	       Uf(I) =GRZ(I,J,K)
	    ENDDO
	    CALL FOUR67(jb1,jq,jp,jsl,ll1,k2,k7,Qi,jndx,Uf,Vf)
	    DO I=1,NROW
	       GRZ(I,J,K) = Vf(I)
	    ENDDO
	   ENDDO
      ENDDO
C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE (  k,j,i,Uf,Vf,jp,jsl,ll1,k2,k4 )
      DO K=1,NROW
	   DO I=1,NROW
	    DO J=1,NROW
	       Uf(J) =GRZ(I,J,K)
	    ENDDO
	    CALL FOUR67(jb1,jq,jp,jsl,ll1,k2,k7,Qi,jndx,Uf,Vf)
	    DO J=1,NROW
	       GRZ(I,J,K) = Vf(J)
	    ENDDO
	   ENDDO
      ENDDO
C.... Open_MP
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+PRIVATE ( k,j,i,Uf,Vf,jp,jsl,ll1,k2,k4 )
      DO J=1,NROW
	   DO I=1,NROW
	    DO K=1,NROW
	       Uf(K) =GRZ(I,J,K)
	    ENDDO
	    CALL FOUR67(jb1,jq,jp,jsl,ll1,k2,k7,Qi,jndx,Uf,Vf)
	    DO K=1,NROW
	       GRZ(I,J,K) = Vf(K)/8.
	    ENDDO
	   ENDDO
      ENDDO
	write (*,*) ' z is done'
      RETURN
      END
C---------------------------------------------
C                       Write current data to  disk/tape
C                        for PMstartM - multiple masses
      SUBROUTINE WriteData(Vscale,AexpV,Wtotal)
C----------------------------------------------
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMinitial.h'	  

c      COMMON / BUFF /	XPt(Nmaxpart),YPt(Nmaxpart),ZPt(Nmaxpart),
c     +			VXt(Nmaxpart),VYt(Nmaxpart),VZt(Nmaxpart)
       COMMON / KINENRG/ SKINE,SX,SY,SZ,SX2,SY2,SZ2
C
      XMAX  = FLOAT(NGRID) + 1.
      XSHF  = FLOAT(NGRID)
      SKINE =0.
      Wtotal =0.
       jstart =1

        write(30,*) ' XXX'
        write(30,'(10g12.5)') (XPt(i),i=1,lspecies(Nspecies),5000)
        write(30,*) ' VVV'
        write(30,'(10g12.5)') (VXt(i),i=1,lspecies(Nspecies),5000)
        write(30,*) ' XXX'
        write(30,'(10g12.5)') (YPt(i),i=1,lspecies(Nspecies),5000)
        write(30,*) ' VVV'
        write(30,'(10g12.5)') (VYt(i),i=1,lspecies(Nspecies),5000)
        write(30,*) ' XXX'
        write(30,'(10g12.5)') (ZPt(i),i=1,lspecies(Nspecies),5000)
        write(30,*) ' VVV'
        write(30,'(10g12.5)') (VZt(i),i=1,lspecies(Nspecies),5000)

       Do j =1,Nspecies
          vvx =0.
          vvy =0.
          vvz =0.
          vx2 =0.
          vy2 =0.
          vz2 =0.
          jend =lspecies(j)
          W     =wspecies(j)
          If(jend.gt.jstart)Then
          Do i=jstart,jend
             vvx = vvx + VXt(i)
             vvy = vvy + VYt(i)
             vvz = vvz + VZt(i)
             vx2 = vx2 + VXt(i)**2
             vy2 = vy2 + VYt(i)**2
             vz2 = vz2 + VZt(i)**2
              SKINE    = SKINE + W*(VXt(i)**2 +VYt(i)**2 +VZt(i)**2)
              Wtotal =Wtotal +w 
          EndDo
           npp =max(jend-jstart+1,1)
           v2 = sqrt((vx2 +vy2+vz2)/npp)/AexpV*Vscale
           vvx =vvx/npp/AexpV*Vscale
           vvy =vvy/npp/AexpV*Vscale
           vvz =vvz/npp/AexpV*Vscale
           write (*,350) vvx,vvy,vvz,v2,jend-jstart+1,j
 350       format('       V(km/s)=',4g11.3, ' Npart=',i6,' Species=',i3)
          jstart =  jend+1
          Endif 
       EndDo
      Ibuff =0
      KROW =0
      Do i=1,lspecies(Nspecies)
C                                Periodical boundary conditions
	       IF(XPt(i).GT.XMAX)	      XPt(i)=XPt(i)-XSHF
	       IF(XPt(i).LE.1.)	            XPt(i)=XPt(i)+XSHF
	       IF(YPt(i).GT.XMAX)	      YPt(i)=YPt(i)-XSHF
	       IF(YPt(i).LE.1.)	            YPt(i)=YPt(i)+XSHF
	       IF(ZPt(i).GT.XMAX)	      ZPt(i)=ZPt(i)-XSHF
	       IF(ZPt(i).LE.1.)	            ZPt(i)=ZPt(i)+XSHF
             Ibuff              = Ibuff +1
             XPAR(Ibuff) = XPt(i)
             YPAR(Ibuff) = YPt(i)
             ZPAR(Ibuff) = ZPt(i)
             VX(Ibuff)      = VXt(i)
             VY(Ibuff)      = VYt(i)
             VZ(Ibuff)      = VZt(i)
             If(  XPAR(Ibuff).lt.1..or.YPAR(Ibuff).lt.1.
     &      .or.ZPAR(Ibuff).lt.1.)write (*,*) ' Error!!!',i, 
     &             XPAR(Ibuff),YPAR(Ibuff),ZPAR(Ibuff)
              If(  XPAR(Ibuff).ge.XMAX.or.YPAR(Ibuff).ge.XMAX
     &      .or.ZPAR(Ibuff).ge.XMAX)write (*,*) ' Error!!!',i, 
     &             XPAR(Ibuff),YPAR(Ibuff),ZPAR(Ibuff)
             If(  XPAR(Ibuff).ge.XMAX)Then
                XPAR(Ibuff) = XPAR(Ibuff) -1.e-4
                write (*,500) XPAR(Ibuff)
             endif 
             If(  YPAR(Ibuff).ge.XMAX)Then
                YPAR(Ibuff) = YPAR(Ibuff) -1.e-4
                write (*,500) YPAR(Ibuff)
             endif 
             If(  ZPAR(Ibuff).ge.XMAX)Then
                ZPAR(Ibuff) = ZPAR(Ibuff) -1.e-4
                write (*,500) ZPAR(Ibuff)
             endif 
 500            format(' fixed boundary. It is now=',g14.7)
            If(Ibuff.ge.NPAGE)Then
                KROW = KROW +1
                write (*,*) ' Write page=',KROW,' i=',i,Ibuff
                CALL WRIROW(KROW,1)
                Ibuff     =0
             EndIf
          EndDo
          If(Ibuff.ne.0)Then    ! write last incomplete page
                KROW = KROW +1
                write (*,*) ' Write page=',KROW,' i_buff=',Ibuff
                CALL WRIROW(KROW,1)
                Ibuff     =0
          EndIf
      RETURN
      END
C------------------------------------------
c                      Rearange species so that the highest resolution
c                                                      is fisrt, lowest resolution - last.
c                                                     Print statistics of spicies
      Subroutine Species_Change(Wtotal)
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMinitial.h'	  

           write (*,300)Nspecies
           write (16,300)Nspecies
 300       format('--- New number of species=',i3,/
     &          3x,'Specie',2x,
     +          'Nparticles Ntotal  Weight/part Weight_Total ')
 310       format(T5,i3,T13,i8,T22,i8,T30,2g10.3)
           Wtotal =0.
           Do i=1,Nspecies
              If(i.eq.1)Then
              Wtotal =wspecies(i)*lspecies(i)
              write (*,310)i,lspecies(i),lspecies(i),wspecies(i),Wtotal
              write (16,310)i,lspecies(i),lspecies(i),wspecies(i),Wtotal
              Else
                 ww =wspecies(i)* (lspecies(i)-lspecies(i-1))
              write (*,310) i,lspecies(i)-lspecies(i-1),
     &                           lspecies(i),wspecies(i),ww
              write (16,310) i,lspecies(i)-lspecies(i-1),
     &                           lspecies(i),wspecies(i),ww
                  Wtotal =Wtotal + ww
              Endif 
           Enddo 
           write (*,*) ' Total weight =',Wtotal

       Return
       End

C---------------------------------------------------
C                                  Read  current data from disk/tape,
C                                  Open files
C                                  Nrecl is the number of values in a record
C                                  Npage is the number of particles in a page
      SUBROUTINE RDTAPE
C---------------------------------------------------
      INCLUDE 'PMparameters.h'
C                                     Open file on a tape
      OPEN(UNIT=9,FILE='PMcrd.DAT',
     +                FORM='UNFORMATTED', STATUS='UNKNOWN')
C                                 Read control information
C                                 and see whether it has proper
C                                 structure
      READ  (9,err=10,end=10) HEADER,
     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     +                  NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     +                   ,Ocurv,extras
      WRITE (*,100) HEADER,
     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                  EKIN,EKIN1,EKIN2,
     +                  NROWC,NGRID,NRECL,Om0,Oml0,hubble,
     +                  Ocurv
      WRITE (16,100) HEADER,
     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                  EKIN,EKIN1,EKIN2,
     +                  NROWC,NGRID,NRECL,Om0,Oml0,hubble,
     +                  Ocurv
100   FORMAT(1X,'Header=>',A45,/
     +            1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     +            1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     +            1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I6,/
     +            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,' Hubble=',f7.3,/
     +            1x,' Omega_curvature=',F7.3)
      IF(NROW.NE.NROWC) THEN
         WRITE (*,*)
     +            ' NROW in PARAMETER and in TAPE-FILE are different'
         write (*,*) ' your configuration may be seriously wrong'
      ENDIF
      IF(NGRID.NE.NGRIDC) THEN
         WRITE (*,*)
     +           ' NGRID in PARAMETER and in TAPE-FILE are different:'
         write (*,*) ' Ngrid=',NGRID,' NgridC=',NGRIDC
         write (*,*) ' your configuration may be seriously wrong'
      ENDIF
C                                         Open work file on disk
 10   NBYTE = NRECL*4
      NACCES= NBYTE / nbyteword
c      NACCES= NACCES/4

      OPEN(UNIT=21,FILE='PMcrs0.DAT',ACCESS='DIRECT',
     +	               STATUS='UNKNOWN',RECL=NACCES)

 
      REWIND 9
      RETURN
      END
C---------------------------------------------
C                       Write current data to  disk/tape
C
      SUBROUTINE WRTAPE
C----------------------------------------------
      INCLUDE 'PMparameters.h'
C                                       write header and control data
      WRITE (9) HEADER,
     +           AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +           TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     +           NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5,
     +           Ocurv,extras
      REWIND 9
      RETURN
      END
C--------------------------------------------------
C                             Write current PAGE of particles (x,v) to disk
C                             NRECL - length of ROW block in words
      SUBROUTINE WRIROW(IROW,Ifile)
C--------------------------------------------------
      INCLUDE 'PMparameters.h'
        WRITE (20+Ifile,REC=IROW) RECDAT
      RETURN
      END
cC--------------------------------------------------
cC                             Set Weights of particles for
cC                             fast access
c      SUBROUTINE SetWeight
cC--------------------------------------------------
c      INCLUDE 'PMparameters.h'
c      If(Nspecies.eq.0)Then  ! old  constant weights
c         N_particles =lspecies(1)
c         Do i=1,N_particles
c            iWeight(i) =PARTW
c         EndDo
c      Else
c         N_particles =lspecies(Nspecies)
c         jstart =1
c         Do j=1,Nspecies
c            jend =lspecies(j)
c            Do k=jstart ,jend
c               iWeight(k) =wspecies(j)
c            EndDo
c            jstart =jend
c         EndDo
c      EndIf
c      write (*,*) ' Set Weights for ',N_particles,' particles'
c      RETURN
c      END
C--------------------------------------------
C                            Read current PAGE of particles  from disk
C                            NRECL - length of ROW block in words
      SUBROUTINE GETROW(IROW,Ifile)
C--------------------------------------------
      INCLUDE 'PMparameters.h'
         READ  (20+Ifile,REC=IROW) RECDAT
      RETURN
      END
C------------------------------------------------
C				                                       random number generator
      FUNCTION RANDd(M)
C------------------------------------------------
      DATA LC,AM,KI,K1,K2,K3,K4,L1,L2,L3,L4
     +	/453815927,2147483648.,2147483647,536870912,131072,256,
     +	 16777216,4,16384,8388608,128/
      ML=M/K1*K1
      M1=(M-ML)*L1
      ML=M/K2*K2
      M2=(M-ML)*L2
      ML=M/K3*K3
      M3=(M-ML)*L3
      ML=M/K4*K4
      M4=(M-ML)*L4
      M5=KI-M
      IF(M1.GE.M5)M1=M1-KI-1
      ML=M+M1
      M5=KI-ML
      IF(M2.GE.M5)M2=M2-KI-1
      ML=ML+M2
      M5=KI-ML
      IF(M3.GE.M5)M3=M3-KI-1
      ML=ML+M3
      M5=KI-ML
      IF(M4.GE.M5)M4=M4-KI-1
      ML=ML+M4
      M5=KI-ML
      IF(LC.GE.M5)ML=ML-KI-1
      M=ML+LC
      RANDd=M/AM
      RETURN
      END
C--------------------------------------
C				normal random numbers
      FUNCTION GAUSS(M)
C--------------------------------------
      X=0.
      DO  I=1,5
         X=X+RANDd(M)
      EndDo
      X2   =1.5491933*(X-2.5)
      GAUSS=X2*(1.-0.01*(3.-X2*X2))
      RETURN
      END
C--------------------------------------
C		normal random numbers
      FUNCTION GAUSS3()
C                        Uses ranlux for homogenous rand numbers   
C                        Uses Box-Muller inversion for gaussian
C                 Excellent quality and speed
c  N=       100000000  GAUSS3+ranlux ------
c  sigma= 0.99997     mean= 0.18692E-03
c  sigma Frac(>sig) Frac(<-sig)   True       n(>)    n(<)
c  1.00 0.1587     0.1586     0.1587     15869095 15857724
c  1.50 0.6682E-01 0.6680E-01 0.6681E-01  6681652  6679768
c  2.00 0.2277E-01 0.2273E-01 0.2275E-01  2276651  2273197
c  2.50 0.6222E-02 0.6206E-02 0.6210E-02   622190   620610
c  3.00 0.1356E-02 0.1347E-02 0.1350E-02   135628   134674
c  4.00 0.3242E-04 0.3145E-04 0.3170E-04     3242     3145
c
C--------------------------------------
      DIMENSION RanNum(2) 
      DATA iFlag/0/
      SAVE gSet
      If(iFlag.eq.0)Then
 1       CALL ranlux(RanNum,2)
         x1 = 2.*RanNum(1)-1.
         x2 = 2.*RanNum(2)-1.
         R =x1**2+x2**2
         If(R.ge.1.)GoTo 1
         F =sqrt(-2.*LOG(R)/R)
         gSet    = x1*F
         GAUSS3 = x2*F
         iFlag = 1
      Else
         GAUSS3 = gSet
         iFlag = 0
      EndIf 
      RETURN
      END


c-------------------------------------------------------------------- 
            
!   LUXURY LEVELS.
!  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
!           and Zaman, very long period, but fails many tests.
!  level 1  (p=48): considerable improvement in quality over level 0,
!           now passes the gap test, but still fails spectral test.
!  level 2  (p=97): passes all known tests, but theoretically still
!           defective.
!  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
!           correlations have very small chance of being observed.
!  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
c---------------------------------------------
      SUBROUTINE ranlux(rvec, lenv)
c                    returns a vector RVEC of LEN   
c                   32-bit random floating point numbers between 
c                   zero (not included) and one (also not incl.).
c           Next call to ranlux gives next LEN of random numbers
c---------------------------------------------
      INTEGER  lenv
      REAL        rvec(lenv)
      INTEGER  iseeds(24)
      INCLUDE 'luxury.h'
      DATA ndskip/ 0, 24, 73, 199, 365 /
      DATA i24,j24,luxlev/24,10,lxdflt/
      DATA notyet/.true./
      DATA in24,kount,mkount,carry/0,0,0,0./


!  NOTYET is .TRUE. if no initialization has been performed yet.
!              Default Initialization by Multiplicative Congruential

      IF (notyet) THEN
        notyet = .false.
        jseed = jsdflt
        inseed = jseed
        WRITE (6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ', jseed
        luxlev = lxdflt
        nskip = ndskip(luxlev)
        lp = nskip + 24
        in24 = 0
        kount = 0
        mkount = 0
      WRITE (6,'(A,I2,A,I4)') ' RANLUX DEFAULT LUXURY LEVEL= ', luxlev,   
     &                            '    p =', lp
      twom24 = 1.
        DO i = 1, 24
          twom24 = twom24 * 0.5
          k = jseed / 53668
          jseed = 40014 * (jseed-k*53668) - k * 12211
          IF (jseed.LT.0) jseed = jseed + icons
          iseeds(i) = MOD(jseed,itwo24)
        END DO
        twom12 = twom24 * 4096.
        DO i = 1, 24
          seeds(i) = REAL(iseeds(i)) * twom24
          next(i) = i - 1
        END DO
        next(1) = 24
        i24 = 24
        j24 = 10
        carry = 0.
        IF (seeds(24).EQ.0.) carry = twom24
      END IF

!          The Generator proper: "Subtract-with-borrow",
!          as proposed by Marsaglia and Zaman,
!          Florida State University, March, 1989

      DO ivec = 1, lenv
        uni = seeds(j24) - seeds(i24) - carry
        IF (uni.LT.0.) THEN
          uni = uni + 1.0
          carry = twom24
        ELSE
          carry = 0.
        END IF
        seeds(i24) = uni
        i24 = next(i24)
        j24 = next(j24)
        rvec(ivec) = uni
      !  small numbers (with less than 12 "significant" bits) are "padded".
        IF (uni.LT.twom12) THEN
          rvec(ivec) = rvec(ivec) + twom24 * seeds(j24)
      !        and zero is forbidden in case someone takes a logarithm
          IF (rvec(ivec).EQ.0.) rvec(ivec) = twom24 * twom24
        END IF
      !        Skipping to luxury.  As proposed by Martin Luscher.
        in24 = in24 + 1
        IF (in24.EQ.24) THEN
          in24 = 0
          kount = kount + nskip
          DO isk = 1, nskip
            uni = seeds(j24) - seeds(i24) - carry
            IF (uni.LT.0.) THEN
              uni = uni + 1.0
              carry = twom24
            ELSE
              carry = 0.
            END IF
            seeds(i24) = uni
            i24 = next(i24)
            j24 = next(j24)
          END DO
        END IF
      END DO
      kount = kount + lenv
      IF (kount.GE.igiga) THEN
        mkount = mkount + 1
        kount = kount - igiga
      END IF
      RETURN

      END ! SUBROUTINE ranlux

c---------------------------------------------
c                    Subroutine to initialize from one or three integers
      SUBROUTINE rluxgo(lux, ins, k1, k2)
c                initializes the generator from  
c               one 32-bit integer INT and sets Luxury Level LUX
c               which is integer between zero and MAXLEV, or if
c               LUX .GT. 24, it sets p=LUX directly.  K1 and K2
c               should be set to zero unless restarting at a break
c               point given by output of RLUXAT
c---------------------------------------------
      INCLUDE 'luxury.h'
      INTEGER  iseeds(24)
      
      IF (lux.LT.0) THEN
        luxlev = lxdflt
      ELSE IF (lux.LE.maxlev) THEN
        luxlev = lux
      ELSE IF (lux.LT.24.OR.lux.GT.2000) THEN
        luxlev = maxlev
        WRITE (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ', lux
      ELSE
        luxlev = lux
        DO ilx = 0, maxlev
          IF (lux.EQ.ndskip(ilx)+24) luxlev = ilx
        END DO
      END IF
      IF (luxlev.LE.maxlev) THEN
        nskip = ndskip(luxlev)
        WRITE (6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',
     &              luxlev,'     P=', nskip + 24
      ELSE
        nskip = luxlev - 24
        WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:', luxlev
      END IF
      in24 = 0
      IF (ins.LT.0) WRITE (6,'(A)')
     &        ' Illegal initialization by RLUXGO, negative input seed'
      IF (ins.GT.0) THEN
        jseed = ins
        WRITE (6,'(A,3I12)')' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
     &                                    jseed, k1, k2
      ELSE
        jseed = jsdflt
       WRITE (6,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
      END IF
      inseed = jseed
      notyet = .false.
      twom24 = 1.
      DO i = 1, 24
        twom24 = twom24 * 0.5
        k = jseed / 53668
        jseed = 40014 * (jseed-k*53668) - k * 12211
        IF (jseed.LT.0) jseed = jseed + icons
        iseeds(i) = MOD(jseed,itwo24)
      END DO
      twom12 = twom24 * 4096.
      DO i = 1, 24
        seeds(i) = REAL(iseeds(i)) * twom24
        next(i) = i - 1
      END DO
      next(1) = 24
      i24 = 24
      j24 = 10
      carry = 0.
      IF (seeds(24).EQ.0.) carry = twom24
      !        If restarting at a break point, skip K1 + IGIGA*K2
      !        Note that this is the number of numbers delivered to
      !        the user PLUS the number skipped (if luxury .GT. 0).
      kount = k1
      mkount = k2
      IF (k1+k2.NE.0) THEN
        DO iouter = 1, k2 + 1
          inner = igiga
          IF (iouter.EQ.k2+1) inner = k1
          DO isk = 1, inner
            uni = seeds(j24) - seeds(i24) - carry
            IF (uni.LT.0.) THEN
              uni = uni + 1.0
              carry = twom24
            ELSE
              carry = 0.
            END IF
            seeds(i24) = uni
            i24 = next(i24)
            j24 = next(j24)
          END DO
        END DO
      !         Get the right value of IN24 by direct calculation
        in24 = MOD(kount,nskip+24)
        IF (mkount.GT.0) THEN
          izip = MOD(igiga, nskip+24)
          izip2 = mkount * izip + in24
          in24 = MOD(izip2, nskip+24)
        END IF
      !       Now IN24 had better be between zero and 23 inclusive
        IF (in24.GT.23) THEN
          WRITE (6,'(A/A,3I11,A,I5)') 
     &       '  Error in RESTARTING with RLUXGO:', '  The values', ins, 
     &                k1, k2, ' cannot occur at luxury level', luxlev
          in24 = 0
        END IF
      END IF
      RETURN
      
      END ! SUBROUTINE rluxgo


C---------------------------------- Read in variables      
      REAL FUNCTION INPUT(text)
C------------------------------------------------
      Character text*(*)
          write (*,'(A,$)')text
          read (*,*) x
          INPUT =x
      Return
      End
C--------------------------------------   Fourier Transform
      SUBROUTINE SETF67(JBC1,JQ1,jbc,jp,jsl,ll1,k2,k3,k4,k7,
     &                  Qi,jndx,Uf,Vf)
c     -----------------------------------------------------
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMinitial.h'	  

      dimension Uf(marr) , Vf(marr)
      dimension Qi(mf67) , jndx(mf67)

      PI=DATAN(1.D0)*4.
       JBC=JBC1
       JQ=JQ1
       IF(JBC.LT.3) GO TO 101
       JQ=JQ-1
  101      K3=2**JQ
       K7=K3/2
       N5=K3/4
       I=1
       JNDX(I)=N5
       QI(I)=0.5*SQRT(2.)
       K=I
       I=I+1
  102     IL=I
       IF(I.EQ.K7) GO TO 104
  103     K1=JNDX(K)/2
       JNDX(I)=K1
       QI(I)=SIN(PI*FLOAT(K1)/FLOAT(K3))
       K1=K7-K1
       I=I+1
       JNDX(I)=K1
       QI(I)=SIN(PI*FLOAT(K1)/FLOAT(K3))
       K=K+1
       I=I+1
       IF(K.EQ.IL) GO TO 102
       GO TO 103
  104     RETURN
       END
C-----------------------------------------
       SUBROUTINE TFOLD(IS,L,ZZZ,k2)
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMinitial.h'	  


       DIMENSION ZZZ(MARR)
       IH2=K2/2-1
       DO 100 I=IS,IH2
       I1=I+L
       I2=K2-I+L
       A=ZZZ(I1)
       ZZZ(I1)=A-ZZZ(I2)
       ZZZ(I2)=A+ZZZ(I2)
  100     CONTINUE
       RETURN
       END
C------------------------------------
       SUBROUTINE NEG(I1,I3,I2,jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMinitial.h'	  
       dimension Uf(marr) , Vf(marr)
       dimension Qi(mf67) , jndx(mf67)
 
       DO 100 K=I1,I3,I2
       Vf(K)=-Vf(K)
  100     CONTINUE
       RETURN
       END
C-------------------------------------
       SUBROUTINE REVNEG(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMinitial.h'	  
       dimension Uf(marr) , Vf(marr)
       dimension Qi(mf67) , jndx(mf67)

       integer  jbc , jp, jsl , ll1 , k2 , k3 , k4 , k7

       DO 100 I=1,K7
       J=K3+1+I
       K=K4+1-I
       A=Vf(J)
       Vf(J)=-Vf(K)
       Vf(K)=-A
  100     CONTINUE
       RETURN
       END
C-------------------------------------
       SUBROUTINE ZEERO(L,jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMinitial.h'	  
       dimension Uf(marr) , Vf(marr)
       dimension Qi(mf67) , jndx(mf67)

       integer  jbc , jp, jsl , ll1 , k2 , k3 , k4 , k7

       DO 100 I=1,K2
       LI=L+I
       Uf(LI-1)=0.0
  100     CONTINUE
       RETURN
       END
C------------------------------------------
       SUBROUTINE TFOLD1(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMinitial.h'	  
       dimension Uf(marr) , Vf(marr)
       dimension Qi(mf67) , jndx(mf67)

       integer  jbc , jp, jsl , ll1 , k2 , k3 , k4 , k7

       II=K2-1
       DO 100 I=1,II
       I1=JSL+I
       I2=LL1-I
       A=Uf(I1)
       Uf(I1)=A+Uf(I2)
       Uf(I2)=A-Uf(I2)
  100     CONTINUE
       RETURN
       END
C-------------------------------------
       SUBROUTINE KFOLD(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMinitial.h'	  
       dimension Uf(marr) , Vf(marr)
       dimension Qi(mf67) , jndx(mf67)

       integer  jbc , jp, jsl , ll1 , k2 , k3 , k4 , k7

       JS1=K2
       I=1
       J5=JSL+K2
       IS1=JSL
       IC1=LL1
       JS1=JS1/2
       IF(JS1.NE.1) GO TO 200
       K1=JNDX(I)
       SN=QI(I)
       IS0=IS1
       IS1=IS1+JS1
       IC0=IC1
       IC1=IC1+JS1
       ODD1=SN*(Uf(IC1)-Uf(IS1))
       Vf(K1+1)=Uf(IC0)+ODD1
       K3MK1=K3-K1
       Vf(K3MK1+1)=Uf(IC0)-ODD1
       IF(JBC.LT.3) GO TO 110
       ODD2=SN*(Uf(IC1)+Uf(IS1))
       K3PK1=K3+K1
       Vf(K3PK1+1)=Uf(IS0)+ODD2
       K4MK1=K4-K1
       Vf(K4MK1+1)=-Uf(IS0)+ODD2
  110     RETURN
  200     SN=QI(I)
       IS1=IS1+JS1
       IC1=IC1+JS1
       J3=IS1+JS1
  210     IS0=IS1-JS1
       IC0=IC1-JS1
       ODD1=SN*(Uf(IC1)-Uf(IS1))
       ODD2=SN*(Uf(IC1)+Uf(IS1))
       Uf(IC1)=Uf(IC0)-ODD1
       Uf(IS1)=-Uf(IS0)+ODD2
       Uf(IC0)=Uf(IC0)+ODD1
       Uf(IS0)=Uf(IS0)+ODD2
       IS1=IS1+1
       IC1=IC1+1
       IF(IS1.NE.J3) GO TO 210
       I=I+1
  300     IS1=JSL
       IC1=LL1
       JS1=JS1/2
       IF(JS1.EQ.1) GO TO 400
  310     SN=QI(I)
       I=I+1
       CS=QI(I)
       IS1=IS1+JS1
       IC1=IC1+JS1
       J3=IS1+JS1
  320     IS0=IS1-JS1
       IC0=IC1-JS1
       ODD1=CS*Uf(IC1)-SN*Uf(IS1)
       ODD2=SN*Uf(IC1)+CS*Uf(IS1)
       Uf(IC1)=Uf(IC0)-ODD1
       Uf(IC0)=Uf(IC0)+ODD1
       Uf(IS1)=-Uf(IS0)+ODD2
       Uf(IS0)=Uf(IS0)+ODD2
       IS1=IS1+1
       IC1=IC1+1
       IF(IS1.NE.J3) GO TO 320
       IS1=IS1+JS1
       IC1=IC1+JS1
       J3=IS1+JS1
  330     IS0=IS1-JS1
       IC0=IC1-JS1
       ODD1=SN*Uf(IC1)-CS*Uf(IS1)
       ODD2=CS*Uf(IC1)+SN*Uf(IS1)
       Uf(IC1)=Uf(IC0)-ODD1
       Uf(IC0)=Uf(IC0)+ODD1
       Uf(IS1)=-Uf(IS0)+ODD2
       Uf(IS0)=Uf(IS0)+ODD2
       IS1=IS1+1
       IC1=IC1+1
       IF(IS1.NE.J3) GO TO 330
       I=I+1
       IF(IS1.EQ.J5) GO TO 300
       GO TO 310
  400     K1=JNDX(I)
       SN=QI(I)
       I=I+1
       CS=QI(I)
       IS0=IS1
       IS1=IS1+JS1
       IC0=IC1
       IC1=IC1+JS1
       ODD1=CS*Uf(IC1)-SN*Uf(IS1)
       Vf(K1+1)=Uf(IC0)+ODD1
       K3MK1=K3-K1
       Vf(K3MK1+1)=Uf(IC0)-ODD1
       IF(JBC.LT.3) GO TO 410
       ODD2=SN*Uf(IC1)+CS*Uf(IS1)
       K3PK1=K3+K1
       Vf(K3PK1+1)=Uf(IS0)+ODD2
       K4MK1=K4-K1
       Vf(K4MK1+1)=-Uf(IS0)+ODD2
  410     IS1=IS1+1
       IC1=IC1+1
       K1=JNDX(I)
       IS0=IS1
       IS1=IS1+JS1
       IC0=IC1
       IC1=IC1+JS1
       ODD1=SN*Uf(IC1)-CS*Uf(IS1)
       Vf(K1+1)=Uf(IC0)+ODD1
       K3MK1=K3-K1
       Vf(K3MK1+1)=Uf(IC0)-ODD1
       IF(JBC.LT.3) GO TO 420
       ODD2=CS*Uf(IC1)+SN*Uf(IS1)
       K3PK1=K3+K1
       Vf(K3PK1+1)=Uf(IS0)+ODD2
       K4MK1=K4-K1
       Vf(K4MK1+1)=-Uf(IS0)+ODD2
  420     IS1=IS1+1
       IC1=IC1+1
       I=I+1
       IF(IS1.NE.J5) GO TO 400
       RETURN
       END

c      ------------------------------------------------------
       SUBROUTINE FOUR67(JBC1,JQ1,jp1,jsl1,ll11,k21,k71,
     &                   Qi,jndx,Uf,Vf)
c      ------------------------------------------------------
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMinitial.h'	  
       dimension Uf(marr) , Vf(marr)
       dimension Qi(mf67) , jndx(mf67)

       integer  jbc , jp, jsl , ll1 , k2 , k3 , k4 , k7

       JBC=JBC1
       JQ=JQ1
       jp = jp1
       jsl = jsl1
       ll1 = ll11
       k2 = k21
       k7 = k71
       A5=0.5*SQRT(2.0)
       K4=2**JQ
       K3=K4
       GO TO (103,103,101,102),JBC
  101     Uf(1)=Uf(1)/2.0
       Uf(K3+1)=Uf(1)
       K2=K3

       CALL TFOLD(0,1,Uf,k2)

       K3=K3/2
       JQ=JQ-1
       GO TO 103
  102     K3=K3/2
       JQ=JQ-1
  103     N5=K3/4
       K7=K3/2
       N11=3*K7
       K31=K3+1
       GO TO(300,400,500,600),JBC
  300     Uf(K31)=0.0
       Uf(1)=0.0
       K2=K3
       DO 301 I=2,JQ

       CALL TFOLD(1,1,Uf,k2)

  301     K2=K2/2
       Vf(K7+1)=Uf(2)
       JF=N5
       JSL=1
       DO 302 JP=2,JQ
       LL1=K2+1

       CALL ZEERO(1,jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)
       CALL KFOLD(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

       I1=3*JF+1
       I2=4*JF
       I3=I1+(K2/2-1)*I2

       CALL NEG(I1,I3,I2,jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

       K2=K2+K2
  302     JF=JF/2
       RETURN
  400     Uf(1)=Uf(1)/2.0
       Uf(K31)=Uf(K31)/2.0
       K2=K3
       DO 401 I=2,JQ

       CALL TFOLD(0,K31-K2,Uf,k2)

  401     K2=K2/2
       LL1=K31-K2
       A=Uf(LL1)+Uf(LL1+2)
       Vf(1)=-A-Uf(LL1+1)
       Vf(K31)=-A+Uf(LL1+1)
       Vf(K7+1)=Uf(LL1)-Uf(LL1+2)
       DO 402 JP=2,JQ
       JSL=K31-K2
       LL1=JSL-K2
       idum = jsl

       CALL ZEERO(idum,jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)
       CALL KFOLD(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

  402     K2=K2+K2

       CALL NEG(1,K31,2,jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

       RETURN
  500     K2=K3
       L2=K4
       DO 501 JP=2,JQ

       CALL TFOLD(1,1,Uf,k2)
       CALL TFOLD(0,L2-K2+1,Uf,k2)

  501     K2=K2/2
       LL1=L2-K2+1
       A=Uf(LL1)+Uf(LL1+2)
       Vf(K7+1)=2.0*(-Uf(LL1)+Uf(LL1+2))
       Vf(1)=2.0*(A+Uf(LL1+1))
       Vf(K31)=2.0*(A-Uf(LL1+1))
       Vf(N11+1)=2.0*Uf(2)
       DO 502 JP=2,JQ
       Uf(K2+1)=2.0*Uf(K2+1)
       JSL=K2+1

       CALL TFOLD1(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

       LL1=LL1-K2
       Uf(LL1)=-2.0*Uf(LL1)

       CALL KFOLD(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

  502     K2=K2+K2
       Vf(1)=Vf(1)*A5
       Vf(K31)=Vf(K31)*A5
       RETURN
  600     Uf(1)=Uf(1)*A5
       Uf(K31)=Uf(K31)*A5
       K2=K3
       DO 601 JP=2,JQ

       CALL TFOLD(0,K31-K2,Uf,k2)
       CALL TFOLD(1,K31,Uf,k2)

  601     K2=K2/2
       LL1=K31-K2
       A=Uf(LL1)+Uf(LL1+2)
       Vf(1)=2.0*(A+Uf(LL1+1))
       Vf(K31)=2.0*(A-Uf(LL1+1))
       Vf(K7+1)=2.0*(-Uf(LL1)+Uf(LL1+2))
       Vf(N11+1)=2.0*Uf(K3+2)
       DO 602 JP=2,JQ
       JSL=K31+K2
       Uf(JSL)=2.0*Uf(JSL)

       CALL TFOLD1(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

       LL1=LL1-K2
       Uf(LL1)=-2.0*Uf(LL1)

       CALL KFOLD(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

  602     K2=K2+K2

       CALL REVNEG(jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

       idum = k3

       CALL NEG(2,idum,2,jbc,jp,jsl,ll1,k2,k3,k4,k7,Qi,jndx,Uf,Vf)

       K2=K4

       CALL TFOLD(1,1,Vf,k2)

       Vf(K4+1)=0.0

       RETURN
       END

