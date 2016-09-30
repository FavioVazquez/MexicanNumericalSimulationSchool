C ______________________________	   POWER SPECTRUM
C		                      25 april 1990
C		                     R. Kates (Potsdam), A. Klypin (NMSU)
C	     NROW = number of particles per row
C	     NGRID= number of cells on a side (normal grid)
       INCLUDE 'PMparameters.h'
       INCLUDE 'PMmesh.h'
       REAL     INPUT
C...................................................................
C			Read data and open files
      OPEN(17,FILE='Spectrum.DAT',STATUS='UNKNOWN')
      CALL RDTAPE

      write (*,*) ' RDTAPE is done'
      WRITE (17,100) HEADER,
     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,wspecies(1),
     +                  EKIN,
     +                  NROWC,NGRID,NRECL,Om0,Oml0,hubble,
     +                  Ocurv,lspecies(1),Nspecies
100   FORMAT(1X,'Header=>',A45,/
     +            1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     +            1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',E12.3,/
     +            1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I9,/
     +            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,' Hubble=',f7.3,/
     +            1x,' Omega_curvature=',F7.3,/
     +            1x,' Number of part.1specie=',i10,' Nspecies=',i3)
      Box = extras(100)
      If(Box.le.1.e-3)Then
          Box =INPUT(' Enter box size in Mpc/h   =')
      EndIf 
C
         CALL DENSIT    !   Find density
            write (*,*) ' Density is ok'
c      Call DensCR
         CALL DENTES(DELR) ! Statistics
         CALL POWER(Box)   ! Power spectrum

      STOP
      END

C-----------------------------------------------------
C		                                                  Compute mean density and rms
      SUBROUTINE DENTES(DELR)
C-----------------------------------------------------
       INCLUDE 'PMparameters.h'
       INCLUDE 'PMmesh.h'
      Dimension x(44),amass(44),imass(44)
      real*8 SUM,SUM2

      Data x/  0.,  0.2, 0.4, 0.6, 0.7,  0.8,  0.9,  1.,  1.1,   1.2,
     .         1.3, 1.5, 1.75, 2., 2.25, 2.5,  3.0,  3.5,  4.,    5.,
     .         7.,  10., 12., 15., 20.,  25.,  30.,  40.,  50.,  70.,
     .       100., 120.,150.,200.,300., 500., 700.,1000.,1200.,1500.,
     .       2.e3, 3.e3,5.e3,1.e4/
      Do i=1,44
         amass(i) =0.
         imass(i) =0
      EndDo
      
      SUM = 0.
      SUM2= 0.
      Nn  = 0
      Am  = 0.
      DO K=1,NGRID
      DO J=1,NGRID
      DO I=1,NGRID
       SUM = SUM + FI(I,J,K)
       SUM2= SUM2+ FI(I,J,K)**2
       Do l=44,1,-1
          Dens =FI(I,J,K) +1.
          If(Dens .GT. x(l))Then
             Nn       =Nn +1
             Am       =Am +Dens
             imass(l) =imass(l) +1
             amass(l) =amass(l) +Dens
             GoTo 10
          EndIf
       EndDo
 10    Continue
      ENDDO
      ENDDO
      ENDDO
      Total =float(NGRID)**3
      DENM = SUM/Total
      DELR = DSQRT(SUM2/Total-DENM**2)
      DensAv=PARTW*NROW**3/float(NGRID)**3
      DensOne= float(NGRID)**3/NROW**3
         WRITE (*,150)  DELR,DENM,DensAv,DensOne
         WRITE (17,150) DELR,DENM,DensAv,DensOne
 150     format(20('-'),' Density Distribution',20('-'),/
     &   '              Density is in units of average density',
     &                       ' in the Box',/20x,
     &                     ' RMS Delta Rho/Rho   =',G11.4,/20x,
     &                     ' Mean Delta Rho/Rho  =',G11.4,/20x,
     &                     ' Average density     =',G11.4,/20x,
     &                     ' One particle density=',G11.4,/
     &' bin    Rho     Rho   Sum(rho_i)/dRho/TotMass        Sum',/
     &'      <center>< left > Normalized density     Number of cells ',/
     &'      < bin  ><border>   Distribution         With this density')  
      Do k=1,44-1
         drho =x(k+1)-x(k)
         rhom =x(k)+drho/2.
         Sn   =drho*Total
            If(Sn .lt.1.e-7)Sn=1.
         write (*,20)  k,rhom,x(k),amass(k)/Sn,imass(k)
         write (17,20) k,rhom,x(k),amass(k)/Sn,imass(k)
 20      Format(i3,2f8.1,2x,G12.3,12X,i12)
      EndDo
         k =44
         write (*,20)  k,x(k),x(k),amass(k),imass(k)
         write (17,20) k,x(k),x(k),amass(k),imass(k)
      RETURN
      END

C----------------------------------------------------
C                                                     power spectrum   P(k)
C                                                     for given density field FI
      SUBROUTINE POWER(Box)
C----------------------------------------------------
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMmesh.h'
      PARAMETER ( Pi=3.141592653, NPOWER =NGRID)
      DIMENSION   DPOWER(NPOWER),iVolume(NPOWER)
      DIMENSION   dk(NPOWER)
      REAL*8       SIGMA,DPOWER

      NSPEC = NGRID/2
      If(Box.le.0..or.Box.gt.1.e+4)Then
         write (*,*) ' Wrong Box size=',Box,' STOP'
         Return
      EndIf
      SIGMA =0.
      Scale = 2.*Pi/Box
      DO M =1,NPOWER
         iVolume(M) =0
         DPOWER(M)  =0.
         dk(M) =0.
      ENDDO
C                                               FFT of density field
      write (*,*) '     Do fft of rho along x-axis'
      IFOUR = INT(ALOG(FLOAT(NGRID))/ALOG(2.)+0.5)
      FT    = (2.*NGRID)**3*(1.*NGRID)**3
      CALL SETF67(3,IFOUR)
      DO K=1,NGRID
         DO J=1,NGRID
            DO I=1,NGRID
               Zf(I) =FI(I,J,K)
            ENDDO
            CALL FOUR67(3,IFOUR)
            DO I=1,NGRID
               FI(I,J,K) = Yf(I)
            ENDDO
         ENDDO
      ENDDO
      write (*,*) '     ....y'
      DO K=1,NGRID
         DO I=1,NGRID
            DO J=1,NGRID
               Zf(J) =FI(I,J,K)
            ENDDO
            CALL FOUR67(3,IFOUR)
            DO J=1,NGRID
               FI(I,J,K) = Yf(J)
            ENDDO
         ENDDO
      ENDDO
      write (*,*) '     ....z'
      DO J=1,NGRID
         IF(J.LE.(NSPEC+1)) THEN
            MJ =(J-1)**2
         ELSE
            MJ =(J-1-NSPEC)**2
         ENDIF
         DO I=1,NGRID
            IF(I.LE.(NSPEC+1)) THEN
               MI =(I-1)**2
            ELSE
               MI =(I-1-NSPEC)**2
            ENDIF
            DO K=1,NGRID
               Zf(K) =FI(I,J,K)
            ENDDO
            CALL FOUR67(3,IFOUR)
            DO K=1,NGRID
               AMP = Yf(K)**2/FT
               IF((I.EQ.1.AND.J.EQ.1).AND.K.EQ.1)THEN
                  FI(I,J,K) = 0.
               ELSE
C                                          AMP is DeltaRho^2
                  SIGMA = SIGMA +AMP
C                                          Find Abs(K)
                  IF(K.LE.(NSPEC+1)) THEN
                     RK =SQRT( FLOAT((K-1)**2     +MI+MJ) )
                  ELSE
                     RK =SQRT( FLOAT((K-1-NSPEC)**2+MI+MJ) )
                  ENDIF
C                                          devide DeltaRho^2 by Abs(K)
                  FI(I,J,K) = AMP / RK
C                                          Power Spectrum
                  MP =INT(RK +0.5)
                  IF(MP.LE.NPOWER)THEN
                     DPOWER(MP) =DPOWER(MP) +AMP
                     iVolume(MP) =iVolume(MP) +1
                     dk(MP)      =dk(MP)  + RK
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
 
      WRITE (*,100) SQRT(SIGMA),Scale*NROW/2.
      WRITE (17,100) SQRT(SIGMA),Scale*NROW/2.
100   FORMAT(5x,'RMS DRho/Rho=',g10.3,' k_Nyquist for Particles=',f6.2,
     &        /4X,'n',2X,'k(h/Mpc)',1X,
     &       'Harmonics ','Power Spectrum(h^3)')
110   FORMAT(1x,I4,F8.3,I8,3X,5G12.4)
      DO I=1,NPOWER
         IF(iVolume(I).GT.0)THEN
            wk  = dk(I)/iVolume(I)*Scale
            Pk  = Box**3*DPOWER(I)/iVolume(I)
             
            WRITE (*,110)  I,wk,iVolume(I),Pk
            WRITE (17,110) I,wk,iVolume(I),Pk
         ENDIF
      ENDDO
      RETURN
      END
C---------------------------------------
C                                                       Density Field    
      SUBROUTINE DENSIT
C---------------------------------------
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMmesh.h'
      	XN   =FLOAT(NGRID)+1.-1.E-7
	      YN   =FLOAT(NGRID)
C				       Subtract mean density
      DO M3=1,NGRID
       DO M2=1,NGRID
	      DO M1=1,NGRID
	        FI(M1,M2,M3) = -1.
	      END DO
       END DO
      END DO

      W     = (FLOAT(NGRID)/FLOAT(NROW))**3
      PARTW = W  
      WPAR = wspecies(1)
      N =0
      xmin =1.e+5
      xmax =-1.e+5
      ymin =1.e+5
      ymax =-1.e+5
      zmin =1.e+5
      zmax =-1.e+5

C				       Loop over particles
       If(Nspecies.eq.0)Then ! old case: 1 complete set
         Npages =NROW          ! Number of pages
         N_in_last=NPAGE       ! Number of particles in last page
      Else                             ! multiple masses
         N_particles =lspecies(Nspecies)   ! Total number of particles
         Npages = (N_particles -1)/NPAGE +1
         N_in_last=N_particles -NPAGE*(Npages-1)
      EndIf     
      write (*,*) Nspecies,
     + lspecies(Nspecies),wspecies(Nspecies)
      Nsp_current =1
      WPAR     = wspecies(1) !*8 *8. 
      Lspec_next=lspecies(1)+1
      DO IROW=1,Npages
            In_page =NPAGE
            If(IROW.eq.Npages)In_page =N_in_last
            If(IROW/50*50.eq.IROW)
     +        write (*,*) ' Read page=',IROW,' N=',N,
     +             ' Npages=',Npages,In_page
            iL = NPAGE*(IROW-1)
C				       Loop over particles
         READ  (21,REC=IROW) RECDAT
	   DO   IN=1,In_page
            N   = N +1
            If(N.ge.Lspec_next)Then   ! next specie
               write (*,*) ' Another specie. Current=',Nsp_current
               write (*,*) '     Particle=',N,WPAR,
     &                         ' lspecies=',lspecies(Nsp_current),
     &                         ' next=',lspecies(Nsp_current+1)
               Nsp_current=Nsp_current+1
               Lspec_next =lspecies(Nsp_current)+1
               WPAR     = wspecies(Nsp_current)
            EndIf
        X=(XPAR(IN)-1.)*NGRID/FLOAT(NGRIDC) +1.
        Y=(YPAR(IN)-1.)*NGRID/FLOAT(NGRIDC) +1.
        Z=(ZPAR(IN)-1.)*NGRID/FLOAT(NGRIDC) +1.
           xmin =MIN(xmin,X)
           ymin =MIN(ymin,Y)
           zmin =MIN(zmin,Z)
           xmax =MAX(xmax,X)
           ymax =MAX(ymax,Y)
           zmax =MAX(zmax,Z)
	        I=INT(X)
	        J=INT(Y)
	        K=INT(Z)
            If(I.le.0)write (*,*) ' X:',X,Y,Z,I,' Irow=',IROW,IN,ifile
            If(J.le.0)write (*,*) ' Y:',X,Y,Z,I,' Irow=',IROW,IN,ifile
            If(K.le.0)write (*,*) ' Z:',X,Y,Z,I,' Irow=',IROW,IN,ifile
	        D1=X-FLOAT(I)
	        D2=Y-FLOAT(J)
	        D3=Z-FLOAT(K)
	        T1=1.-D1
	        T2=1.-D2
	        T3=1.-D3
	        T2W =T2*WPAR
	        D2W =D2*WPAR
	        I1=I+1
	           IF(I1.GT.NGRID)I1=1
	        J1=J+1
	           IF(J1.GT.NGRID)J1=1
	        K1=K+1
	           IF(K1.GT.NGRID)K1=1
		      FI(I ,J ,K ) =FI(I ,J ,K ) +T3*T1*T2W
		      FI(I1,J ,K ) =FI(I1,J ,K ) +T3*D1*T2W
		      FI(I ,J1,K ) =FI(I ,J1,K ) +T3*T1*D2W
		      FI(I1,J1,K ) =FI(I1,J1,K ) +T3*D1*D2W
    
		      FI(I ,J ,K1) =FI(I ,J ,K1) +D3*T1*T2W
		      FI(I1,J ,K1) =FI(I1,J ,K1) +D3*D1*T2W
		      FI(I ,J1,K1) =FI(I ,J1,K1) +D3*T1*D2W
		      FI(I1,J1,K1) =FI(I1,J1,K1) +D3*D1*D2W
	      ENDDO
        ENDDO 

      write (*,*) ' Number of particles read=', N
      write (*,*) ' X:',xmin,xmax
      write (*,*) ' Y:',ymin,ymax
      write (*,*) ' Z:',zmin,zmax
        RETURN
      END

C------------------------------------
c      Subroutine densCR
c      INCLUDE 'PMparameters.h'
c      INCLUDE 'PMmesh.h'
c      	XN   =FLOAT(NGRID)+1.-1.E-7
c	      YN   =FLOAT(NGRID)

c      SUMM = 0.
c      WD       =0.
c      open(40,file='dens_cr.dat',status='unknown')
c      read(40,*) FI
c      Do k=1,NROW
c      Do j=1,NROW
c      Do i=1,NROW
c         SUMM =SUMM +FI(i,j,k)
c         WD      =WD +FI(i,j,k)**2
c      EndDo
c      EndDo
c      EndDo
c      sigma = sqrt(MAX(1.d-10,WD/NROW**3))
c      write (*,*) ' Average density=',SUMM/NROW**3
c      write (*,*) ' RMS density=',sigma
c       IFOUR = INT(ALOG(FLOAT(NROW))/ALOG(2.)+0.5)
c       Write (*,*) ' Ifour =',IFOUR,' Nrow=',NROW

c       Return
c       End
