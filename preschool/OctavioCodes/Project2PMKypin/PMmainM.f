C-------------------------------------------------
C                                                     Serial PM code
C	     16 November 1997  Anatoly Klypin (aklypin@nmsu.edu)
C	                           Astronomy Department, NMSU
C-------------------------------------------------
C	       NROW = number of particles in one dimension
C	       NGRID= number of cells        in one dimension
C          AEXPN= expansion parameter = 1/(1+z)
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMmesh.h'

C			                                     Read data and open files
      OPEN(16,FILE='Result.DAT',STATUS='UNKNOWN')
      CALL RDTAPE
      CALL SetWeight
         write (*,*) ' RDTAPE is done'
      WRITE (*,'(A,$)') ' Enter number of steps to go => '
      READ (*,*) NUMBS		         ! Make this number of steps
      IF(AEXPN.GE.1.)THEN         ! change this if you need to run
	      WRITE (*,*) ' Cannot run over a=1' !  beyond a=1
	      STOP
      ENDIF

C					                                     Main loop
      DO 200 IRN=1,NUMBS
         CALL DENSIT            ! Define density
            write (*,*)         ' Density is ok'
         CALL POTENT            ! Define potential
            WRITE (*,*)         ' POTENT IS DONE'
         CALL MOVE              ! Move particles
            WRITE (*,*)         ' MOVE   IS DONE'
         CALL ADDTIME           ! Advance time 
      IF(AEXPN.GE.1.) GOTO 300  ! Do again if a < 1
200   CONTINUE

300   CALL WRTAPE               ! Write data to the disk
      STOP
      END
C--------------------------------------------------
C	                                               Advance Aexpn, Istep, tIntg,...
C	                                              AEXPN = currnet expansion parameter
C                                                   ISTEP = current step
C                                                  ASTEP = step in the expansion parameter			    
      SUBROUTINE ADDTIME
C--------------------------------------------------
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMmesh.h'
      COMMON /ENERG/	 ENKIN,ENPOT

      ISTEP = ISTEP + 1
      AEXPN = AEXPN + ASTEP
C		                     Energy conservation
      IF(ISTEP.EQ.1)THEN
        EKIN1 = EKIN
        EKIN2 = 0.
        EKIN  = ENKIN
        AU0   = AEXP0*ENPOT
        AEU0  = AEXP0*ENPOT + AEXP0*(EKIN+EKIN1)/2.
        TINTG = 0.
        WRITE (*,40) ISTEP,AEXPN,EKIN,ENPOT,AU0,AEU0
        WRITE (16,40) ISTEP,AEXPN,EKIN,ENPOT,AU0,AEU0
40      FORMAT('**** STEP=',I3,' A=',F10.4,' E KIN=',E12.4,
     .                  ' E POT=',E12.4,/'      AU0,AEU0=',2E12.4)
      ELSE
        EKIN2 = EKIN1
        EKIN1 = EKIN
        EKIN  = ENKIN
        TINTG = TINTG +ASTEP*(EKIN1 +(EKIN -2.*EKIN1 +EKIN2)/24.)
        ERROR = ((AEXPN-ASTEP)*((EKIN+EKIN1)/2.+ENPOT)-AEU0+TINTG)/
     .                  ((AEXPN-ASTEP)*ENPOT)*100.
        WRITE (*,50)  ISTEP,AEXPN,ERROR,EKIN,ENPOT,TINTG
        WRITE (16,50) ISTEP,AEXPN,ERROR,EKIN,ENPOT,TINTG
      ENDIF
50    FORMAT('*** ',I4,' A=',F7.4,' Error(%)=',f7.2,
     .                  ' Ekin=',E11.3,' Epot=',E11.3,' Intg=',E11.3)
      RETURN
      END
C------------------------------------------------------------
C
C     Advance each particle:	    dA	  AEXPN     dA
C	   by one step	     I______._______I_______.______I	  ->  A
C			    i-1     .	    i	    .	  i+1	step
C				    .	  {Fi}	    .	   .
C				  { vx }  { x }     .	   .
C				  { vy }  { y }     ^	   .
C				    ._______________^	   .
C					    .		   ^
C					    .______________^
C
C			  0.5
C		  dP = - A     * Grad(Fi) * dA ; A =AEXPN
C			  i		i
C				   3/2
C		  dX =	 P(new)/A	* dA ; A      =AEXPN+dA/2
C				 i+1/2		i+1/2
C
      SUBROUTINE MOVE
C------------------------------------------------
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMmesh.h'
      COMMON /ENERG/	ENKIN,ENPOT
C                                    PCONST = factor to change velocities
C                                   XCONST = factor to change coordinates
C	                               Note: 0.5 is needed in Pconst because
C		                          Fi(i+1)-Fi(i-1) is used as gradient
      PCONST = - SQRT(AEXPN/(Om0+Oml0*AEXPN**3+Ocurv*AEXPN))*ASTEP*0.5
      Ahalf  =   AEXPN+ASTEP/2.
      XCONST =   ASTEP/
     .           SQRT(Ahalf*(Om0+Oml0*Ahalf**3+Ocurv*Ahalf))/Ahalf
	   SVEL = 0.                      ! counter for \Sum(v_i**2)
	   SPHI = 0.                      ! counter for \Sum(phi_i)
	   XN   = FLOAT(NGRID)+1.-1.E-8   ! N+1
	   YN   = FLOAT(NGRID)            ! N
      ifile      =1
      If(Nspecies.eq.0)Then ! old case: 1 complete set
         Npages =NROW          ! Number of pages
         N_in_last=NPAGE       ! Number of particles in last page
      Else                             ! multiple masses
         N_particles =lspecies(Nspecies)   ! Total number of particles
         Npages = (N_particles -1)/NPAGE +1
         N_in_last=N_particles -NPAGE*(Npages-1)
      EndIf
c      write (*,*) ' Pages=',Npages,' Species=',Nspecies
c      write (*,*) ' N_in_last=',N_in_last

      DO  IROW=1, Npages         ! Loop over particle pages
            In_page =NPAGE
            If(IROW.eq.Npages)In_page =N_in_last
c            write (*,*)' Move page=',IROW,' file=',ifile,' N=',In_page
            iL = NPAGE*(IROW-1)
         CALL GETROW(IROW,ifile) ! read in a page of particles
	    DO  IN=1, In_page          ! Loop over particles
  	       X=XPAR(IN)
	       Y=YPAR(IN)
	       Z=ZPAR(IN)
	       VVX=VX(IN)
	       VVY=VY(IN)
	       VVZ=VZ(IN)
                Ipart =IN+iL                     ! current particle number
                WPAR =iWeight(Ipart)   ! particles weight
 	       I=INT(X)
	       J=INT(Y)
	       K=INT(Z)
	       D1=X-FLOAT(I)
	       D2=Y-FLOAT(J)
	       D3=Z-FLOAT(K)
	       T1=1.-D1
	       T2=1.-D2
	       T3=1.-D3
	       T2T1 =T2*T1
	       T2D1 =T2*D1
	       D2T1 =D2*T1
	       D2D1 =D2*D1
	       I1=I+1
	          IF(I1.GT.NGRID)I1=1
	       J1=J+1
	          IF(J1.GT.NGRID)J1=1
	       K1=K+1
	          IF(K1.GT.NGRID)K1=1
	       K2=K+2
	          IF(K2.GT.NGRID)K2=K2-NGRID
	       K3=K-1
	          IF(K3.LT.1    )K3=NGRID
	       F111 =FI(I ,J ,K )  !  Read potential to Fij vars
	       F211 =FI(I1,J ,K )
	       F121 =FI(I ,J1,K )
	       F221 =FI(I1,J1,K )
   
	       F112 =FI(I ,J ,K1)
	       F212 =FI(I1,J ,K1)
	       F122 =FI(I ,J1,K1)
	       F222 =FI(I1,J1,K1)
   
	       F113 =FI(I ,J ,K2)
	       F213 =FI(I1,J ,K2)
	       F123 =FI(I ,J1,K2)
	       F223 =FI(I1,J1,K2)
   
	       F110 =FI(I ,J ,K3)
	       F210 =FI(I1,J ,K3)
	       F120 =FI(I ,J1,K3)
	       F220 =FI(I1,J1,K3)
   
	       I2=I+2
	          IF(I2.GT.NGRID)I2=I2-NGRID
	       J2=J+2
	          IF(J2.GT.NGRID)J2=J2-NGRID
	       F311 =FI(I2,J ,K )
	       F321 =FI(I2,J1,K )
	       F131 =FI(I ,J2,K )
	       F231 =FI(I1,J2,K )
   
	       F312 =FI(I2,J ,K1)
	       F322 =FI(I2,J1,K1)
	       F132 =FI(I ,J2,K1)
	       F232 =FI(I1,J2,K1)
   
	       I0=I-1
	          IF(I0.LT.1)I0=NGRID
	       J0=J-1
	          IF(J0.LT.1)J0=NGRID
	       F011 =FI(I0,J ,K )
	       F021 =FI(I0,J1,K )
	       F101 =FI(I ,J0,K )
	       F201 =FI(I1,J0,K )
   
	       F012 =FI(I0,J ,K1)
	       F022 =FI(I0,J1,K1)
	       F102 =FI(I ,J0,K1)
	       F202 =FI(I1,J0,K1)
C				 Find {2*gradient} in nods
	          GX111 =F211 -F011
	          GX211 =F311 -F111
	          GX121 =F221 -F021
	          GX221 =F321 -F121
       
	          GX112 =F212 -F012
	          GX212 =F312 -F112
	          GX122 =F222 -F022
	          GX222 =F322 -F122
       
	          GY111 =F121 -F101
	          GY211 =F221 -F201
	          GY121 =F131 -F111
	          GY221 =F231 -F211
       
	          GY112 =F122 -F102
	          GY212 =F222 -F202
	          GY122 =F132 -F112
	          GY222 =F232 -F212
       
	          GZ111 =F112 -F110
	          GZ211 =F212 -F210
	          GZ121 =F122 -F120
	          GZ221 =F222 -F220
       
	          GZ112 =F113 -F111
	          GZ212 =F213 -F211
	          GZ122 =F123 -F121
	          GZ222 =F223 -F221
C				 Interpolate to the point
      GX=PCONST*(T3*(
     .		    T2T1*GX111+T2D1*GX211 +D2T1*GX121+D2D1*GX221 )+
     .		 D3*(
     .		    T2T1*GX112+T2D1*GX212 +D2T1*GX122+D2D1*GX222 ))

      GY=PCONST*(T3*(
     .		    T2T1*GY111+T2D1*GY211 +D2T1*GY121+D2D1*GY221 )+
     .		 D3*(
     .		    T2T1*GY112+T2D1*GY212 +D2T1*GY122+D2D1*GY222 ))

      GZ=PCONST*(T3*(
     .		    T2T1*GZ111+T2D1*GZ211 +D2T1*GZ121+D2D1*GZ221 )+
     .		 D3*(
     .		    T2T1*GZ112+T2D1*GZ212 +D2T1*GZ122+D2D1*GZ222 ))

C				 Find potential of the point
      FP=	 T3*(
     .		    T2T1*F111+T2D1*F211 +D2T1*F121+D2D1*F221 )+
     .		 D3*(
     .		    T2T1*F112+T2D1*F212 +D2T1*F122+D2D1*F222 )
	    SPHI = SPHI + FP*WPAR

          If(ipart.le.3)write (70,80) ipart,IROW,X,Y,Z,VVX,VVY,VVZ,
     &                     WPART,GX,GY,GZ,PCONST,XCONST
 80       format(' i=',2i7,' x=',3g11.4,' v=',3g11.4,/3x,f7.3,' g=',
     &                 3g11.4,' const=',2g11.4)

	    VVX =VVX+GX             ! Move points
	    VVY =VVY+GY
	    VVZ =VVZ+GZ
	    X	=X  +VVX*XCONST
	    Y	=Y  +VVY*XCONST
	    Z	=Z  +VVZ*XCONST
		  IF(X.LT.1.)X=X+YN     ! Periodical conditions
		  IF(X.GE.XN)X=X-YN
		  IF(Y.LT.1.)Y=Y+YN
		  IF(Y.GE.XN)Y=Y-YN
		  IF(Z.LT.1.)Z=Z+YN
		  IF(Z.GE.XN)Z=Z-YN
C			                                            Kinetic energy
	    SVEL=SVEL+(VVX**2+VVY**2+VVZ**2)*WPAR
C			                                        
	     XPAR(IN)=X            ! Write new coordinates
	     YPAR(IN)=Y
	     ZPAR(IN)=Z
	     VX(IN)=VVX
	     VY(IN)=VVY
	     VZ(IN)=VVZ
      ENDDO
C			   Write row to disk
       CALL WRIROW(IROW,ifile)
      ENDDO
C			   Set energies:
C			   Kin energy now at A(i+1/2)
C			   Pot energy		at A(i)
      ENKIN = SVEL / 2. / (AEXPN+ASTEP/2.)**2
      ENPOT = SPHI /2.
      RETURN
      END
C-------------------------------------------------
C	    Find potential on Grid FI:	DENSITY    ->	POTENTIAL
C
C		   O 1		    ^ - Fourier component
C		   |
C	     1	   |-4	 1	^      ^	2Pi
C	     O-----O-----O     Fi    =	Rho	/ (2cos(---  (i-1))+
C		   |	   i,j		i,j	Ngrid
C		   |
C		   O 1			  2Pi
C				       2cos(---  (j-1))-4)
C		       ^			Ngrid
C		       Fi	= 1 (?)
C			 11
C		   2
C		NABLA  Fi = 3/2  /A * (Rho - <Rho>) ;
C		   X
C			      <Rho> = 1

      SUBROUTINE POTENT
C---------------------------------------------
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMmesh.h'
      DIMENSION       GREENC(NGRID)

C			Set  a coefficient in Poisson eq.
C	       		  (2*NGRID)**3 - from FFT
      TRFI=1.5/AEXPN/(2.*NGRID)**3*Om0
C			Set Green function components
      P16 = 2.*3.14159265
      NGRID2=NGRID/2+2
      DO I=1,NGRID
	       XX =P16*(I-1.)/NGRID
	       GREENC(I) =2.*COS(XX)
	       IF(I.GE.NGRID2) GREENC(I) =-GREENC(I)
      END DO
C			Ngrid = 2**IQ
		       IQ = INT(ALOG(FLOAT(NGRID))/ALOG(2.)+0.5)
      Ndex=1
      IB1  =3
C					  ALONG X-DIRECTION
      CALL SETF67(IB1,IQ)
100   CONTINUE
      DO K=1,NGRID
	    DO J=1,NGRID
		  DO I=1,NGRID
			Zf(I) =FI(I,J,K)
		  END DO

		  CALL FOUR67(IB1,IQ)

		  DO I=1,NGRID
			FI(I,J,K) =Yf(I)
		  END DO
	    END DO
C					  ALONG Y-DIRECTION
	    DO I=1,NGRID
		  DO J=1,NGRID
			Zf(J) =FI(I,J,K)
		  END DO

		  CALL FOUR67(IB1,IQ)

		  DO J=1,NGRID
			FI(I,J,K) =Yf(J)
		  END DO
	    END DO
      END DO
c					  EXIT IF IT IS THE SECOND LOOP
      IF(Ndex.EQ.2) RETURN
C					  ALONG Z-DIRECTION

      DO J=1,NGRID
	    IB1 =3
	    DO I=1,NGRID
		  DO K=1,NGRID
			Zf(K) =FI(I,J,K)
		  END DO

		  CALL FOUR67(IB1,IQ)

		  DO K=1,NGRID
			FI(I,J,K) =Yf(K)
		  END DO
	    END DO
C					  BACK IN Z
	    A3 = GREENC(J) -6.
	    IB1= 4
	    DO I=1,NGRID
		  A2 = GREENC(I) + A3
		  DO K=1,NGRID
			A1 =A2 +GREENC(K)
			IF(ABS(A1).LT.1.E-4) A1=1.
			Zf(K) =FI(I,J,K)*TRFI/A1
		  END DO

		  CALL FOUR67(IB1,IQ)

		  DO K=1,NGRID
			FI(I,J,K) =Yf(K)
		  END DO
	    END DO
      END DO
      Ndex =2
      GOTO 100

      END

C------------------------------------------------------
      SUBROUTINE DENSIT
C------------------------------------------------------
      INCLUDE 'PMparameters.h'
      INCLUDE 'PMmesh.h'
      REAL*8      summass, summ2
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

      ifile      =1
      If(Nspecies.eq.0)Then ! old case: 1 complete set
         Npages =NROW          ! Number of pages
         N_in_last=NPAGE       ! Number of particles in last page
      Else                             ! multiple masses
         N_particles =lspecies(Nspecies)   ! Total number of particles
         Npages = (N_particles -1)/NPAGE +1
         N_in_last=N_particles -NPAGE*(Npages-1)
      EndIf
c      write (*,*) ' Pages=',Npages,' Species=',Nspecies
c      write (*,*) ' N_in_last=',N_in_last

      DO  IROW=1, Npages         ! Loop over particle pages
            In_page =NPAGE
            If(IROW.eq.Npages)In_page =N_in_last
c            write (*,*)' Move page=',IROW,' file=',ifile,' N=',In_page
            iL = NPAGE*(IROW-1)
         CALL GETROW(IROW,ifile) ! read in a page of particles
	    DO  IN=1, In_page          ! Loop over particles
                Ipart =IN+iL                     ! current particle number
                WPAR =iWeight(Ipart)   ! particles weight
	       X=XPAR(IN)
	       Y=YPAR(IN)
	       Z=ZPAR(IN)
	       I=INT(X)
	       J=INT(Y)
	       K=INT(Z)
c                 If(I.le.0)write (*,*) ' X:',X,Y,Z,I,' Irow=',IROW,IN,ifile
c                 If(J.le.0)write (*,*) ' Y:',X,Y,Z,I,' Irow=',IROW,IN,ifile
c                 If(K.le.0)write (*,*) ' Z:',X,Y,Z,I,' Irow=',IROW,IN,ifile
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
             IF(Ipart.eq.1)Then 
                write (70,90) ipart, WPAR,X,Y,Z,
     &                              i,j,k,T1,T2,T3
 90             format(' i=',i6,' w=',f8.2,3f8.3,3i3,/
     &                        12x,3f8.3)
c                STOP
             EndIf

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
      summass =0.
      summ2     =0.
      do k=1,NGRID
         Do j=1,NGRID
            Do i=1,NGRID
               dens =FI(i,j,k)
               summass =summass +dens
               summ2     =summ2     +dens**2
            Enddo 
            Enddo 
            Enddo 
            write (*,*) ' Summ density=',summass
            write (*,*) ' Sum2 density=',summ2

      RETURN
      END
