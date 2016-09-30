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
      READ  (9,err=20,end=20) HEADER,
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
     +          1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     +            1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     +            1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I9,/
     +            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,' Hubble=',f7.3,/
     +            1x,' Omega_curvature=',F7.3)
      IF(NROW.NE.NROWC) THEN
         WRITE (*,*)
     +            ' NROW in PARAMETER and in TAPE-FILE are different'
      ENDIF
      IF(NGRID.NE.NGRIDC) THEN
         WRITE (*,*)
     +           ' NGRID in PARAMETER and in TAPE-FILE are different:'
         write (*,*) ' Ngrid=',NGRID,' NgridC=',NGRIDC
      ENDIF
C                                         Open work file on disk
      write (*,*)  ' start openning PMcrs0'
 10   NBYTE = NRECL*4
c 10   NBYTE = NRECL*2
      NACCES= NBYTE / nbyteword

c      NACCES= NACCES/4

      OPEN(UNIT=21,FILE='PMcrs0.DAT',ACCESS='DIRECT',
     +	               STATUS='UNKNOWN',RECL=NACCES)


c      OPEN(UNIT=21,FILE='PMcrs0.DAT',ACCESS='DIRECT',
c     +                 STATUS='UNKNOWN',RECL=1)



      REWIND 9
       write (*,*)  ' done openning PMcrs0'
      RETURN
 20   write (*,*) ' Error reading the header file: Abort'
      stop
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
C--------------------------------------------------
C                             Set Weights of particles for
C                             fast access
      SUBROUTINE SetWeight
C--------------------------------------------------
      INCLUDE 'PMparameters.h'
      If(Nspecies.eq.0)Then  ! old  constant weights
         N_particles =lspecies(1)
         Do i=1,N_particles
            iWeight(i) =PARTW
         EndDo
      Else
         N_particles =lspecies(Nspecies)
         jstart =1
         Do j=1,Nspecies
            jend =lspecies(j)
            Do k=jstart ,jend
               iWeight(k) =wspecies(j)
            EndDo
            jstart =jend
         EndDo
      EndIf
      write (*,*) ' Set Weights for ',N_particles,' particles'
      RETURN
      END
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

C---------------------------------- Read in variables      
      REAL FUNCTION INPUT(text)
C------------------------------------------------
      Character text*(*)
          write (*,'(A,$)')text
          read (*,*) x
          INPUT =x
      Return
      End
C--------------------------------------   Fourier Transorm
      SUBROUTINE SETF67(IBC1,IQ1)
C----------------------------------------------------
      INCLUDE 'PMparameters.h'
      PI=DATAN(1.D0)*4.
       IBC=IBC1
       IQ=IQ1
       IF(IBC.LT.3) GO TO 101
       IQ=IQ-1
  101      N3=2**IQ
       N7=N3/2
       N5=N3/4
       I=1
       INDEX(I)=N5
       SI(I)=0.5*SQRT(2.)
       K=I
       I=I+1
  102     IL=I
       IF(I.EQ.N7) GO TO 104
  103     K1=INDEX(K)/2
       INDEX(I)=K1
       SI(I)=SIN(PI*FLOAT(K1)/FLOAT(N3))
       K1=N7-K1
       I=I+1
       INDEX(I)=K1
       SI(I)=SIN(PI*FLOAT(K1)/FLOAT(N3))
       K=K+1
       I=I+1
       IF(K.EQ.IL) GO TO 102
       GO TO 103
  104     RETURN
       END
C-----------------------------------------
       SUBROUTINE TFOLD(IS,L,ZZZ)
       INCLUDE 'PMparameters.h'
       DIMENSION ZZZ(NARR)
       IH2=N2/2-1
       DO 100 I=IS,IH2
       I1=I+L
       I2=N2-I+L
       A=ZZZ(I1)
       ZZZ(I1)=A-ZZZ(I2)
       ZZZ(I2)=A+ZZZ(I2)
  100     CONTINUE
       RETURN
       END
C------------------------------------
       SUBROUTINE NEG(I1,I3,I2)
       INCLUDE 'PMparameters.h'
       DO 100 K=I1,I3,I2
       Yf(K)=-Yf(K)
  100     CONTINUE
       RETURN
       END
C-------------------------------------
       SUBROUTINE REVNEG
       INCLUDE 'PMparameters.h'
       DO 100 I=1,N7
       J=N3+1+I
       K=N4+1-I
       A=Yf(J)
       Yf(J)=-Yf(K)
       Yf(K)=-A
  100     CONTINUE
       RETURN
       END
C-------------------------------------
       SUBROUTINE ZEERO(L)
       INCLUDE 'PMparameters.h'
       DO 100 I=1,N2
       LI=L+I
       Zf(LI-1)=0.0
  100     CONTINUE
       RETURN
       END
C------------------------------------------
       SUBROUTINE TFOLD1
       INCLUDE 'PMparameters.h'
       II=N2-1
       DO 100 I=1,II
       I1=ISL+I
       I2=L1-I
       A=Zf(I1)
       Zf(I1)=A+Zf(I2)
       Zf(I2)=A-Zf(I2)
  100     CONTINUE
      RETURN
       END
C-------------------------------------
      SUBROUTINE KFOLD
       INCLUDE 'PMparameters.h'
       JS1=N2
       I=1
       J5=ISL+N2
       IS1=ISL
       IC1=L1
       JS1=JS1/2
       IF(JS1.NE.1) GO TO 200
       K1=INDEX(I)
       SN=SI(I)
       IS0=IS1
       IS1=IS1+JS1
       IC0=IC1
       IC1=IC1+JS1
       ODD1=SN*(Zf(IC1)-Zf(IS1))
       Yf(K1+1)=Zf(IC0)+ODD1
       N3MK1=N3-K1
       Yf(N3MK1+1)=Zf(IC0)-ODD1
       IF(IBC.LT.3) GO TO 110
       ODD2=SN*(Zf(IC1)+Zf(IS1))
       N3PK1=N3+K1
       Yf(N3PK1+1)=Zf(IS0)+ODD2
       N4MK1=N4-K1
       Yf(N4MK1+1)=-Zf(IS0)+ODD2
  110     RETURN
  200     SN=SI(I)
       IS1=IS1+JS1
       IC1=IC1+JS1
       J3=IS1+JS1
  210     IS0=IS1-JS1
       IC0=IC1-JS1
       ODD1=SN*(Zf(IC1)-Zf(IS1))
       ODD2=SN*(Zf(IC1)+Zf(IS1))
       Zf(IC1)=Zf(IC0)-ODD1
       Zf(IS1)=-Zf(IS0)+ODD2
       Zf(IC0)=Zf(IC0)+ODD1
       Zf(IS0)=Zf(IS0)+ODD2
       IS1=IS1+1
       IC1=IC1+1
       IF(IS1.NE.J3) GO TO 210
       I=I+1
  300     IS1=ISL
       IC1=L1
       JS1=JS1/2
       IF(JS1.EQ.1) GO TO 400
  310     SN=SI(I)
       I=I+1
       CS=SI(I)
       IS1=IS1+JS1
       IC1=IC1+JS1
       J3=IS1+JS1
  320     IS0=IS1-JS1
       IC0=IC1-JS1
       ODD1=CS*Zf(IC1)-SN*Zf(IS1)
       ODD2=SN*Zf(IC1)+CS*Zf(IS1)
       Zf(IC1)=Zf(IC0)-ODD1
       Zf(IC0)=Zf(IC0)+ODD1
       Zf(IS1)=-Zf(IS0)+ODD2
       Zf(IS0)=Zf(IS0)+ODD2
       IS1=IS1+1
       IC1=IC1+1
       IF(IS1.NE.J3) GO TO 320
       IS1=IS1+JS1
       IC1=IC1+JS1
       J3=IS1+JS1
  330     IS0=IS1-JS1
       IC0=IC1-JS1
       ODD1=SN*Zf(IC1)-CS*Zf(IS1)
       ODD2=CS*Zf(IC1)+SN*Zf(IS1)
       Zf(IC1)=Zf(IC0)-ODD1
       Zf(IC0)=Zf(IC0)+ODD1
       Zf(IS1)=-Zf(IS0)+ODD2
       Zf(IS0)=Zf(IS0)+ODD2
       IS1=IS1+1
       IC1=IC1+1
       IF(IS1.NE.J3) GO TO 330
       I=I+1
       IF(IS1.EQ.J5) GO TO 300
       GO TO 310
  400     K1=INDEX(I)
       SN=SI(I)
       I=I+1
       CS=SI(I)
       IS0=IS1
       IS1=IS1+JS1
       IC0=IC1
       IC1=IC1+JS1
       ODD1=CS*Zf(IC1)-SN*Zf(IS1)
       Yf(K1+1)=Zf(IC0)+ODD1
       N3MK1=N3-K1
       Yf(N3MK1+1)=Zf(IC0)-ODD1
       IF(IBC.LT.3) GO TO 410
       ODD2=SN*Zf(IC1)+CS*Zf(IS1)
       N3PK1=N3+K1
       Yf(N3PK1+1)=Zf(IS0)+ODD2
       N4MK1=N4-K1
       Yf(N4MK1+1)=-Zf(IS0)+ODD2
  410     IS1=IS1+1
       IC1=IC1+1
       K1=INDEX(I)
       IS0=IS1
       IS1=IS1+JS1
       IC0=IC1
       IC1=IC1+JS1
       ODD1=SN*Zf(IC1)-CS*Zf(IS1)
       Yf(K1+1)=Zf(IC0)+ODD1
       N3MK1=N3-K1
       Yf(N3MK1+1)=Zf(IC0)-ODD1
       IF(IBC.LT.3) GO TO 420
       ODD2=CS*Zf(IC1)+SN*Zf(IS1)
       N3PK1=N3+K1
       Yf(N3PK1+1)=Zf(IS0)+ODD2
       N4MK1=N4-K1
       Yf(N4MK1+1)=-Zf(IS0)+ODD2
  420     IS1=IS1+1
       IC1=IC1+1
       I=I+1
       IF(IS1.NE.J5) GO TO 400
       RETURN
       END
C-----------------------------------------
       SUBROUTINE FOUR67(IBC1,IQ1)
       INCLUDE 'PMparameters.h'
       IBC=IBC1
       IQ=IQ1
       A5=0.5*SQRT(2.0)
       N4=2**IQ
       N3=N4
       GO TO (103,103,101,102),IBC
  101     Zf(1)=Zf(1)/2.0
       Zf(N3+1)=Zf(1)
       N2=N3
       CALL TFOLD(0,1,Zf)
       N3=N3/2
       IQ=IQ-1
       GO TO 103
  102     N3=N3/2
       IQ=IQ-1
  103     N5=N3/4
       N7=N3/2
       N11=3*N7
       N31=N3+1
       GO TO(300,400,500,600),IBC
  300     Zf(N31)=0.0
       Zf(1)=0.0
       N2=N3
       DO 301 I=2,IQ
       CALL TFOLD(1,1,Zf)
  301     N2=N2/2
       Yf(N7+1)=Zf(2)
       JF=N5
       ISL=1
       DO 302 IP=2,IQ
       L1=N2+1
       CALL ZEERO(1)
       CALL KFOLD
       I1=3*JF+1
       I2=4*JF
       I3=I1+(N2/2-1)*I2
       CALL NEG(I1,I3,I2)
       N2=N2+N2
  302     JF=JF/2
       RETURN
  400     Zf(1)=Zf(1)/2.0
       Zf(N31)=Zf(N31)/2.0
       N2=N3
       DO 401 I=2,IQ
       CALL TFOLD(0,N31-N2,Zf)
  401     N2=N2/2
       L1=N31-N2
       A=Zf(L1)+Zf(L1+2)
       Yf(1)=-A-Zf(L1+1)
       Yf(N31)=-A+Zf(L1+1)
       Yf(N7+1)=Zf(L1)-Zf(L1+2)
       DO 402 IP=2,IQ
       ISL=N31-N2
       L1=ISL-N2
       CALL ZEERO(ISL)
       CALL KFOLD
  402     N2=N2+N2
       CALL NEG(1,N31,2)
       RETURN
  500     N2=N3
       L2=N4
       DO 501 IP=2,IQ
       CALL TFOLD(1,1,Zf)
       CALL TFOLD(0,L2-N2+1,Zf)
  501     N2=N2/2
       L1=L2-N2+1
       A=Zf(L1)+Zf(L1+2)
       Yf(N7+1)=2.0*(-Zf(L1)+Zf(L1+2))
       Yf(1)=2.0*(A+Zf(L1+1))
       Yf(N31)=2.0*(A-Zf(L1+1))
       Yf(N11+1)=2.0*Zf(2)
       DO 502 IP=2,IQ
       Zf(N2+1)=2.0*Zf(N2+1)
       ISL=N2+1
       CALL TFOLD1
       L1=L1-N2
       Zf(L1)=-2.0*Zf(L1)
       CALL KFOLD
  502     N2=N2+N2
       Yf(1)=Yf(1)*A5
       Yf(N31)=Yf(N31)*A5
       RETURN
  600     Zf(1)=Zf(1)*A5
       Zf(N31)=Zf(N31)*A5
       N2=N3
       DO 601 IP=2,IQ
       CALL TFOLD(0,N31-N2,Zf)
       CALL TFOLD(1,N31,Zf)
  601     N2=N2/2
       L1=N31-N2
       A=Zf(L1)+Zf(L1+2)
       Yf(1)=2.0*(A+Zf(L1+1))
       Yf(N31)=2.0*(A-Zf(L1+1))
       Yf(N7+1)=2.0*(-Zf(L1)+Zf(L1+2))
       Yf(N11+1)=2.0*Zf(N3+2)
       DO 602 IP=2,IQ
       ISL=N31+N2
       Zf(ISL)=2.0*Zf(ISL)
       CALL TFOLD1
       L1=L1-N2
       Zf(L1)=-2.0*Zf(L1)
       CALL KFOLD
  602     N2=N2+N2
       CALL REVNEG
       CALL NEG(2,N3,2)
       N2=N4
       CALL TFOLD(1,1,Yf)
       Yf(N4+1)=0.0
       RETURN
       END
