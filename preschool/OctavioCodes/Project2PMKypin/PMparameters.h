c      IMPLICIT REAL*8 (D)
      PARAMETER (NROW =128)      ! maximum number of particles in 1D
      PARAMETER (NGRID =256)      ! zero-level mesh in 1D
                                  ! must be equal to ng in ART
      PARAMETER (Nmaxpart = 256**3+5e6)  ! max number of particles for this run
                                  ! Nmaxpart must be less than NROW**3
      PARAMETER (nbyteword = 1)   ! defines length of direct-access:1or4

      PARAMETER (NPAGE = NROW**2)    ! number of particles in a record
      PARAMETER (NMAX=NGRID/2)
      PARAMETER (NRECL= NPAGE*6)  ! number of words in a record
c      PARAMETER (NRECL= NPAGE)  ! number of words in a record
      PARAMETER (NARR =NGRID+1)   !      need in FFT
      PARAMETER (NF67=NGRID/2)       !      need in FFT
c      PARAMETER (NARR =MAX(2*NROW+1,NGRID+1))   !      need in FFT
c      PARAMETER (NF67=MAX(NROW,NGRID/2))

      COMMON / CONTROL/ AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     +                  NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     +                  ,Ocurv,extras(100)
      COMMON / HEADDR/  HEADER
      CHARACTER*45      HEADER
      COMMON /FOURAR/Zf(NARR),Yf(NARR)
      COMMON /F67COM/
     +                 IBC,      IP,       ISL,     L1,     N2,
     +                 N3,       N4,        N7,
     +                 SI(NF67),    INDEX(NF67)
      COMMON / ROW /	XPAR(NPAGE),YPAR(NPAGE),ZPAR(NPAGE),
     +			VX(NPAGE),VY(NPAGE),VZ(NPAGE)
      COMMON / BWEIG/ iWeight(Nmaxpart)
      REAL  iWeight
      DIMENSION     RECDAT(NRECL),wspecies(10),lspecies(10)
      EQUIVALENCE    (RECDAT(1),XPAR(1))
     +                               ,(wspecies(1),extras(1)),
     +                               (lspecies(1),extras(11))








