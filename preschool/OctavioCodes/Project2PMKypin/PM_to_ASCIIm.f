C ______________________________	 Convert PM format to ASCII 
C		             November 1997         A. Klypin (aklypin@nmsu.edu) 
C		                     
C          Scales internal PM coordinates and velocities 
C          to  coordinates in Mpc/h and km/s
C          In PM the coordinates are in the range 1 - (NGRID+1)
C                     velocities are P = a_expansion*V_pec/(x_0H_0)
C                     where     x_0 = comoving cell_size=Box/Ngrid
C                                    H_0 = Hubble at z=0
C                   
C	     NROW = number of particles in 1D
C	     NGRID= number of cells        in 1D
       INCLUDE 'PMparameters.h'
       REAL     INPUT
       logical inside
       Character  FileASCII*50,File1*120,File2*120
C...................................................................
C			Read data and open files
      WRITE (*,'(A,$)') 'Enter Name for output ASCII file = '
      READ  (*,'(A)') FileASCII
      OPEN(17,FILE=FileASCII,STATUS='UNKNOWN')
      WRITE (*,'(A,$)')'Enter snapshot (0-for standard files)= '
      READ  (*,*) moment
      If(moment.le.0)Then
         File1 ='PMcrd.DAT'
         File2 ='PMcrs0.DAT'
      Else
         write(File1,'(a,i5.5,a)')'PMcrd_',moment,'.DAT'
         write(File2,'(a,i5.5,a)')'PMcrs0_',moment,'.DAT'
      EndIf 
      write(*,*) ' Open files:'
      write (*,'(a)')  File1
      write (*,'(a)')  File2
      CALL RDdata(File1,File2)

      write (*,*) ' RDTAPE is done'
      Fraction =INPUT(' Enter fraction of particles     =')
      If(Fraction.le.0.)Stop
      Ifraction = INT(1./Fraction+0.01)
      Radius =INPUT(' Enter radius of sphere or zero for all in box  =')
      Rad2 = 0.
      If(Radius.gt.0.)Then    
         write (*,'(A,$)') ' Enter Center: x,y,z     ='
         read  (*,*) xc,yc,zc
         xl =xc -Radius
         xr =xc +Radius
         yl =yc -Radius 
         yr =yc +Radius 
         zl =zc -Radius  
         zr =zc +Radius  
         Rad2 = Radius**2

      EndIf
      xrand     = 0.
         ncount = 0
      write (*,*) ' Every ',Ifraction,' particle will be selected'
      WRITE (17,100) HEADER,
     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                  EKIN,
     +                  NROWC,NGRID,NRECL,Om0,Oml0,hubble,
     +                  xc,yc,zc,Fraction
100   FORMAT(1X,'Header=>',A45,/
     +          1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     +           1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',E12.3,/
     +           1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I6,/
     +           1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,' Hubble=',f7.3,/
     +           1x,' Center=',3F7.3,' Fraction =',F7.3)

      If(extras(100).ne.0.)Then
         Box = extras(100)
      else
        Box =INPUT(' Enter box size in Mpc/h   =')
      endif 
      write (*,*) ' Box size is =  ',Box,'Mpc/h'
      BoxV =Box*100.    ! Box size in km/s
      ScaleV = BoxV/AEXPN/NGRID  ! scale factor for Velocities
      ScaleC = Box/NGRID         ! scale factor for Coordinates

      xmax =-1.e+9
      xmin = 1.e+9
      vmax =-1.e+9
      vmin = 1.e+9
      Icount =0
      ifile      =1
         N_particles =lspecies(Nspecies)   ! Total number of particles
         Npages = (N_particles -1)/NPAGE +1
         N_in_last=N_particles -NPAGE*(Npages-1)
      write (*,*) ' Pages=',Npages,' Species=',Nspecies
      write (*,*) ' N_in_last=',N_in_last
C				       Loop over particles
       DO IROW=1,Npages
            In_page =NPAGE
            If(IROW.eq.Npages)In_page =N_in_last
c            If(IROW.eq.2)
             write (*,'(" Read page=",i4," min/max=",10g11.3)')
     &             IROW,xmax,xmin,vmax,vmin
            iL = NPAGE*(IROW-1)
	      CALL GETROW(IROW,ifile)
	      DO IN=1,In_page
	        X  =ScaleC* (XPAR(IN)-1.)
	        Y  =ScaleC* (YPAR(IN)-1.)
	        Z  =ScaleC* (ZPAR(IN)-1.)
             Vxs=ScaleV* VX(IN)
             Vys=ScaleV* VY(IN)
             Vzs=ScaleV* VZ(IN)
             Icount =Icount +1

             xmax =MAX(xmax,X,Y,Z)
             xmin =MIN(xmin,X,Y,Z)
             vmax =MAX(vmax,Vxs,Vys,Vzs)
             vmin =MIN(vmin,Vxs,Vys,Vzs)
                Ipart =IN+iL             ! current particle number
                Do i =1,Nspecies
                   if(Ipart.le.lspecies(i))Then
                      W =wspecies(i)
                      goto 50
                   endif 
                EndDo 
 50        If(Fraction.lt.0.999)xrand =RANDd(Nseed) ! random fraction
           inside =.true.
           If(Radius.gt.0.)Then
              if(X.lt.xl.or.X.gt.xr)inside =.false.
              if(Y.lt.yl.or.Y.gt.yr)inside =.false.
              if(Z.lt.zl.or.Z.gt.zr)inside =.false.
           EndIf
              If(xrand.le.Fraction.and.inside)Then
                 rr = (X-xc)**2+(Y-yc)**2+(Z-zc)**2
                 If(rr.le.Rad2)ncount =ncount +1
c                 write (17,'(3f10.5,f8.3)') X,Y,Z,W
                         write (17,200) X,Y,Z,Vxs,Vys,Vzs,W
 200              Format(3F10.5,3F9.2,f8.4)
              EndIf  
c              If(X.lt.0..or.Y.lt.0..or.Z.lt.0.)
c     &              write(70,210)Ipart,X,Y,Z,Vxs,Vys,Vzs,W
c              If(X.gt.Box.or.Y.gt.Box.or.Z.gt.Box)
c     &              write(70,210)Ipart,X,Y,Z,Vxs,Vys,Vzs,W
              Vv =sqrt(Vxs**2+Vys**2+Vzs**2+0.001)
              If(VV.gt.10000.)
     &              write(70,210)Ipart,X,Y,Z,Vxs,Vys,Vzs,W
              
 210          format(i8,3G11.3,3G11.3,f8.3)
	      ENDDO
        ENDDO 

      write (*,*) ' Scaled Coordinates were written to:',
     &            FileASCII
      write (*,*) ' Number of particles          =',Icount
      IF(Radius.gt.0.)
     &write (*,*) ' Number of particles inside Rad=',ncount
      write (*,*) ' Min/Max of coordinates(Mpc/h)=',xmin,xmax
      write (*,*) ' Min/Max of velocities (km/s) =',vmin,vmax

      END


C---------------------------------------------------
C                                  Read  current data from disk/tape,
C                                  Open files
C                                  Nrecl is the number of values in a record
C                                  Npage is the number of particles in a page
      SUBROUTINE RDdata(File1,File2)
C---------------------------------------------------
      INCLUDE 'PMparameters.h'
      CHARACTER  File1*(*),File2*(*)
C                                     Open file on a tape
      OPEN(UNIT=9,FILE=File1,
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
         write (*,*)  ' NROW,NGRID (PMparamters.h)=',NROW,NGRID
         write (*,*)  ' NROW,NGRID (PMcrd.DAT)    =',NROWC,NGRIDC
      ENDIF
      IF(NGRID.NE.NGRIDC) THEN
         WRITE (*,*)
     +           ' NGRID in PARAMETER and in TAPE-FILE are different:'
         write (*,*) ' Ngrid=',NGRID,' NgridC=',NGRIDC
      ENDIF
C                                         Open work file on disk
 10   NBYTE = NRECL*4
      NACCES= NBYTE / nbyteword

      OPEN(UNIT=21,FILE=File2,ACCESS='DIRECT',
     +	               STATUS='UNKNOWN',RECL=NACCES)


      REWIND 9
      RETURN
 20   write (*,*) ' Error reading the header file: Abort'
      stop
      END

