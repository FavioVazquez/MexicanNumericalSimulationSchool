      PARAMETER (Nmax       =150000)   ! maximum number of particles
      PARAMETER (Nsnapshots =   10)   ! maximum number of snapshots
      PARAMETER (pi         =  3.14159265)   
      PARAMETER (pi4        =  4.*pi)   
      REAL*8           Coords,dt,time,epsilon
      CHARACTER*80     FileNAMES(Nsnapshots)
      COMMON /MAINDATA/Coords(10,Nmax)
      COMMON /CONTROL/ dt,time,Np,iStep,epsilon,RunTime,
     &                 dPrintTime
      COMMON /AUXILIARY/ eps2,Rmaxinit,Rmininit,rs      
c              Coords: 1-3 = x,y,z
c                      4-6 = vx,vy,vz
c                      7-9 = gx,gy,gz
c                      10  = mass
      DATA FileNAMES/'Snap.1.dat','Snap.2.dat','Snap.3.dat',
     &               'Snap.4.dat','Snap.5.dat','Snap.6.dat',
     &               'Snap.7.dat','Snap.8.dat','Snap.9.dat',
     &               'Snap.10.dat'/
