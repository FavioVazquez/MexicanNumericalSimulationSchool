c-------------------------------------------------------------------- 
c               Find EnergY of the system
      SUBROUTINE Energies(Ekin,Epot)
c-------------------------------------------------------------------- 
      INCLUDE 'nbody1.h'
      REAL*8   acc,rr,Ekin,Epot

      Ekin =0.
      Epot =0.

      Do i=1,Np         ! kinetic energy
            Ekin = Ekin +(Coords(4,i)**2+
     &                    Coords(5,i)**2+
     &                    Coords(6,i)**2)*Coords(10,i)
      EndDo 
      Ekin = Ekin/2.

      Do i=1,Np-1       ! sum contributions pair-wise
         Do j=i+1,Np
            rr = sqrt(
     &           (Coords(1,i)-Coords(1,j))**2+
     &           (Coords(2,i)-Coords(2,j))**2+
     &           (Coords(3,i)-Coords(3,j))**2+ eps2 ) 
            Epot = Epot - Coords(10,i)*Coords(10,j)/rr
         EndDo 
      EndDo 
      Return
      End
c-------------------------------------------------------------------- 
c               write snapshot of the system to file
      SUBROUTINE WriteSnapshot(yc)
c-------------------------------------------------------------------- 
      INCLUDE 'nbody1.h'
      REAL*8 yc(3)
                 ! construct file name using time
                 ! split:  time = i1.i2
      i1   = INT(time)
      i2   = MIN(INT((time-i1)*1000.+0.5),999)
      write(FileNAMES(1),'(a,i3.3,".",i3.3,a)')'DATA/Snap',i1,i2,'.dat'
      open(10,file=FileNAMES(1))
         write(10,'(i6,2g11.4,i8,3g12.4))')
     &      Np,dt,time,iStep,epsilon,RunTime,
     &                 PrintTime
         Do i=1,Np
            write(10,20) (Coords(k,i)-yc(k),k=1,3),
     &                   (Coords(k,i),k=4,10)
         EndDo 
      close(10)
 20   format(10g13.5)
      Return
      End
c-------------------------------------------------------------------- 
c               save current state of the system
      SUBROUTINE SaveMoment
c-------------------------------------------------------------------- 
      INCLUDE 'nbody1.h'

      close(1)
      Open(1,file='particlesA.dat',status='unknown', action='write')
      Open(30,file='DATA/particlesA.dat')

         write(1,*)'',
     &      Np,dt,time,iStep,epsilon,RunTime,
     &                 PrintTime
         write(30,10)
     &      Np,dt,time,iStep,epsilon,RunTime,
     &                 PrintTime
         Do i=1,Np
            write(1,*)'', (Coords(k,i),k=1,10)
            write(30,20) (Coords(k,i),k=1,10)
         EndDo 


 10   format(i6,2g11.4,i8,3g12.4)
 20   format(10g13.5)
      Return
      End
c-------------------------------------------------------------------- 
c               
      SUBROUTINE Analyze(indx,yc)
c-------------------------------------------------------------------- 
      PARAMETER (Nbin =18)  ! number of bins for mass profile
      INCLUDE 'nbody1.h'
      REAL*8 xc(3),vc(3),yc(3),wc(3),ycn(3),wcn(3),Rmax,Raverage,
     &      aMass,aMass2
      DIMENSION Radii(Nbin),aMrad(Nbin)
      DIMENSION dDen(-100:Nbin),rDen(-100:Nbin),nDen(-100:Nbin),
     &          sDen(-100:Nbin),uDen(-100:Nbin)
      DATA Radii/0.05,0.10,0.16,0.2,0.25,0.3,0.35,0.4,0.5,
     &           0.6,0.80,1.0,2.0,4.0,
     &           5.0,10.0,15.0,20.0/

         i1   = INT(time)
         i2   = MIN(INT((time-i1)*1000.+0.5),999)
      write(FileNAMES(1),'(a,i3.3,".",i3.3,a)')'DATA/Prof',i1,i2,'.dat'
      open(9,file=FileNAMES(1))
      Rtoolarge =Radii(Nbin)*10.
      drlog     = 0.12            ! step in log for dens profile
      Do k=1,3
         xc(k) =0.
         vc(k) =0.
         yc(k) =0.
         wc(k) =0.
      EndDo
      If(indx.eq.0)Then
         write(8,120)(Radii(i),i=1,Nbin)
         indx = 1
 120   format(T48,20f8.2)
      EndIf 

      ! find center of mass and average velocity
         Raverage =0.
         Rmax     =0.
         aMass2   =0.
         aMass    =0.
         Do i=1,Np         
            aM    = Coords(10,i)
            aMass = aMass + aM
            rr    = sqrt(Coords(1,i)**2+Coords(2,i)**2+Coords(3,i)**2)
            Do k=1,3
                xc(k) = xc(k) + Coords(k,i)*aM
                vc(k) = vc(k) + Coords(k+3,i)*aM
            EndDo 
            If(rr.lt.Rtoolarge)Then
               aMass2= aMass2 + aM
               Do k=1,3
                  yc(k) = yc(k) + Coords(k,i)*aM
                  wc(k) = wc(k) + Coords(k+3,i)*aM
               EndDo 
            EndIf 
         EndDo
         Do k=1,3       ! normalize results
            xc(k) =xc(k)/aMass
            vc(k) =vc(k)/aMass
            yc(k) =yc(k)/aMass2
            wc(k) =wc(k)/aMass2
         EndDo

      ! find center of the densiest part of the system
         Dens0 = Amass2/Rtoolarge**3
      Do iter=1,150
         Rtoolarge =Rtoolarge/1.1
         Raverage =0.
         aMass2   =0.
         inside   =0
         Do k=1,3
            ycn(k) =0.
            wcn(k) =0.
         EndDo

         Do i=1,Np         
            aM    = Coords(10,i)
            dx    = Coords(1,i) -yc(1)
            dy    = Coords(2,i) -yc(2)
            dz    = Coords(3,i) -yc(3)
            rr    = sqrt(dx**2+dy**2+dz**2)
            If(rr.lt.Rtoolarge)Then
               aMass2= aMass2 + aM
               inside= inside +1
               Do k=1,3
                  ycn(k) = ycn(k) + Coords(k,i)*aM
                  wcn(k) = wcn(k) + Coords(k+3,i)*aM
               EndDo 
            EndIf 
         EndDo 
         Do k=1,3       ! normalize results
            yc(k) =ycn(k)/aMass2
            wc(k) =wcn(k)/aMass2
         EndDo
         Dens1 = Amass2/Rtoolarge**3
         If(Dens1.lt.Dens0.or.inside.lt.Np/3)goto300
         Dens0 = Dens1
      EndDo 
      
 300  Do i=1,Nbin
         aMrad(i) =0.
      EndDo 
      Rtoolarge =Radii(Nbin)
      aMass2   =0.
      Raverage =0.
      Rmax     =0.
      do i=-100,Nbin
         dDen(i) =0.
         rDen(i) =0.
         nDen(i) =0
         sDen(i) =0.
         uDen(i) =0.
      EndDo 

      Do i=1,Np      ! get statistics relative to the center      
         aM      = Coords(10,i)
         dx      = Coords(1,i) -yc(1)
         dy      = Coords(2,i) -yc(2)
         dz      = Coords(3,i) -yc(3)
         rr      = sqrt(dx**2+dy**2+dz**2)
         vv      = Coords(4,i)**2+Coords(5,i)**2+Coords(6,i)**2
         vr      = (Coords(4,i)*dx+Coords(5,i)*dy+Coords(6,i)*dz)
     &                                             /max(rr,1.e-10)
         ind     = INT(log10(rr)/drlog+1.e5)-1e5
         ind     = max(min(ind,Nbin),-100)
         dDen(ind) = dDen(ind) + aM
         rDen(ind) = rDen(ind) + rr*aM
         nDen(ind) = nDen(ind) + 1
         sDen(ind) = sDen(ind) + vr*aM
         uDen(ind) = uDen(ind) + vv*aM
         
         If(rr.lt.Rtoolarge)Then
            aMass2 =aMass2 + aM
            Raverage = Raverage + rr**2*aM
            do j=1,Nbin
               if(rr.lt.Radii(j))Then
                  aMrad(j) =aMrad(j) +aM
                  goto 10
               EndIf 
            EndDo 
 10         Rmax =max(Rmax,rr)
         EndIf 
      EndDo
      Do j=2,Nbin   ! make mass cumulative
         aMrad(j) =aMrad(j) + aMrad(j-1)
      EndDo

       Raverage = sqrt(Raverage/aMass2)   ! normalize results
       write(8,101)time,iStep,aMass,aMass2,Rmax,Raverage,
     &    (aMrad(i),i=1,Nbin)
 100   format(a,f8.4,a,i8,a,2f8.3,a,2f8.3,20f8.4)
 101   format(f8.4,i8,2f8.3,2f8.3,20f8.4)
                            
 110   format(i4,2g12.4,3x,2g12.4,3x,g12.3,i5,4x,4f8.3)
                   !  write density profile
       write(9,100)' t=',time,' Step=',iStep,
     &                    ' Mass=',aMass,aMass2,
     &                   ' Rmax/Raverage=',Rmax,Raverage
       write(9,120)
       aMass2 =0.
       Do i=-100,Nbin
          If(nDen(i).gt.0)Then
             aM = dDen(i)
             aMass2 = aMass2 +aM
             rr = rDen(i)/aM
             vr = sDen(i)/aM
             v2 = sqrt(max( (uDen(i)/aM-vr**2),1.e-10))
             rout = 10.**((i+1)*drlog)
             rin  = 10.**(i*drlog)
             vol= 4.*pi/3.*(rout**3-rin**3)
             den= aM/vol
             dn = nDen(i)/vol
             write(9,110)i,rr,rin,den,dn,aM/nDen(i),nDen(i),
     &                   vr,v2,aMass2
          EndIf 
       EndDo 
       close(9)
      Return
      End
