!-------------------------------------------------------------------- 
!                    Run simple N-body simulations
!                    using leap-frog integration scheme
!
!                    Logic: - read parameters, coordinates, velocities, 
!                             masses from file 'particles.dat' 
!                           - run for RunTime
!                           - save the data for continuation of the run 
!-------------------------------------------------------------------- 
      INCLUDE 'nbody1.h'
      REAL*8    Ekin,Epot,Etot,yc(3)

                          ! open files
      Open(1,file='particlesA.dat',status='old', action='read') !form='unformatted')
      Open(2,file='DATA/trajectories.dat')
      Open(3,file='DATA/energies.dat',position='append')
      Open(8,file='DATA/analysis.dat',position='append')

                          ! read data
      Read(1,*) Np,dt,time,iStep,epsilon,RunTime,PrintTime
      RunTime =0.10
      write(*,*) ' Number of particles=',Np,' Time=',time
      write(*,*) ' Time step          =',dt,' Step=',iStep
      write(*,*) ' Force softening    =',epsilon
      write(*,*) ' Run Time           =',RunTime
      Do i=1,Np
      read(1,*)(Coords(k,i),k=1,10)
!            read(30,20) (Coords(k,i),k=1,10)
      rewind(1)
      EndDo 
         E0 =-0.88451

                          ! set parameters
      Nsteps     = INT(RunTime/dt)  ! number of steps for this run
      dPrintTime = 0.10
      PrintTime  = time +dPrintTime ! next output time
      Ntraj      = 5                ! save this number of trajectories
      eps2       = epsilon**2
      iSnapshot  = 1                ! current snapshot
      ienergy    = 50              ! save energies this often
      indx       = 0 
	  
	  write(*,*) ' Number of steps=',Nsteps,' Time=',time
      write(*,*) ' Time step          =',dt,' Step=',iStep
      write(*,*) ' next output time    =',PrintTime
      
	  
      Do i=1,Nsteps                 ! main loop
         Call GetAccelerations
         Call MoveParticles
         time = time + dt
         iStep= iStep+ 1
	      write(*,*) 'time and step', time,iStep  
		 	 
         If(mod(i,ienergy).eq.0)Then
         Call WriteTrajectories(Ntraj)
            Call Energies(Ekin,Epot)
            Etot = Ekin+Epot
            write(*,20) time,iStep,Ekin,Epot,(Etot/E0-1)*100.,Ekin/Etot
            write(3,20) time,iStep,Ekin,Epot,Etot,Ekin/Etot
!            If(ABS(Etot/E0-1)*100..gt.0.1)write (*,*)' too large Error'
            Call Analyze(indx,yc)
         EndIf 
         If(time.ge.PrintTime)Then
            Call WriteSnapshot(yc)
            PrintTime  = PrintTime + dPrintTime
         EndIf 
      EndDo 
      Call SaveMoment
 20   format(f8.3,i8,6g13.5)

      Stop
      End
!-------------------------------------------------------------------- 
!               Find Acceleations for all particles
      SUBROUTINE GetAccelerations
!-------------------------------------------------------------------- 
      INCLUDE 'nbody1.h'
      REAL*8   acc,rr

      Do i=1,Np         ! set acceleration counters to zero
         Do k=7,9
            Coords(k,i) =0.
         EndDo 
      EndDo 
       Do i=1,Np-1       ! sum contributions pair-wise
         Do j=i+1,Np
            rr = (sqrt(
     &           (Coords(1,i)-Coords(1,j))**2+
     &           (Coords(2,i)-Coords(2,j))**2+
     &           (Coords(3,i)-Coords(3,j))**2+ eps2 ) )**3
            Do k=1,3
               acc           =  (Coords(k,i)-Coords(k,j))/rr
               Coords(k+6,i) =   Coords(k+6,i) -Coords(10,j)*acc
               Coords(k+6,j) =   Coords(k+6,j) +Coords(10,i)*acc
            EndDo 
         EndDo 
      EndDo 
      Return
      End
!-------------------------------------------------------------------- 
!               
      SUBROUTINE MoveParticles
!-------------------------------------------------------------------- 
      INCLUDE 'nbody1.h'

      Do i=1,Np         
         Do k=1,3
            Coords(k+3,i) = Coords(k+3,i) + Coords(k+6,i)*dt
            Coords(k  ,i) = Coords(k  ,i) + Coords(k+3,i)*dt
         EndDo 
      EndDo 
      Return
      End
!-------------------------------------------------------------------- 
!               write Ntraj trajectories into file 2
      SUBROUTINE WriteTrajectories(Ntraj)
!-------------------------------------------------------------------- 
      INCLUDE 'nbody1.h'

      write(2,10)time,iStep,((Coords(i,j),i=1,3),j=1,Ntraj)

 10   format(g12.5,i8,20(2x,3f9.4))
      Return
      End
