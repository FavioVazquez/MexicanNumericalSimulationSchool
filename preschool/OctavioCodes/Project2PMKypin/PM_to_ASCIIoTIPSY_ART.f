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
       INCLUDE 'PMparART.h'
c       REAL     INPUT
       Character  FileASCII*50
       character sntime*6
c       character sntime*5  !para simulaciones con HART
       character tail*4
       character tailasc*4
       character tailtip*4
       character head1*6
       character head2*7
       dimension xpart(10000000,11),ypart(10000000,11),zpart(10000000,11
     +),vxpart(10000000,11),vypart(10000000,11),vzpart(10000000,11),nums
     +(11),x(100000000),y(100000000),z(100000000)!,vx(100000000),vy(10000
c     +0000),vz(100000000),
      dimension vxx(100000000),vyy(100000000),vzz(100000000),
     &lspecies1(10)
      Box  =1.0 !(box size Mpc/h)                                                                                                               
      scrip=1 !cada cuantas particulas se van a escribir en el archivo de salida tipsy
      BoxV =Box*100.    ! Box size in km/s                                                                                                 
      tail='.DAT'
      tailasc='.asc'
      tailtip='.tip'
      head1='PMcrda'
      head2='PMcrs0a'
      xmax =-1.e+9
      xmin = 1.e+9
      vmax =-1.e+9
      vmin = 1.e+9


C...................................................................
C			Read data and open files
c      WRITE (*,'(A,$)') 'Enter File Name for ASCII data = '
c      READ  (*,'(A)') FileASCII

      write(*,*)'Write snapshot time in scale factor units, e.p. 0.6000'
      read(*,*)sntime

c       Narg = IARGC ()
c       write(*,*)narg
cc       call getarg(2,sntime)
c       sntime='0.6920'
c      OPEN(17,FILE=FileASCII,STATUS='UNKNOWN')
      open(3,file=head1//sntime//tail,form='unformatted',status='unknown
     &')
c      open(3,file='PMcrd.DAT_ini',form='unformatted',status='unknown
c     &')
      read      (3,err=10,end=10) HEADER,
     &              AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &              TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &              NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5
     &              ,Ocurv,extras
      ScaleV = BoxV/AEXPN/NGRID  ! scale factor for Velocities                                                                            \
                                                        
      ScaleC = Box*aexpn/NGRID         ! scale factor for Coordinates 
c      WRITE (17,100) HEADER,
c     +                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
c     +                  EKIN,
c     +                  NROWC,NGRID,NRECL,Om0,Oml0,hubble,
c     +                  Ocurv
 100     FORMAT(1X,'Header=>',A45,/
     +            1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/                                                                  
     +            1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',E12.3,/
     +            1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I6,/
     +            1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,' Hubble=',f7.3,/
     +            1x,' Omega_curvature=',F7.3)

 10   npaget    = nrow**2    ! # of particles in a page
      nrecli    = npaget * 6 ! length of particle row in words
      nbyte  = nrecli * 4
      nbyteword=1!  defines length of direct-access record, 4 para HART, 1 para el resto
      nacces = nbyte / nbyteword
      xn=float(ngrid)+1.-1.E-7
      yn=float(ngrid)
c      nrowc=1024
      open ( 1 , file =head2//sntime//tail, access = 'direct',
     &           status = 'unknown', recl = nacces      )
c       open ( 1 , file ='PMcrs0.DAT_ini', access = 'direct',
c     &           status = 'unknown', recl = nacces      )

c         lspecies(0)=0
c         do i=1,nspecies
c          lspecies1(i)=(lspecies(i)-lspecies(i-1))
c         enddo
c         do i=1,nspecies
c           if (lspecies1(i).lt.0.0) then
c           lspecies(i)=0
c           else
c           if (lspecies1(i-1).lt.0.0) then
c           lspecies(i)=lspecies1(i)+lspecies(i-2)
c           else
c           lspecies(i)=lspecies1(i)+lspecies(i-1)
c           endif 
c          endif
c         enddo

         N_particles =lspecies(Nspecies)   ! Total number of particles                                                                    
         do i=1,nspecies
          write(*,*)'Particulas de la especie',i,' = ',lspecies(i),'weig
     &ht =', wspecies(i)
         enddo
         Npages = (N_particles -1)/NPAGE +1
         N_in_last=N_particles -NPAGE*(Npages-1)
      write (*,*) ' Pages=',Npages,' Species=',Nspecies
      write (*,*) ' Nparticles=',N_particles

      DO  IROW=1, Npages         ! Loop over particle pages                                                                                
            In_page =NPAGE
            If(IROW.eq.Npages)In_page =N_in_last
c            write (*,*)' Read page=',IROW,' file=',ifile,' N=',In_page                                                                    
            iL = NPAGE*(IROW-1)
         CALL GETROW(IROW,1) ! read in a page of particles                                                                                 
         DO  IN=1, In_page          ! Loop over particles                                                                                  
                ip =IN+iL                     ! current particle number                                                                    
c                  WPAR =iWeight(ip)   ! particles weight                                                                                  
c check rounding errors (REAL*8 <=> REAL*4)                                                                                                
                x(ip) = xpar(in)
c                IF(x(ip).LT.1.) x(ip) = x(ip) + yn
c                IF(x(ip).GT.xn) x(ip) = x(ip) - yn
                y(ip) = ypar(in)
c                IF(y(ip).LT.1.) y(ip) = y(ip) + yn
c                IF(y(ip).GT.xn) y(ip) = y(ip) - yn
                z(ip) = zpar(in)
c                IF(z(ip).LT.1.) z(ip) = z(ip) + yn
c                IF(z(ip).GT.xn) z(ip) = z(ip) - yn
             If(  x(ip).ge.xn)Then
                x(ip) = x(ip) -1.e-4
                write (*,*) x(ip),ip,1
             endif
             If(  y(ip).ge.xn)Then
                y(ip) = y(ip) -1.e-4
                write (*,*) y(ip),ip,2
             endif
             If(  z(ip).ge.xn)Then
                z(ip) = z(ip) -1.e-4
                write (*,*) z(ip),ip,3
             endif

                vxx(ip) = vx(in)
                vyy(ip) = vy(in)
                vzz(ip) = vz(in)
c           write(*,*)vx(in),vy(in),vz(in)
           X(ip)  =ScaleC* (X(Ip)-1.)                                                                                             
           Y(ip)  =ScaleC* (Y(Ip)-1.)
           Z(ip)  =ScaleC* (Z(Ip)-1.)          
           Vxx(ip)=ScaleV* VXx(Ip)
           Vyy(ip)=ScaleV* VYy(Ip)
           Vzz(ip)=ScaleV* VZz(Ip) 

c           write(17,*)ip,x(ip),y(ip),z(ip),vx(ip),vy(ip),vz(ip)

           xmax =MAX(xmax,X(ip),Y(ip),Z(ip))
           xmin =MIN(xmin,X(ip),Y(ip),Z(ip))
           vmax =MAX(vmax,Vxx(ip),Vyy(ip),Vzz(ip))
           vmin =MIN(vmin,Vxx(ip),Vyy(ip),Vzz(ip))

         Enddo
      Enddo

      hfact=hubble
c      halomass1 = halomass*hfact
      diskmass1 = diskmass/hfact
      pdmass=diskmass1/diskpart
      qhmass=0.
      do i=1,nspecies
        if(i.eq.1) then
       qhmass=qhmass+(lspecies(i)-diskpart)
        else
       qhmass=qhmass+lspecies(i)*2*(i-1)
        endif
      enddo
      phmass=halomass1/qhmass
         xmean=0.
         ymean=0.
         zmean=0.
      do i=1,100000
         xmean=xmean+x(i)
         ymean=ymean+y(i)
         zmean=zmean+z(i)
      enddo
         xmean=xmean*1000./100000
         ymean=ymean*1000./100000
         zmean=zmean*1000./100000
       write(*,*)xmean/hubble,ymean/hubble,zmean/hubble
c      do i=diskpart,diskpart+100000
c         xmean=xmean+x(i)
c         ymean=ymean+y(i)
c         zmean=zmean+z(i)
c      enddo
c         xmean=xmean*1000./100000
c         ymean=ymean*1000./100000
c         zmean=zmean*1000./100000
c        write(*,*)xmean/hubble,ymean/hubble,zmean/hubble


         lspecies(0)=0
         do i=1,nspecies
          lspecies1(i)=(lspecies(i)-lspecies(i-1))
         enddo
         do i=1,nspecies
           if (lspecies1(i).lt.0.0) then
           lspecies(i)=lspecies(i-1)
           else
            lspecies(i)=lspecies1(i)+lspecies(i-1)
           endif
         enddo

         N_particles =lspecies(Nspecies)   ! Total number of particles                                                                    
         ip=n_particles
         do i=1,nspecies
          write(*,*)'Particulas de la especie',i,' = ',lspecies(i),'weig
     &ht =', wspecies(i)
         enddo
       nconta=0
       do i=1,n_particles,scrip
           nconta=nconta+1
       enddo
c      write(*,*)'To ascii, 0, to tipsy 1'
c      read(*,*)nopti
      nopti=0
      If(nopti.eq.1) then 
      write(*,*)'Número de partícules',lspecies(nspecies)
      write(*,*)'Número estrelles',diskpart
      OPEN(17,FILE=head2//sntime//tailtip,STATUS='UNKNOWN')
      write(17,*)nconta,0,int(diskpart/scrip)
      write(17,*)3
      write(17,*)0
       ncont=0      
      do j=1,lspecies(1)-diskpart,scrip
         write(17,*)pdmass/(3.17e15)
       ncont=ncont+1
      enddo
       write(*,*)ncont
c      ncont1=0
      do j=2,nspecies
      do i=1,(lspecies(j)-lspecies(j-1)),scrip
c         ncont1=ncont1+1
c         if(ncont1.lt.(n_particles-diskpart)) then
         write(17,*)pdmass*(2**(j-1))/(3.17e15)
         ncont=ncont+1
c         else
c         go to 1223
c         endif
         enddo
      enddo
c 1223 continue
      write(*,*)ncont
      do i=1,diskpart,scrip
         write(17,*)pdmass/(3.17e15)
         ncont=ncont+1
      enddo
         write(*,*)ncont
      if (ncont.lt.int(ip/scrip)) then
         do i=1,int((ip/scrip-ncont))
            write(17,*)pdmass*(2**(nspecies-1))/(3.17e15)
          enddo
      endif
      write(*,*)'Número de partícules',ncont
c      enddo
c        write(*,*)ncont
      
      do i=1+diskpart,ip,scrip
c         x(i)=(x(i)*1000.-xmean)/hubble
c         y(i)=(y(i)*1000.-ymean)/hubble
c         z(i)=(z(i)*1000.-zmean)/hubble
      
         write(17,*)(x(i)*1000.-xmean)/hubble/2.87e4
c           xmax =MAX(xmax,X(i),Y(i),Z(i))
c           xmin =MIN(xmin,X(i),Y(i),Z(i))
           ncont=ncont+1
      enddo
         write(*,*)ncont
      do i=1,diskpart,scrip
         write(17,*)(x(i)*1000.-xmean)/hubble/2.87e4
         ncont=ncont+1
      enddo
         write(*,*)ncont
      do i=1+diskpart,ip,scrip
         write(17,*)(y(i)*1000.-ymean)/hubble/2.87e4
         ncont=ncont+1
      enddo
      write(*,*)ncont
      do i=1,diskpart,scrip
         write(17,*)(y(i)*1000.-ymean)/hubble/2.87e4
         ncont=ncont+1
      enddo
      write(*,*)ncont
      do i=1+diskpart,ip,scrip
         write(17,*)(z(i)*1000.-zmean)/hubble/2.87e4
         ncont=ncont+1
      enddo
       write(*,*)ncont
      do i=1,diskpart,scrip
          write(17,*)(z(i)*1000.-zmean)/hubble/2.87e4
           ncont=ncont+1
      enddo
      write(*,*)ncont
      do i=1+diskpart,ip,scrip
         write(17,*)vxx(i)/691.0d0
           vmax =MAX(vmax,Vxx(i),vyy(i),vzz(i))
           vmin =MIN(vmin,Vxx(i),Vyy(i),Vzz(i))
       ncont=ncont+1
      enddo
      write(*,*)ncont
      do i=1,diskpart,scrip
          write(17,*)vxx(i)/691.0d0
           ncont=ncont+1
      enddo
      write(*,*)ncont
      do i=1+diskpart,ip,scrip
         write(17,*)vyy(i)/691.0d0
      ncont=ncont+1
      enddo
      write(*,*)ncont
      do i=1,diskpart,scrip
          write(17,*)vyy(i)/691.0d0
           ncont=ncont+1
      enddo
      write(*,*)ncont
      do i=1+diskpart,ip,scrip
         write(17,*)vzz(i)/691.0d0
       ncont=ncont+1
      enddo
       write(*,*)ncont
      do i=1,diskpart,scrip
          write(17,*)vzz(i)/691.0d0
           ncont=ncont+1
      enddo
      write(*,*)ncont
      close(17)
      else
      OPEN(17,FILE=head2//sntime//tailasc,STATUS='UNKNOWN')
      do i=1,diskpart
         x(i)=(x(i)*1000.-xmean)/hubble
         y(i)=(y(i)*1000.-ymean)/hubble
         z(i)=(z(i)*1000.-zmean)/hubble
         write(17,101)x(i),y(i),z(i),vxx(i),vyy(i),vzz(i),
     &pdmass,1
      enddo
      do i=diskpart+1,lspecies(1)
         x(i)=(x(i)*1000.-xmean)/hubble
         y(i)=(y(i)*1000.-ymean)/hubble
         z(i)=(z(i)*1000.-zmean)/hubble
         write(17,101)x(i),y(i),z(i),vxx(i),vyy(i),vzz(i),
     &pdmass,2
      enddo
      do j=2,nspecies
         do i=lspecies(j-1)+1,lspecies(j)
         x(i)=(x(i)*1000.-xmean)/hubble
         y(i)=(y(i)*1000.-ymean)/hubble
         z(i)=(z(i)*1000.-zmean)/hubble
      write(17,101)x(i),y(i),z(i),vxx(i),vyy(i),vzz(i),pdmass
     &*(2**(j-1)),2
         enddo
      enddo
      close(17)
      endif

      write (*,*) ' Scaled Coordinates were written to:',
     &            FileASCII
      write (*,*) ' Min/Max of coordinates(Mpc/h)=',xmin,xmax
      write (*,*) ' Min/Max of velocities (km/s) =',vmin,vmax
101     FORMAT(F11.6,1x,F11.6,1x,F11.6,1x,F11.6,1x,F11.6,1x,F11.6,1x,
     +F9.0,1X,I1)
      END


      subroutine GetRow ( irow , ifile )
c     ----------------------------------
c
c     input:  irow - number of record
c             ifile - number of i/o unit (0-99)
c     nrecl - length of ROW block in words
c
      include 'PMparART.h'
      read  ( ifile , rec = irow ) recdat
      return
      end


