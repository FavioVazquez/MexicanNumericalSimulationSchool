PROGRAM pi_mc

   IMPLICIT NONE

   INTEGER, PARAMETER :: num_int = 1000000


   INTEGER :: i, ncirc=0;
   REAL(8):: pi=0, x, y
   REAL(8):: r = 1.0
   REAL(8):: intentos

!$OMP PARALLEL DO PRIVATE(x,y) REDUCTION(+:ncirc)

   DO i=1,num_int
      x = rand()
      y = rand()
      IF( (x*x + y*y) <= r*r) THEN
         ncirc = ncirc + 1
      END IF
   END DO
!$OMP END PARALLEL DO

   intentos = num_int
   pi = 4.0 * (ncirc/intentos)

   PRINT*,pi

END PROGRAM pi_mc
