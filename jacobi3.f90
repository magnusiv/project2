PROGRAM jacobi3

IMPLICIT NONE

INTEGER, PARAMETER :: dp = KIND(1.0D0)
INTEGER, PARAMETER :: n = 1000
REAL(dp),PARAMETER :: rm = 27.0_dp
INTEGER, PARAMETER :: max_rot = n**3
REAL(dp),PARAMETER :: step = rm/FLOAT(n)
REAL(dp),DIMENSION(n,n) :: a,E
REAL(dp),DIMENSION(n) :: d
REAL(dp),DIMENSION(4) :: omega
REAL(dp),PARAMETER :: threshold = 10.0_dp**(-6.0_dp)
INTEGER :: i,p,q,j,l,rotations
REAL(dp) :: sum,offdiag1,offdiag2,t,s,c,tau,gamma
CHARACTER(len=1024) :: FMT = '(I2.2)'
CHARACTER(len=1024) :: stri

omega = [0.01_dp,0.5_dp,1.0_dp,5.0_dp]

DO l =1,4
   DO i = 1, n ! Zeroing
      d(i) = 0.0_dp
      DO j = 1, n
         a(i,j) = 0.0_dp
         E(i,j) = 0.0_dp
      END DO
   END DO

   DO i = 1, n ! Initialize arrays
      IF (i .LT. n) THEN
         a(i,i+1) = -1.0_dp/step**2.0_dp
         a(i+1,i) = a(i,i+1)
         a(i,i)   = 2.0_dp/step**2.0_dp + (i*omega(l)*step)**2.0_dp +1.0_dp/(i*step)
      ELSE
         a(i,i)   = 2.0_dp/step**2.0_dp + (i*omega(l)*step)**2.0_dp +1.0_dp/(i*step)
      END IF
      d(i)   = a(i,i)
      E(i,i) = 1.0_dp
   END DO
   
   rotations = 0
   DO WHILE (rotations .LT. max_rot) ! Rotations loop
      
      sum = 0.0_dp
      DO p = 1,n-1 ! sum off-diag
         DO q = p+1,n
            sum = sum + ABS(a(p,q))**2.0_dp
         END DO
      END DO
      sum = SQRT(2.0_dp*sum)
      !WRITE(*,*) sum
   
      IF (sum .LT. threshold) THEN
         WRITE(stri,FMT) 10+l
         OPEN(UNIT=10+l,FILE='sol'//TRIM(stri)//'.dat',STATUS='new')
         rotations = MINLOC(d, DIM=1) ! dummy index, location of minimum
         DO j = 1,n
            WRITE(10+l,*) E(j,rotations)
         END DO
         CLOSE(10+l)
         rotations = max_rot
      END IF

      sum = 0.0_dp
      DO i = 1,n ! Finding max element
         DO j = i+1,n
            offdiag1 = ABS(a(i,j)) ! Element
            IF (offdiag1 .GT. sum) THEN
               sum = offdiag1
               p = i
               q = j
            END IF
         END DO
      END DO

         tau      = ( a(q,q)-a(p,p) )/( 2.0_dp*a(p,q) )
         t        = 1.0_dp/( ABS(tau)+SQRT(tau**2.0_dp+1) )
         IF (tau .LT. 0.0_dp) t = -t
         c        = 1.0_dp/SQRT(1.0_dp+t**2.0_dp)
         s        = t*c
         gamma    = s/(1.0_dp+c)
         offdiag1 = t*a(p,q)
         a(p,p)   = a(p,p) - offdiag1
         d(p)     = a(p,p)
         a(q,q)   = a(q,q) + offdiag1
         d(q)     = a(q,q)
         a(p,q)   = 0.0_dp
         ! Main rotation algorithm:
         DO j=1,p-1
            offdiag1 = a(j,p)
            offdiag2 = a(j,q)
            a(j,p)   = c*offdiag1 - s*offdiag2
            a(j,q)   = c*offdiag2 + s*offdiag1
         END DO
         DO j=p+1,q-1
            offdiag1 = a(p,j)
            offdiag2 = a(j,q)
            a(p,j)   = c*offdiag1 - s*offdiag2
            a(j,q)   = c*offdiag2 + s*offdiag1
         END DO
         DO j=q+1,n
            offdiag1 = a(p,j)
            offdiag2 = a(q,j)
            a(p,j)   = c*offdiag1 - s*offdiag2
            a(q,j)   = c*offdiag2 + s*offdiag1
         END DO
         DO j=1,n
            offdiag1 = E(j,p)
            offdiag2 = E(j,q)
            E(j,p)   = c*offdiag1 - s*offdiag2
            E(j,q)   = c*offdiag2 + s*offdiag1
         END DO
         rotations = rotations + 1
      END DO
END DO

END PROGRAM jacobi3
