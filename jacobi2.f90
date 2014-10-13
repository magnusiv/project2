PROGRAM jacobi2

IMPLICIT NONE

INTEGER, PARAMETER :: dp = KIND(1.0D0)
INTEGER, PARAMETER :: n = 80
REAL(dp),PARAMETER :: rm = 4.0_dp
INTEGER, PARAMETER :: max_rot = n**3
REAL(dp),PARAMETER :: step = rm/FLOAT(n)
REAL(dp),DIMENSION(n,n) :: a
REAL(dp),DIMENSION(n) :: d
REAL(dp),PARAMETER :: threshold = 10.0_dp**(-6.0_dp)
INTEGER :: i,p,q,j,l,rotations
REAL(dp) :: sum,offdiag1,offdiag2,t,s,c,tau,gamma,omega


OPEN(UNIT=11,FILE='sol2.dat',STATUS='new')
OPEN(UNIT=12,FILE='omegas.dat',STATUS='new')

DO l =1,2*n
   IF (l .LT. 1.5_dp*n) THEN
      omega = LOG(1.0_dp+FLOAT(l)/100.0_dp )
   ELSE
      omega = LOG(2.2_dp)+FLOAT(i)/10.0_dp-(3.0_dp*l)/200.0_dp
   END IF
   WRITE(12,*) 1.0_dp/omega

   DO i = 1, n ! Zeroing
      d(i) = 0.0_dp
      DO j = 1, n
         a(i,j) = 0.0_dp
      END DO
   END DO

   DO i = 1, n ! Initialize arrays
      IF (i .LT. n) THEN
         a(i,i+1) = -1.0_dp/step**2.0_dp
         a(i+1,i) = a(i,i+1)
         a(i,i)   = 2.0_dp/step**2.0_dp + (i*omega*step)**2.0_dp +1.0_dp/(i*step)
      ELSE
         a(i,i)   = 2.0_dp/step**2.0_dp + (i*omega*step)**2.0_dp +1.0_dp/(i*step)
      END IF
      d(i)   = a(i,i)
   END DO
   
   rotations = 0
   DO WHILE (rotations .LT. max_rot)

      sum = 0.0_dp
      DO p = 1,n-1 ! sum off-diag
         DO q = p+1,n
            sum = sum + ABS(a(p,q))**2.0_dp
         END DO
      END DO
      sum = SQRT(2.0_dp*sum)
      !WRITE(*,*) sum
   
      IF (sum .LT. threshold) THEN
         rotations = max_rot
         WRITE(11,*) MINVAL(d)
      END IF

      sum = 0.0_dp
      DO i = 1,n
         DO j = i+1,n
            offdiag1 = ABS(a(i,j))
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

         rotations = rotations + 1


      END DO
END DO

CLOSE(11)
CLOSE(12)

END PROGRAM jacobi2
