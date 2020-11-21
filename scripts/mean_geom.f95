PROGRAM GEOM_MEAN

      ! DECLARATION
      IMPLICIT NONE
      REAL, DIMENSION(10) :: X = (/ 3, 6, 2, 1, 8, 9, 5, 4, 1, 7 /)
      REAL :: R = 1
      ! CALCULUS
      INTEGER :: I
      DO I = 1, 10
            R = R * X(I)
      END DO

      WRITE(*,*) "PRODUCT IS ", R

      R = R ** (1.0/10.0)

      WRITE(*,*) "GEOMETRIC MEAN IS", R

      DO I = 1, 10
            R = R + X(I)
      END DO

      WRITE(*,*) "THE SUM IS ", R

      WRITE(*,*) "THE MEAN IS ", R/10.0

END PROGRAM GEOM_MEAN

