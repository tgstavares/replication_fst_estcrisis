module normal
  use prec
contains

  FUNCTION rnorm() RESULT( fn_val )

    !   Generate a random normal deviate using the polar method.
    !   Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
    !              normal variables', Siam Rev., vol.6, 260-264, 1964.

    IMPLICIT NONE
    REAL  :: fn_val

    ! Local variables

    REAL            :: u, sum
    REAL, SAVE      :: v, sln
    LOGICAL, SAVE   :: second = .FALSE.
    REAL, PARAMETER :: one = 1.0, vsmall = TINY( one )

    IF (second) THEN
       ! If second, use the second random number generated on last call

       second = .false.
       fn_val = v*sln

    ELSE
       ! First call; generate a pair of random normals

       second = .true.
       DO
          CALL RANDOM_NUMBER( u )
          CALL RANDOM_NUMBER( v )
          u = SCALE( u, 1 ) - one
          v = SCALE( v, 1 ) - one
          sum = u*u + v*v + vsmall         ! vsmall added to prevent LOG(zero) / zero
          IF(sum < one) EXIT
       END DO
       sln = SQRT(- SCALE( LOG(sum), 1 ) / sum)
       fn_val = u*sln
    END IF

    RETURN
  END FUNCTION rnorm




  !!Subroutine for standard normal cdf
  subroutine normal_01_cdf ( x, cdf )
    !
    !*******************************************************************************
    !
    !! NORMAL_01_CDF evaluates the Normal 01 CDF.
    !
    !
    !  Reference:
    !
    !    A G Adams,
    !    Areas Under the Normal Curve,
    !    Algorithm 39,
    !    Computer j.,
    !    Volume 12, pages 197-198, 1969.
    !
    !  Modified:
    !
    !    10 February 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real X, the argument of the CDF.
    !
    !    Output, real CDF, the value of the CDF.
    !
    implicit none
    !
    real(dp), parameter :: a1 = 0.398942280444E+00
    real(dp), parameter :: a2 = 0.399903438504E+00
    real(dp), parameter :: a3 = 5.75885480458E+00
    real(dp), parameter :: a4 = 29.8213557808E+00
    real(dp), parameter :: a5 = 2.62433121679E+00
    real(dp), parameter :: a6 = 48.6959930692E+00
    real(dp), parameter :: a7 = 5.92885724438E+00
    real(dp), parameter :: b0 = 0.398942280385E+00
    real(dp), parameter :: b1 = 3.8052E-08
    real(dp), parameter :: b2 = 1.00000615302E+00
    real(dp), parameter :: b3 = 3.98064794E-04
    real(dp), parameter :: b4 = 1.98615381364E+00
    real(dp), parameter :: b5 = 0.151679116635E+00
    real(dp), parameter :: b6 = 5.29330324926E+00
    real(dp), parameter :: b7 = 4.8385912808E+00
    real(dp), parameter :: b8 = 15.1508972451E+00
    real(dp), parameter :: b9 = 0.742380924027E+00
    real(dp), parameter :: b10 = 30.789933034E+00
    real(dp), parameter :: b11 = 3.99019417011E+00
    real(dp) cdf
    real(dp) q
    real(dp) x
    real(dp) y
    !
    !  |X| <= 1.28.
    !
    if ( abs ( x ) <= 1.28 ) then

       y = 0.5E+00 * x**2

       q = 0.5E+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
            + a6 / ( y + a7 ) ) ) )
       !
       !  1.28 < |X| <= 12.7
       !
    else if ( abs ( x ) <= 12.7E+00 ) then

       y = 0.5E+00 * x**2

       q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
            + b2 / ( abs ( x ) + b3 &
            + b4 / ( abs ( x ) - b5 &
            + b6 / ( abs ( x ) + b7 &
            - b8 / ( abs ( x ) + b9 &
            + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
       !
       !  12.7 < |X|
       !
    else

       q = 0.0E+00

    end if
    !
    !  Take account of negative X.
    !
    if ( x < 0.0E+00 ) then
       cdf = q
    else
       cdf = 1.0E+00 - q
    end if

    return
  end subroutine normal_01_cdf

  subroutine normal_01_pdf ( x, pdf )

    !*****************************************************************************80
    !
    !! NORMAL_01_PDF evaluates the Normal 01 PDF.
    !
    !  Discussion:
    !
    !    The Normal 01 PDF is also called the "Standard Normal" PDF, or
    !    the Normal PDF with 0 mean and standard deviation 1.
    !
    !    PDF(X) = exp ( - 0.5 * X^2 ) / sqrt ( 2 * PI )
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    04 December 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument of the PDF.
    !
    !    Output, real ( kind = 8 ) PDF, the value of the PDF.
    !
    implicit none

    real ( kind = 8 ) pdf
    real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
    real ( kind = 8 ) x

    pdf = exp ( -0.5D+00 * x * x ) / sqrt ( 2.0D+00 * r8_pi )

    return
  end subroutine normal_01_pdf

end module normal
