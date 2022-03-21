module linear_op
  ! for now: 	matrix mult .x.
  ! 		matrix trans tr()
  !		matrix inverser inv()
  interface operator (.x.)
     module procedure matrix_product
  end interface operator (.x.)

contains
  function matrix_product(a,b) result (w)
    use prec
    implicit none
    real(dp), dimension(:,:), intent(in) :: a
    real(dp), dimension(:,:), intent(in) :: b
    real(dp), dimension(size(a,1),size(b,2)) :: w
    w = matmul(a,b)
  end function matrix_product

  function tr(A)
    use prec
    real(dp), dimension(:,:), intent(in):: A
    real(dp), dimension(size(A,2),size(A,1)) :: tr
    tr = transpose(A)
  end function tr
  
  ! function inv(A)
  !   ! Returns the inverse of a matrix calculated by finding the LU
  !   ! decomposition.  Depends on LAPACK.
  !   use prec
  !   real(dp), dimension(:,:), intent(in) :: A
  !   real(dp), dimension(size(A,1),size(A,2)) :: inv
  !   real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
  !   integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  !   integer :: n, info

  !   ! External procedures defined in LAPACK
  !   external DGETRF
  !   external DGETRI

  !   ! Store A in inv to prevent it from being overwritten by LAPACK
  !   inv = A
  !   n = size(A,1)

  !   ! DGETRF computes an LU factorization of a general M-by-N matrix A
  !   ! using partial pivoting with row interchanges.
  !   call DGETRF(n, n, inv, n, ipiv, info)

  !   if (info /= 0) then
  !      stop 'Matrix is numerically singular!'
  !   end if

  !   ! DGETRI computes the inverse of a matrix using the LU factorization
  !   ! computed by DGETRF.
  !   call DGETRI(n, inv, n, ipiv, work, n, info)

  !   if (info /= 0) then
  !      stop 'Matrix inversion failed!'
  !   end if
  ! end function inv

end module linear_op
