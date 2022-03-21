module uteis
  use prec
contains
  real(dp) function indic(true)
    implicit none
    logical, intent(in)::true
    indic = 0.0d0
    if(true.eqv..true.) indic = 1.0d0
  end function indic
  
  !!Subroutine to create equaly-spaced grids
  subroutine grid(ng,lb,ub,vector)
    ! lb - lower bound
    ! ub - upper bound
    ! ng - number of grids
    ! vector - output vector
    implicit none
    integer, intent(in)                     :: ng
    real(dp), intent(in)                     :: lb, ub
    real(dp), intent(out),dimension(1:ng)    :: vector
    integer                                 :: i
    if(ng.eq.1) then
       vector(1) = (ub+lb)*0.5d0
    else
       vector(1) = lb
       do i=2,ng
          vector(i) = vector(i-1) + (ub-lb)/(real(ng,dp)-1.0d0)
       enddo
       vector(ng)= ub
    end if
  end subroutine grid

  subroutine position(x,v,N,yloc1,yloc2)
    ! x (in) - value to be located
    ! v (in) - vector
    ! N (in) - vector dimension
    ! yloc1 (out) - location with v(yloc(1)) <= x
    ! yloc2 (out) - location with x <= v(yloc(2))
    implicit none
    integer,intent(in)                  ::N
    real(dp),intent(in),dimension(1:N)   ::v
    real(dp),intent(in)                  ::x
    integer,intent(out)                 ::yloc1
    integer,intent(out)                 ::yloc2
    integer                             ::i
    if(x<=v(1))then
       yloc1 = 1
       yloc2 = 1
    elseif(x>=v(N))then
       yloc1 = N
       yloc2 = N
    else
       i=1;do while(v(i)<=x);i=i+1;enddo
       if(v(i-1)==x)then
          yloc1 = i-1
          yloc2 = i-1
       else
          yloc1 = i-1
          yloc2 = i
       endif
    endif
  end subroutine position

  integer function posit(n,xv,x)
    implicit none
    integer n,res(1)
    real(dp) xv(n),x
    res=minloc(abs(xv-x))
    posit=res(1)
  end function posit

end module uteis
