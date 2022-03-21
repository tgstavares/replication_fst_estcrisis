module quadratic_equation
  use prec
contains
  subroutine fresolvente(a,b,c,x1,x2)
    !
    ! Solves:
    !    ax**2 + bx + c = 0
    ! => (x - x1)*(x - x2) = 0
    !
    ! erro = 1 if complex
    ! erro = 0 if real
    !
    implicit none
    real(dp), intent(in) ::a,b,c
    real(dp), intent(out)::x1,x2
    integer erro

    if(b**2 - 4.0d0*a*c.lt.0.0d0)then
       erro = 1
       print*,'Error: complex root in fresolvent'
    else
       x1 = (-b + sqrt(b**2 - 4.0d0*a*c))/(2.0d0*a)
       x2 = (-b - sqrt(b**2 - 4.0d0*a*c))/(2.0d0*a)
       erro = 0
    end if
  end subroutine fresolvente
end module quadratic_equation
