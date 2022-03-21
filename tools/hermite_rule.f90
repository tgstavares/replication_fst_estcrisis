program hermite_rule
  use gqrul_int
  implicit none
  integer, parameter::dp=selected_real_kind(15,100)
  integer n,i
  real(dp), allocatable:: x(:),w(:)

  open(1,file='hermite.dat')
  do n=1,25
     allocate (x(n))
     allocate (w(n))
     call d_gqrul(n,x,w,iweigh=4) ! Hermite
     w = w / sum(w)
     do i=1,n
        write(1,'(I5,2(F20.14))')n,x(i),w(i)
     end do
     deallocate(x)
     deallocate(w)  
  end do
  close(1)
end program hermite_rule
  
