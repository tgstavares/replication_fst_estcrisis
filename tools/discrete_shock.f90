module discrete_shock
  use prec
  use normal
  use constants
contains
  
  subroutine tauchen86(rho,sig,rl,nsd,p,g)
    ! Tauchen 1986 method to discretize AR(1) processes
    !       y(t) = rho y(t-1) + sigma e(t),   e(t)~N(0,1)
    ! Inputs:
    !  rho - auto-correlation coefficient
    !  sig - standard deviation of the inovation
    !  rl - dimension of discretization .ge. 2
    !  nsd - state space with nsd uncond stdev dist around the mean
    ! Outputs:
    !  p - Rl*Rl Markov chain
    !  g - state space
    implicit none
    integer,intent(in)                      :: rl
    real(dp),intent(in)                      :: rho,sig,nsd
    real(dp),intent(out),dimension(1:Rl,1:Rl):: p
    real(dp),intent(out),dimension(1:Rl)     :: g
    integer                                 :: i,j
    real(dp)                                 :: x, y, sige, stp
    sige = sig/sqrt(1.0d0-rho**2)
    g(1) = -nsd*sige
    g(rl) =  nsd*sige
    stp = nsd*2.0d0*sige/(real(rl,dp)-1.0d0)
    do i=2,rl
       g(i) = g(i-1) + stp
    end do
    

    do i=1,Rl
       call normal_01_cdf((g(1)-rho*g(i)+stp*0.5d0)/sig,x); p(i,1) = x
       p(i,Rl) = 1 - p(i,1)
       do j=2,Rl-1
          call normal_01_cdf((g(j)-rho*g(i)+stp*0.5d0)/sig,x);
          call normal_01_cdf((g(j)-rho*g(i)-stp*0.5d0)/sig,y);
          p(i,j) =  x - y
          p(i,Rl) = p(i,Rl) - p(i,j)
       enddo
    enddo
  end subroutine tauchen86

  subroutine tauchen_hussey91(mu,rho,sig,p,z)
    ! Tauchen and Hussey 1991 quadrature rule to discretize AR(1) processes
    !      y(t) = (1-rho)mu + rho y(t-1) + sig e(t),   e(t)~N(0,1)
    ! Inputs:
    !  mu - unconditional mean
    !  rho - auto-correlation coefficient
    !  sig - standard deviation of inovation
    ! Outputs
    !  p - transition matrix
    !  z - state space
    implicit none
    real(dp), intent(in) :: mu,rho,sig
    real(dp), dimension(:,:), intent(out)::p
    real(dp), dimension(:), intent(out)::z
    integer  j,i,n,n_file(325) ! 1 + 2 + ... + 25
    real(dp) x_file(325),w_file(325)
    real(dp), allocatable::x(:),w(:)
    logical existe
    n = size(z);
    if(n.gt.25) stop 'Error at tauchen_hussey91 - reduce n'
    inquire(file='tools/hermite.dat',exist=existe); if(existe.eqv..false.) stop 'Error at tauchen_hussey91 - no hermite.dat'
    open(1,file='tools/hermite.dat')
    do i=1,325
       read(1,'(I5,2(F20.14))')n_file(i),x_file(i),w_file(i)
    end do
    allocate (x(n)); allocate (w(n))
    j = 1
    do i=1,325
       if(n_file(i).eq.n) then
          x(j) = x_file(i)
          w(j) = w_file(i)
          j = j + 1
       end if
    end do
    do j=1,n
       z(j) = sqrt(2.0d0)*sig*x(j) + mu
    end do
    do i=1,n
       do j=1,n
          p(i,j) = w(j) * (1.0d0/sqrt(pi)) * &
               exp(-(z(j)-((1-rho)*mu+rho*z(i)))**2/(2.0d0*sig**2)) / &
               exp(-(z(j)-mu)**2/(2.0d0*sig**2))
       end do
       p(i,:) = p(i,:)/sum(p(i,:))
    end do
  end subroutine tauchen_hussey91

    
  SUBROUTINE rouwenhorst95(rho,sig,p,y)
    ! Rowenhurst method to approximate univariate AR(1) process by Markov chain
    !      y(t) = rho y(t-1) + sig e(t),   e(t)~N(0,1)
    ! INPUTS: rho - serial correlation coefficient,
    !         sigma - coefficient of variation
    ! OUTPUT: P is an n-by-n matrix of Markov transition probabilities
    !         y is an n-by-1 vector of symmetric and evenly-spaced Markov state space
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: rho, sig 
    REAL(dp), DIMENSION(:,:), INTENT(OUT) :: p
    REAL(dp), DIMENSION(:), INTENT(OUT) :: y
    REAL(dp) :: ybar,q,sigma
    INTEGER :: n
    sigma = sig / sqrt(1 - rho**2)
    n=size(y)
    IF (size(p,dim=1)/=n .or. size(p,dim=2)/=n) THEN
       PRINT '(a,i3,a,i3)', 'rouwenhorst: p must be a square matrix of size ',n,' x ',n
       STOP 'program terminated by rouwenhorst'
    END IF
    ybar=sigma*sqrt(real(n-1,dp))
    q=(1+rho)/2
    CALL rhmat(p)
    CALL grids(y,-ybar,ybar,1.0_dp)
  CONTAINS
    RECURSIVE SUBROUTINE rhmat(p)
      IMPLICIT NONE
      REAL(dp), DIMENSION(:,:), INTENT(OUT) :: p
      REAL(dp), DIMENSION(size(p,dim=1)-1,size(p,dim=2)-1) :: p1
      INTEGER :: h
      h=size(p,dim=1)
      IF (size(p,dim=2)/=h) STOP 'P must be a square matrix'
      IF (h<2) STOP 'P must be at least 2-by-2 matrix'
      IF (h==2) THEN
         p=reshape((/q,1-q,1-q,q/),(/2,2/))
      ELSE
         CALL rhmat(p1)
         p=0
         p(1:h-1,1:h-1)=q*p1
         p(1:h-1,2:h)=(1-q)*p1+p(1:h-1,2:h)
         p(2:h,1:h-1)=(1-q)*p1+p(2:h,1:h-1)
         p(2:h,2:h)=q*p1+p(2:h,2:h)
         p(2:h-1,:)=p(2:h-1,:)/2
      END IF
    END SUBROUTINE rhmat
  END SUBROUTINE rouwenhorst95

  SUBROUTINE grids(x,xmin,xmax,s)
    ! Purpose: Generate grid x on [xmin,xmax] using spacing parameter s set as follows:
    ! s=1linear spacing
    ! s>1left skewed grid spacing with power s
    ! 0<s<1right skewed grid spacing with power s
    ! s<0geometric spacing with distances changing by a factor -s^(1/(n-1)), (>1 grow, <1 shrink)
    ! s=-1logarithmic spacing with distances changing by a factor (xmax-xmin+1)^(1/(n-1))
    ! s=0logarithmic spacing with distances changing by a factor (xmax/xmin)^(1/(n-1)), only if xmax,xmin>0
    IMPLICIT NONE
    REAL(dp), DIMENSION(:), INTENT(OUT) :: x
    REAL(dp), INTENT(IN) :: xmin,xmax,s
    REAL(dp) :: c ! growth rate of grid subintervals for logarithmic spacing
    INTEGER :: n,i
    n=size(x)
    FORALL(i=1:n) x(i)=(i-1)/real(n-1,dp)
    IF (s>0.0_dp) THEN
       x=x**s*(xmax-xmin)+xmin
       IF (s==1.0_dp) THEN
          !PRINT '(a,i8,a,f6.3,a,f6.3,a)', 'Using ',n,' equally spaced grid points over domain [',xmin,',',xmax,']'
       ELSE
          !PRINT '(a,i8,a,f6.3,a,f6.3,a,f6.3,a)', 'Using ',n,' skewed spaced grid points with power ',s,' over domain [',xmin,',',xmax,']'
       END IF
    ELSE
       IF (s==-1.0_dp) THEN
          c=xmax-xmin+1
          !ELSEIF (s==0.0_dp) THEN
          !IF (xmin>0.0_dp) THEN
          !c=xmax/xmin
          !ELSE
          !STOP 'grid: can not use logarithmic spacing for nonpositive values'
          !END IF
       ELSE
          c=-s
       END IF
       !PRINT '(a,i8,a,f6.3,a,f6.3,a,f6.3,a)', 'Using ',n,' logarithmically spaced grid points with growth rate ',c,' over domain [',xmin,',',xmax,']'
       x=((xmax-xmin)/(c-1))*(c**x)-((xmax-c*xmin)/(c-1))
    END IF
  END SUBROUTINE grids

end module discrete_shock
