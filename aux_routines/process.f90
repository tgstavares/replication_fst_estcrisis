module process
  use prec
contains

  subroutine get_eqspace(n,min,max,x)
    implicit none
    integer, intent(in)::n
    real(dp), intent(in)::min,max
    real(dp), intent(out)::x(n)
    integer i
    real(dp) step
    if(n.eq.1)then
       x(1) = (max+min)/2d0       
    else
       step = (max-min)/dble(n-1)
       x(1) = min
       do i=2,n-1
          x(i) = x(i-1) + step
       end do
       x(n) = max
    end if
  end subroutine get_eqspace

  subroutine get_markov(nv,min,max,rho,sigma,prbv,v)
    use normal
    implicit none
    integer, intent(in)::nv
    real(dp), intent(in)::min,max,rho,sigma
    real(dp), intent(out)::prbv(nv,nv),v(nv)
    integer i,j
    real(dp) step,x,y
    if(nv.eq.1)then
       v(1) = (min+max)/2d0
       prbv(1,1) = 1d0
    else
       ! GRID
       call get_eqspace(nv,min,max,v)
       ! TRANSITION MATRIX
       step = v(2)-v(1)
       do i=1,nv
          call normal_01_cdf( (v(1)-rho*v(i) + step*0.5d0) / sigma,x); prbv(i,1) = x
          do j=2,nv-1
             call normal_01_cdf((v(j)-rho*v(i) + step*0.5d0) / sigma,x)
             call normal_01_cdf((v(j)-rho*v(i) - step*0.5d0) / sigma,y)
             prbv(i,j) = x-y
          end do
          call normal_01_cdf( (v(nv)-rho*v(i) - step*0.5d0) / sigma,y); prbv(i,nv) = 1-y
          prbv(i,:) = prbv(i,:) / sum(prbv(i,:))
       end do
    end if
  end subroutine get_markov

  subroutine get_markovforts(n,mu1,mu2,sigma1,sigma2,rho1,rho2,zmin,zmax,prb1,prb2,t)
    use normal
    implicit none
    integer, intent(in)::n
    real(dp), intent(in)::mu1,mu2,sigma1,sigma2,rho1,rho2,zmin,zmax	! mu1,mu2 ARE UNCONDITIONAL MEANS
    real(dp), intent(out)::prb1(n,n),prb2(n,n),t(n)
    real(dp) stepf,stepb,cdff,cdfb
    integer i,j
    ! CONSTRUCT GRID
    if(n.le.2) then
       stop 'ERROR - in the choice of markov chains; put something larger than 2'
    elseif(mod(n,2).eq.0)then
       stop 'ERROR - in the choice of markov chains; put something that is odd'
    else
       ! FIND BOUNDS
       t(1) = zmin
       t(n) = zmax
       t(n/2+1) = 0.5d0*(mu1+mu2)
       do i=1,n/2
          t(n/2+1 + i) = t(n/2+1) + y( (dble(i)/dble(n/2)) * x(abs(t(n) - t(n/2+1))) )
          t(n/2+1 - i) = t(n/2+1) - y( (dble(i)/dble(n/2)) * x(abs(t(1) - t(n/2+1))) )
       end do       
    end if
    ! GET PROBABILITIES
    do i=1,n
       stepf = t(2) - t(1)
       call normal_01_cdf(((t(1)+stepf*0.5d0) - ((1d0-rho1)*mu1 + rho1*t(i)))/sigma1,cdff)
       prb1(i,1) = cdff
       do j=2,n-1
          stepf = t(j+1) - t(j)
          stepb = t(j) - t(j-1)
          call normal_01_cdf(((t(j)+stepf*0.5d0) - ((1d0-rho1)*mu1 + rho1*t(i)))/sigma1,cdff)
          call normal_01_cdf(((t(j)-stepb*0.5d0) - ((1d0-rho1)*mu1 + rho1*t(i)))/sigma1,cdfb)
          prb1(i,j) = cdff-cdfb
       end do
       stepb = t(n) - t(n-1)
       call normal_01_cdf(((t(n)-stepb*0.5d0) - ((1d0-rho1)*mu1 + rho1*t(i)))/sigma1,cdfb)
       prb1(i,n) = 1d0 - cdfb
       stepf = t(2) - t(1)
       call normal_01_cdf(((t(1)+stepf*0.5d0) - ((1d0-rho2)*mu2 + rho2*t(i)))/sigma2,cdff)
       prb2(i,1) = cdff
       do j=2,n-1
          stepf = t(j+1) - t(j)
          stepb = t(j) - t(j-1)
          call normal_01_cdf(((t(j)+stepf*0.5d0) - ((1d0-rho2)*mu2 + rho2*t(i)))/sigma2,cdff)
          call normal_01_cdf(((t(j)-stepb*0.5d0) - ((1d0-rho2)*mu2 + rho2*t(i)))/sigma2,cdfb)
          prb2(i,j) = cdff-cdfb
       end do
       stepb = t(n) - t(n-1)
       call normal_01_cdf(((t(n)-stepb*0.5d0) - ((1d0-rho2)*mu2 + rho2*t(i)))/sigma2,cdfb)
       prb2(i,n) = 1d0 - cdfb
    end do
  contains
    real(dp) function x(y)
      implicit none
      real(dp) y
      x = y
      !x = y**(1d0/1d0)
      !x = log(y)
    end function x
    real(dp) function y(x)
      implicit none
      real(dp) x
      y = x
      !y = x**1d0
      !y = exp(x)
    end function y
  end subroutine get_markovforts

  subroutine teste_markov(tsimul,na,a,prba,sigma,rho,mu)
    implicit none
    integer, intent(in)::tsimul,na
    real(dp), intent(in)::a(na),prba(na,na),sigma,rho,mu
    integer flag,s1,s3,io0,io1,t,n,ia0,ia1,i
    integer, dimension(:),allocatable::s2
    real(dp) shocks(tsimul),apath(tsimul),fcumul(na)
    real stdev,average,average2,averagee,averagee2
    call system_clock(COUNT=s3)
    call random_seed(size=s1); allocate(s2(s1));
    do i=1,s1
       s2(i) = s3 + i
    end do
    call random_seed(put=s2)
    call random_number(shocks)
    ia0 = na / 2 + 1
    do t=1,tsimul
       apath(t) = a(ia0)
       fcumul(1) = prba(ia0,1)
       if(shocks(t).le.fcumul(1))then
          ia1 = 1
       else
          do i=2,na
             fcumul(i) = fcumul(i-1) + prba(ia0,i)
             if(fcumul(i).ge.shocks(t))then
                ia1 = i
                exit
             end if
          end do
       end if
       ia0 = ia1
    end do
    call sd(real(apath(1:tsimul)),tsimul,0,stdev)
    call mean(real(apath(1:tsimul)),tsimul,0,average)
    average2 = real(sum(apath(1:tsimul))/dble(tsimul))
    call mean(real(exp(apath(1:tsimul))),tsimul,0,averagee)
    averagee2 = real(sum(exp(apath(1:tsimul)))/dble(tsimul))
    print*,"Test Markov chains: ",tsimul
    write(*,'(10a15)')'nmean','nmean','mean','nstd','std','nmeanexp','nmeanexp','meanexp'
    write(*,'(10f15.4)')average,average2,mu,stdev,sigma/sqrt(1d0-rho**2),averagee,averagee2,exp(mu+0.5d0*sigma**2/(1d0-rho**2))
  end subroutine teste_markov

  subroutine get_markov_vpath(nv,v,rho,sigma,prbv)
    use normal
    implicit none
    integer, intent(in)::nv
    real(dp), intent(in)::v(nv),rho,sigma
    real(dp), intent(out)::prbv(nv,nv)
    integer i,j
    real(dp) x,y
    ! transition matrix
    do i=1,nv
       call normal_01_cdf( (v(1)-rho*v(i)   + (v(2)-v(1)  )*0.5d0) / sigma,x); prbv(i,1) = x
       do j=2,nv-1          
          call normal_01_cdf((v(j)-rho*v(i) - (v(j)-v(j-1))*0.5d0) / sigma,y)
          call normal_01_cdf((v(j)-rho*v(i) + (v(j+1)-v(j))*0.5d0) / sigma,x)
          prbv(i,j) = x-y
       end do
       call normal_01_cdf( (v(nv)-rho*v(i) - (v(nv)-v(nv-1))*0.5d0) / sigma,y); prbv(i,nv) = 1-y
       prbv(i,:) = prbv(i,:) / sum(prbv(i,:))
    end do
  end subroutine get_markov_vpath

  subroutine get_pareto(n,slow,supp,shape,sgrid,sprob)
    implicit none
    integer, intent(in)::n
    real(dp), intent(in)::slow,supp,shape
    real(dp), intent(out)::sgrid(n),sprob(n)
    integer i
    if(n.eq.1)then
       sgrid(1) = 0.5d0*(log(supp)+log(slow))
       sgrid = exp(sgrid)
       sprob(1) = 1d0
    else
       do i=1,n
          sgrid(i) = log(slow) + dble(i-1)*(log(supp)-log(slow))/dble(n-1)
       end do
       sgrid = exp(sgrid)
       do i=1,n-1
          sprob(i) = (slow/sgrid(i))**shape-(slow/sgrid(i+1))**shape
       end do
       sprob(n) = (slow/sgrid(n))**shape-(slow/(2d0*sgrid(n)-sgrid(n-1)))**shape
       sprob = sprob/(1d0-(slow/(2d0*sgrid(n)-sgrid(n-1)))**shape)
    end if
  end subroutine get_pareto

  subroutine get_transigtion_sz(ns,nz,sigmaz,rho_z,szprob)
  end subroutine get_transigtion_sz

end module process
