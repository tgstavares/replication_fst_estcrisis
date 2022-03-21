program montecarlo
  use prec; use parameters; use constants; use process_alt; use uteis
  implicit none

  real(dp) pk,c2,sigma_e
  integer kfact
  real(dp) k(nk),o(no),prbo(no,no)

  integer, parameter:: ng=15,n=n_simul*t_remain
  integer i1,i2,i3,i1max,i2max,i3max
  real(dp) &
       grid_pk(ng),min_pk,max_pk, &
       grid_c2(ng),min_c2,max_c2, &
       grid_sigma_e(ng),min_sigma_e,max_sigma_e,le(n,2)
  real(dp) mle(ng,ng,ng),mlemax
  integer i
  real(dp) data(n,9),likelihood

  ! interpolation
  real(dp) &
       ttk(nk+kk),  tto(no+ko), &
       wwork(nk*no + max(2*kk*(nk+1),2*ko*(no+1))), &
       bbcoef(nk,no) , &
       wwork2s(3*max(kk,ko)+ko)
  real(dp) db2val
  integer iflag

  ! simplex
  real(dp) lb(3),ub(3),guess(3),sol(3),fsol,reqmin,step(3)
  integer konvge,kcount,icount,numres,ifault

  open(1,file='data/export.txt',position="rewind")
  !open(1,file='data/data_trunc.txt',position="rewind")
  do i=1,n     
     read(1,*) data(i,:)
  end do
  close(1)

  !call tauchen86(rho_o,sigma_o,no,4d0,prbo,o)
  call geto(rho_o,sigma_o,no,omin,omax,prbo,o); o = o + mu_o
  call getk(nk,kmin,kmax,delta,k,kfact)

  min_pk = 0.85d0; max_pk = 1.15d0
  min_sigma_e = 0.04d0; max_sigma_e = 0.05d0
  min_c2 = 0.05d0; max_c2 = 0.5d0

  call grid(ng,min_pk,max_pk,grid_pk)
  call grid(ng,min_sigma_e,max_sigma_e,grid_sigma_e)
  call grid(ng,min_c2,max_c2,grid_c2)

  mlemax = -infty
  do i1=ng/2+1,ng/2+1
     do i3=ng/2+1-3,ng/2+1
        do i2=1,ng
           pk = grid_pk(i1)
           c2 = grid_c2(i2)
           sigma_e = grid_sigma_e(i3)
           call get_mle(pk,c2,sigma_e,likelihood,le)
           print*,pk,sigma_e,c2,"     Log Likehood = ",likelihood

           ! open(1,file='data/le.txt',position="rewind")
           ! do i=1,n
           !    write(1,*) le(i)
           ! end do
           ! close(1)
           ! stop

           mle(i1,i2,i3) = likelihood
           if(likelihood.gt.mlemax)then
              mlemax=likelihood
              i1max = i1
              i2max = i2
              i3max = i3
           end if
        end do
     end do
  end do
  call get_mle(grid_pk(i1max),grid_c2(i2max),grid_sigma_e(i3max),likelihood,le)
  open(1,file='data/le.txt',position="rewind")
  do i=1,n
     !if(le(i,2).le.0.5d0)
     write(1,*) le(i,:)
  end do
  close(1)  
  
  print*,"Max achieved in       :",grid_pk(i1max),grid_sigma_e(i3max),grid_c2(i2max)
  print*,"With log likelihood of:",mlemax
  call get_mle(1d0,0.5d0,0.05d0,likelihood,le)
  print*,"Log likelihood at parm:",likelihood
  open(1,file='data/le_exact.txt',position="rewind")
  do i=1,n
     !if(le(i,2).le.0.5d0)
     write(1,*) le(i,:)
  end do
  close(1)  

contains
  subroutine get_mle(pk,c2,sigma_e,mle,le)
    use process_alt; use normal; 
    implicit none
    real(dp), intent(in) ::pk,c2,sigma_e
    real(dp), intent(out)::mle,le(n,2)
    integer i
    real(dp) vg(nk,no,ne),expec(nk,no),rhs,lhs,likelihood
    ! interpolation
    real(dp) &
         ttk(nk+kk),  tto(no+ko), &
         wwork(nk*no + max(2*kk*(nk+1),2*ko*(no+1))), &
         bbcoef(nk,no) , &
         wwork2s(3*max(kk,ko)+ko)
    real(dp) db2val,xx
    integer iflag
    call compvg(pk,c2,k,o,vg)
    call get_equil(pk,c2,sigma_e,k,kfact,o,prbo,vg,expec)
    iflag = 0; call db2ink(k,nk,o,no,expec,nk,kk,ko,ttk,tto,bbcoef,wwork,iflag)
    likelihood = 0d0
    do i=1,n
       rhs = beta*db2val(data(i,8),data(i,3),1,0,ttk,tto,nk,no,kk,ko,bbcoef,wwork2s)
       
       lhs = (pk+c2*(data(i,8)-data(i,6)*(1d0-delta))/data(i,6))
       le(i,1) = log(rhs) - log(lhs)
       !le(i) = rnorm()*(0.05d0)
       !print*,le(i),data(i,9),data(i,8)-data(i,6)*(1d0-delta)
       if(data(i,8)-data(i,6)*(1d0-delta).gt.1d-4)then
          call normal_01_pdf((le(i,1))/sigma_e,xx)
          likelihood = likelihood + log(xx/sigma_e)
          le(i,2) = 0d0
       else
         call normal_01_cdf((le(i,1))/sigma_e,xx)
         likelihood = likelihood + log((1d0-xx))
         le(i,2) = 1d0
       end if
    end do
    mle = likelihood    
  end subroutine get_mle

  real(dp) function fn(x)
    implicit none
    real(dp) x(3),pk,ce,sigma_e,mle
    pk = x(1)
    c2 = x(2)
    sigma_e = x(3)
    call get_mle(pk,c2,sigma_e,mle,le)
    fn = mle
  end function fn

  subroutine get_equil(pk,c2,sigma_e,k,kfact,o,prbo,vguess,expec)
    use omp_lib; use discrete_shock
    implicit none
    integer, intent(in) ::kfact
    real(dp), intent(in) ::pk,c2,sigma_e,k(nk),o(no),prbo(no,no),vguess(nk,no,ne)
    real(dp), intent(out)::expec(nk,no)
    ! main variables
    real(dp) e(ne),prbe(ne,ne),vg(nk,no,ne),v(nk,no,ne)
    ! Auxiliary variables
    real(dp) dif,vcurr,vmax,prcurr,mu_e
    real(dp) pprcurr(nk,no),iinv(nk,nk,ne)
    integer ik,io,ie,iter,j,i,imax
    real(8) t1,t2
    t1 = omp_get_wtime()
    vg = vguess
    mu_e  = 0d0!-sigma_e**2/2d0
    !call gete(mu_e,rho_e,sigma_e,ne,prbe,e)
    call tauchen86(rho_e,sigma_e,ne,4d0,prbe,e)!; e = e - sigma_e**2/2d0
    do ik=1,nk
       do io=1,no
          pprcurr(ik,io)=profit(k(ik),o(io))
       end do
       do ie=1,ne
          do i=1,nk
             iinv(i,ik,ie) = invadj(pk,c2,k(ik),k(i))*exp(e(ie))
          end do
       end do
    end do
    ! Iterations
    iter = 1; dif=1.0; !print*,'Running......'
    do while(dif.gt.erro2)
       do io=1,no
          do ik=1,nk
             expec(ik,io) = 0d0
             do j=1,ne
                expec(ik,io) = expec(ik,io) + prbe(1,j)*sum(vg(ik,:,j)*prbo(io,:))
             end do
          end do
       end do
       do ie=1,ne
          do io=1,no
             imax = kfact+1
             do ik=1,nk
                vmax  = -infty
                do i=max(ik-kfact,imax),nk
                   vcurr = pprcurr(ik,io) - iinv(i,ik,ie) + beta*expec(i,io)
                   if(vcurr-vmax.gt.erro1)then
                      vmax = vcurr; imax = i
                   else
                      exit
                   end if
                end do
                v(ik,io,ie) = vmax
             end do
          end do
       end do
       dif = maxval(abs(vg-v))
       vg = v
       !write(*,*) dif,iter
       !iter = iter + 1
    end do
    t2 = omp_get_wtime()
    !print*,'Equilibrium time: ',t2-t1
  end subroutine get_equil
end program montecarlo
  
