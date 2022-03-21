module equilibrium
  use prec; use parameters; use globals
contains

  subroutine get_equil(par,m,showoutput,getinterpcoeff)
    use values
    implicit none
    type(model), intent(inout):: m
    type(param), intent(in):: par
    character(len=1), intent(in), optional::showoutput,getinterpcoeff
    integer flag_showoutput,flag_getinterpcoeff,iflag
    integer i,ik0,ip0,iz0,ia0,ix0,ik1,ip1
    real(dp) difv,cashflow
    integer  ilinx,iliny,ier
    real(dp) bcxmin(1),bcxmax(1),bcymin(1),bcymax(1),wk(5*max(nk,np))    

    ! PRINT OR SKIP OUTPUT
    flag_showoutput = 1
    if(present(showoutput))then
       if(showoutput(1:1).eq.'n'.or.showoutput(1:1).eq.'N') flag_showoutput = 0
    end if
    ! GENERATE INTERP 
    flag_getinterpcoeff = 0
    if(present(getinterpcoeff))then
       if(getinterpcoeff(1:1).eq.'y'.or.getinterpcoeff(1:1).eq.'Y') flag_getinterpcoeff = 1
    end if

    if((flag_showoutput.eq.1)) print*,""    
    ! DO VALUE FUNCTION ITERERATIONS
    difv = 1d0
    do while(difv.gt.erro4)
       m%iter = m%iter + 1
       !if(m%hflag.gt.0) m%hflag = m%hflag - 1

       ! HOWARD ACCELERATION
       i = 0
       !if(m%iter.gt.1.and.naccel.gt.0.and.m%hflag.eq.0.and.m%jflag.eq.0)then
       !if(m%iter.gt.1.and.naccel.gt.0.and.m%jflag.eq.0.and.m%iter.lt.11)then
       if(m%iter.gt.3.and.naccel.gt.0.and.m%jflag.eq.0.and.m%iter.lt.11)then
          do i=1,naccel

             !$acc parallel loop collapse(5) private(ix0,ia0,iz0,ip1,ik1)
             !$omp parallel do collapse(5) private(ix0,ia0,iz0,ip1,ik1)
             do ix0=1,nx
                do ia0=1,na
                   do iz0=1,nz
                      do ip1=1,np
                         do ik1=1,nk
                            call gget_ev(ik1,ip1,iz0,ia0,ix0,m%pz,m%pa,m%px,m%vg,m%evg(ik1,ip1,iz0,ia0,ix0))
                         end do
                      end do
                   end do
                end do
             end do
             !$omp end parallel do

             !$acc parallel loop collapse(5) private(ix0,ia0,iz0,ip0,ik0,ik1,ip1,cashflow)
             !$omp parallel do collapse(5) private(ix0,ia0,iz0,ip0,ik0,ik1,ip1,cashflow)
             do ix0=1,nx
                do ia0=1,na
                   do iz0=1,nz
                      do ip0=1,np
                         do ik0=1,nk
                            ik1 = m%ikp(ik0+nk*(ip0-1+np*(iz0-1+nz*(ia0-1+na*(ix0-1)))))
                            ip1 = m%ipp(ik0+nk*(ip0-1+np*(iz0-1+nz*(ia0-1+na*(ix0-1)))))
                            call cashflowr(par,m%r,m%tau,m%taudiv,m%k(ik0),m%p(ip0),m%z(iz0),m%a(ia0),m%k(ik1),m%p(ip1),cashflow)
                            m%vg(ik0+nk*(ip0-1+np*(iz0-1+nz*(ia0-1+na*(ix0-1))))) = cashflow + beta*m%evg(ik1,ip1,iz0,ia0,ix0)
                         end do
                      end do
                   end do
                end do
             end do
             !$omp end parallel do

          end do
       end if

       ! COMPUTE EXPECTED FUTURE VALUE

       !$acc parallel loop collapse(5) private(ix0,ia0,iz0,ip1,ik1)
       !$omp parallel do collapse(5) private(ix0,ia0,iz0,ip1,ik1)
       do ix0=1,nx
          do ia0=1,na
             do iz0=1,nz
                do ip1=1,np
                   do ik1=1,nk
                      call gget_ev(ik1,ip1,iz0,ia0,ix0,m%pz,m%pa,m%px,m%vg,m%evg(ik1,ip1,iz0,ia0,ix0))
                   end do
                end do
             end do
          end do
       end do
       !$omp end parallel do

       ! COMPUTE CURRENT VALUES OF ADJUSTMENT AND INACTION

       !$acc parallel loop collapse(5)
       !$omp parallel do collapse(5) private(ix0,ia0,iz0,ip0,ik0)
       do ix0=1,nx
          do ia0=1,na
             do iz0=1,nz
                do ip0=1,np
                   do ik0=1,nk
                      call get_vc(par,m%sparam,m%r,m%tau,m%taudiv,m%k(ik0),m%p(ip0),m%z(iz0),m%a(ia0),m%k,m%p,m%evg(:,:,iz0,ia0,ix0), &
                           m%vc( ik0+nk*(ip0-1+np*(iz0-1+nz*(ia0-1+na*(ix0-1))))), &
                           m%ikp(ik0+nk*(ip0-1+np*(iz0-1+nz*(ia0-1+na*(ix0-1))))), &
                           m%ipp(ik0+nk*(ip0-1+np*(iz0-1+nz*(ia0-1+na*(ix0-1))))))
                   end do
                end do
             end do
          end do
       end do
       !$omp end parallel do

       difv = sum(abs(m%vc-m%vg)/abs(m%vg+1))/dble(nk*np*nz*na)       
       m%vg = m%vc

       if(difv.le.erro4*8d0)then
          m%jflag = 1
       else
          m%jflag = 0
       end if
       if((flag_showoutput.eq.1)) &
            write(*,'(3i5,f17.7,12f17.2)') &
            m%iter,i,m%hflag,difv, &
            m%k(minval(m%ikp)),m%k(1),m%k(maxval(m%ikp)),m%k(nk), &
            m%p(minval(m%ipp)),m%p(1),m%p(maxval(m%ipp)),m%p(np), &
            minval(m%vg)
       if(m%iter.ge.nvfi) exit
       !if(mod(m%iter,9).eq.0.and.m%iter.gt.2.and.m%hflag.eq.0) then
       !   m%hflag = 25
       !end if
    end do

    ! COMPUTE COEFF OF INTERPOLATION
    if(flag_getinterpcoeff.eq.1)then
       do ix0=1,nx
          do ia0=1,na
             do iz0=1,nz
                m%coefevg(1,1,:,:,iz0,ia0,ix0) = m%evg(    :,:,iz0,ia0,ix0)
                call bcspline(m%k,nk,m%p,np, m%coefevg(:,:,:,:,iz0,ia0,ix0),nk,0,bcxmin,0,bcxmax,0,bcymin,0,bcymax,wk,5*max(nk,np),ilinx,iliny,ier)
             end do
          end do
       end do
    end if

    if((flag_showoutput.eq.1)) print*,""
  end subroutine get_equil
end module equilibrium
