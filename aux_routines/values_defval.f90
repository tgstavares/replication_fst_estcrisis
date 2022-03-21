module values
  use prec; use parameters; use globals
contains

  ! ROUTINES FOR NON VECTORIZED ARRAYS
  subroutine cashflowr(par,r,tau,taudiv,k0,p0,z0,a0,k1,p1,qb,cf)
    implicit none
    type(param), intent(in)::par
    real(dp), intent(in) ::tau,r,taudiv,k0,p0,z0,a0,k1,p1,qb
    real(dp), intent(out)::cf
    real(dp) adj,adjs,inv0
    !$acc routine
    inv0 = k1-k0*(1d0-par%delta)
    adj  = 0.5d0 * par%a * (inv0 * k0**(-1))**2 * k0!  + par%gamma*k0
    adjs = 0d0
    if(abs(inv0).lt.1d-7)	adj  = 0d0
    if(inv0.lt.-1d-7) 		adjs =-par%ps*inv0
    ! CAHSFLOW WITH TAXES
    cf   = (1d0-taudiv)*((1d0-tau)*exp(z0+a0) * k0**theta + tau*par%delta*k0 - inv0 - adj - adjs + &
         (1d0+((1d0+r)/qb-1d0)*(1d0-tau))**(-1) * p1*k1 - p0*k0) - kappa
         !qb*(1d0+r*(1d0-tau))**(-1) * p1*k1 - p0*k0) - kappa
  end subroutine cashflowr

  ! subroutine colconstr(par,sparam,k0,z0,a0,k1,plimit)
  !   implicit none
  !   type(param), intent(in):: par
  !   real(dp), intent(in) ::sparam,k0,z0,a0,k1
  !   real(dp), intent(out)::plimit
  !   !$acc routine
  !   ! UNCOMMENT THE RELEVANT TYPE
  !   !borrowconstraint = sparam								! ASSET BASED COLLATERAL CONSTRAINT
  !   !borrowconstraint = sparam*profit(k0,z0,a0)/k1					! PROFIT BASED COLLATERAL CONSTRAINT
  !   !borrowconstraint = max((1d0+r)*(sparam - (phi/xi)*profit(k0,z0,a0)/k1),pmin)	! JERMAN AND QUADRINI 2012
  !   !plimit = sparam*max(sparam2*exp(z0+a0) * k0**theta * k1**(-1),sparam1)
  !   plimit = sparam*max(par%sparam2*exp(rhoz*z0+sigmaz**2*0.5d0 + rhoa*a0+sigmaa**2*0.5d0) * k0**theta * k1**(-1),par%sparam1)
  ! end subroutine colconstr

  subroutine get_vc(par,sparam,r,beta,tau,taudiv,k0,p0,z0,a0,k,p,ev,qb,v,qopt,ikp,ipp)
    implicit none
    type(param), intent(in):: par
    real(dp), intent(in) ::sparam,r,beta,tau,taudiv,k0,p0,z0,a0,k(nk),p(np),ev(nk,np),qb(nk,np)
    integer,  intent(out)::ikp,ipp
    real(dp), intent(out)::v,qopt
    integer ik1,ip1
    real(dp) vc,cf
    !$acc routine
    ikp = 1
    ipp = 1
    v   =-100000d0
    qopt= 1d0/(1d0+r)
    do ik1=1,nk
       do ip1=1,np
          call cashflowr(par,r,tau,taudiv,k0,p0,z0,a0,k(ik1),p(ip1),qb(ik1,ip1),cf)
          if(cf.ge.dexit)then
             vc = cf + beta*ev(ik1,ip1)
             if(vc.gt.v+0.000001d0)then             
                v   = vc
                qopt= qb(ik1,ip1)/(1d0+r)
                ikp = ik1
                ipp = ip1
             end if
          end if          
       end do
    end do
  end subroutine get_vc

  subroutine get_vc_nn(par,sparam,r,beta,tau,taudiv,nnk,nnp,k0,p0,z0,a0,k,p,coefevg,coefqbg,v,qopt,kp,pp)
    implicit none
    type(param), intent(in):: par
    integer,  intent(in) ::nnk,nnp
    real(dp), intent(in) ::sparam,r,beta,tau,taudiv,k0,p0,z0,a0,k(nk),p(np),coefevg(4,4,nk,np),coefqbg(4,4,nk,np)
    real(dp), intent(out)::v,qopt,kp,pp
    integer ik1,ip1
    real(dp) vc,cf,plimit,expec,qb1,k1,p1
    integer ier,iselect(6)
    iselect = 0
    iselect(1) = 1
    v   =-1000d0
    do ik1=1,nnk
       do ip1=1,nnp
          k1 = k(1) + dble(ik1-1) * (k(nk)-k(1)) / dble(nnk-1)
          p1 = p(1) + dble(ip1-1) * (p(np)-p(1)) / dble(nnp-1)
          call bcspeval(k1,p1,iselect,qb1,k,nk,p,np,0,0,coefqbg,nk,ier)
          call cashflowr(par,r,tau,taudiv,k0,p0,z0,a0,k1,p1,qb1,cf)
          if(cf.ge.dexit)then
             call bcspeval(k1,p1,iselect,expec,k,nk,p,np,0,0,coefevg,nk,ier)
             vc = cf + beta*expec
             if(vc.gt.v+0.00001d0)then
                v    = vc
                qopt = qb1
                kp   = k1
                pp   = p1
             end if
          end if
       end do
    end do
  end subroutine get_vc_nn

  subroutine get_kp_bp(par,m,iz0,ia0,ix0,k0,p0,v,kp,pp,qq,nnk,nnp)
    implicit none
    type(param), intent(in):: par
    type(model), intent(in):: m
    integer, intent(in)::iz0,ia0,ix0
    real(dp), intent(in)::k0,p0
    real(dp), intent(out)::v,kp,pp,qq
    integer, intent(in), optional:: nnk,nnp
    integer ik1,ip1
    real(dp) kpa,kpi,ppa,ppi,va,vi,z0,a0
    a0 = m%a(ia0)
    z0 = m%z(iz0)
    if(present(nnk).and.present(nnp))then
       call get_vc_nn(par,m%sparam,m%r,m%beta,m%tau,m%taudiv,nnk,nnp,k0,p0,z0,a0,m%k,m%p,m%coefevg(:,:,:,:,iz0,ia0,ix0),m%coefqbg(:,:,:,:,iz0,ia0,ix0),v,qq,kp,pp)
    else
       call get_vc(par,m%sparam,m%r,m%beta,m%tau,m%taudiv,k0,p0,z0,a0,m%k,m%p,m%evg(:,:,iz0,ia0,ix0),m%qbg(:,:,iz0,ia0,ix0),v,qq,ik1,ip1)
       kp = m%k(ik1)
       pp = m%p(ip1)
       !print*,ik1,ip1,iz0,qq
    end if

  end subroutine get_kp_bp

  ! ROUTINES FOR VECTORIZED ARRAYS
  subroutine gget_ev(delta,dwlk,k1,p1,z,a,ik1,ip1,iz0,ia0,ix0,prbz,prba,prbx,vg,ev,qb)
    implicit none
    real(dp), intent(in) ::delta,dwlk,k1,p1,z(nz),a(na)
    integer,  intent(in) ::ik1,ip1,iz0,ia0,ix0
    real(dp), intent(in) ::prbz(nz*nz),prba(na*na),prbx(nx*nx),vg(nk*np*nz*na*nx)
    real(dp), intent(out)::ev,qb
    integer iz1,ia1,ix1
    real(dp) vc,prob
    !$acc routine
    qb = 0d0
    ev = 0d0
    do ix1=1,nx
       do ia1=1,na
          do iz1=1,nz
             vc   = vg(ik1+nk*(ip1-1+np*(iz1-1+nz*(ia1-1+na*(ix1-1)))))
             prob = prba(ia0+na*(ia1-1))*prbz(iz0+nz*(iz1-1+nz*(ix0-1)))*prbx(ix0+nx*(ix1-1))
             if(vc.gt.0d0)then
                ev = ev + prob * vc                
                qb = qb + prob
             else
                qb = qb + prob * min(1d0,dwlk/p1)
                !qb = qb + prob*min(1d0,dwlk*(exp(z(iz1)+a(ia1))*k1**theta+(1d0-delta)*k1)/(k1*p1))
             end if
          end do
       end do
    end do
  end subroutine gget_ev

  
  ! ROUTINES FOR VECTORIZED ARRAYS
  subroutine gget_probdef(delta,dwlk,k1,p1,z,a,ik1,ip1,iz0,ia0,ix0,prbz,prba,prbx,vg,qb)
    implicit none
    real(dp), intent(in) ::delta,dwlk,k1,p1,z(nz),a(na)
    integer,  intent(in) ::ik1,ip1,iz0,ia0,ix0
    real(dp), intent(in) ::prbz(nz*nz),prba(na*na),prbx(nx*nx),vg(nk*np*nz*na*nx)
    real(dp), intent(out)::qb
    integer iz1,ia1,ix1
    real(dp) vc,prob
    !$acc routine
    qb = 0d0
    do ix1=1,nx
       do ia1=1,na
          do iz1=1,nz
             vc   = vg(ik1+nk*(ip1-1+np*(iz1-1+nz*(ia1-1+na*(ix1-1)))))
             prob = prba(ia0+na*(ia1-1))*prbz(iz0+nz*(iz1-1+nz*(ix0-1)))*prbx(ix0+nx*(ix1-1))
             if(vc.gt.0d0)then
                qb = qb + prob
             end if
          end do
       end do
    end do
  end subroutine gget_probdef

end module values
