module values
  use prec; use parameters; use globals
contains

  ! ROUTINES FOR NON VECTORIZED ARRAYS
  subroutine cf1f(par,tau,k0,p0,z0,a0,cf1)
    implicit none
    type(param), intent(in)::par
    real(dp), intent(in) ::tau,k0,p0,z0,a0
    real(dp), intent(out)::cf1
    !$acc routine
    cf1 = (1d0-tau)*exp(z0+a0)*k0**theta + tau*par%delta*k0 + k0*(1d0-par%delta) - p0*k0
  end subroutine cf1f
  
  subroutine cf2f(par,tau,k0,k1,p1,cf2)
    implicit none
    type(param), intent(in)::par
    real(dp), intent(in) ::tau,k0,k1,p1
    real(dp), intent(out)::cf2
    real(dp) adj,adjs,inv0
    !$acc routine
    inv0 = k1-k0*(1d0-par%delta)
    adj  = 0.5d0 * par%a * (inv0 * k0**(-1))**2 * k0  + par%gamma*k0
    adjs = 0d0
    if(abs(inv0).lt.1d-7)	adj  = 0d0
    if(inv0.lt.-1d-7) 		adjs =-par%ps*inv0
    ! CAHSFLOW WITH TAXES
    if(p1.ge.0d0)then
       cf2 = -k1 -adj -adjs + p1*k1/(1d0+r*(1d0-tau)) - mukappa
    else
       cf2 = -k1 -adj -adjs + p1*k1/(1d0+r-0.01d0) - mukappa
    end if
  end subroutine cf2f

  subroutine cfef(tau,k1,p1,cfe)
    implicit none
    real(dp), intent(in) ::tau,k1,p1
    real(dp), intent(out)::cfe
    cfe = -k1 + p1*k1/(1d0+r*(1d0-tau)) + cet !-ce*k1
    if(cfe.lt.0d0) cfe = cfe*(1d0+ceq)    
  end subroutine cfef
  
  subroutine colconstr(par,sparam,k0,z0,a0,k1,plimit)
    implicit none
    type(param), intent(in):: par
    real(dp), intent(in) ::sparam,k0,z0,a0,k1
    real(dp), intent(out)::plimit
    !$acc routine
    plimit = sparam*max(par%sparam2*exp(rhoz*z0+sigmaz**2*0.5d0 + rhoa*a0+sigmaa**2*0.5d0) * k1**theta * k1**(-1),par%sparam1)
  end subroutine colconstr

  subroutine get_vc(par,sparam,tau,taudiv,k0,p0,z0,a0,k,p,ev,v,ikp,ipp)
    implicit none
    type(param), intent(in):: par
    real(dp), intent(in) ::sparam,tau,taudiv,k0,p0,z0,a0,k(nk),p(np),ev(nk,np)
    integer,  intent(out)::ikp,ipp
    real(dp), intent(out)::v
    integer ik1,ip1
    real(dp) vc,cf,cf1,cf2,plimit
    !$acc routine
    ikp = 1
    ipp = 1
    v   =-1000d0
    call cf1f(par,tau,k0,p0,z0,a0,cf1)
    do ik1=1,nk
       do ip1=1,np
          call colconstr(par,sparam,k0,z0,a0,k(ik1),plimit)
          call cf2f(par,tau,k0,k(ik1),p(ip1),cf2)
          cf = cf1 + cf2
          if(cf.ge.0d0.and.p(ip1).le.plimit)then
             vc = cf + beta*ev(ik1,ip1)
             if(vc.gt.v+0.00001d0)then             
                v   = vc
                ikp = ik1
                ipp = ip1
             end if
          end if
          if(p(ip1).gt.plimit) exit
       end do
    end do
  end subroutine get_vc

  subroutine get_vc_nn(par,sparam,tau,taudiv,nnk,nnp,k0,p0,z0,a0,k,p,coefevg,v,kp,pp)
    implicit none
    type(param), intent(in):: par
    integer,  intent(in) ::nnk,nnp
    real(dp), intent(in) ::sparam,tau,taudiv,k0,p0,z0,a0,k(nk),p(np),coefevg(4,4,nk,np)
    real(dp), intent(out)::v,kp,pp
    integer ik1,ip1
    real(dp) vc,cf1,cf2,cf,plimit,expec,k1,p1
    integer ier,iselect(6)
    iselect = 0
    iselect(1) = 1
    v   =-1000d0
    call cf1f(par,tau,k0,p0,z0,a0,cf1)
    do ik1=1,nnk
       do ip1=1,nnp
          k1 = k(1) + dble(ik1-1) * (k(nk)-k(1)) / dble(nnk-1)
          p1 = p(1) + dble(ip1-1) * (p(np)-p(1)) / dble(nnp-1)
          call colconstr(par,sparam,k0,z0,a0,k1,plimit)          
          call cf2f(par,tau,k0,k1,p1,cf2)
          cf = cf1 + cf2
          if(cf.ge.0d0.and.p1.le.plimit)then
             call bcspeval(k1,p1,iselect,expec,k,nk,p,np,0,0,coefevg,nk,ier)
             vc = cf + beta*expec
             if(vc.gt.v+0.00001d0)then
                v  = vc
                kp = k1
                pp = p1
             end if
          end if
          if(p1.gt.plimit) exit
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
       call get_vc_nn(par,m%sparam,m%tau,m%taudiv,nnk,nnp,k0,p0,z0,a0,m%k,m%p,m%coefevg(:,:,:,:,iz0,ia0,ix0),v,kp,pp)
    else
       call get_vc(par,m%sparam,m%tau,m%taudiv,k0,p0,z0,a0,m%k,m%p,m%evg(:,:,iz0,ia0,ix0),v,ik1,ip1)
       kp = m%k(ik1)
       pp = m%p(ip1)
    end if
    qq = 1d0/(1d0+r)
  end subroutine get_kp_bp

  ! ROUTINES FOR VECTORIZED ARRAYS
  subroutine gget_ev(ik1,ip1,iz0,ia0,ix0,prbz,prba,prbx,vg,vx,ev)
    implicit none
    integer,  intent(in) ::ik1,ip1,iz0,ia0,ix0
    real(dp), intent(in) ::prbz(nz*nz),prba(na*na),prbx(nx*nx),vg(nk*np*nz*na*nx),vx(nk*np*nz*na)
    real(dp), intent(out)::ev
    integer iz1,ia1,ix1
    !$acc routine
    ev = 0d0
    do ix1=1,nx
       do ia1=1,na
          do iz1=1,nz
             ev = ev + &
                  prba(ia0+na*(ia1-1))*prbz(iz0+nz*(iz1-1+nz*(ix0-1)))*prbx(ix0+nx*(ix1-1)) * &

                  !max(vg(ik1+nk*(ip1-1+np*(iz1-1+nz*(ia1-1+na*(ix1-1))))),0d0)
                  
                  max( &
                  vg(ik1+nk*(ip1-1+np*(iz1-1+nz*(ia1-1+na*(ix1-1))))),	&
                  vx(ik1+nk*(ip1-1+np*(iz1-1+nz*(ia1-1))))		&
                  )
             
          end do
       end do
    end do
  end subroutine gget_ev

  
  subroutine get_ve(par,sparam,tau,taudiv,k0,p0,z0,a0,k,p,ev,v,ikp,ipp)
    implicit none
    type(param), intent(in):: par
    real(dp), intent(in) ::sparam,tau,taudiv,k0,p0,z0,a0,k(nk),p(np),ev(nk,np)
    integer,  intent(out)::ikp,ipp
    real(dp), intent(out)::v
    integer ik1,ip1
    real(dp) vc,cfe,plimit
    !$acc routine
    ikp = 1
    ipp = 1
    v   =-1000d0
    
    do ik1=1,nk
       do ip1=1,np
          call colconstr(par,sparam,0d0,z0,a0,k(ik1),plimit)
          call cfef(tau,k(ik1),p(ip1),cfe)          
          if(p(ip1).le.plimit)then
             vc = cfe + beta*ev(ik1,ip1)
             if(vc.gt.v+0.00001d0)then             
                v   = vc
                ikp = ik1
                ipp = ip1
             end if
          end if
          if(p(ip1).gt.plimit) exit
       end do
    end do
  end subroutine get_ve

  
  subroutine get_kpe_bpe(par,m,iz0,ia0,ix0,k0,p0,v,kp,pp,qq,nnk,nnp)
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
       stop 'FATAL ERROR - change code to allow for this'
       !call get_vc_nn(par,m%sparam,m%tau,m%taudiv,nnk,nnp,k0,p0,z0,a0,m%k,m%p,m%coefevg(:,:,:,:,iz0,ia0,ix0),v,kp,pp)
    else
       call get_ve(par,m%sparam,m%tau,m%taudiv,k0,p0,z0,a0,m%k,m%p,m%evg(:,:,iz0,ia0,ix0),v,ik1,ip1)
       kp = m%k(ik1)
       pp = m%p(ip1)
    end if
    qq = 1d0/(1d0+r)
  end subroutine get_kpe_bpe


end module values
