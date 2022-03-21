program gplot
#ifdef _OPENMP
  use omp_lib
#endif
#ifdef _OPENACC
  use openacc
#endif
  use prec; use parameters; use globals; use process; use equilibrium; use values
  implicit none

  type(model) m
  type(param) par
  integer i,j,ik0,ik1,ip0,ip1,iz0,iz1,ia0,ix0,iaux(1),nnk,nnp,iflag_searchwide
  real(dp) e_k,e_b
  real(dp) k0,k1,p0,p1,z0,z1,a0,a1,inv0,val0,div0,qq,plimit,pdef
  integer ier,iselect(6)
  iselect = 0
  iselect(1) = 1

  print*,""
  print*,"Initiate program ..."
  print*,""

#ifdef _OPENMP
  t1 = omp_get_wtime()
#else
  call cpu_time(t1)
#endif

  ! ALLOCATE ARRAYS FROM GLOBALS.F90

  ! ALLOCATE VARIABLES FROM GLOBALS
  allocate(m%k(nk),m%p(np),m%z(nz),m%a(na),m%xsigmaz(nx))
  allocate(m%prbz(nz,nz,nx),m%prba(na,na),m%prbx(nx,nx))
  allocate(m%evg(nk,np,nz,na,nx),m%qbg(nk,np,nz,na,nx))
  allocate(m%ikp(nk*np*nz*na*nx),m%ipp(nk*np*nz*na*nx))
  allocate(m%vg(nk*np*nz*na*nx),m%vc(nk*np*nz*na*nx),m%pz(nz*nz*nx),m%pa(na*na),m%px(nx*nx))
  allocate(m%coefevg(4,4,nk,np,nz,na,nx),m%coefqbg(4,4,nk,np,nz,na,nx))

  ! RECOVER EQUILIBRIUM
  m%prbx(1,:) = (/1.0d0-uncfreq, uncfreq/); m%prbx(1,:) = m%prbx(1,:)/sum(m%prbx(1,:)) 	! NOT SURE EXACTLY WHY WE NEED THESE TWO LINES OF CODE. 
  m%prbx(2,:) = (/1.0d0-uncpers, uncpers/); m%prbx(2,:) = m%prbx(2,:)/sum(m%prbx(2,:)) 	! IF THEY'RE NOT INCLUDED THE ACC VERSION BREAKS DOWN
  open(1,file='functions/k.txt',	position="rewind");read(1,*)m%k		;close(1)
  open(1,file='functions/p.txt',	position="rewind");read(1,*)m%p		;close(1)
  open(1,file='functions/a.txt',	position="rewind");read(1,*)m%a		;close(1)
  open(1,file='functions/z.txt',	position="rewind");read(1,*)m%z		;close(1)
  open(1,file='functions/xsigmaz.txt',	position="rewind");read(1,*)m%xsigmaz	;close(1)
  open(1,file='functions/prba.txt',	position="rewind");read(1,*)m%prba	;close(1)
  open(1,file='functions/prbz.txt',	position="rewind");read(1,*)m%prbz	;close(1)
  open(1,file='functions/prbx.txt',	position="rewind");read(1,*)m%prbx	;close(1)
  open(1,file='functions/pa.txt',	position="rewind");read(1,*)m%pa	;close(1)
  open(1,file='functions/pz.txt',	position="rewind");read(1,*)m%pz	;close(1)
  open(1,file='functions/px.txt',	position="rewind");read(1,*)m%px	;close(1)
  open(1,file='functions/sparam.txt',	position="rewind");read(1,*)m%sparam	;close(1)
  open(1,file='functions/tau.txt',	position="rewind");read(1,*)m%tau	;close(1)
  open(1,file='functions/taudiv.txt',	position="rewind");read(1,*)m%taudiv	;close(1)
  open(1,file='functions/r.txt',	position="rewind");read(1,*)m%r		;close(1)
  open(1,file='functions/beta.txt',	position="rewind");read(1,*)m%beta	;close(1)
  open(1,file='functions/dwl.txt',	position="rewind");read(1,*)m%dwl	;close(1)
  open(1,file='functions/vg.txt',	position="rewind");read(1,*)m%vg	;close(1)
  open(1,file='functions/ikp.txt',	position="rewind");read(1,*)m%ikp	;close(1)
  open(1,file='functions/ipp.txt',	position="rewind");read(1,*)m%ipp	;close(1)
  open(1,file='functions/iter.txt',	position="rewind");read(1,*)m%iter	;close(1)
  open(1,file='functions/jflag.txt',	position="rewind");read(1,*)m%jflag	;close(1)
  open(1,file='functions/hflag.txt',	position="rewind");read(1,*)m%hflag	;close(1)
  open(1,file='functions/delta.txt',	position="rewind");read(1,*)par%delta	;close(1)
  open(1,file='functions/gamma.txt',	position="rewind");read(1,*)par%gamma	;close(1)
  open(1,file='functions/aa.txt',	position="rewind");read(1,*)par%a	;close(1)
  open(1,file='functions/ps.txt',	position="rewind");read(1,*)par%ps	;close(1)
  open(1,file='functions/sparam1.txt',	position="rewind");read(1,*)par%sparam1	;close(1)
  open(1,file='functions/sparam2.txt',	position="rewind");read(1,*)par%sparam2	;close(1)
  print*,"Recover equilibrium ..."
  print*,""
  call get_equil(par,m,showoutput='y',getinterpcoeff='y')

  ! WANT A LARGER SEARCH?
  iflag_searchwide = 0
  nnk = nk*1
  nnp = np*1

  ! WRITE POLICY AND PRICE FUNCTIONS
  print*,"Write policies for capital, investment, investment rate, and hiring ..."
  print*,""

  open(1,file='data/e_k.txt',position="rewind");read(1,*)e_k;close(1)
  open(1,file='data/e_b.txt',position="rewind");read(1,*)e_b;close(1)

  open(1,file='plot_data/Policyfnk.txt',position="rewind")
  open(2,file='plot_data/Policyfnz.txt',position="rewind")
  open(3,file='plot_data/Pricedebtfnp.txt',position="rewind")

  ! FIX SOME A AND X
  ix0 = 1
  ia0 = na/2+1
  a0 = m%a(ia0)
  do i=1,2

     ! HIGH AND LOW LEVERAGE FOR POLICY FUNCTIONS
     if(i.eq.1)then
        iaux = minloc(abs(m%p  - e_b*0.90)); ip0  = iaux(1)
     elseif(i.eq.2)then
        iaux = minloc(abs(m%p  - e_b*1.10)); ip0  = iaux(1)
     end if
     p0 = m%p(ip0)
     !p0 = m%p(15)
     print*,p0

     ! FIX ON SOME Z - CHANGE K
     iz0 = nz/2+1
     z0 = m%z(iz0)
     do ik0=1,nk
        k0 = m%k(ik0)
        if(iflag_searchwide.eq.1)then
           call get_kp_bp(par,m,iz0,ia0,ix0,k0,p0,val0,k1,p1,qq,nnk,nnp)
        else
           call get_kp_bp(par,m,iz0,ia0,ix0,k0,p0,val0,k1,p1,qq)
        end if
        inv0 = k1 - (1d0-par%delta)*k0
        call cashflowr(par,m%r,m%tau,m%taudiv,k0,p0,z0,a0,k1,p1,qq,div0)
        !call colconstr(par,m%sparam,k0,z0,a0,k1,plimit)
        write(1,'(15f20.6)')k0,p0,z0,k1,p1,inv0,inv0/k0,p1,div0,val0,a0,plimit
     end do
     write(1,'(//)')

     ! FIX SOME K - CHANGE Z
     iaux = minloc(abs(m%k  - e_k*0.70d0));  ik0  = iaux(1)
     k0 = m%k(ik0)
     do iz0=1,nz
        z0 = m%z(iz0)
        if(iflag_searchwide.eq.1)then
           call get_kp_bp(par,m,iz0,ia0,ix0,k0,p0,val0,k1,p1,qq,nnk,nnp)
        else
           call get_kp_bp(par,m,iz0,ia0,ix0,k0,p0,val0,k1,p1,qq)
        end if
        inv0 = k1 - (1d0-par%delta)*k0
        call cashflowr(par,m%r,m%tau,m%taudiv,k0,p0,z0,a0,k1,p1,qq,div0)
        write(2,'(15f20.6)')k0,p0,z0,k1,p1,inv0,inv0/k0,p1,div0,val0,a0
     end do
     write(2,'(//)')

     ! HIGH AND LOW Z FOR PRICE FUNCTIONS - FIX K1 AND CHANGE P1
     if(i.eq.1)then
        iz0 = nz/2+1+2
     elseif(i.eq.2)then
        iz0 = nz/2+1-2
     end if
     iaux = minloc(abs(m%k  - e_k*1.0d0));  ik0  = iaux(1)
     k0 = m%k(ik0)
     k1 = k0
     z0 = m%z(iz0)

     print*,'look',i,e_k,k0
     
     if(iflag_searchwide.eq.1)then
        !call get_kp_bp(par,m,iz0,ia0,ix0,k0,p0,val0,k1,p1,qq,nnk,nnp)
        iselect = 0
        iselect(1) = 1
        do j=1,nnp
           p0 = m%p(1) + dble(j-1)*(m%p(np)-m%p(1))/dble(nnp-1)
           call bcspeval(m%k(ik0),p0,iselect,qq,m%k,nk,m%p,np,0,0,m%coefqbg(:,:,:,:,iz0,ia0,ix0),nk,ier)
           !print*,iz0,ia0,ix0,p0,qq
           qq = qq/(1d0+m%r)
           write(3,'(15f20.6)')k0,p0,z0,p0,qq,qq*p0*k0,a0
        end do
        !stop
     else
        do ip0=1,np        
           p0 = m%p(ip0)

           !call get_kp_bp(par,m,iz0,ia0,ix0,k0,p0,val0,k1,p1,qq)
           qq = m%qbg(ik0,ip0,iz0,ia0,ix0)
           call gget_probdef(par%delta,m%beta,m%dwl,m%k(ik0),m%p(ip0),m%z,m%a,ik0,ip0,iz0,ia0,ix0,m%pz,m%pa,m%px,m%vg,pdef)
           !call gget_probdef(par%delta,m%dwlk,m%k(ik0),m%p(ip0),m%z,m%a,ik0,ip0,iz0,ia0,ix0,m%pz,m%pa,m%px,m%vg,pdef)

           qq = qq/(1d0+m%r)
           write(3,'(15f20.6)')k0,p0,z0,p0,qq,qq*p0*k0,a0,pdef
        end do
     end if

     write(3,'(//)')
  end do
  close(1)
  close(2)
  close(3)

  ! CONLCUDE PROGRAM
#ifdef _OPENMP
  t2 = omp_get_wtime()
#else
  call cpu_time(t2)
#endif
  print*,'Program execution time: ',t2-t1
  print*,""


  ! DEALLOCATE ARRAYS
  deallocate(m%k,m%p,m%z,m%a,m%xsigmaz)
  deallocate(m%prbz,m%prba,m%prbx)
  deallocate(m%evg)
  deallocate(m%ikp,m%ipp)
  deallocate(m%vg,m%vc,m%pz,m%pa,m%px)
  deallocate(m%coefevg)

end program gplot
