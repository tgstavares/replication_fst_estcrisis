program greek_crisis_exp
#ifdef _OPENMP
  use omp_lib
#endif
#ifdef _OPENACC
  use openacc
#endif
  use prec; use parameters; use globals; use process; use equilibrium; use simulation
  implicit none

  type(model) m,m_2,m_1,m0,m1,m2,m3,m4
  type(param) par
  integer, parameter:: tsimul = 114, nsimul = 50000, nvars = 24
  integer n,t,tdiscard,tsimul1,iloc(1),tt0,tt1,iflag_turnashocks,iflag_turnxshocks,iflag_creditshock,iflag_saved
  real(dp) ee,sparam,tau,taudiv
  integer,  allocatable:: iapath(:),ixpath(:)
  real(dp), allocatable:: shocksz(:,:),shocksa(:),shocksx(:),data(:,:,:),initconds(:,:)
  logical, allocatable:: mask(:,:),mask2(:)
  integer, parameter:: nstats = 35
  real statistics(nstats)
  real(dp) meanIK,d0IK,d1IK,d2IK,d3IK,d4IK
  real(dp) meanLP,d0LP,d1LP,d2LP,d3LP,d4LP
  real(dp) meanLP2,d0LP2,d1LP2,d2LP2,d3LP2,d4LP2
  real(dp) meanAV,d0AV,d1AV,d2AV,d3AV,d4AV
  real(dp) meanAP,d0AP,d1AP,d2AP,d3AP,d4AP
  real(dp) meanAI,d0AI,d1AI,d2AI,d3AI,d4AI
  real(dp) meanAD,d0AD,d1AD,d2AD,d3AD,d4AD
  real scale_param

#ifdef _OPENMP
  t1 = omp_get_wtime()
#else
  call cpu_time(t1)
#endif

  ! CHECK NUMBER OF OMP PROCS
#ifdef _OPENMP
  print*,"OMP Procs initiated:   ",OMP_GET_NUM_PROCS()
  print*,"OMP Threads initiated: ",OMP_GET_MAX_THREADS()
  print*,""
#endif

  ! ALLOCATE ARRAYS FOR SIMULATIONS
  allocate(shocksz(nsimul,tsimul),shocksa(tsimul),shocksx(tsimul))
  allocate(iapath(tsimul),ixpath(tsimul))
  allocate(data(nsimul,tsimul,nvars),initconds(nsimul,6))

  ! ALLOCATE VARIABLES FROM GLOBALS
  allocate(m%k(nk),m%p(np),m%z(nz),m%a(na),m%xsigmaz(nx))
  allocate(m%prbz(nz,nz,nx),m%prba(na,na),m%prbx(nx,nx))
  allocate(m%evg(nk,np,nz,na,nx))
  allocate(m%ikp(nk*np*nz*na*nx),m%ipp(nk*np*nz*na*nx))
  allocate(m%vg(nk*np*nz*na*nx),m%vc(nk*np*nz*na*nx),m%pz(nz*nz*nx),m%pa(na*na),m%px(nx*nx))
  allocate(m%coefevg(4,4,nk,np,nz,na,nx))

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
  call get_equil(par,m,showoutput='y')

  print*,"Simulate Greek crisis ..."
  print*,""

  ! PROGRAM OPTIONS
  iflag_turnashocks = 1
  iflag_turnxshocks = 1
  iflag_creditshock = 1
  iflag_saved       = 1
  sparam = m%sparam
  tau    = m%tau
  taudiv = m%taudiv
  ! GENERATE SHOCK SERIES FOR FIRMS AND SET A=0
  call get_simulshocks(nsimul,tsimul,shocksz,shocksa,shocksx,randomseed='n')

  iloc = minloc(abs(log(1.00d0)-m%a))
  iapath = iloc(1)
  if(iflag_turnashocks.eq.1)then

     iloc = minloc(abs( 0.15d0 - m%a)); iapath(103) = iloc(1)
     iloc = minloc(abs( 0.13d0 - m%a)); iapath(104) = iloc(1)
     iloc = minloc(abs( 0.07d0 - m%a)); iapath(105) = iloc(1)
     iloc = minloc(abs( 0.04d0 - m%a)); iapath(106) = iloc(1)
     iloc = minloc(abs( 0.11d0 - m%a)); iapath(107) = iloc(1)
     iloc = minloc(abs( 0.05d0 - m%a)); iapath(108) = iloc(1)
     iloc = minloc(abs( 0.01d0 - m%a)); iapath(109) = iloc(1)
     iloc = minloc(abs(-0.11d0 - m%a)); iapath(110) = iloc(1)
     iloc = minloc(abs(-0.18d0 - m%a)); iapath(111) = iloc(1)
     iloc = minloc(abs(-0.24d0 - m%a)); iapath(112) = iloc(1)
     iloc = minloc(abs(-0.08d0 - m%a)); iapath(113) = iloc(1)
     iloc = minloc(abs(-0.01d0 - m%a)); iapath(114) = iloc(1)

  end if

  ixpath = 1
  if(iflag_turnxshocks.eq.1)then
     ixpath(107) = 2
     ixpath(108) = 1
     ixpath(109) = 1
     ixpath(110) = 2
     ixpath(111) = 2
     ixpath(112) = 1
     ixpath(113) = 1
     ixpath(114) = 1
  end if

  initconds = 0d0
  ! GENERATE THE SIMULATION FOR THE SPECIFIC SEQUENCE OF SHOCKS
  if(iflag_creditshock.eq.0)then
     tt1 = 5
     call get_simulations(par,m,nsimul,tsimul-tt1,nvars,shocksz(:,1:tsimul-tt1),shocksa(1:tsimul-tt1),shocksx(1:tsimul-tt1),data(:,1:tsimul-tt1,:),iapath=iapath(1:tsimul-tt1),ixpath=ixpath(1:tsimul-tt1))
     do tt0 = 110,114
        ! GENERATE INITIAL CONTIONS TO FEED IN SIMULATIONS
        initconds(:,1) = data(:,tt0-1,6) 	! k0
        initconds(:,2) = data(:,tt0-1,7) 	! p0
        initconds(:,3) = data(:,tt0-1,5) 	! z_1
        initconds(:,4) = data(:,tt0-1,9) 	! t0
        initconds(:,5) = data(:,tt0-1,3) 	! k_1
        initconds(:,6) = data(:,tt0-1,10)	! inv_1
        call get_simulations(par,m,nsimul,1,nvars,shocksz(:,tt0),shocksa(tt0),shocksx(tt0),data(:,tt0,:),initconds=initconds,iapath=iapath(tt0),ixpath=ixpath(tt0),replace_exit='n')
        !print*,count(data(:,tt0,9).le.0.1d0)
     end do
  else
     tt1 = 5+2
     call get_simulations(par,m,nsimul,tsimul-tt1,nvars,shocksz(:,1:tsimul-tt1),shocksa(1:tsimul-tt1),shocksx(1:tsimul-tt1),data(:,1:tsimul-tt1,:),iapath=iapath(1:tsimul-tt1),ixpath=ixpath(1:tsimul-tt1))

     m%iter   = 0
     m%hflag  = 0
     m%jflag  = 0
     m_2= m
     m_1= m
     m0 = m
     m1 = m
     m2 = m
     m3 = m
     m4 = m

     ! CHOOSE SPECIFIC SEQUENCE
     ee = 0.0d0	! TO SIMULATE MORE PRONOUNCED CREDIT SUPPLY SHOCKS
     
     ! ! SPARAM UNDER BASELINE
     ! m_2%sparam = sparam*(1d0 - 0.001d0*(1d0+ee))
     ! m_1%sparam = sparam*(1d0 - 0.013d0*(1d0+ee))
     ! m0%sparam = sparam*(1d0 - 0.0220d0*(1d0+ee))
     ! m1%sparam = sparam*(1d0 - 0.043d0*(1d0+ee))
     ! m2%sparam = sparam*(1d0 - 0.068d0*(1d0+ee))
     ! m3%sparam = sparam*(1d0 - 0.089d0*(1d0+ee))
     ! m4%sparam = sparam*(1d0 - 0.095d0*(1d0+ee))     

     ! SPARAM UNDER PROVISION
     m_2%sparam = sparam*(1d0 + 0.00d0*(1d0+ee))
     m_1%sparam = sparam*(1d0 - 0.013d0*(1d0+ee))
     m0%sparam  = sparam*(1d0 - 0.030d0*(1d0+ee))
     m1%sparam  = sparam*(1d0 - 0.070d0*(1d0+ee))    
     m2%sparam  = sparam*(1d0 - 0.136d0*(1d0+ee))
     m3%sparam  = sparam*(1d0 - 0.180d0*(1d0+ee))
     m4%sparam  = sparam*(1d0 - 0.217d0*(1d0+ee))
     
     ! ! SPARAM UNDER NPL
     ! m_2%sparam = sparam*(1d0 + 0.015d0*(1d0+ee))
     ! m_1%sparam = sparam*(1d0 - 0.015d0*(1d0+ee))
     ! m0%sparam  = sparam*(1d0 - 0.060d0*(1d0+ee))
     ! m1%sparam  = sparam*(1d0 - 0.146d0*(1d0+ee))    
     ! m2%sparam  = sparam*(1d0 - 0.290d0*(1d0+ee))
     ! m3%sparam  = sparam*(1d0 - 0.350d0*(1d0+ee))
     ! m4%sparam  = sparam*(1d0 - 0.400d0*(1d0+ee))

     ! TAU
     m_2%tau= 0.25d0
     m_1%tau= 0.35d0
     m0%tau = 0.24d0
     m1%tau = 0.20d0
     m2%tau = 0.20d0
     m3%tau = 0.26d0
     m4%tau = 0.26d0
     
     ! ! TAUDIV
     ! m0%taudiv = 0.05d0
     ! m1%taudiv = 0.08d0
     ! m2%taudiv = 0.10d0
     ! m3%taudiv = 0.15d0
     ! m4%taudiv = 0.15d0
          
     do tt0 = 108,114
        select case (tt0)

        case (108)
           if(iflag_saved.eq.1)then
              open(1,file='functions/vg_2.txt',		position="rewind");read(1,*)m_2%vg	;close(1)
              open(1,file='functions/ikp_2.txt',	position="rewind");read(1,*)m_2%ikp	;close(1)
              open(1,file='functions/ipp_2.txt',	position="rewind");read(1,*)m_2%ipp	;close(1)
              !open(1,file='functions/iter0.txt',	position="rewind");read(1,*)m0%iter	;close(1)
              !open(1,file='functions/hflag0.txt',	position="rewind");read(1,*)m0%hflag	;close(1)
              !open(1,file='functions/jflag0.txt',	position="rewind");read(1,*)m0%jflag	;close(1)
           else
              !m0%iter   = 0
              !m0%hflag  = 0
              !m0%jflag  = 0
           end if
           write(*,'(a40,2f12.4,i10)')'##### Getting simulation for sparam: ',m_2%sparam,m_2%tau,tt0
           call get_equil(par,m_2,showoutput='y')
           open(1,file='functions/sparam__2.txt',position="rewind");write(1,*)m_2%sparam;close(1)
           open(1,file='functions/vg_2.txt',	position="rewind");write(1,*)m_2%vg	;close(1)
           open(1,file='functions/ikp_2.txt',	position="rewind");write(1,*)m_2%ikp	;close(1)
           open(1,file='functions/ipp_2.txt',	position="rewind");write(1,*)m_2%ipp	;close(1)
           open(1,file='functions/iter_2.txt',	position="rewind");write(1,*)m_2%iter	;close(1)
           open(1,file='functions/hflag_2.txt',	position="rewind");write(1,*)m_2%hflag	;close(1)
           open(1,file='functions/jflag_2.txt',	position="rewind");write(1,*)m_2%jflag	;close(1)

        case (109)
           if(iflag_saved.eq.1)then
              open(1,file='functions/vg_1.txt',		position="rewind");read(1,*)m_1%vg	;close(1)
              open(1,file='functions/ikp_1.txt',	position="rewind");read(1,*)m_1%ikp	;close(1)
              open(1,file='functions/ipp_1.txt',	position="rewind");read(1,*)m_1%ipp	;close(1)
              !open(1,file='functions/iter0.txt',	position="rewind");read(1,*)m0%iter	;close(1)
              !open(1,file='functions/hflag0.txt',	position="rewind");read(1,*)m0%hflag	;close(1)
              !open(1,file='functions/jflag0.txt',	position="rewind");read(1,*)m0%jflag	;close(1)
           else
              !m0%iter   = 0
              !m0%hflag  = 0
              !m0%jflag  = 0
           end if
           write(*,'(a40,2f12.4,i10)')'##### Getting simulation for sparam: ',m_1%sparam,m_1%tau,tt0
           call get_equil(par,m_1,showoutput='y')
           open(1,file='functions/sparam__1.txt',position="rewind");write(1,*)m_1%sparam;close(1)
           open(1,file='functions/vg_1.txt',	position="rewind");write(1,*)m_1%vg	;close(1)
           open(1,file='functions/ikp_1.txt',	position="rewind");write(1,*)m_1%ikp	;close(1)
           open(1,file='functions/ipp_1.txt',	position="rewind");write(1,*)m_1%ipp	;close(1)
           open(1,file='functions/iter_1.txt',	position="rewind");write(1,*)m_1%iter	;close(1)
           open(1,file='functions/hflag_1.txt',	position="rewind");write(1,*)m_1%hflag	;close(1)
           open(1,file='functions/jflag_1.txt',	position="rewind");write(1,*)m_1%jflag	;close(1)

        case (110)
           if(iflag_saved.eq.1)then
              open(1,file='functions/vg0.txt',		position="rewind");read(1,*)m0%vg	;close(1)
              open(1,file='functions/ikp0.txt',		position="rewind");read(1,*)m0%ikp	;close(1)
              open(1,file='functions/ipp0.txt',		position="rewind");read(1,*)m0%ipp	;close(1)
              !open(1,file='functions/iter0.txt',	position="rewind");read(1,*)m0%iter	;close(1)
              !open(1,file='functions/hflag0.txt',	position="rewind");read(1,*)m0%hflag	;close(1)
              !open(1,file='functions/jflag0.txt',	position="rewind");read(1,*)m0%jflag	;close(1)
           else
              !m0%iter   = 0
              !m0%hflag  = 0
              !m0%jflag  = 0
           end if
           write(*,'(a40,2f12.4,i10)')'##### Getting simulation for sparam: ',m0%sparam,m0%tau,tt0
           call get_equil(par,m0,showoutput='y')
           open(1,file='functions/sparam_0.txt',position="rewind");write(1,*)m0%sparam	;close(1)
           open(1,file='functions/vg0.txt',	position="rewind");write(1,*)m0%vg	;close(1)
           open(1,file='functions/ikp0.txt',	position="rewind");write(1,*)m0%ikp	;close(1)
           open(1,file='functions/ipp0.txt',	position="rewind");write(1,*)m0%ipp	;close(1)
           open(1,file='functions/iter0.txt',	position="rewind");write(1,*)m0%iter	;close(1)
           open(1,file='functions/hflag0.txt',	position="rewind");write(1,*)m0%hflag	;close(1)
           open(1,file='functions/jflag0.txt',	position="rewind");write(1,*)m0%jflag	;close(1)

        case(111)
           if(iflag_saved.eq.1)then
              open(1,file='functions/vg1.txt',		position="rewind");read(1,*)m1%vg	;close(1)
              open(1,file='functions/ikp1.txt',		position="rewind");read(1,*)m1%ikp	;close(1)
              open(1,file='functions/ipp1.txt',		position="rewind");read(1,*)m1%ipp	;close(1)
              !open(1,file='functions/iter1.txt',	position="rewind");read(1,*)m1%iter	;close(1)
              !open(1,file='functions/hflag1.txt',	position="rewind");read(1,*)m1%hflag	;close(1)
              !open(1,file='functions/jflag1.txt',	position="rewind");read(1,*)m1%jflag	;close(1)
           else
              !m1%iter   = 0
              !m1%hflag  = 0
              !m1%jflag  = 0
           end if           
           write(*,'(a40,2f12.4,i10)')'##### Getting simulation for sparam: ',m1%sparam,m1%tau,tt0
           call get_equil(par,m1,showoutput='y')
           open(1,file='functions/sparam_1.txt',position="rewind");write(1,*)m1%sparam	;close(1)
           open(1,file='functions/vg1.txt',	position="rewind");write(1,*)m1%vg	;close(1)
           open(1,file='functions/ikp1.txt',	position="rewind");write(1,*)m1%ikp	;close(1)
           open(1,file='functions/ipp1.txt',	position="rewind");write(1,*)m1%ipp	;close(1)
           open(1,file='functions/iter1.txt',	position="rewind");write(1,*)m1%iter	;close(1)
           open(1,file='functions/hflag1.txt',	position="rewind");write(1,*)m1%hflag	;close(1)
           open(1,file='functions/jflag1.txt',	position="rewind");write(1,*)m1%jflag	;close(1)           

        case(112)
           if(iflag_saved.eq.1)then
              open(1,file='functions/vg2.txt',		position="rewind");read(1,*)m2%vg	;close(1)
              open(1,file='functions/ikp2.txt',		position="rewind");read(1,*)m2%ikp	;close(1)
              open(1,file='functions/ipp2.txt',		position="rewind");read(1,*)m2%ipp	;close(1)
              !open(1,file='functions/iter2.txt',	position="rewind");read(1,*)m2%iter	;close(1)
              !open(1,file='functions/hflag2.txt',	position="rewind");read(1,*)m2%hflag	;close(1)
              !open(1,file='functions/jflag2.txt',	position="rewind");read(1,*)m2%jflag	;close(1)
           else
              !m2%iter   = 0
              !m2%hflag  = 0
              !m2%jflag  = 0
           end if           
           write(*,'(a40,2f12.4,i10)')'##### Getting simulation for sparam: ',m2%sparam,m2%tau,tt0
           call get_equil(par,m2,showoutput='y')
           open(1,file='functions/sparam_2.txt',position="rewind");write(1,*)m2%sparam	;close(1)
           open(1,file='functions/vg2.txt',	position="rewind");write(1,*)m2%vg	;close(1)
           open(1,file='functions/ikp2.txt',	position="rewind");write(1,*)m2%ikp	;close(1)
           open(1,file='functions/ipp2.txt',	position="rewind");write(1,*)m2%ipp	;close(1)
           open(1,file='functions/iter2.txt',	position="rewind");write(1,*)m2%iter	;close(1)
           open(1,file='functions/hflag2.txt',	position="rewind");write(1,*)m2%hflag	;close(1)
           open(1,file='functions/jflag2.txt',	position="rewind");write(1,*)m2%jflag	;close(1)

        case(113)
           if(iflag_saved.eq.1)then
              open(1,file='functions/vg3.txt',		position="rewind");read(1,*)m3%vg	;close(1)
              open(1,file='functions/ikp3.txt',		position="rewind");read(1,*)m3%ikp	;close(1)
              open(1,file='functions/ipp3.txt',		position="rewind");read(1,*)m3%ipp	;close(1)
              !open(1,file='functions/iter3.txt',	position="rewind");read(1,*)m3%iter	;close(1)
              !open(1,file='functions/hflag3.txt',	position="rewind");read(1,*)m3%hflag	;close(1)
              !open(1,file='functions/jflag3.txt',	position="rewind");read(1,*)m3%jflag	;close(1)
           else
              !m3%iter   = 0
              !m3%hflag  = 0
              !m3%jflag  = 0
           end if           
           write(*,'(a40,2f12.4,i10)')'##### Getting simulation for sparam: ',m3%sparam,m3%tau,tt0
           call get_equil(par,m3,showoutput='y')
           open(1,file='functions/sparam_3.txt',position="rewind");write(1,*)m3%sparam	;close(1)
           open(1,file='functions/vg3.txt',	position="rewind");write(1,*)m3%vg	;close(1)
           open(1,file='functions/ikp3.txt',	position="rewind");write(1,*)m3%ikp	;close(1)
           open(1,file='functions/ipp3.txt',	position="rewind");write(1,*)m3%ipp	;close(1)
           open(1,file='functions/iter3.txt',	position="rewind");write(1,*)m3%iter	;close(1)
           open(1,file='functions/hflag3.txt',	position="rewind");write(1,*)m3%hflag	;close(1)
           open(1,file='functions/jflag3.txt',	position="rewind");write(1,*)m3%jflag	;close(1)

        case(114)
           if(iflag_saved.eq.1)then
              open(1,file='functions/vg4.txt',		position="rewind");read(1,*)m4%vg	;close(1)
              open(1,file='functions/ikp4.txt',		position="rewind");read(1,*)m4%ikp	;close(1)
              open(1,file='functions/ipp4.txt',		position="rewind");read(1,*)m4%ipp	;close(1)
              !open(1,file='functions/iter4.txt',	position="rewind");read(1,*)m4%iter	;close(1)
              !open(1,file='functions/hflag4.txt',	position="rewind");read(1,*)m4%hflag	;close(1)
              !open(1,file='functions/jflag4.txt',	position="rewind");read(1,*)m4%jflag	;close(1)
           else
              !m4%iter   = 0
              !m4%hflag  = 0
              !m4%jflag  = 0
           end if           
           write(*,'(a40,2f12.4,i10)')'##### Getting simulation for sparam: ',m4%sparam,m4%tau,tt0
           call get_equil(par,m4,showoutput='y')
           open(1,file='functions/sparam_4.txt',position="rewind");write(1,*)m4%sparam	;close(1)
           open(1,file='functions/vg4.txt',	position="rewind");write(1,*)m4%vg	;close(1)
           open(1,file='functions/ikp4.txt',	position="rewind");write(1,*)m4%ikp	;close(1)
           open(1,file='functions/ipp4.txt',	position="rewind");write(1,*)m4%ipp	;close(1)
           open(1,file='functions/iter4.txt',	position="rewind");write(1,*)m4%iter	;close(1)
           open(1,file='functions/hflag4.txt',	position="rewind");write(1,*)m4%hflag	;close(1)
           open(1,file='functions/jflag4.txt',	position="rewind");write(1,*)m4%jflag	;close(1)
        end select

        ! GENERATE INITIAL CONTIONS TO FEED IN SIMULATIONS
        initconds(:,1) = data(:,tt0-1,6) ! k0
        initconds(:,2) = data(:,tt0-1,7) ! p0
        initconds(:,3) = data(:,tt0-1,5) ! z_1
        initconds(:,4) = data(:,tt0-1,9) ! t0
        initconds(:,5) = data(:,tt0-1,3) ! k_1
        initconds(:,6) = data(:,tt0-1,10)! inv_1

        select case (tt0)
        case (108)
           call get_simulations(par,m_2,nsimul,1,nvars,shocksz(:,tt0),shocksa(tt0),shocksx(tt0),data(:,tt0,:),initconds=initconds,iapath=iapath(tt0),ixpath=ixpath(tt0),replace_exit='n')
        case (109)
           call get_simulations(par,m_1,nsimul,1,nvars,shocksz(:,tt0),shocksa(tt0),shocksx(tt0),data(:,tt0,:),initconds=initconds,iapath=iapath(tt0),ixpath=ixpath(tt0),replace_exit='n')
        case (110)
           call get_simulations(par,m0,nsimul,1,nvars,shocksz(:,tt0),shocksa(tt0),shocksx(tt0),data(:,tt0,:),initconds=initconds,iapath=iapath(tt0),ixpath=ixpath(tt0),replace_exit='n')
        case (111)
           call get_simulations(par,m1,nsimul,1,nvars,shocksz(:,tt0),shocksa(tt0),shocksx(tt0),data(:,tt0,:),initconds=initconds,iapath=iapath(tt0),ixpath=ixpath(tt0),replace_exit='n')
        case (112)
           call get_simulations(par,m2,nsimul,1,nvars,shocksz(:,tt0),shocksa(tt0),shocksx(tt0),data(:,tt0,:),initconds=initconds,iapath=iapath(tt0),ixpath=ixpath(tt0),replace_exit='n')
        case (113)
           call get_simulations(par,m3,nsimul,1,nvars,shocksz(:,tt0),shocksa(tt0),shocksx(tt0),data(:,tt0,:),initconds=initconds,iapath=iapath(tt0),ixpath=ixpath(tt0),replace_exit='n')
        case (114)
           call get_simulations(par,m4,nsimul,1,nvars,shocksz(:,tt0),shocksa(tt0),shocksx(tt0),data(:,tt0,:),initconds=initconds,iapath=iapath(tt0),ixpath=ixpath(tt0),replace_exit='n')
        end select

     end do
  end if

  ! SHOW SEQUENCE OF STATS FOR SIMULATION
  write(*,'(1a4,1a8,22a11)')'t','obs','a','x','e_ik','m_ik','a_ik', &
       'e_bkp','m_bkp','q66_bkp','a_bkp', &
       'ina','spike','e_k', &
       'a_salesk','sd_mpk','a_sales','a_inv','a_b'
  open(1,file='plot_data/Results_cal_estimation.txt',position="rewind")
  scale_param = 1000000.0
  do t=101,tsimul
     n = 3
     tdiscard = t; tsimul1 = tdiscard + 0
     allocate(mask(nsimul,1))
     mask = data(:,tdiscard:tsimul1,9).gt.0.5d0
     call get_stats(data(:,tdiscard:tsimul1,:),nstats,statistics,showoutput='n',cond=mask)     
     write(*,'(1i4,1i8,22f11.4)')t,count(mask),m%a(iapath(t)),m%xsigmaz(ixpath(t)), &
          statistics(1),statistics(2),statistics(3), &
          statistics(4),statistics(5),statistics(30),statistics(6), &
          statistics(9),statistics(12),statistics(23), &
          statistics(25),statistics(26), &
          statistics(24)/scale_param,statistics(29)/scale_param, &
          statistics(27)/scale_param
     write(1,'(1i4,1i8,22f11.4)')t,count(mask),m%a(iapath(t)),m%xsigmaz(ixpath(t)), &
          statistics(1),statistics(2),statistics(3), &
          statistics(4),statistics(5),statistics(30),statistics(6), &
          statistics(9),statistics(12),statistics(23), &
          statistics(25),statistics(26), &
          statistics(24)/scale_param,statistics(29)/scale_param, &
          statistics(27)/scale_param
     deallocate(mask)
  end do
  close(1)

  ! SELECT WHICH CASE YOU WANT TO SAVE TO FILE

  !open(3,file='plot_data/Distributions_cal_estimation_trash.txt',position="rewind")
  !open(3,file='plot_data/Distributions_cal_estimation_PU_nonexporters.txt',position="rewind")
  !open(3,file='plot_data/Distributions_cal_estimation_PUT_nonexporters.txt',position="rewind")
  !open(3,file='plot_data/Distributions_cal_estimation_PUTebkacc_nonexporters.txt',position="rewind")
  open(3,file='plot_data/Distributions_cal_estimation_PUTebkprov_nonexporters.txt',position="rewind")
  !open(3,file='plot_data/Distributions_cal_estimation_PUTebknpl.txt',position="rewind")  

  do t=101,tsimul ! SEE LEGACY CODES FOR CONSISTENCY
     do n=1,nsimul
        write(3,'(2i10,21f16.5)')t,n,m%sparam, &
             1.0d0-0d0,data(n,t,19),data(n,t,12), &
             data(n,t,23),data(n,t,21),data(n,t,5),data(n,t,3),data(n,t,4),data(n,t,18), &
             data(n,t,6),data(n,t,7)
     end do
  end do
  write(3,'(//)')
  close(3)

  ! VERIFY DIFFERENCES BETWEEN NORMAL AND CRISIS PERIOD
  print*,""  
  meanIK = sum(data(:,103:107,12))/dble(nsimul*5)
  meanLP = sum(data(:,103:107, 7))/dble(nsimul*5)
  meanLP2= sum(data(:,103:107,16))/sum(data(:,103:107,6))
  meanAV = sum(data(:,103:107,20))/dble(5)
  meanAP = sum(data(:,103:107,18))/dble(5)
  meanAI = sum(data(:,103:107,10))/dble(5)
  meanAD = sum(data(:,103:107,16))/dble(5)

  allocate(mask2(nsimul))
  mask2 = data(:,110,9).gt.0.5d0
  d0IK = sum(data(:,110,12),mask2)/dble(count(mask2)) - meanIK
  d0LP = sum(data(:,110, 7),mask2)/dble(count(mask2))/meanLP - 1d0
  d0LP2= sum(data(:,110,16),mask2)/sum(data(:,110,6),mask2)/meanLP2 - 1d0
  d0AV = sum(data(:,110,20),mask2)/meanAV
  d0AP = sum(data(:,110,18),mask2)/meanAP
  d0AI = sum(data(:,110,10),mask2)/meanAI
  d0AD = sum(data(:,110,16),mask2)/meanAD
  deallocate(mask2)
  allocate(mask2(nsimul))
  mask2 = data(:,111,9).gt.0.5d0  
  d1IK = sum(data(:,111,12),mask2)/dble(count(mask2)) - meanIK
  d1LP = sum(data(:,111, 7),mask2)/dble(count(mask2))/meanLP - 1d0
  d1LP2= sum(data(:,111,16),mask2)/sum(data(:,111,6),mask2)/meanLP2 - 1d0
  d1AV = sum(data(:,111,20),mask2)/meanAV
  d1AP = sum(data(:,111,18),mask2)/meanAP
  d1AI = sum(data(:,111,10),mask2)/meanAI
  d1AD = sum(data(:,111,16),mask2)/meanAD
  deallocate(mask2)
  allocate(mask2(nsimul))
  mask2 = data(:,112,9).gt.0.5d0  
  d2IK = sum(data(:,112,12),mask2)/dble(count(mask2)) - meanIK
  d2LP = sum(data(:,112, 7),mask2)/dble(count(mask2))/meanLP - 1d0
  d2LP2= sum(data(:,112,16),mask2)/sum(data(:,112,6),mask2)/meanLP2 - 1d0
  d2AV = sum(data(:,112,20),mask2)/meanAV
  d2AP = sum(data(:,112,18),mask2)/meanAP
  d2AI = sum(data(:,112,10),mask2)/meanAI
  d2AD = sum(data(:,112,16),mask2)/meanAD
  deallocate(mask2)
  allocate(mask2(nsimul))
  mask2 = data(:,113,9).gt.0.5d0  
  d3IK = sum(data(:,113,12),mask2)/dble(count(mask2)) - meanIK
  d3LP = sum(data(:,113, 7),mask2)/dble(count(mask2))/meanLP - 1d0
  d3LP2= sum(data(:,113,16),mask2)/sum(data(:,113,6),mask2)/meanLP2 - 1d0
  d3AV = sum(data(:,113,20),mask2)/meanAV
  d3AP = sum(data(:,113,18),mask2)/meanAP
  d3AI = sum(data(:,113,10),mask2)/meanAI
  d3AD = sum(data(:,113,16),mask2)/meanAD
  deallocate(mask2)
  allocate(mask2(nsimul))
  mask2 = data(:,114,9).gt.0.5d0  
  d4IK = sum(data(:,114,12),mask2)/dble(count(mask2)) - meanIK
  d4LP = sum(data(:,114, 7),mask2)/dble(count(mask2))/meanLP - 1d0
  d4LP2= sum(data(:,114,16),mask2)/sum(data(:,114,6),mask2)/meanLP2 - 1d0
  d4AV = sum(data(:,114,20),mask2)/meanAV
  d4AP = sum(data(:,114,18),mask2)/meanAP
  d4AI = sum(data(:,114,10),mask2)/meanAI
  d4AD = sum(data(:,114,16),mask2)/meanAD
  deallocate(mask2)

  write(*,'(a15,6f18.4)')'Diff mean inv rate',meanik,d0ik,d1ik,d2ik,d3ik,d4ik
  write(*,'(a15,6f18.4)')'Diff mean leverage',meanlp,d0lp,d1lp,d2lp,d3lp,d4lp
  write(*,'(a15,6f18.4)')'Diff agg  leverage',meanlp2,d0lp2,d1lp2,d2lp2,d3lp2,d4lp2
  write(*,'(a15,6f18.4)')'Diff agg value    ',meanav,d0av,d1av,d2av,d3av,d4av
  write(*,'(a15,6f18.4)')'Diff agg profit   ',meanap,d0ap,d1ap,d2ap,d3ap,d4ap
  write(*,'(a15,6f18.4)')'Diff agg invest   ',meanai,d0ai,d1ai,d2ai,d3ai,d4ai
  write(*,'(a15,6f18.4)')'Diff agg debt     ',meanad,d0ad,d1ad,d2ad,d3ad,d4ad
  print*,""

  ! CHECK COMPUTATION TIME
  ! CONCLUDE PROGRAM
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

  deallocate(iapath,ixpath)
  deallocate(shocksz,shocksa,shocksx)
  deallocate(data,initconds)
end program greek_crisis_exp
