program main
#ifdef _OPENMP
  use omp_lib
#endif
#ifdef _OPENACC
  use openacc
#endif
  use prec; use parameters; use globals; use process; use equilibrium
  implicit none

  type(model) m
  type(param) par
  integer iflag_newguess,iflag_pregrida,iflag_whichsample
  integer ik0,ip0

  print*,""
  print*,"Initiate program ..."
  print*,""

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

  ! PROGRAM OPTIONS
  iflag_newguess = 1
  iflag_pregrida = 1
  iflag_whichsample = 1

  ! SET SEARCH PARAM
  par%delta   = 0.1338d0
  par%gamma   = 0.00002d0
  par%a       = 0.532d0
  par%ps      = 0.02400d0
  par%sparam1 = 0.609d0
  par%sparam2 = 1.8751d0
  
  m%sparam    = 1d0
  m%tau       = 0.32d0
  m%taudiv    = 0.0d0

  ! ALLOCATE ARRAYS FROM GLOBALS.F90
  allocate(m%k(nk),m%p(np),m%z(nz),m%a(na),m%xsigmaz(nx))
  allocate(m%prbz(nz,nz,nx),m%prba(na,na),m%prbx(nx,nx))
  allocate(m%evg(nk,np,nz,na,nx))
  allocate(m%ikp(nk*np*nz*na*nx),m%ipp(nk*np*nz*na*nx))
  allocate(m%vg(nk*np*nz*na*nx),m%vc(nk*np*nz*na*nx),m%vx(nk*np*nz*na),m%pz(nz*nz*nx),m%pa(na*na),m%px(nx*nx))
  allocate(m%coefevg(4,4,nk,np,nz,na,nx))

  ! SET UP GRIDS
  do ik0 = nk,1,-1
     m%k(ik0) = kmax*(1d0-par%delta)**(dble(nk-ik0)/dble(nk_divs))
  end do
  do ip0 = 1,np
     m%p(ip0) = pmin + dble(ip0-1)*(pmax - pmin)/dble(np-1)
  end do
  if(nx.eq.1)then
     m%xsigmaz(1) = sigmaz
     m%prbx(1,1)  = 1d0
     zmax	  = nsd*sigmaz/sqrt(1d0-rhoz**2d0)
     zmin         =-zmax
     call get_markov(nz,zmin,zmax,rhoz,sigmaz,m%prbz,m%z)
     m%z = m%z + mu1
  else
     m%xsigmaz   = (/sigmaz , sjump*sigmaz/)
     m%prbx(1,:) = (/1.0d0-uncfreq, uncfreq/); m%prbx(1,:) = m%prbx(1,:)/sum(m%prbx(1,:))
     m%prbx(2,:) = (/1.0d0-uncpers, uncpers/); m%prbx(2,:) = m%prbx(2,:)/sum(m%prbx(2,:))
     zmax	  = nsd*sjump*sigmaz/sqrt(1d0-rhoz**2d0)
     zmin         =-zmax
     call get_markovforts(nz,mu1,mu2,m%xsigmaz(1),m%xsigmaz(2),rhoz,rhoz,zmin,zmax,m%prbz(:,:,1),m%prbz(:,:,2),m%z)
     ! TEST PROBABILITY CHAINS
     call teste_markov(5000000,nz,m%z,m%prbz(:,:,1),m%xsigmaz(1),rhoz,mu1)
     call teste_markov(5000000,nz,m%z,m%prbz(:,:,2),m%xsigmaz(2),rhoz,mu2)  
  end if
  
  if(iflag_pregrida.eq.1)then
     if(na.ne.11) stop 'Fatal error: please check grids for agg shock'
     ! THESE NUMBERS REFLECT THE GRID OF OBSERVED AGGREGATE PROFITABILITY SHOCKS DURING THE CRISIS
     if(iflag_whichsample.eq.1)then
        ! WHOLE SAMPLE
        
        m%a(1 ) = log(0.777d0)
        m%a(2 ) = log(0.794d0)
        m%a(3 ) = log(0.834d0)
        m%a(4 ) = log(0.885d0)
        m%a(5 ) = log(0.934d0)
        m%a(6 ) = log(0.955d0)
        m%a(7 ) = log(1.000d0)
        m%a(8 ) = log(1.010d0)
        m%a(9 ) = log(1.029d0)
        m%a(10) = log(1.041d0)
        m%a(11) = log(1.110d0)


     end if
     call get_markov_vpath(na,m%a,rhoa,sigmaa,m%prba)
     amin = m%a(1)
     amax = m%a(na)
  else
     amax = 3d0*(sigmaa/sqrt(1d0-rhoa**2))
     amin =-amax
     call get_markov(na,amin,amax,rhoa,sigmaa,m%prba,m%a)
  end if

  ! PACK GRIDS FOR SHOCKS
  m%pz = pack(m%prbz,1.eq.1)
  m%pa = pack(m%prba,1.eq.1)
  m%px = pack(m%prbx,1.eq.1)

  ! LOAD GUESSES FOR VFI
  if(iflag_newguess.eq.1)then
     m%vg     = 1000d0
     m%iter   = 0
     m%jflag  = 0
     m%hflag  = 0
  else
     open(1,file='functions/sparam.txt',position="rewind");read(1,*)m%sparam	;close(1)
     open(1,file='functions/tau.txt',	position="rewind");read(1,*)m%tau	;close(1)
     open(1,file='functions/taudiv.txt',position="rewind");read(1,*)m%taudiv	;close(1)
     open(1,file='functions/vg.txt',	position="rewind");read(1,*)m%vg	;close(1)
     open(1,file='functions/ikp.txt',	position="rewind");read(1,*)m%ikp	;close(1)
     open(1,file='functions/ipp.txt',	position="rewind");read(1,*)m%ipp	;close(1)
     open(1,file='functions/iter.txt',	position="rewind");read(1,*)m%iter	;close(1)
     open(1,file='functions/jflag.txt',	position="rewind");read(1,*)m%jflag	;close(1)
     open(1,file='functions/hflag.txt',	position="rewind");read(1,*)m%hflag	;close(1)
  end if

  ! RUN VFI
  write(*,'(a25,3f15.3)')"--> pmin ,pmax  : ",m%p(1) ,m%p(np)
  write(*,'(a25,2f15.3)')"--> kmin ,kmax  : ",m%k(1) ,m%k(nk)
  write(*,'(a25,2f15.3)')"--> zmin ,zmax  : ",m%z(1) ,m%z(nz)
  write(*,'(a25,2f15.3)')"--> amin ,amax  : ",m%a(1), m%a(na)
  print*,""
  print*,"Value function iteration ..."
  print*,""

  call get_equil(par,m,showoutput='y')

  ! END PROGRAM AND CHECK RUN TIME
#ifdef _OPENMP
  t2 = omp_get_wtime()
#else
  call cpu_time(t2)
#endif
  print*,'Program execution time: ',t2-t1
  print*,""

  ! SAVE RESULTS
  open(1,file='functions/k.txt',	position="rewind");write(1,*)m%k	;close(1)
  open(1,file='functions/p.txt',	position="rewind");write(1,*)m%p	;close(1)
  open(1,file='functions/a.txt',	position="rewind");write(1,*)m%a	;close(1)
  open(1,file='functions/z.txt',	position="rewind");write(1,*)m%z	;close(1)
  open(1,file='functions/xsigmaz.txt',	position="rewind");write(1,*)m%xsigmaz	;close(1)
  open(1,file='functions/prba.txt',	position="rewind");write(1,*)m%prba	;close(1)
  open(1,file='functions/prbz.txt',	position="rewind");write(1,*)m%prbz	;close(1)
  open(1,file='functions/prbx.txt',	position="rewind");write(1,*)m%prbx	;close(1)
  open(1,file='functions/pa.txt',	position="rewind");write(1,*)m%pa	;close(1)
  open(1,file='functions/pz.txt',	position="rewind");write(1,*)m%pz	;close(1)
  open(1,file='functions/px.txt',	position="rewind");write(1,*)m%px	;close(1)
  open(1,file='functions/sparam.txt',	position="rewind");write(1,*)m%sparam	;close(1)
  open(1,file='functions/tau.txt',	position="rewind");write(1,*)m%tau	;close(1)
  open(1,file='functions/taudiv.txt',	position="rewind");write(1,*)m%taudiv	;close(1)
  open(1,file='functions/vg.txt',	position="rewind");write(1,*)m%vg	;close(1)
  open(1,file='functions/ikp.txt',	position="rewind");write(1,*)m%ikp	;close(1)
  open(1,file='functions/ipp.txt',	position="rewind");write(1,*)m%ipp	;close(1)
  open(1,file='functions/iter.txt',	position="rewind");write(1,*)m%iter	;close(1)
  open(1,file='functions/jflag.txt',	position="rewind");write(1,*)m%jflag	;close(1)
  open(1,file='functions/hflag.txt',	position="rewind");write(1,*)m%hflag	;close(1)
  open(1,file='functions/delta.txt',	position="rewind");write(1,*)par%delta	;close(1)
  open(1,file='functions/gamma.txt',	position="rewind");write(1,*)par%gamma	;close(1)
  open(1,file='functions/aa.txt',	position="rewind");write(1,*)par%a	;close(1)
  open(1,file='functions/ps.txt',	position="rewind");write(1,*)par%ps	;close(1)
  open(1,file='functions/sparam1.txt',	position="rewind");write(1,*)par%sparam1;close(1)
  open(1,file='functions/sparam2.txt',	position="rewind");write(1,*)par%sparam2;close(1)

  ! DEALLOCATE ARRAYS
  deallocate(m%k,m%p,m%z,m%a,m%xsigmaz)
  deallocate(m%prbz,m%prba,m%prbx)
  deallocate(m%evg)
  deallocate(m%ikp,m%ipp)
  deallocate(m%vg,m%vc,m%vx,m%pz,m%pa,m%px)
  deallocate(m%coefevg)

end program main
