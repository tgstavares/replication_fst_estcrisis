program stats
#ifdef _OPENMP
  use omp_lib
#endif
#ifdef _OPENACC
  use openacc
#endif
  use prec; use parameters; use globals; use process; use equilibrium; use simulation
  implicit none

  type(model) m
  type(param) par
  integer, parameter:: tsimul = 107, tdiscard = 103, nsimul = 20000, nfirms=1000, nentrants=40, nvars = 27
  real(dp), allocatable:: data(:,:,:),shocksz(:,:),shocksa(:),shocksx(:),initconds(:,:)
  integer,  allocatable:: iapath(:),ixpath(:)
  logical, allocatable:: mask(:,:)
  integer, parameter:: nstats = 35
  real statistics(nstats),eb,ev,ek,ebk,eqq,epr
  integer n,t,auxmax(2),auxmin(2),nmax,tmax,nmin,nminv,tmin,iflag_asindata,iloc(1),tt1,tt0

  print*,""
  print*,"Initiate program ..."
  print*,""

#ifdef _OPENMP
  t1 = omp_get_wtime()
#else
  call cpu_time(t1)
#endif

  ! PROGRAM OPTIONS
  iflag_asindata = 1

  ! ALLOCATE ARRAYS FOR SIMULATIONS
  allocate(shocksz(nsimul,tsimul),shocksa(tsimul),shocksx(tsimul))
  allocate(iapath(tsimul),ixpath(tsimul))
  allocate(data(nsimul,tsimul,nvars),initconds(nsimul,8))

  ! ALLOCATE VARIABLES FROM GLOBALS
  allocate(m%k(nk),m%p(np),m%z(nz),m%a(na),m%xsigmaz(nx))
  allocate(m%prbz(nz,nz,nx),m%prba(na,na),m%prbx(nx,nx))
  allocate(m%evg(nk,np,nz,na,nx))
  allocate(m%ikp(nk*np*nz*na*nx),m%ipp(nk*np*nz*na*nx))
  allocate(m%vg(nk*np*nz*na*nx),m%vc(nk*np*nz*na*nx),m%vx(nk*np*nz*na),m%pz(nz*nz*nx),m%pa(na*na),m%px(nx*nx))
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

  ! SIMULATIONS - 1ST SHOCKS; 2ND SIMULATIONS; 3RD STATS (MAY USE A MASK)
  print*,"-> Generate simulation for N firms, T periods, discard: ",nsimul,tsimul,tdiscard
  print*,""
  call get_simulshocks(nsimul,tsimul,shocksz,shocksa,shocksx,randomseed='n')

  if(iflag_asindata.eq.1)then
     iloc = minloc(abs(log(1.0d0) - m%a))
     iapath = iloc(1)
     ixpath = 1
     iloc = minloc(abs(log(1.041d0) - m%a)); iapath(103) = iloc(1)
     iloc = minloc(abs(log(1.029d0) - m%a)); iapath(104) = iloc(1)
     iloc = minloc(abs(log(0.954d0) - m%a)); iapath(105) = iloc(1)
     iloc = minloc(abs(log(0.965d0) - m%a)); iapath(106) = iloc(1)
     iloc = minloc(abs(log(1.010d0) - m%a)); iapath(107) = iloc(1)
     tt1 = 5
     call get_simulations(par,m,nfirms,nsimul,tsimul-tt1,nvars, &
          shocksz(:,1:tsimul-tt1), &
          shocksa(1:tsimul-tt1), &
          shocksx(1:tsimul-tt1), &
          data(:,1:tsimul-tt1,:),iapath=iapath(1:tsimul-tt1),ixpath=ixpath(1:tsimul-tt1),replace_exit='y',endog_entry='y')
     do tt0=103,107
        initconds(:,1) = data(:,tt0-1,6) 	! k0
        initconds(:,2) = data(:,tt0-1,7) 	! p0
        initconds(:,3) = data(:,tt0-1,5) 	! z_1
        initconds(:,4) = data(:,tt0-1,9) 	! t0
        initconds(:,5) = data(:,tt0-1,3) 	! k_1
        initconds(:,6) = data(:,tt0-1,10)	! inv_1
        initconds(:,7) = data(:,tt0-1,2)	! t
        initconds(:,8) = data(:,tt0-1,27)	! mask-cond
        call get_simulations(par,m,nfirms,nsimul,1,nvars,shocksz(:,tt0),shocksa(tt0),shocksx(tt0),data(:,tt0,:),&
             initconds=initconds,iapath=iapath(tt0),ixpath=ixpath(tt0),replace_exit='y',endog_entry='y',massentrants=nentrants)
        !call get_simulations(par,m,nfirms,nsimul,1,nvars,shocksz(:,tt0),shocksa(tt0),shocksx(tt0),data(:,tt0,:),&
        !     initconds=initconds,iapath=iapath(tt0),ixpath=ixpath(tt0),replace_exit='y',endog_entry='y')
     end do
  else
     call get_simulations(par,m,nfirms,nsimul,tsimul,nvars,shocksz,shocksa,shocksx,data)
  end if

  !stop
  
  ! do n=1,nsimul
  !    do t=1,tsimul
  !       write(*,'(3f10.0,3f15.4,2f10.0,3f15.4)')data(n,t,1),data(n,t,2),data(n,t,9),data(n,t,3),data(n,t,6),data(n,t,20),data(n,t,25),data(n,t,26),data(n,t,3),data(n,t,6)
  !    end do
  !    stop
  ! end do
  
  allocate(mask(nsimul,tsimul-tdiscard))
  mask = data(:,tdiscard:tsimul,9).gt.1.5d0
  call get_stats(data(:,tdiscard:tsimul,:),nstats,statistics,showoutput='y',cond=mask)
  open(1,file='data/e_k.txt',position="rewind");write(1,*)statistics(23);close(1)
  open(1,file='data/e_b.txt',position="rewind");write(1,*)statistics(31);close(1)

  ! ADDITIONAL STATS - MOSTLY IMPORTANT MAX AND MIN CAHSFLOW
  print*,''
  auxmax = maxloc(data(:,tdiscard:tsimul,19),mask)
  nmax = auxmax(1)
  tmax = tdiscard + auxmax(2) - 1
  write(*,'(a40,10f12.4)')'max (CF,prft,I,CF/prft,CF/I): ',data(nmax,tmax,19), &
       data(nmax,tmax,18),data(nmax,tmax,10), &
       data(nmax,tmax,19)/data(nmax,tmax,18), &
       data(nmax,tmax,19)/data(nmax,tmax,10)
  auxmin = minloc(data(:,tdiscard:tsimul,19),mask)
  nmin = auxmin(1)
  tmin = tdiscard + auxmin(2) - 1
  write(*,'(a40,10f12.4)')'min (CF,prft,I,CF/prft,CF/I): ',data(nmin,tmin,19), &
       data(nmin,tmin,18),data(nmin,tmin,10), &
       data(nmin,tmin,19)/data(nmin,tmin,18), &
       data(nmin,tmin,19)/data(nmin,tmin,10)
  !print*,maxval(data(:,tdiscard:tsimul,19)),minval(data(:,tdiscard:tsimul,19))

  ! ADDITIONAL STATS - CHECK EXIT FREQUENCY
  print*,''
  write(*,'(a15,i12)')'-> neg div:',count(pack(data(:,tdiscard:tsimul,19),mask).lt.0d0)
  write(*,'(a15,2i12,f12.3)')'-> neg val:',count(pack(data(:,tdiscard:tsimul,20),mask).lt.0d0),count(mask),dble(count(pack(data(:,tdiscard:tsimul,20),mask).lt.0d0))/dble(count(mask))

  ! ADDTIONAL STATS - CHECK BOUNDS AT SIMULATION; MAX MIN OF VARIABLES; AVERAGES
  ek = sum(data(:,tdiscard:tsimul,3 ),mask)/dble(count(mask))
  eb = sum(data(:,tdiscard:tsimul,4 ),mask)/dble(count(mask))
  ebk= sum(data(:,tdiscard:tsimul,15),mask)/dble(count(mask))
  eqq= sum(data(:,tdiscard:tsimul,17),mask)/dble(count(mask))
  ev = sum(data(:,tdiscard:tsimul,20),mask)/dble(count(mask))
  epr= sum(data(:,tdiscard:tsimul,18),mask)/dble(count(mask))
  write(*,'(a15,5f12.4)')'max min k :',maxval(data(:,tdiscard:tsimul,3),mask),minval(data(:,tdiscard:tsimul,3),mask),m%k(nk),m%k(1)
  write(*,'(a15,5f12.4)')'#maxmin k :',dble(count(pack(data(:,tdiscard:tsimul,3),mask).ge.m%k(nk)-erro5))/dble(count(mask)),dble(count(pack(data(:,tdiscard:tsimul,3),mask).le.m%k(1)+erro5))/dble(count(mask))
  write(*,'(a15,5f12.4)')'max min l :',maxval(data(:,tdiscard:tsimul,4),mask),minval(data(:,tdiscard:tsimul,4),mask),m%p(np),m%p(1)
  write(*,'(a15,5f12.4)')'#maxmin l :',dble(count(pack(data(:,tdiscard:tsimul,4),mask).ge.m%p(np)-erro5))/dble(count(mask)),dble(count(pack(data(:,tdiscard:tsimul,4),mask).le.m%p(1)+erro5))/dble(count(mask))
  write(*,'(a15,5f12.4)')'max min p :',maxval(data(:,tdiscard:tsimul,15),mask),minval(data(:,tdiscard:tsimul,15),mask)
  write(*,'(a15,5f12.4)')'max min qq:',maxval(data(:,tdiscard:tsimul,17),mask),minval(data(:,tdiscard:tsimul,17),mask)
  write(*,'(a15,5f12.4)')'max min v :',maxval(data(:,tdiscard:tsimul,20),mask),minval(data(:,tdiscard:tsimul,20),mask)
  write(*,'(a15,5f12.4)')'E(k)      :',ek
  write(*,'(a15,5f12.4)')'E(v)      :',ev
  write(*,'(a15,5f12.4)')'E(l)      :',eb
  write(*,'(a15,5f12.4)')'E(p)      :',ebk
  write(*,'(a15,5f12.4)')'E(qq)     :',eqq
  write(*,'(a15,5f12.4)')'E(pr)     :',epr
  print*,''

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
  deallocate(m%vg,m%vc,m%vx,m%pz,m%pa,m%px)
  deallocate(m%coefevg)

  deallocate(shocksz,shocksa,shocksx)
  deallocate(iapath,ixpath)
  deallocate(data)

end program stats
