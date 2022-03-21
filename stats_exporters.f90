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
  integer, parameter:: tsimul = 107, tdiscard = 103, nsimul = 50000, nvars = 24
  real(dp), allocatable:: data(:,:,:),shocksz(:,:),shocksa(:),shocksx(:)
  integer,  allocatable:: iapath(:),ixpath(:)
  integer, parameter:: nstats = 35
  real statistics(nstats),eb,ev,ek,ebk,eqq,epr
  integer n,t,auxmax(2),auxmin(2),nmax,tmax,nmin,nminv,tmin,iflag_asindata,iflag_savedata,iflag_saveshocks,iflag_useshocks,iloc(1)

  print*,""
  print*,"Initiate program ..."
  print*,""

#ifdef _OPENMP
  t1 = omp_get_wtime()
#else
  call cpu_time(t1)
#endif

  ! PROGRAM OPTIONS
  iflag_asindata   = 1
  iflag_useshocks  = 1
  iflag_saveshocks = 0
  iflag_savedata   = 0

  ! ALLOCATE ARRAYS FOR SIMULATIONS
  allocate(shocksz(nsimul,tsimul),shocksa(tsimul),shocksx(tsimul))
  allocate(iapath(tsimul),ixpath(tsimul))
  allocate(data(nsimul,tsimul,nvars))  

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

  ! SIMULATIONS - 1ST SHOCKS; 2ND SIMULATIONS; 3RD STATS (MAY USE A MASK)
  print*,"-> Generate simulation for N firms, T periods, discard: ",nsimul,tsimul,tdiscard
  print*,""
  if(iflag_useshocks.eq.0)then
     call get_simulshocks(nsimul,tsimul,shocksz,shocksa,shocksx,randomseed='n')
  else
     open(1,file='functions/shocksz_exporters.dat',position="rewind");read(1,*)shocksz;close(1)
     open(1,file='functions/shocksa_exporters.dat',position="rewind");read(1,*)shocksa;close(1)
     open(1,file='functions/shocksx_exporters.dat',position="rewind");read(1,*)shocksx;close(1)
  end if

  if(iflag_asindata.eq.1)then
     if(iflag_useshocks.eq.0)then
        iloc = minloc(abs(0.00d0 - m%a)); iapath = iloc(1)
        ixpath = 1

        ! iloc = minloc(abs( 0.03d0 - m%a)); iapath(100) = iloc(1)
        ! iloc = minloc(abs( 0.06d0 - m%a)); iapath(101) = iloc(1)
        ! iloc = minloc(abs( 0.09d0 - m%a)); iapath(102) = iloc(1)

        iloc = minloc(abs(-0.06d0 - m%a)); iapath(103) = iloc(1)
        iloc = minloc(abs(-0.08d0 - m%a)); iapath(104) = iloc(1)
        iloc = minloc(abs(-0.15d0 - m%a)); iapath(105) = iloc(1)
        iloc = minloc(abs( 0.02d0 - m%a)); iapath(106) = iloc(1)
        iloc = minloc(abs(-0.05d0 - m%a)); iapath(107) = iloc(1)
     else
        open(1,file='functions/iapath_exporters.dat',position="rewind");read(1,*)iapath;close(1)
        open(1,file='functions/ixpath_exporters.dat',position="rewind");read(1,*)ixpath;close(1)
     end if
     call get_simulations(par,m,nsimul,tsimul,nvars,shocksz,shocksa,shocksx,data,iapath=iapath,ixpath=ixpath)
  else
     call get_simulations(par,m,nsimul,tsimul,nvars,shocksz,shocksa,shocksx,data)
  end if

  if(iflag_saveshocks.eq.1)then
     open(1,file='functions/shocksz_exporters.dat',position="rewind");write(1,*)shocksz;close(1)
     open(1,file='functions/shocksa_exporters.dat',position="rewind");write(1,*)shocksa;close(1)
     open(1,file='functions/shocksx_exporters.dat',position="rewind");write(1,*)shocksx;close(1)
     open(1,file='functions/iapath_exporters.dat',position="rewind");write(1,*)iapath;close(1)
     open(1,file='functions/ixpath_exporters.dat',position="rewind");write(1,*)ixpath;close(1)
  end if

  ! n=1
  ! do t=1,tsimul
  !    write(*,'(10f15.4)') &
  !         data(n,t,1),data(n,t,2),data(n,t,5), &
  !         data(n,t,21),data(n,t,22),data(n,t,23)
  ! end do
  ! stop

  call get_stats(data(:,tdiscard:tsimul,:),nstats,statistics,showoutput='y')
  open(1,file='data/e_k.txt',position="rewind");write(1,*)statistics(23);close(1)
  open(1,file='data/e_b.txt',position="rewind");write(1,*)statistics(31);close(1)

  ! ADDITIONAL STATS - MOSTLY IMPORTANT MAX AND MIN CAHSFLOW
  print*,''
  auxmax = maxloc(data(:,tdiscard:tsimul,19))
  nmax = auxmax(1)
  tmax = tdiscard + auxmax(2) - 1
  write(*,'(a40,10f12.4)')'max (CF,prft,I,CF/prft,CF/I): ',data(nmax,tmax,19), &
       data(nmax,tmax,18),data(nmax,tmax,10), &
       data(nmax,tmax,19)/data(nmax,tmax,18), &
       data(nmax,tmax,19)/data(nmax,tmax,10)
  auxmin = minloc(data(:,tdiscard:tsimul,19))
  nmin = auxmin(1)
  tmin = tdiscard + auxmin(2) - 1
  write(*,'(a40,10f12.4)')'min (CF,prft,I,CF/prft,CF/I): ',data(nmin,tmin,19), &
       data(nmin,tmin,18),data(nmin,tmin,10), &
       data(nmin,tmin,19)/data(nmin,tmin,18), &
       data(nmin,tmin,19)/data(nmin,tmin,10)
  !print*,maxval(data(:,tdiscard:tsimul,19)),minval(data(:,tdiscard:tsimul,19))

  ! ADDITIONAL STATS - CHECK EXIT FREQUENCY
  print*,''
  write(*,'(a15,i12)')'-> neg div:',count(data(:,tdiscard:tsimul,19).lt.0.0d0)
  write(*,'(a15,2i12,f12.3)')'-> neg val:',count(data(:,tdiscard:tsimul,9 ).lt.1.5d0),nsimul*(tsimul-tdiscard+1),dble(count(data(:,tdiscard:tsimul,9 ).lt.1.5d0))/dble(nsimul*(tsimul-tdiscard+1))

  ! ADDTIONAL STATS - CHECK BOUNDS AT SIMULATION; MAX MIN OF VARIABLES; AVERAGES
  ek = sum(data(:,tdiscard:tsimul,3))/dble(nsimul*(tsimul-tdiscard+1))
  eb = sum(data(:,tdiscard:tsimul,4))/dble(nsimul*(tsimul-tdiscard+1))
  ebk= sum(data(:,tdiscard:tsimul,15))/dble(nsimul*(tsimul-tdiscard+1))
  eqq= sum(data(:,tdiscard:tsimul,17))/dble(nsimul*(tsimul-tdiscard+1))
  ev = sum(data(:,tdiscard:tsimul,20))/dble(nsimul*(tsimul-tdiscard+1))
  epr= sum(data(:,tdiscard:tsimul,18))/dble(nsimul*(tsimul-tdiscard+1))
  write(*,'(a15,5f12.4)')'max min k :',maxval(data(:,tdiscard:tsimul,3)),minval(data(:,tdiscard:tsimul,3)),m%k(nk),m%k(1)
  write(*,'(a15,5f12.4)')'#maxmin k :',dble(count(data(:,tdiscard:tsimul,3).ge.m%k(nk)-erro5))/dble(nsimul*(tsimul-tdiscard+1)),dble(count(data(:,tdiscard:tsimul,3).le.m%k(1)+erro5))/dble(nsimul*(tsimul-tdiscard+1))
  write(*,'(a15,5f12.4)')'max min l :',maxval(data(:,tdiscard:tsimul,4)),minval(data(:,tdiscard:tsimul,4)),m%p(np),m%p(1)
  write(*,'(a15,5f12.4)')'#maxmin l :',dble(count(data(:,tdiscard:tsimul,4).ge.m%p(np)-erro5))/dble(nsimul*(tsimul-tdiscard+1)),dble(count(data(:,tdiscard:tsimul,4).le.m%p(1)+erro5))/dble(nsimul*(tsimul-tdiscard+1))
  write(*,'(a15,5f12.4)')'max min p :',maxval(data(:,tdiscard:tsimul,15)),minval(data(:,tdiscard:tsimul,15))
  write(*,'(a15,5f12.4)')'max min qq:',maxval(data(:,tdiscard:tsimul,17)),minval(data(:,tdiscard:tsimul,17))
  write(*,'(a15,5f12.4)')'max min v :',maxval(data(:,tdiscard:tsimul,20)),minval(data(:,tdiscard:tsimul,20))
  write(*,'(a15,5f12.4)')'E(k)      :',ek
  write(*,'(a15,5f12.4)')'E(v)      :',ev
  write(*,'(a15,5f12.4)')'E(l)      :',eb
  write(*,'(a15,5f12.4)')'E(p)      :',ebk
  write(*,'(a15,5f12.4)')'E(qq)     :',eqq
  write(*,'(a15,5f12.4)')'E(pr)     :',epr

  if(iflag_savedata.eq.1)then
     open(1,file='data/simulations.txt',position="rewind")
     do n=1,nsimul
        do t=tdiscard-20,tsimul
           write(1,'(10f15.6)') data(n,t,1),data(n,t,2),data(n,t,3),data(n,t,4),data(n,t,5),data(n,t,21)
        end do
     end do
  end if

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

  deallocate(shocksz,shocksa,shocksx)
  deallocate(iapath,ixpath)
  deallocate(data)

end program stats
