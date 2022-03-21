program smm_estimation_par
  use prec; use parameters; use globals; use process; use equilibrium; use simulation
  use globals_montecarlo; use sort_tiago
  implicit none
  include 'mpif.h'

  type(model) m
  type(param) par
  integer, parameter:: tsimul = 107, tdiscard = 103, nsimul = 50000, nvars = 24, niiter = 0, nmomentsall = 35
  real(dp), allocatable:: data(:,:,:),shocksz(:,:),shocksa(:),shocksx(:)
  integer,  allocatable:: iapath(:),ixpath(:)
  integer i,j,idl,iaa,ips,is1,is2,iloc(1),ik0,ip0
  real(dp) moments_data(nmoments,1),wmatrix(nmoments,nmoments),wmatrix4(4,4),identity(nmoments,nmoments),identity3(4,4),&
       qdist(ndl*naa*nps*ns1*ns2),qdist2(ndl*naa*nps*ns1*ns2),qdist3(ndl*naa*nps*ns1*ns2),qdist4(ndl*naa*nps*ns1*ns2)
  integer, dimension(ndl*naa*nps*ns1*ns2)::idx,xdl,xgm,xaa,xps,xs1,xs2
  real moments(nmomentsall)

  !    MPI
  integer status(MPI_STATUS_SIZE)
  integer numtasks, rankk, ierr, resultlen, ndst,tag,resto
  integer tempo1,tempo2,clock_rate,clock_max
  character (len=20) :: name

  ! CHECK MPI  
  call MPI_INIT(ierr)
  call SYSTEM_CLOCK(tempo1,clock_rate,clock_max)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rankk, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
  call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)

  if(rankk.eq.0) print*,"" !MPI
  if(rankk.eq.0) print*,"Prepared to run with:" !MPI
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call SYSTEM_CLOCK(tempo1,clock_rate,clock_max)
  do i=1,numtasks
     if(rankk.eq.i-1) write(*,'(a10,a30,a10,i5,a10,i5)')'computer ',name,'rank ',rankk,'of ',numtasks
     !if(rankk.eq.i-1) print*,name,rankk,numtasks
  end do
  !write(*,'(a10,a30,a10,i5,a10,i5)')'computer ',name,'rank ',rankk,'of ',numtasks
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  if(rankk.eq.0) print*,"" !MPI

  ! TRY TO CALIBRATE AVAILABLE CORES TO SPECIFIC PARALLEL PROGRAM REQUIREMENTS
  call MPI_BARRIER(MPI_COMM_WORLD, ierr) ! MPI
  if(mod(ndl*naa*nps*ns1*ns2,numtasks).eq.0)then
     ! equal split of jobs on processes
     ndst = max(ndl*naa*nps*ns1*ns2/numtasks,1)
  else
     ndst = max(floor(real(ndl*naa*nps*ns1*ns2)/real(numtasks-1)),1)
     resto= ndl*naa*nps*ns1*ns2-ndst*(numtasks-1)
     if(resto.gt.ndst) then
        if(numtasks.gt.ndst)then
           stop "Erro - numtasks.gt.ndst - adjsut numtasks"
        end if
     end if
  end if

  !call MPI_BARRIER(MPI_COMM_WORLD, ierr) ! MPI
  if(rankk.eq.0)then
     print*,""
     write(*,'(a60,2i5)')"Number of iters per thread / total number of threads:",numtasks,ndst
     print*,""
     !stop
  end if

  ! BARRIER HERE!!
  call MPI_BARRIER(MPI_COMM_WORLD, ierr) ! MPI

  ! SET COMMON ELEMENTS
  allocate(shocksz(nsimul,tsimul),shocksa(tsimul),shocksx(tsimul))
  allocate(iapath(tsimul),ixpath(tsimul))
  allocate(data(nsimul,tsimul,nvars))
  allocate(m%k(nk),m%p(np),m%z(nz),m%a(na),m%xsigmaz(nx))
  allocate(m%prbz(nz,nz,nx),m%prba(na,na),m%prbx(nx,nx))
  allocate(m%evg(nk,np,nz,na,nx))
  allocate(m%ikp(nk*np*nz*na*nx),m%ipp(nk*np*nz*na*nx))
  allocate(m%vg(nk*np*nz*na*nx),m%vc(nk*np*nz*na*nx),m%pz(nz*nz*nx),m%pa(na*na),m%px(nx*nx))
  allocate(m%coefevg(4,4,nk,np,nz,na,nx))

  ! ! SET UP GRIDS
  ! do ip0 = 1,np
  !    m%p(ip0) = pmin + dble(ip0-1)*(pmax - pmin)/dble(np-1)
  ! end do
  ! if(nx.eq.1)then
  !    m%xsigmaz(1) = sigmaz
  !    m%prbx(1,1)  = 1d0
  !    zmax	  = nsd*sigmaz/sqrt(1d0-rhoz**2d0)
  !    zmin         =-zmax
  !    call get_markov(nz,zmin,zmax,rhoz,sigmaz,m%prbz,m%z)
  ! else
  !    m%xsigmaz   = (/sigmaz , sjump*sigmaz/)
  !    m%prbx(1,:) = (/1.0d0-uncfreq, uncfreq/)!; m%prbx(1,:) = m%prbx(1,:)/sum(m%prbx(1,:))
  !    m%prbx(2,:) = (/1.0d0-uncpers, uncpers/)!; m%prbx(2,:) = m%prbx(2,:)/sum(m%prbx(2,:))
  !    zmax	  = nsd*sjump*sigmaz/sqrt(1d0-rhoz**2d0)
  !    zmin         =-zmax
  !    call get_markovforts(nz,mu1,mu2,m%xsigmaz(1),m%xsigmaz(2),rhoz,rhoz,zmin,zmax,m%prbz(:,:,1),m%prbz(:,:,2),m%z)
  ! end if

  ! if(na.ne.11) stop 'Fatal error: please check grids for agg shock'
  ! ! THESE NUMBERS REFLECT THE GRID OF OBSERVED AGGREGATE PROFITABILITY SHOCKS DURING THE CRISIS
  ! m%a(1 ) =-0.18d0
  ! m%a(2 ) =-0.16d0
  ! m%a(3 ) =-0.11d0
  ! m%a(4 ) =-0.05d0
  ! m%a(5 ) = 0.00d0
  ! m%a(6 ) = 0.02d0
  ! m%a(7 ) = 0.03d0
  ! m%a(8 ) = 0.08d0
  ! m%a(9 ) = 0.10d0
  ! m%a(10) = 0.11d0
  ! m%a(11) = 0.18d0

  ! call get_markov_vpath(na,m%a,rhoa,sigmaa,m%prba)
  ! amin = m%a(1)
  ! amax = m%a(na)

  ! ! PACK GRIDS FOR SHOCKS
  ! m%pz = pack(m%prbz,1.eq.1)
  ! m%pa = pack(m%prba,1.eq.1)
  ! m%px = pack(m%prbx,1.eq.1)

  ! ! LOAD GUESSES FOR VFI
  ! m%sparam = 1d0
  ! m%tau    = 0.31d0
  ! m%taudiv = 0.00d0
  ! m%vg     = 1000d0
  ! m%iter   = 0
  ! m%jflag  = 0
  ! m%hflag  = 0

  if(rankk.eq.0)then
     m%prbx(1,:) = (/1.0d0-uncfreq, uncfreq/); m%prbx(1,:) = m%prbx(1,:)/sum(m%prbx(1,:)) 	! NOT SURE EXACTLY WHY WE NEED THESE TWO LINES OF CODE. 
     m%prbx(2,:) = (/1.0d0-uncpers, uncpers/); m%prbx(2,:) = m%prbx(2,:)/sum(m%prbx(2,:)) 	! IF THEY'RE NOT INCLUDED THE ACC VERSION BREAKS DOWN
     open(1,file='functions/k.txt',		position="rewind");read(1,*)m%k		;close(1)
     open(1,file='functions/p.txt',		position="rewind");read(1,*)m%p		;close(1)
     open(1,file='functions/a.txt',		position="rewind");read(1,*)m%a		;close(1)
     open(1,file='functions/z.txt',		position="rewind");read(1,*)m%z		;close(1)
     open(1,file='functions/xsigmaz.txt',	position="rewind");read(1,*)m%xsigmaz	;close(1)
     open(1,file='functions/prba.txt',		position="rewind");read(1,*)m%prba	;close(1)
     open(1,file='functions/prbz.txt',		position="rewind");read(1,*)m%prbz	;close(1)
     open(1,file='functions/prbx.txt',		position="rewind");read(1,*)m%prbx	;close(1)
     open(1,file='functions/pa.txt',		position="rewind");read(1,*)m%pa	;close(1)
     open(1,file='functions/pz.txt',		position="rewind");read(1,*)m%pz	;close(1)
     open(1,file='functions/px.txt',		position="rewind");read(1,*)m%px	;close(1)
     open(1,file='functions/sparam.txt',	position="rewind");read(1,*)m%sparam	;close(1)
     open(1,file='functions/tau.txt',		position="rewind");read(1,*)m%tau	;close(1)
     open(1,file='functions/taudiv.txt',	position="rewind");read(1,*)m%taudiv	;close(1)
     open(1,file='functions/vg.txt',		position="rewind");read(1,*)m%vg	;close(1)
     open(1,file='functions/ikp.txt',		position="rewind");read(1,*)m%ikp	;close(1)
     open(1,file='functions/ipp.txt',		position="rewind");read(1,*)m%ipp	;close(1)
     open(1,file='functions/iter.txt',		position="rewind");read(1,*)m%iter	;close(1)
     open(1,file='functions/jflag.txt',		position="rewind");read(1,*)m%jflag	;close(1)
     open(1,file='functions/hflag.txt',		position="rewind");read(1,*)m%hflag	;close(1)
     open(1,file='functions/delta.txt',		position="rewind");read(1,*)par%delta	;close(1)
     open(1,file='functions/gamma.txt',		position="rewind");read(1,*)par%gamma	;close(1)
     open(1,file='functions/aa.txt',		position="rewind");read(1,*)par%a	;close(1)
     open(1,file='functions/ps.txt',		position="rewind");read(1,*)par%ps	;close(1)
     open(1,file='functions/sparam1.txt',	position="rewind");read(1,*)par%sparam1	;close(1)
     open(1,file='functions/sparam2.txt',	position="rewind");read(1,*)par%sparam2	;close(1)

     open(1,file='functions/shocksz.dat',position="rewind");read(1,*)shocksz;close(1)
     open(1,file='functions/shocksa.dat',position="rewind");read(1,*)shocksa;close(1)
     open(1,file='functions/shocksx.dat',position="rewind");read(1,*)shocksx;close(1)
     open(1,file='functions/iapath.dat',position="rewind");read(1,*)iapath;close(1)
     open(1,file='functions/ixpath.dat',position="rewind");read(1,*)ixpath;close(1)
  end if
  call MPI_BARRIER(MPI_COMM_WORLD, ierr) ! MPI  
  call MPI_BCAST(m%k,nk,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%p,np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%a,na,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%z,nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%xsigmaz,nx,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%prba,na*na,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%prbz,nz*nz*nx,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%prbx,nx*nx,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%pa,na*na,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%pz,nz*nz*nx,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%px,nx*nx,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%sparam,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%tau,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%taudiv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%vg,nk*np*nz*na*nx,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%ikp,nk*np*nz*na*nx,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%ipp,nk*np*nz*na*nx,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%iter,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%jflag,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(m%hflag,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(par%delta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(par%gamma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(par%a,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(par%ps,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(par%sparam1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(par%sparam2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(shocksz,nsimul*tsimul,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(shocksa,tsimul,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(shocksx,tsimul,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(iapath,tsimul,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ixpath,tsimul,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)

  ! call MPI_BARRIER(MPI_COMM_WORLD, ierr) ! MPI
  ! print*,m%prbx(2,1),m%ikp(10),par%a,iapath(103),rankk
  ! call MPI_BARRIER(MPI_COMM_WORLD, ierr) ! MPI
  ! if(rankk.eq.0) stop

  ! BARRIER HERE!!
  call MPI_BARRIER(MPI_COMM_WORLD, ierr) ! MPI

  ! MOMENTS USED:
  moments_data(1,1) = 0.115539211719137d0 	! Average investment
  moments_data(2,1) = 0.195777879761358d0	! Inaction at 1pc
  moments_data(3,1) = 0.166590178981184d0 	! Positive spike 20pc ik    
  moments_data(4,1) = 0.061220743460303d0 	! Negative invest 1pc ik
  moments_data(5,1) = 0.166431092880951d0	! Stdev of invest rate (trim with 1%)
  moments_data(6,1) = 0.590442311583295d0	! Average leveraeg
  moments_data(7,1) = 0.0122315903053963d0	! Slope coefficient

  ! WEIGHTING MATRIX:
  wmatrix(1,:) = (/132.49976989236214,9.621792258716074,-39.04848615173803,22.620860720160156,-38.153122086794085,0.36648876433001704,1.520158091539538/)
  wmatrix(2,:) = (/9.621792258716075,5.780141695537245,-1.321127377807546,2.5478324319199954,-2.8577694433972205,-0.12974007172222177,0.023165335647234594/)
  wmatrix(3,:) = (/-39.04848615173804,-1.321127377807545,18.98727667103684,-5.047201139249203,8.060736862504331,-0.865040656042619,-0.3901496847594708/)
  wmatrix(4,:) = (/22.620860720160156,2.5478324319199945,-5.047201139249202,16.935333335412203,-7.478189364937309,-0.18463637718518008,0.31199055985147445/)
  wmatrix(5,:) = (/-38.15312208679409,-2.85776944339722,8.06073686250434,-7.478189364937315,21.12201413205612,-0.22397048074930503,0.008334433634992497/)
  wmatrix(6,:) = (/0.3664887643300148,-0.12974007172222193,-0.8650406560426185,-0.18463637718518047,-0.2239704807493042,7.424914146355202,0.3730340411106479/)
  wmatrix(7,:) = (/1.5201580915395365,0.023165335647234508,-0.39014968475947054,0.3119905598514742,0.008334433634992977,0.3730340411106479,8.934291834950383/)

  wmatrix4(1,:) = (/ 5.0774,    1.4994,    0.9005,   -0.0856/)
  wmatrix4(2,:) = (/ 1.4994,    7.4012,    1.5934,   -3.2001/)
  wmatrix4(3,:) = (/ 0.9005,    1.5934,   13.0647,   -0.9718/)
  wmatrix4(4,:) = (/-0.0856,   -3.2001,   -0.9718,   10.1110/)

  identity = 0d0
  do i=1,nmoments
     identity(i,i) = 1d0
  end do
  identity3 = 0d0
  do i=1,4
     identity3(i,i) = 1d0
  end do

  ! GRID FOR THE SEARCH

  smm%min_dl = 0.0900d0;	smm%max_dl = smm%min_dl+0.04d0
  smm%min_aa = 0.3000d0;	smm%max_aa = smm%min_aa+0.50d0
  smm%min_ps = 0.0010d0;	smm%max_ps = smm%min_ps+0.05d0
  smm%min_s1 = 0.5000d0;	smm%max_s1 = smm%min_s1+0.30d0
  smm%min_s2 = 0.5000d0;	smm%max_s2 = smm%min_s2+1.50d0

  call get_eqspace(ndl,smm%min_dl,smm%max_dl,smm%grid_dl)
  !call get_eqspace(ngm,smm%min_gm,smm%max_gm,smm%grid_gm)
  call get_eqspace(naa,smm%min_aa,smm%max_aa,smm%grid_aa)
  call get_eqspace(nps,smm%min_ps,smm%max_ps,smm%grid_ps)
  call get_eqspace(ns1,smm%min_s1,smm%max_s1,smm%grid_s1)
  call get_eqspace(ns2,smm%min_s2,smm%max_s2,smm%grid_s2)

  ! ! SHOCK SEQUENCE
  ! call get_simulshocks(nsimul,tsimul,shocksz,shocksa,shocksx,randomseed='n')
  ! iloc = minloc(abs(0.00d0 - m%a)); iapath = iloc(1)
  ! ixpath = 1
  ! iloc = minloc(abs(0.11d0 - m%a)); iapath(103) = iloc(1)
  ! iloc = minloc(abs(0.10d0 - m%a)); iapath(104) = iloc(1)
  ! iloc = minloc(abs(0.02d0 - m%a)); iapath(105) = iloc(1)
  ! iloc = minloc(abs(0.03d0 - m%a)); iapath(106) = iloc(1)
  ! iloc = minloc(abs(0.08d0 - m%a)); iapath(107) = iloc(1)

  ! BARRIER HERE!!
  call MPI_BARRIER(MPI_COMM_WORLD, ierr) ! MPI

  ! KEEP TRACK OF INDEXES
  do is2=1,ns2
     do is1=1,ns1
        do ips=1,nps
           do iaa=1,naa
              do idl=1,ndl
                 j = idl + ndl*(iaa-1 + naa*(ips-1 + nps*(is1-1 + ns1*(is2-1))))
                 xdl(j) = idl
                 !xgm(j) = igm
                 xaa(j) = iaa
                 xps(j) = ips
                 xs1(j) = is1
                 xs2(j) = is2
              end do
           end do
        end do
     end do
  end do

  ! BARRIER HERE!!
  call MPI_BARRIER(MPI_COMM_WORLD, ierr) ! MPI

  if(rankk.eq.0)then
     write(*,'(a40,a2,7f8.4,a2,a7,a2)')'','|',moments_data(:,1),'|'
     write(*,'(5a8,a2,7a8,a2,4a8,a4,a2,a8,a2,4a9,a2,a5)')&
          'dlt','a','ps','s1','s2','|', &
          'e.i','ina','p-spk','n-i','sd.i','e.l','slp','|', &
          'dist1','ecd1','dist2','ecld2','itr','|', &
          'sd.l','|', &
          'smaxk','smink','smaxp','sminp','|','rankk'
  end if

  ! BARRIER HERE!!
  call MPI_BARRIER(MPI_COMM_WORLD, ierr) ! MPI
  do j=1+ndst*rankk,ndst*(rankk+1) ! MPI
     if(j.le.ndl*naa*nps*ns1*ns2)then

        par%delta   = smm%grid_dl(xdl(j))
        par%a       = smm%grid_aa(xaa(j))
        par%ps      = smm%grid_ps(xps(j))
        par%sparam1 = smm%grid_s1(xs1(j))
        par%sparam2 = smm%grid_s2(xs2(j))                    

        do ik0 = nk,1,-1
           m%k(ik0) = kmax*(1d0-par%delta)**(dble(nk-ik0)/dble(nk_divs))
        end do
        call get_smm(par,qdist(j),qdist2(j),qdist3(j),qdist4(j),moments)

        write(*,'(5f8.4,a2,7f8.4,a2,4f8.4,i4,a2,f8.4,a2,4f9.5,a2,i5)') &
             par%delta,par%a,par%ps,par%sparam1,par%sparam2,'|', &
             moments(1),moments(9),moments(10),moments(8),moments(13),moments(4),moments(32),'|', &
             qdist(j),qdist2(j),qdist3(j),qdist4(j),m%iter-niiter,'|', &
             moments(16),'|', &
             dble(count(data(:,tdiscard:tsimul,3).ge.m%k(nk)-erro5))/dble(nsimul*(tsimul-tdiscard+1)), &
             dble(count(data(:,tdiscard:tsimul,3).le.m%k( 1)+erro5))/dble(nsimul*(tsimul-tdiscard+1)), &
             dble(count(data(:,tdiscard:tsimul,4).ge.m%p(np)-erro5))/dble(nsimul*(tsimul-tdiscard+1)), &
             dble(count(data(:,tdiscard:tsimul,4).le.m%p( 1)+erro5))/dble(nsimul*(tsimul-tdiscard+1)), &
             '|',rankk
     end if
  end do

  ! BARRIER HERE AND A TAG
  tag = 1234 ! MPI
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  if(numtasks.gt.1)then
     do i=1,numtasks-2
        if(rankk .eq. i)then
           Call MPI_Send( qdist(1+i*ndst:i*ndst+ndst),ndst,MPI_DOUBLE_PRECISION,0,tag, MPI_COMM_WORLD, ierr)
           Call MPI_Send(qdist2(1+i*ndst:i*ndst+ndst),ndst,MPI_DOUBLE_PRECISION,0,tag, MPI_COMM_WORLD, ierr)
           Call MPI_Send(qdist3(1+i*ndst:i*ndst+ndst),ndst,MPI_DOUBLE_PRECISION,0,tag, MPI_COMM_WORLD, ierr)
           Call MPI_Send(qdist4(1+i*ndst:i*ndst+ndst),ndst,MPI_DOUBLE_PRECISION,0,tag, MPI_COMM_WORLD, ierr)
        endif
        if(rankk .eq. 0)then
           Call MPI_Recv( qdist(1+i*ndst:i*ndst+ndst),ndst,MPI_DOUBLE_PRECISION,i,tag, MPI_COMM_WORLD, status,ierr)
           Call MPI_Recv(qdist2(1+i*ndst:i*ndst+ndst),ndst,MPI_DOUBLE_PRECISION,i,tag, MPI_COMM_WORLD, status,ierr)
           Call MPI_Recv(qdist3(1+i*ndst:i*ndst+ndst),ndst,MPI_DOUBLE_PRECISION,i,tag, MPI_COMM_WORLD, status,ierr)
           Call MPI_Recv(qdist4(1+i*ndst:i*ndst+ndst),ndst,MPI_DOUBLE_PRECISION,i,tag, MPI_COMM_WORLD, status,ierr)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     end do
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     i=numtasks-1
     j=ndl*naa*nps*ns1*ns2-1+i*ndst + 1
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     if(rankk .eq. i)then
        Call MPI_Send( qdist(1+i*ndst:ndl*naa*nps*ns1*ns2),ndst,MPI_DOUBLE_PRECISION,0,tag, MPI_COMM_WORLD, ierr)
        Call MPI_Send(qdist2(1+i*ndst:ndl*naa*nps*ns1*ns2),ndst,MPI_DOUBLE_PRECISION,0,tag, MPI_COMM_WORLD, ierr)
        Call MPI_Send(qdist3(1+i*ndst:ndl*naa*nps*ns1*ns2),ndst,MPI_DOUBLE_PRECISION,0,tag, MPI_COMM_WORLD, ierr)
        Call MPI_Send(qdist4(1+i*ndst:ndl*naa*nps*ns1*ns2),ndst,MPI_DOUBLE_PRECISION,0,tag, MPI_COMM_WORLD, ierr)
     endif
     if(rankk .eq. 0)then
        Call MPI_Recv( qdist(1+i*ndst:ndl*naa*nps*ns1*ns2),ndst,MPI_DOUBLE_PRECISION,i,tag, MPI_COMM_WORLD, status,ierr)
        Call MPI_Recv(qdist2(1+i*ndst:ndl*naa*nps*ns1*ns2),ndst,MPI_DOUBLE_PRECISION,i,tag, MPI_COMM_WORLD, status,ierr)
        Call MPI_Recv(qdist3(1+i*ndst:ndl*naa*nps*ns1*ns2),ndst,MPI_DOUBLE_PRECISION,i,tag, MPI_COMM_WORLD, status,ierr)
        Call MPI_Recv(qdist4(1+i*ndst:ndl*naa*nps*ns1*ns2),ndst,MPI_DOUBLE_PRECISION,i,tag, MPI_COMM_WORLD, status,ierr)
     endif
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end if

  ! BARRIER HERE
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  if(rankk.eq.0)then


     print*,''
     print*,'Sort results on Eculidian 2'
     call qsort_dble(ndl*naa*nps*ns1*ns2,qdist4,idx)
     do j=1,ndl*naa*nps*ns1*ns2
        write(*,'(5f10.5,a2,f8.4)') &
             smm%grid_dl(xdl(idx(j))),&
             smm%grid_aa(xaa(idx(j))),&
             smm%grid_ps(xps(idx(j))),&
             smm%grid_s1(xs1(idx(j))),&
             smm%grid_s2(xs2(idx(j))),'|',&
             qdist4(idx(j))
     end do

     print*,''
     print*,'Sort results on Distance 2'
     call qsort_dble(ndl*naa*nps*ns1*ns2,qdist3,idx)
     do j=1,ndl*naa*nps*ns1*ns2
        write(*,'(5f10.5,a2,f8.4)') &
             smm%grid_dl(xdl(idx(j))),&
             smm%grid_aa(xaa(idx(j))),&
             smm%grid_ps(xps(idx(j))),&
             smm%grid_s1(xs1(idx(j))),&
             smm%grid_s2(xs2(idx(j))),'|',&
             qdist3(idx(j))
     end do

     print*,''
     print*,'Sort results on Eculidian'
     call qsort_dble(ndl*naa*nps*ns1*ns2,qdist2,idx)
     open(1,file='data/res_ecl.txt',position="rewind")
     do j=1,ndl*naa*nps*ns1*ns2
        write(*,'(5f10.5,a2,f8.4)') &
             smm%grid_dl(xdl(idx(j))),&
             smm%grid_aa(xaa(idx(j))),&
             smm%grid_ps(xps(idx(j))),&
             smm%grid_s1(xs1(idx(j))),&
             smm%grid_s2(xs2(idx(j))),'|',&
             qdist2(idx(j))
        write(1,'(5f10.5,a2,f8.4)') &
             smm%grid_dl(xdl(idx(j))),&
             smm%grid_aa(xaa(idx(j))),&
             smm%grid_ps(xps(idx(j))),&
             smm%grid_s1(xs1(idx(j))),&
             smm%grid_s2(xs2(idx(j))),'|',&
             qdist2(idx(j))
     end do
     close(1)

     print*,''
     print*,'Sort results on Distance'
     call qsort_dble(ndl*naa*nps*ns1*ns2,qdist,idx)
     open(1,file='data/res_dst.txt',position="rewind")
     do j=1,ndl*naa*nps*ns1*ns2
        write(*,'(5f10.5,a2,f8.4)') &
             smm%grid_dl(xdl(idx(j))),&
             smm%grid_aa(xaa(idx(j))),&
             smm%grid_ps(xps(idx(j))),&
             smm%grid_s1(xs1(idx(j))),&
             smm%grid_s2(xs2(idx(j))),'|',&
             qdist(idx(j))
        write(1,'(5f10.5,a2,f8.4)') &
             smm%grid_dl(xdl(idx(j))),&
             smm%grid_aa(xaa(idx(j))),&
             smm%grid_ps(xps(idx(j))),&
             smm%grid_s1(xs1(idx(j))),&
             smm%grid_s2(xs2(idx(j))),'|',&
             qdist(idx(j))
     end do
     close(1)
  end if

  deallocate(shocksz,shocksa,shocksx)
  deallocate(iapath,ixpath)
  deallocate(data)
  deallocate(m%k,m%p,m%z,m%a,m%xsigmaz)
  deallocate(m%prbz,m%prba,m%prbx)
  deallocate(m%evg)
  deallocate(m%ikp,m%ipp)
  deallocate(m%vg,m%vc,m%pz,m%pa,m%px)
  deallocate(m%coefevg)

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  if(rankk.eq.0) then
     call SYSTEM_CLOCK(tempo2,clock_rate,clock_max)
     print*,''
     print*,'Total execution time',real(tempo2-tempo1)/real(clock_rate)
  end if

  call MPI_FINALIZE(ierr)

contains

  subroutine get_smm(par,dist,dist2,dist3,dist4,simul_moments)
    use linear_op
    implicit none
    type(param), intent(in):: par
    real(dp), intent(out)::dist,dist2,dist3,dist4
    real, intent(out):: simul_moments(nmomentsall)
    real(dp) distmat(1,1)
    real(dp) moments_model(nmoments,1),selmom_data(4,1),selmom_model(4,1)
    integer, parameter:: nstats = nmomentsall
    real stats(nstats)

    m%iter = niiter
    call get_equil(par,m,showoutput='n')
    if(m%iter.gt.nvfi) print*,"NO CONVERGENCE"
    call get_simulations(par,m,nsimul,tsimul,nvars,shocksz,shocksa,shocksx,data,iapath=iapath,ixpath=ixpath)
    call get_stats(data(:,tdiscard:tsimul,:),nstats,stats,showoutput='n')

    moments_model(1,1) = real(stats(1 ),8)	! Average investment
    moments_model(2,1) = real(stats(9 ),8)	! Inaction at 1pc
    moments_model(3,1) = real(stats(10),8)	! Positive spike 20pc ik
    moments_model(4,1) = real(stats(8 ),8)	! Negative invest 1pc ik    
    moments_model(5,1) = real(stats(13),8)	! Stdev of invest rate (trim with 1%)
    moments_model(6,1) = real(stats(4 ),8)	! Average leverage
    moments_model(7,1) = real(stats(32),8)	! Slope coefficient
    !moments_model(4,1) = real(stats(11),8)
    !moments_model(5,1) = real(stats(7 ),8)    
    distmat = tr(moments_model-moments_data).x.(wmatrix).x.(moments_model-moments_data)
    dist    = distmat(1,1)
    distmat = tr(moments_model-moments_data).x.(identity).x.(moments_model-moments_data)
    dist2   = sqrt(distmat(1,1)/dble(nmoments))

    selmom_data(1,1) = moments_data(2,1)	! Inaction at 1pc
    selmom_data(2,1) = moments_data(3,1)	! Positive spike 20pc ik
    selmom_data(3,1) = moments_data(4,1)	! Negative invest 1pc ik
    selmom_data(4,1) = moments_data(5,1)	! Stdev of invest rate (trim with 1%)
    selmom_model(1,1) = moments_model(2,1)
    selmom_model(2,1) = moments_model(3,1)
    selmom_model(3,1) = moments_model(4,1)
    selmom_model(4,1) = moments_model(5,1)
    distmat = tr(selmom_model-selmom_data).x.(wmatrix4).x.(selmom_model-selmom_data)
    dist3 = distmat(1,1)
    distmat = tr(selmom_model-selmom_data).x.(identity3).x.(selmom_model-selmom_data)
    dist4 = sqrt(distmat(1,1)/dble(4))

    simul_moments = 0.0
    simul_moments(1:nstats) = stats(1:nstats)
  end subroutine get_smm

end program smm_estimation_par





