program smm_estimation
  use prec; use parameters; use globals; use process; use equilibrium; use simulation
  use globals_montecarlo; use sort_tiago
  implicit none

  type(model) m
  type(param) par
  integer, parameter:: tsimul = 107, tdiscard = 103, nsimul = 50000, nvars = 24, niiter = 0, nmomentsall = 35
  real(dp), allocatable:: data(:,:,:),shocksz(:,:),shocksa(:),shocksx(:)
  integer,  allocatable:: iapath(:),ixpath(:)
  integer i,j,idl,iaa,ips,is1,is2,iloc(1),ik0,iflag_generateshocks
  real(dp) moments_data(nmoments,1),wmatrix(nmoments,nmoments),wmatrix4(4,4),identity(nmoments,nmoments),identity3(4,4),&
       qdist(ndl*naa*nps*ns1*ns2),qdist2(ndl*naa*nps*ns1*ns2),qdist3(ndl*naa*nps*ns1*ns2),qdist4(ndl*naa*nps*ns1*ns2)
  integer, dimension(ndl*naa*nps*ns1*ns2)::idx,xdl,xaa,xps,xs1,xs2
  real moments(nmomentsall),dchange

  iflag_generateshocks = 0
  
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

  ! MOMENTS USED:
  moments_data(1,1) = 0.109741186223729d0 	! Average investment
  moments_data(2,1) = 0.152542372881356d0	! Inaction at 1pc
  moments_data(3,1) = 0.144406779661017d0 	! Positive spike 20pc ik    
  moments_data(4,1) = 0.0725423728813559d0 	! Negative invest 1pc ik
  moments_data(5,1) = 0.150687829159272d0	! Stdev of invest rate (trim with 1%)
  moments_data(6,1) = 0.568584205084746d0	! Average leveraeg
  moments_data(7,1) = 0.0217399718163962d0	! Slope coefficient

  !moments_data(4,1) = 0.010d0	! Negative spike 20pc ik
  !moments_data(5,1) = 0.740d0 	! Positive invest 1pc ik

  ! WEIGHTING MATRIX:
  wmatrix(1,:) = (/107.66326032498057,7.556894725289942,-34.83911987039023,17.010281991464176,-25.172265233240417,0.252657894510155,1.5634025357296255/)
  wmatrix(2,:) = (/7.556894725289944,6.443888987743151,-0.8692717321618624,1.8392792299101866,-1.4848568188477818,-0.5808927010453249,-0.2437938283631378/)
  wmatrix(3,:) = (/-34.83911987039023,-0.8692717321618623,19.67977370249942,-3.4349357288762734,3.7227288761116935,-1.1199827456766607,-0.3992655053194529/)
  wmatrix(4,:) = (/17.010281991464176,1.839279229910186,-3.4349357288762725,14.698771516298216,-4.796343170704165,-0.21706481875267247,0.14052023690535345/)
  wmatrix(5,:) = (/-25.17226523324042,-1.4848568188477813,3.72272887611169,-4.796343170704165,21.420208847586906,0.8914139578883276,-0.6487392179242273/)
  wmatrix(6,:) = (/0.2526578945101549,-0.5808927010453248,-1.1199827456766607,-0.21706481875267247,0.8914139578883273,8.420253086989458,0.969906532014356/)
  wmatrix(7,:) = (/1.563402535729626,-0.24379382836313782,-0.39926550531945315,0.14052023690535348,-0.6487392179242273,0.969906532014356,7.041448276098913/)

  wmatrix4(1,:) = (/ 5.8592,    1.5116,    0.6239,    0.3334/)
  wmatrix4(2,:) = (/ 1.5116,    8.2707,    2.0403,   -4.2928/)
  wmatrix4(3,:) = (/ 0.6239,    2.0403,   12.0025,   -0.7946/)
  wmatrix4(4,:) = (/ 0.3334,   -4.2928,   -0.7946,   15.4052/)

  identity = 0d0
  do i=1,nmoments
     identity(i,i) = 1d0
  end do
  identity3 = 0d0
  do i=1,4
     identity3(i,i) = 1d0
  end do

  ! GRID FOR THE SEARCH

  smm%min_dl = 0.10970d0;	smm%max_dl = smm%min_dl
  smm%min_aa = 0.57500d0;	smm%max_aa = smm%min_aa
  smm%min_ps = 0.00950d0;	smm%max_ps = smm%min_ps
  smm%min_s1 = 0.57400d0;	smm%max_s1 = smm%min_s1
  smm%min_s2 = 1.43000d0+0.05d0;	smm%max_s2 = smm%min_s2

  call get_eqspace(ndl,smm%min_dl,smm%max_dl,smm%grid_dl)
  !call get_eqspace(ngm,smm%min_gm,smm%max_gm,smm%grid_gm)
  call get_eqspace(naa,smm%min_aa,smm%max_aa,smm%grid_aa)
  call get_eqspace(nps,smm%min_ps,smm%max_ps,smm%grid_ps)
  call get_eqspace(ns1,smm%min_s1,smm%max_s1,smm%grid_s1)
  call get_eqspace(ns2,smm%min_s2,smm%max_s2,smm%grid_s2)

  if(iflag_generateshocks.eq.1)then
     ! SHOCK SEQUENCE
     call get_simulshocks(nsimul,tsimul,shocksz,shocksa,shocksx,randomseed='n')
     iloc = minloc(abs(0.00d0 - m%a)); iapath = iloc(1)
     ixpath = 1
     iloc = minloc(abs(0.11d0 - m%a)); iapath(103) = iloc(1)
     iloc = minloc(abs(0.10d0 - m%a)); iapath(104) = iloc(1)
     iloc = minloc(abs(0.02d0 - m%a)); iapath(105) = iloc(1)
     iloc = minloc(abs(0.03d0 - m%a)); iapath(106) = iloc(1)
     iloc = minloc(abs(0.08d0 - m%a)); iapath(107) = iloc(1)
  else
     open(1,file='functions/shocksz_exporters.dat',position="rewind");read(1,*)shocksz;close(1)
     open(1,file='functions/shocksa_exporters.dat',position="rewind");read(1,*)shocksa;close(1)
     open(1,file='functions/shocksx_exporters.dat',position="rewind");read(1,*)shocksx;close(1)
     open(1,file='functions/iapath_exporters.dat',position="rewind");read(1,*)iapath;close(1)
     open(1,file='functions/ixpath_exporters.dat',position="rewind");read(1,*)ixpath;close(1)
  end if

  write(*,'(a40,a2,7f8.4,a2,a7,a2)')'','|',moments_data(:,1),'|'
  write(*,'(5a8,a2,7a8,a2,4a8,a4,a2,a8,a2,4a9)')&
       'dlt','a','ps','s1','s2','|', &
       'e.i','ina','p-spk','n-i','sd.i','e.l','slp','|', &
       'dist1','ecd1','dist2','ecld2','itr','|', &
       'sd.l','|', &
       'smaxk','smink','smaxp','sminp'

  !open (11,file='plot_data/Identification.txt',position="rewind")
  do is2=1,ns2
     do is1=1,ns1
        do ips=1,nps
           do iaa=1,naa
              !do igm=1,ngm
              do idl=1,ndl

                 j = idl + ndl*(iaa-1 + naa*(ips-1 + nps*(is1-1 + ns1*(is2-1))))
                 xdl(j) = idl
                 !xgm(j) = igm
                 xaa(j) = iaa
                 xps(j) = ips
                 xs1(j) = is1
                 xs2(j) = is2

                 par%delta   = smm%grid_dl(idl)
                 !par%gamma   = smm%grid_gm(igm)
                 par%a       = smm%grid_aa(iaa)
                 par%ps      = smm%grid_ps(ips)
                 par%sparam1 = smm%grid_s1(is1)
                 par%sparam2 = smm%grid_s2(is2)                    

                 do ik0 = nk,1,-1
                    m%k(ik0) = kmax*(1d0-par%delta)**(dble(nk-ik0)/dble(nk_divs))
                 end do
                 call get_smm(par,qdist(j),qdist2(j),qdist3(j),qdist4(j),moments)

                 write(*,'(5f8.4,a2,7f8.4,a2,4f8.4,i4,a2,f8.4,a2,4f9.5)') &
                      par%delta,par%a,par%ps,par%sparam1,par%sparam2,'|', &
                      moments(1),moments(9),moments(10),moments(8),moments(13),moments(4),moments(32),'|', &
                      qdist(j),qdist2(j),qdist3(j),qdist4(j),m%iter-niiter,'|', &
                      moments(16),'|', &
                      dble(count(data(:,tdiscard:tsimul,3).ge.m%k(nk)-erro5))/dble(nsimul*(tsimul-tdiscard+1)), &
                      dble(count(data(:,tdiscard:tsimul,3).le.m%k( 1)+erro5))/dble(nsimul*(tsimul-tdiscard+1)), &
                      dble(count(data(:,tdiscard:tsimul,4).ge.m%p(np)-erro5))/dble(nsimul*(tsimul-tdiscard+1)), &
                      dble(count(data(:,tdiscard:tsimul,4).le.m%p( 1)+erro5))/dble(nsimul*(tsimul-tdiscard+1))

                 !write(11,'(2f15.6)')par%sparam1,moments(4)
                 !end do
              end do
           end do
        end do
     end do
  end do
  !close(11)

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

  ! This is for the Jacobian
  write(*,'(9f20.9)') qdist(1),smm%grid_dl(1),smm%grid_aa(1),smm%grid_ps(1),smm%grid_s1(1),smm%grid_s2(1)
  write(*,'(9f20.9)') qdist(1),moments(1),moments(9),moments(10),moments(8),moments(13),moments(4),moments(32)
  write(*,'(9f20.9)') &
       qdist(1)-0.056751762, &
       moments(1)-0.110604621, &
       moments(9)-0.161699995, &
       moments(10)-0.184783995, &
       moments(8)-0.068584003, &
       moments(13)-0.109168939, &
       moments(4)-0.569049597, &
       moments(32)-0.025841406
  dchange = real(0.05d0)
  print*,'delta',dchange
  write(*,'(9f20.9)') &
       (qdist(1)-0.056751762)/dchange, &
       (moments(1)-0.110604621)/dchange, &
       (moments(9)-0.161699995)/dchange, &
       (moments(10)-0.184783995)/dchange, &
       (moments(8)-0.068584003)/dchange, &
       (moments(13)-0.109168939)/dchange, &
       (moments(4)-0.569049597)/dchange, &
       (moments(32)-0.025841406)/dchange
  print*,''

  deallocate(shocksz,shocksa,shocksx)
  deallocate(iapath,ixpath)
  deallocate(data)
  deallocate(m%k,m%p,m%z,m%a,m%xsigmaz)
  deallocate(m%prbz,m%prba,m%prbx)
  deallocate(m%evg)
  deallocate(m%ikp,m%ipp)
  deallocate(m%vg,m%vc,m%pz,m%pa,m%px)
  deallocate(m%coefevg)


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

end program smm_estimation





