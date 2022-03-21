module simulation
  use prec; use parameters; use globals
contains

  subroutine get_simulshocks(nsimul,tsimul,shocksz,shocksa,shocksx,randomseed)
    implicit none
    integer, intent(in)::nsimul,tsimul
    real(dp), intent(out)::shocksz(nsimul,tsimul),shocksa(tsimul),shocksx(tsimul)
    character(len=1),intent(in), optional::randomseed
    integer flag,s1,s3,iz0,iz1,t,n,i
    integer, dimension(:),allocatable::s2
    flag = 1
    if(present(randomseed))then
       if(randomseed(1:1).eq.'y'.or.randomseed(1:1).eq.'Y') flag = 0
    end if
    ! JUST GENERATE RANDOM NUMBERS
    s3 = 12503
    if(flag.eq.0) call system_clock(COUNT=s3)
    call random_seed(size=s1); allocate(s2(s1));
    do i=1,s1
       s2(i) = s3 + i
    end do
    call random_seed(put=s2)
    call random_number(shocksz)
    call random_number(shocksa)
    call random_number(shocksx)
  end subroutine get_simulshocks

  subroutine get_simulations(par,m,nsimul,tsimul,nvars,shocksz,shocksa,shocksx,data,initconds,iapath,ixpath,replace_exit)
    use omp_lib; use values
    implicit none
    type(param), intent(in):: par
    type(model), intent(in):: m
    integer, intent(in)::nsimul,tsimul,nvars
    real(dp), intent(in)::shocksz(nsimul,tsimul),shocksa(tsimul),shocksx(tsimul)
    real(dp), intent(in), optional::initconds(nsimul,6)
    integer,  intent(in), optional::iapath(tsimul),ixpath(tsimul)
    character(len=1),intent(in), optional::replace_exit
    real(dp), intent(out)::data(nsimul,tsimul,nvars)
    integer iflag_initconds,iflag_iapath,iflag_ixpath,iflag_replace,n,t,i,iz0,ia0,ix0,iz_1,ia_1,ix_1,iaux(1),ik0,ip0
    real(dp) k0,p0,z_1,a_1,t0,k_1,inv_1,fcumulz(nz),fcumula(na),fcumulx(nx),z0,a0,v0,k1,p1,qq,inv0,div0,nworth,cashflow
    real(dp) ek,ep,ez

    ! CHECK OPTIONS
    if(present(initconds))then
       iflag_initconds = 1
    else
       iflag_initconds = 0
    end if
    if(present(iapath))then
       iflag_iapath = 1
    else
       iflag_iapath = 0
    end if
    if(present(ixpath))then
       iflag_ixpath = 1
    else
       iflag_ixpath = 0
    end if
    if(present(replace_exit))then
       if(replace_exit(1:1).eq.'n'.or.replace_exit(1:1).eq.'N')then
          iflag_replace = 0
       else
          iflag_replace = 1
       end if
    else
       iflag_replace = 1
    end if

    data = 0d0
    do t=1,tsimul

       ! FIRST THE AGGREGATE SHOCKS SEQUENCE
       if(iflag_iapath.eq.1)then
          if(t.gt.1) then
             a_1 = a0
          else
             a_1 = m%a(iapath(1))
          end if
          ia0 = iapath(t)
          a0 = m%a(ia0)
       else
          if(t.gt.1) then
             a_1 = a0
          else
             a_1 = m%a(na/2+1)
          end if
          iaux = minloc(abs(m%a - a_1)); ia_1 = iaux(1)          
          fcumula(1) = m%prba(ia_1,1)
          if(shocksa(t).le.fcumula(1))then
             ia0 = 1
          else
             do i=2,na
                fcumula(i) = fcumula(i-1) + m%prba(ia_1,i)
                if(fcumula(i).gt.shocksa(t))then
                   ia0 = i
                   exit
                end if
             end do
          end if
          a0 = m%a(ia0)
       end if

       if(iflag_ixpath.eq.1)then          
          if(t.gt.1) then
             ix_1 = ix0
          else
             ix_1 = ixpath(1)
          end if
          ix0 = ixpath(t)
       else
          if(t.gt.1) then
             ix_1 = ix0
          else
             ix_1 = 1
          end if
          fcumulx(1) = m%prbx(ix_1,1)
          if(shocksx(t).le.fcumulx(1))then
             ix0 = 1
          else
             do i=2,nx
                fcumulx(i) = fcumulx(i-1) + m%prbx(ix_1,i)
                if(fcumulx(i).gt.shocksx(t))then
                   ix0 = i
                   exit
                end if
             end do
          end if
       end if

       !print*,ix0,ix_1
       !stop
       
       ! SECOND IDIOSYNCRATIC SHOCK AND SIMUL OUTCOME
       !$omp parallel do collapse(1) private(n,k0,p0,z_1,t0,k_1,inv_1,iaux,fcumulz,iz0,i,z0,v0,k1,p1,qq,inv0,div0,cashflow)
       do n=1,nsimul
          ! INITIAL VARIABLES OR CARRIED FROM PREV PERIOD
          if(iflag_initconds.eq.1)then
             k0    = initconds(n,1)
             p0    = initconds(n,2)
             z_1   = initconds(n,3)
             t0    = initconds(n,4) + 1d0
             k_1   = initconds(n,5)
             inv_1 = initconds(n,6)
             if(t0.le.0.5d0) t0 = -1d0
          else
             if(t.gt.1)then
                k0    = data(n,t-1,6)
                p0    = data(n,t-1,7)
                z_1   = data(n,t-1,5)
                t0    = data(n,t-1,9) + 1d0
                k_1   = data(n,t-1,3)
                inv_1 = data(n,t-1,10)
                if(t0.le.0.5d0) t0 = -1d0
             else
                k0    = m%k(nk)
                p0    = m%p(1)
                z_1   = m%z(nz/2+1)
                t0    = 1d0
                k_1   = k0
                inv_1 = par%delta*k0
             end if
          end if

          ! GET REALIZATION OF IDIOSYNCRATIC SHOCK
          iaux = minloc(abs(m%z - z_1)); iz_1 = iaux(1)          
          fcumulz(1) = m%prbz(iz_1,1,ix0)
          if(shocksz(n,t).le.fcumulz(1))then
             iz0 = 1
          else
             do i=2,nz
                fcumulz(i) = fcumulz(i-1) + m%prbz(iz_1,i,ix0)
                if(fcumulz(i).gt.shocksz(n,t))then
                   iz0 = i
                   exit
                end if
             end do
          end if
          z0 = m%z(iz0)

          ! SIMULATION OUTCOMES
          if(t0.le.0.5d0)then
             k1   = k0
             p1   = p0
             qq   = 0d0
             div0 = -1d0
             v0   = -1d0
          else             
             call get_kp_bp(par,m,iz0,ia0,ix0,k0,p0,v0,k1,p1,qq)
             qq = qq*(1d0+m%r)
             inv0 = k1 - (1d0-par%delta)*k0
             call cashflowr(par,m%r,m%tau,m%taudiv,k0,p0,z0,a0,k1,p1,qq,div0)
             qq = qq/(1d0+m%r)
             !print*,div0,k0,p0,z0,k1,p1,qq
          end if
          if(v0.lt.vexit.or.div0.lt.dexit) t0 = -1d0 ! MARK CASE IF EXIT

          !print*,k0,z0,exp(z0+a0) * k0**theta,div0,inv0,p0,p1
          !stop

          ! STORE IN ARRAY
          data(n,t,1)  = dble(n)
          data(n,t,2)  = dble(t)
          data(n,t,3)  = k0
          data(n,t,4)  = p0
          data(n,t,5)  = z0
          data(n,t,6)  = k1
          data(n,t,7)  = p1
          data(n,t,8)  = z_1
          data(n,t,9)  = t0
          data(n,t,10) = inv0
          data(n,t,11) = inv_1          
          data(n,t,12) = inv0/k0
          data(n,t,13) = inv_1/k_1
          data(n,t,14) = (z0+a0)-(z_1+a_1)
          data(n,t,15) = p0*k0
          data(n,t,16) = p1*k1
          data(n,t,17) = qq
          data(n,t,18) = exp(z0+a0) * k0**theta
          data(n,t,19) = div0
          data(n,t,20) = v0
          data(n,t,21) = a0
          data(n,t,22) = a_1
          data(n,t,23) = dble(ix0)
          data(n,t,24) = dble(ix_1)
          data(n,t,25) = qq
          
          !if(n.eq.1.and.t.le.50) write(*,'(2i10,15f15.4)')n,t,t0,a0,z0,k0,k1,p0,p1,inv0,inv0/k0,qq,div0,v0
          !if(n.eq.1.and.t.ge.1.and.t.le.241)write(*,'(21f10.3)')data(n,t,1:20)
       end do
       !$omp end parallel do

       if(iflag_replace.eq.1)then
          ! FOR FIRMS THAT EXIT REPLACE THEM AT INDUSTRY AVERAGES
          ek = sum(data(:,t,3),data(:,t,9).gt.0.0d0)/dble(count(data(:,t,9).gt.0.0d0))
          ep = sum(data(:,t,4),data(:,t,9).gt.0.0d0)/dble(count(data(:,t,9).gt.0.0d0))
          ez = sum(data(:,t,5),data(:,t,9).gt.0.0d0)/dble(count(data(:,t,9).gt.0.0d0))
          iaux = minloc(abs(m%k - ek)); ik0 = iaux(1)
          iaux = minloc(abs(m%p - ep)); ip0 = iaux(1)
          iaux = minloc(abs(m%z - ez)); iz0 = iaux(1)

          !$omp parallel do collapse(1) private(n,k0,p0,z0,z_1,t0,k_1,inv_1,v0,k1,p1,qq,inv0,div0,cashflow)
          do n=1,nsimul
             if(data(n,t,9).le.0.0d0) then
                k0    = m%k(ik0)
                p0    = m%p(ip0)
                z0    = m%z(iz0)
                z_1   = z0
                t0    = dble(1)
                k_1   = k0
                inv_1 = par%delta*k0
                call get_kp_bp(par,m,iz0,ia0,ix0,k0,p0,v0,k1,p1,qq)
                qq    = qq*(1d0+m%r)  
                inv0  = k1 - (1d0-par%delta)*k0
                call cashflowr(par,m%r,m%tau,m%taudiv,k0,p0,z0,a0,k1,p1,qq,div0)
                qq    = qq/(1d0+m%r)
                ! ALLOW FOR A FATAL ERROR
                if(v0.lt.vexit.or.div0.lt.dexit) stop 'FATAL ERROR: Industry collapse by default'
                data(n,t,1)  = dble(n)
                data(n,t,2)  = dble(t)
                data(n,t,3)  = k0
                data(n,t,4)  = p0
                data(n,t,5)  = z0
                data(n,t,6)  = k1
                data(n,t,7)  = p1
                data(n,t,8)  = z_1
                data(n,t,9)  = t0
                data(n,t,10) = inv0
                data(n,t,11) = inv_1          
                data(n,t,12) = inv0/k0
                data(n,t,13) = inv_1/k_1
                data(n,t,14) = (z0+a0)-(z_1+a_1)
                data(n,t,15) = p0*k0
                data(n,t,16) = p1*k1
                data(n,t,17) = qq
                data(n,t,18) = exp(z0+a0) * k0**theta
                data(n,t,19) = div0
                data(n,t,20) = v0
                data(n,t,21) = a0
                data(n,t,22) = a_1
                data(n,t,23) = dble(ix0)
                data(n,t,24) = dble(ix_1)
                data(n,t,25) = qq
             end if
          end do
          !$omp end parallel do
       end if
    end do
    !stop
  end subroutine get_simulations

  subroutine get_stats(data,nstats,stats,showoutput,cond)
    implicit none
    integer, intent(in)::nstats
    real(dp),intent(in)::data(:,:,:)
    real,intent(out)::stats(nstats)
    logical,intent(in), optional::cond(:,:)
    character(len=1), intent(in), optional::showoutput
    logical mask(size(data,dim=1),size(data,dim=2))
    integer flag_showoutput,n
    real e_ik,m_ik,a_ik,p_ik,n_ik,ina,p_spike,n_spike,sd_o,aut_o,sd_ik,aut_ik,corr_iko,corr_ikdo,e_bp,m_bp,a_bp,e_k,sd_bp,e_b
    real corr_k1o,corr_bpo,corr_ikbp,spike,covari,vari,slope
    real a_sales,a_k,a_salesk,sd_mpk,a_b,a_b1,a_i,q66_bp
    real mean_a,mean_z,var_z,agg_profits,agg_capital

    ! TAKE CARE OF OPTIONAL VARIABLES
    flag_showoutput = 1
    if(present(showoutput))then
       if(showoutput(1:1).eq.'n'.or.showoutput(1:1).eq.'No') flag_showoutput = 0
    end if
    if(present(cond))then
       mask = cond
    else
       mask = .true.
    end if

    ! RELEVANT OBSERVATIONS ARE THE ONES IMPLIED BY THE MASK
    n = count(mask)    
    call mean(real(pack(data(:,:,12),mask)),n,0,e_ik)
    call mean(real(pack(data(:,:, 7),mask)),n,0,e_bp)
    call mean(real(pack(data(:,:, 3),mask)),n,0,e_k)
    call mean(real(pack(data(:,:, 4),mask)),n,0,e_b)    
    call mean(real(pack(data(:,:, 5),mask)),n,0,mean_a)    
    call sampp(real(pack(data(:,:,12),mask)),n,0.5,0,m_ik)
    call sampp(real(pack(data(:,:, 7),mask)),n,0.5,0,m_bp)
    call sampp(real(pack(data(:,:, 7),mask)),n,0.66,0,q66_bp)    
    call sd(real(pack(data(:,:,12),mask)),n,0,sd_ik)
    call sd(real(pack(data(:,:, 7),mask)),n,0,sd_bp)
    call sd(real(pack(data(:,:, 5),mask)),n,0,sd_o)
    call sd(real(pack(data(:,:,18)/data(:,:,3),mask)),n,0,sd_mpk)
    call corr(real(pack(data(:,:, 5),mask)),real(    pack(data(:,:,12),mask)) ,n,0,corr_iko)
    call corr(real(pack(data(:,:, 5),mask)),real(log(pack(data(:,:, 6),mask))),n,0,corr_k1o)
    call corr(real(pack(data(:,:, 5),mask)),real(    pack(data(:,:, 7),mask)) ,n,0,corr_bpo)
    call corr(real(pack(data(:,:,12),mask)),real(    pack(data(:,:, 7),mask)) ,n,0,corr_ikbp)
    call corr(real(pack(data(:,:, 5),mask)),real(    pack(data(:,:, 8),mask)) ,n,0,aut_o)
    call corr(real(pack(data(:,:,12),mask)),real(    pack(data(:,:,13),mask)) ,n,0,aut_ik)
    call corr(real(pack(data(:,:,14),mask)),real(    pack(data(:,:,12),mask)) ,n,0,corr_ikdo)
    a_ik    = real(sum(pack(data(:,:,10),mask))/sum(pack(data(:,:,3),mask)))
    a_bp    = real(sum(pack(data(:,:,16),mask))/sum(pack(data(:,:,6),mask)))
    a_sales = real(sum(pack(data(:,:,18),mask)))
    a_k     = real(sum(pack(data(:,:, 3),mask)))
    a_b     = real(sum(pack(data(:,:,15),mask)))
    a_b1    = real(sum(pack(data(:,:,16),mask)))
    a_i     = real(sum(pack(data(:,:,10),mask)))
    p_ik    = real(count(pack(data(:,:,12),mask).gt. 0.01))/real(n)
    n_ik    = real(count(pack(data(:,:,12),mask).lt.-0.01))/real(n)    
    p_spike = real(count(pack(data(:,:,12),mask).gt. 0.2))/real(n)
    n_spike = real(count(pack(data(:,:,12),mask).lt.-0.2))/real(n)
    ina     = real(count(abs(pack(data(:,:,12),mask)).lt.0.01))/real(n)    
    spike   = p_spike+n_spike    
    a_salesk= a_sales / a_k
    call cov(real(pack(data(:,:,7),mask)),real(pack(data(:,:,14),mask)) ,n,0,covari)
    call var(real(pack(data(:,:,7),mask)),n,0,vari)
    slope = covari/vari


    
    ! EXPORT A STATS ARRAY
    stats(1) = e_ik;	  stats(2) = m_ik;	stats(3) = a_ik
    stats(4) = e_bp;	  stats(5) = m_bp;	stats(6) = a_bp
    stats(7) = p_ik;	  stats(8) = n_ik;	stats(9) = ina
    stats(10)= p_spike;	  stats(11)= n_spike;	stats(12)= spike
    stats(13)= sd_ik;	  stats(14)= corr_iko;	stats(15)= corr_k1o
    stats(16)= sd_bp;	  stats(17)= corr_bpo;	stats(18)= corr_ikbp
    stats(19)= sd_o;	  stats(20)= aut_o;	stats(21)= aut_ik
    stats(22)= corr_ikdo; stats(23)= e_k;	stats(24)= a_sales
    stats(25)= a_salesk;  stats(26)= sd_mpk;	stats(27)= a_b
    stats(28)= a_b1;	  stats(29)= a_i;       stats(30)= q66_bp
    stats(31)= e_b;	  stats(32)= slope

    ! PRINT SOME OF THE STATS
    if(flag_showoutput.eq.1)then
       write(*,'(22a11)') &
            'obs','E(ik)','med(ik)','pos(ik)','neg(ik)','inaction', &
            'pos spk','neg spk','sd(o)','aut(o)','sd(ik)', &
            'aut(ik)','cr(ik,o)','cr(ik,do)','sl(do,bp)', &
            'E(bp)','med(bp)','q66(bp)','sd(bp)','E(k)'
       write(*,'(i11,21f11.4)') &
            n,e_ik,m_ik,p_ik,n_ik,ina, &
            p_spike,n_spike,sd_o,aut_o,sd_ik, &
            aut_ik,corr_iko,corr_ikdo, &
            slope,e_bp,m_bp,q66_bp,sd_bp,e_k
    end if
  end subroutine get_stats

  
end module simulation
