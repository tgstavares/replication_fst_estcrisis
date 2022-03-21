module getparams
  use prec
contains

  subroutine read_data_moments(sector,era)
    use globals_montecarlo        
    use csv_module; use iso_fortran_env, only: wp => real64
    implicit none
    integer, intent(in)::sector,era
    integer, parameter:: nera = 2,nsector = 14
    integer j,i

    type(csv_file) :: f
    character(len=30),dimension(:),allocatable :: header
    real(wp),dimension(:),allocatable :: m1
    logical :: status_ok
    integer,dimension(:),allocatable :: itypes
    
    call f%read('data/sectoral-moments-for-dynamic-estimation-v2.csv',header_row=1,status_ok=status_ok)
    call f%get_header(header,status_ok)
    call f%variable_types(itypes,status_ok)

    j = nera*(sector-1)+(era+1)
    write(*,'(a40,3i5)')'Reading moments   in line/sector/era:',j,sector,era

    call f%get(3 ,m1 ,status_ok); smm%moments_data(1 ,1) = m1(j) ! average_inv
    call f%get(4 ,m1 ,status_ok); smm%moments_data(2 ,1) = m1(j) ! pos_spike
    call f%get(5 ,m1 ,status_ok); smm%moments_data(3 ,1) = m1(j) ! neg_spike
    call f%get(6 ,m1 ,status_ok); smm%moments_data(4 ,1) = m1(j) ! inaction_1pc
    call f%get(7 ,m1 ,status_ok); smm%moments_data(5 ,1) = m1(j) ! neg_inv_rt
    call f%get(8 ,m1 ,status_ok); smm%moments_data(6 ,1) = m1(j) ! neg_inv_rt_l1pc
    call f%get(9 ,m1 ,status_ok); smm%moments_data(7 ,1) = m1(j) ! autocorrir
    call f%get(10,m1 ,status_ok); smm%moments_data(8 ,1) = m1(j) ! corr_lnomega_ir
    call f%get(11,m1 ,status_ok); smm%moments_data(9 ,1) = m1(j) ! corr_lnomega_lnk
    call f%get(12,m1 ,status_ok); smm%moments_data(10,1) = m1(j) ! corr_innov_ir
    call f%get(13,m1 ,status_ok); smm%moments_data(11,1) = m1(j) ! struct_b0_const
    call f%get(14,m1 ,status_ok); smm%moments_data(12,1) = m1(j) ! struct_rho
    call f%get(15,m1 ,status_ok); smm%moments_data(13,1) = m1(j) ! invar_sd_innov
    call f%get(16,m1 ,status_ok); smm%moments_data(14,1) = m1(j) ! capital_coeff
    call f%get(17,m1 ,status_ok); smm%moments_data(15,1) = m1(j) ! invar_sd_lnomega    
    call f%get(18,m1 ,status_ok); smm%moments_data(16,1) = m1(j) ! variant_mean_lnomega
    call f%get(19,m1 ,status_ok); smm%moments_data(17,1) = m1(j) ! inaction5pc
    call f%get(20,m1 ,status_ok); smm%moments_data(18,1) = m1(j) ! p25ir
    call f%get(21,m1 ,status_ok); smm%moments_data(19,1) = m1(j) ! medir
    call f%get(22,m1 ,status_ok); smm%moments_data(20,1) = m1(j) ! p75ir
    call f%get(23,m1 ,status_ok); smm%moments_data(21,1) = m1(j) ! sdir
    call f%get(24,m1 ,status_ok); smm%moments_data(22,1) = m1(j) ! min_lnk
    call f%get(25,m1 ,status_ok); smm%moments_data(23,1) = m1(j) ! max_lnk
    call f%get(26,m1 ,status_ok); smm%moments_data(24,1) = m1(j) ! min_lnomega
    call f%get(27,m1 ,status_ok); smm%moments_data(25,1) = m1(j) ! max_lnomega

    call f%destroy()
  end subroutine read_data_moments

  subroutine read_init_conditions_rev(sector,era,year,ninit,datainit)
    use csv_module; use iso_fortran_env, only: wp => real64
    implicit none
    integer, intent(in)::sector,era,year,ninit
    real(dp), intent(out)::datainit(3,ninit)
    type(csv_file) :: f
    character(len=30),dimension(:),allocatable :: header
    real(wp),dimension(:),allocatable :: m1,m2,m3,m0
    logical :: status_ok
    integer,dimension(:),allocatable :: itypes
    integer, parameter:: nera = 2,nsector = 14
    integer j,i,ii
    character(len=80), dimension(nsector*nera):: res_init
    real(dp) lk_lb,lk_ub,o_lb,o_ub


    j = nera*(sector-1)+(era+1)
    write(*,'(a40,3i5)')'Reading initconds in line/sector/era:',j,sector,era

    call f%read('data/sectoral-moments-for-dynamic-estimation-v2.csv',header_row=1,status_ok=status_ok)
    call f%get_header(header,status_ok)
    call f%variable_types(itypes,status_ok)
    call f%get(24 ,m0 ,status_ok); lk_lb = m0(j)
    call f%get(25 ,m0 ,status_ok); lk_ub = m0(j)
    call f%get(26 ,m0 ,status_ok); o_lb = m0(j)
    call f%get(27 ,m0 ,status_ok); o_ub = m0(j)
    call f%destroy()
    
    if(sector.ne.14) stop "ERROR - WRONG SECTOR FOR THIS EXERCISE"

    if(era.eq.0)then
       if(year.eq.2003.or.year.eq.2007)then
          if(year.eq.2003)then
             call f%read('data/new_initcond/pulled-yearly-joint-kolev-year-2003-6buckets.csv', &
                  header_row=1,status_ok=status_ok)
          else
             call f%read( &
                  'data/new_initcond/pulled-yearly-joint-kolev-year-2007-6buckets.csv', &
                  header_row=1,status_ok=status_ok)
          end if
       else
          stop "ERROR - Wrong year for era 0"
       end if
    else
       if(year.eq.2010.or.year.eq.2011.or.year.eq.2012.or.year.eq.2013.or.year.eq.2014)then
          if(year.eq.2010)then
             call f%read('data/new_initcond/pulled-yearly-joint-kolev-year-2010-6buckets.csv', &
                  header_row=1,status_ok=status_ok)
          elseif(year.eq.2011)then
             call f%read('data/new_initcond/pulled-yearly-joint-kolev-year-2011-6buckets.csv', &
                  header_row=1,status_ok=status_ok)
          elseif(year.eq.2012)then
             call f%read('data/new_initcond/pulled-yearly-joint-kolev-year-2012-6buckets.csv', &
                  header_row=1,status_ok=status_ok)
          elseif(year.eq.2013)then
             call f%read('data/new_initcond/pulled-yearly-joint-kolev-year-2013-6buckets.csv', &
                  header_row=1,status_ok=status_ok)
          else
             call f%read('data/new_initcond/pulled-yearly-joint-kolev-year-2014-6buckets.csv', &
                  header_row=1,status_ok=status_ok)
          end if
       else
          stop "ERROR - Wrong year for era 1"
       end if
    end if

    call f%get_header(header,status_ok)
    call f%variable_types(itypes,status_ok)
    call f%get(2 ,m1 ,status_ok);
    call f%get(3 ,m2 ,status_ok);
    call f%get(4 ,m3 ,status_ok); !expected leverage
    !call f%get(5 ,m3 ,status_ok); !predicted leverage
    
    ii = 1
    do i=1,40000
       if(m1(i).lt.lk_ub.and.m1(i).gt.lk_lb)then
          if(m2(i).lt.o_ub.and.m2(i).gt.o_lb)then
             datainit(1,ii) = exp(m1(i))
             datainit(2,ii) = m2(i)
             datainit(3,ii) = m3(i)
             if(ii.eq.ninit) exit
             ii = ii+1
          end if
       end if
    end do    
    if(ii.lt.ninit) stop "ERROR - Too few data points given bounds"
    call f%destroy()
  end subroutine read_init_conditions_rev

end module getparams
