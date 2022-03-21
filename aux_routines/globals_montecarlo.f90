module globals_montecarlo
  use prec; use parameters

  integer, parameter:: ndl=10
  !integer, parameter:: ngm=1
  integer, parameter:: naa=10
  integer, parameter:: nps=10
  integer, parameter:: ns1=10
  integer, parameter:: ns2=10
  integer, parameter:: nmoments=7

  type estimation     
     real(dp) grid_dl(ndl), min_dl, max_dl, true_dl
     !real(dp) grid_gm(ngm), min_gm, max_gm, true_gm
     real(dp) grid_aa(naa), min_aa, max_aa, true_aa
     real(dp) grid_ps(nps), min_ps, max_ps, true_ps
     real(dp) grid_s1(ns1), min_s1, max_s1, true_s1
     real(dp) grid_s2(ns2), min_s2, max_s2, true_s2     
  end type estimation

  type(estimation) smm
end module globals_montecarlo
