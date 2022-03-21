module globals
  use prec

  ! MODEL VARIABLES
  type model
     ! VECTORS
     real(dp), allocatable, dimension(:)	    ::k,p,z,a,xsigmaz
     ! NON-VECTORIZED ARRAYS
     real(dp), allocatable, dimension(:,:)	    ::prba,prbx
     real(dp), allocatable, dimension(:,:,:)	    ::prbz
     real(dp), allocatable, dimension(:,:,:,:,:)    ::evg,qbg
     real(dp), allocatable, dimension(:,:,:,:,:,:,:)::coefevg,coefqbg
     ! VECTORIZED ARRAYS
     integer,  allocatable, dimension(:)	    ::ikp,ipp
     real(dp), allocatable, dimension(:)	    ::vg,vc,pz,pa,px
     ! OTHER VARIABLES
     integer iter,hflag,jflag
     real(dp) sparam,tau,taudiv, r,beta,dwlk
  end type model

  ! PARAMETERS SUBJECTED TO SEARCH
  type param
     real(dp) delta,gamma,a,ps,sparam1,sparam2
  end type param
  
  ! OTHER PARAMETERS TO BE INITIATED DURING THE PROGRAM
  real(dp) t1,t2,amin,amax,zmin,zmax
  
end module globals
