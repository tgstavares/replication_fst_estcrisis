module parameters
  use prec

  ! Grids
  integer, parameter:: &
       ! SIZE OF GRIDS
       nx	= 2, & 
       na       = 11, &
       nz	= 21, & 	! 85 15
       nk	= 91, & 	! 65
       np       = 61, & 	! 79 40 33
       nk_divs  = 2, &
       nvfi     = 350, &
       naccel   = 10
  
  ! Parameters
  real(dp), parameter:: &
       r	= 0.02d0, &		! INTEREST RATE
       beta     = 1d0/(1d0+r), &	! DISCOUNT FACTOR
       sigmaa 	= 0.065d0, &		! VOL AGG SHOCK
       rhoa 	= 0.704d0, &		! PER AGG SHOCK
       sigmaz  	= 0.392d0, &		! VOL IDI SHOCK
       rhoz    	= 0.733d0, &		! PER IDI SHOCK
       theta	= 0.727d0, &		! CURVATURE PROFIT

       kmax     = 6500d0, &
       kmin     = 110d0, &       

       pmax     = 0.82d0, &
       pmin	= 0.0d0, &

       mukappa  = 12d0, &
       auxis    = 1.0d0, &
       ceq      = 0.1d0, &
       cet      = 10d0, &

       sjump    = 0.416d0/sigmaz, &	! JUMP IN VOL DUE TO UNC
       muz      = 0d0, &
       mu1      = muz - 0.5d0 * (sigmaz      )**2 / (1d0-rhoz**2) , &
       mu2      = muz - 0.5d0 * (sjump*sigmaz)**2 / (1d0-rhoz**2), &
       nsd      = 5d0, &

       uncfreq 	= 1d0-0.741d0, &  	! COND PROB OF UNC SHOCK
       uncpers 	= 0.330d0 		! COND UNC SHOCK PERSISTENCE

                     
  ! Control parameters - small and large numbers
  real(dp), parameter:: &
       infty = 1.0d15, &
       erro3 = 1.0d-3, &
       erro4 = 1.0d-5, &
       erro5 = 1.0d-5, &
       erro6 = 1.0d-6, &
       erro7 = 1.0d-7, &
       erro8 = 1.0d-8, &
       erro9 = 1.0d-9, &
       erro10= 1.0d-10, &
       erro12= 1.0d-12

end module parameters
