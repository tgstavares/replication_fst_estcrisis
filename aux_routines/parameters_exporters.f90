module parameters
  use prec

  ! Grids
  integer, parameter:: &
       ! SIZE OF GRIDS
       nx	= 2, & 
       na       = 11, &
       nz	= 21, & 	! 85 15
       nk	= 95, & 	! 65
       np       = 61, & 	! 79 40 33
       nk_divs  = 2, &
       nvfi     = 350, &
       naccel   = 30
  
  ! Parameters
  real(dp), parameter:: &
       r	= 0.025d0, &		! INTEREST RATE
       beta     = 1d0/(1d0+r), &	! DISCOUNT FACTOR
       sigmaa 	= 0.0866873d0, &		! VOL AGG SHOCK
       rhoa 	= 0.3928568d0, &		! PER AGG SHOCK
       sigmaz  	= 0.3692738d0, &		! VOL IDI SHOCK
       rhoz    	= 0.7228288d0, &		! PER IDI SHOCK
       theta	= 0.7181045d0, &		! CURVATURE PROFIT
       vexit    =-0d10, &		! DEFAULT EXIT
       dexit    =-0d10, &		! LIQUIDITY EXIT

       kmax     = 4000d0, & ! 15000
       kmin     = 110d0, &       

       pmax     = 0.80d0, & ! 0.7
       pmin	= 0.40d0, & ! 0,0

       sjump    = 0.4178976d0/sigmaz, &	! JUMP IN VOL DUE TO UNC
       muz      = 0d0, &
       mu1      = muz - 0.5d0 * (sigmaz      )**2 / (1d0-rhoz**2) , &
       mu2      = muz - 0.5d0 * (sjump*sigmaz)**2 / (1d0-rhoz**2), &
       nsd      = 5d0, &

       uncfreq 	= 1d0-0.6921704d0, &  	! COND PROB OF UNC SHOCK
       uncpers 	= 0.4338174d0 		! COND UNC SHOCK PERSISTENCE

                     
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
