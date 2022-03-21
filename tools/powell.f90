module powell_routines
contains

  SUBROUTINE POWELL(FUNC,P,XI,N,NP,FTOL,ITER,FRET)
  !-----------------------------------------------------------
  ! Minimization of a function  FUNC of N variables  (FUNC is
  ! not an argument, it is a fixed function name). Input con-
  ! sists of an initial starting point P  that is a vector of
  ! length N; an initial matrix XI  whose  logical dimensions
  ! are N by N, physical dimensions NP by NP, and whose columns
  ! contain the initial set of directions (usually the N unit
  ! vectors); and FTOL, the fractional tolerance in the func-
  ! tion value such that failure to decrease by more than this
  ! amount on one iteration signals doneness. On output, P is
  ! set to the best point found, XI is the then-current direc-
  ! tion set,  FRET is the returned function value at P,  and
  ! ITER is the number of iterations taken. The routine LINMIN
  ! is used.
  !------------------------------------------------------------
  IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
  PARAMETER(NMAX=2,ITMAX=100)
  DIMENSION P(NP),XI(NP,NP),PT(NMAX),PTT(NMAX),XIT(NMAX)

  interface
     DOUBLE PRECISION PURE FUNCTION FUNC(ARG)
       DOUBLE PRECISION, intent(in) :: ARG(2)
     end function FUNC
  end interface

  FRET=  FUNC(P) 

  DO J=1,N
     PT(J)=P(J)       !Save initial pont
  END DO
  ITER=0
1 ITER=ITER+1
  FP=FRET
  IBIG=0
  DEL=0.0D+0           !Will be the biggest function decrease.
  DO I=1,N           !In each iteration, loop over all directions in the set.
     DO J=1,N         !Copy the direction
        XIT(J)=XI(J,I)
     END DO
     FPTT=FRET
     CALL LINMIN(FUNC,P,XIT,N,FRET)  !Minimize along it.
     IF (DABS(FPTT-FRET).GT.DEL) THEN
        DEL=DABS(FPTT-FRET)
        IBIG=I
     END IF
  END DO
  IF (2.D0*DABS(FP-FRET).LE.FTOL*(DABS(FP)+DABS(FRET))) RETURN !Termination criterion
  IF (ITER.EQ.ITMAX) Then
     !Pause ' Powell exceeding maximum iterations.'
     return
  END IF
  DO J=1,N
     PTT(J)=2.D0*P(J)-PT(J)  !Construct the extrapolated point and the average
     XIT(J)=P(J)-PT(J)       !direction moved. Save the old starting point.
     PT(J)=P(J)
  END DO
  FPTT= FUNC(PTT) ! objective_function_nobind(PTT) !          !Function value at extrapolated point.
  IF (FPTT.GE.FP) GO TO 1   !One reason not to use new direction. 
  T=2.D0*(FP-2.D0*FRET+FPTT)*(FP-FRET-DEL)**2-DEL*(FP-FPTT)**2
  IF (T.GE.0.D0) GO TO 1    !Other reason not to use new direction.
  CALL LINMIN(FUNC,P,XIT,N,FRET) !Move to the minimum of the new direction.
  DO J=1,N                  !and save the new direction
     XI(J,IBIG)=XIT(J)
  END DO
  GO TO 1
END SUBROUTINE POWELL

SUBROUTINE LINMIN(FUNC,P,XI,N,FRET)
  !----------------------------------------------------------
  ! Given an N dimensional point P and a N dimensional direc-
  ! tion XI, moves and resets P to where the function FUNC(P)
  ! takes on a minimum along the direction XI from P, and 
  ! replaces XI by the actual vector displacement that P was
  ! moved. Also returns as FRET the value of FUNC at the
  ! returned location P. This is actually all accomplished by
  ! calling the routines MNBRAK and BRENT.
  !----------------------------------------------------------
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  PARAMETER(NMAX=50,TOL=1.D-4)
  DIMENSION P(N),XI(N)
  COMMON /F1COM/ PCOM(NMAX),XICOM(NMAX),NCOM

  interface
     DOUBLE PRECISION PURE FUNCTION FUNC(ARG)
       DOUBLE PRECISION, intent(in) :: ARG(2)
     end function FUNC
  end interface
  NCOM=N
  DO J=1,N
     PCOM(J)=P(J)
     XICOM(J)=XI(J)
  END DO
  AX=0.D0
  XX=1.D0
  BX=2.D0
  CALL MNBRAK(FUNC,AX,XX,BX,FA,FX,FB)
  FRET=BRENT(FUNC,AX,XX,BX,TOL,XMIN)
  DO J=1,N
     XI(J)=XMIN*XI(J)
     P(J)=P(J)+XI(J)
  END DO
  RETURN
END SUBROUTINE LINMIN

DOUBLE PRECISION FUNCTION F1DIM(FUNC,X)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  PARAMETER(NMAX=50)
  COMMON /F1COM/ PCOM(NMAX),XICOM(NMAX),NCOM
  DIMENSION XT(NMAX)
  !EXTERNAL objective_function_nobind
  ! EXTERNAL FUNC

  interface
     DOUBLE PRECISION PURE FUNCTION FUNC(ARG)
       DOUBLE PRECISION, intent(in) :: ARG(2)
     end function FUNC
  end interface

  DO J=1, NCOM
     XT(J)=PCOM(J)+X*XICOM(J)
  END DO
  F1DIM = FUNC(XT) !objective_function_nobind(XT) !
  RETURN
END FUNCTION F1DIM

SUBROUTINE MNBRAK(FUNC,AX,BX,CX,FA,FB,FC)
  !----------------------------------------------------------------------
  !Given a Function F1DIM(X), and given distinct initial points AX and
  !BX, this routine searches in the downhill direction (defined by the
  !F1DIMtion as evaluated at the initial points) and returns new points
  !AX, BX, CX which bracket a minimum of the Function. Also returned
  !are the Function values at the three points, FA, FB and FC.
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  PARAMETER(GOLD=1.618034,GLIMIT=100.,TINY=1.D-20)
  interface
     DOUBLE PRECISION PURE FUNCTION FUNC(ARG)
       DOUBLE PRECISION, intent(in) :: ARG(2)
     end function FUNC
  end interface

  !The first parameter is the default ratio by which successive intervals
  !are magnified; the second is the maximum magnification allowed for
  !a parabolic-fit step.
  !----------------------------------------------------------------------
  FA=F1DIM(FUNC,AX)
  FB=F1DIM(FUNC,BX)
  IF(FB.GT.FA) THEN
     DUM=AX
     AX=BX
     BX=DUM
     DUM=FB
     FB=FA
     FA=DUM
  ENDIF
  CX=BX+GOLD*(BX-AX)
  FC=F1DIM(FUNC,CX)
1 IF(FB.GE.FC) THEN
     R=(BX-AX)*(FB-FC)
     Q=(BX-CX)*(FB-FA)
     U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
     ULIM=BX+GLIMIT*(CX-BX)
     IF((BX-U)*(U-CX).GT.0) THEN
        FU=F1DIM(FUNC,U)
        IF(FU.LT.FC) THEN
           AX=BX
           FA=FB
           BX=U
           FB=FU
           GOTO 1
        ELSE IF(FU.GT.FB) THEN
           CX=U
           FC=FU
           GOTO 1
        ENDIF
        U=CX+GOLD*(CX-BX)
        FU=F1DIM(FUNC,U)
     ELSE IF((CX-U)*(U-ULIM).GT.0) THEN
        FU=F1DIM(FUNC,U)
        IF(FU.LT.FC) THEN
           BX=CX
           CX=U
           U=CX+GOLD*(CX-BX)
           FB=FC
           FC=FU
           FU=F1DIM(FUNC,U)
        ENDIF
     ELSE IF((U-ULIM)*(ULIM-CX).GE.0) THEN
        U=ULIM
        FU=F1DIM(FUNC,U)
     ELSE
        U=CX+GOLD*(CX-BX)
        FU=F1DIM(FUNC,U)
     ENDIF
     AX=BX
     BX=CX
     CX=U
     FA=FB
     FB=FC
     FC=FU
     GOTO 1
  ENDIF
  RETURN
END SUBROUTINE MNBRAK

DOUBLE PRECISION FUNCTION BRENT(FUNC,AX,BX,CX,TOL,XMIN)
  !-------------------------------------------------------------------
  !Given a function F1DIM, and a bracketing triplet of abscissas
  !AX,BX,CX (such that BX is between AX and CX, and F(BX) is less 
  !than both F(AX) and F(CX)), this routine isolates the minimum 
  !to a fractional precision of about TOL using Brent's method.
  !The abscissa of the minimum is returned in XMIN, and the minimum
  !function value is returned as BRENT, the returned function value.
  !-------------------------------------------------------------------
  !INTEGER,PARAMETER :: ITMAX=100 ! The maximum number of iterations.
  !DOUBLE PRECISION,PARAMETER :: CGOLD=0.3819660D0,ZEPS=1.D-10 ! CGOLD is the golden ratio, ZEPS is a small number that protects
  !PARAMETER(ITMAX=100,CGOLD=.3819660,ZEPS=1.D-10)
  !Maximum allowed number of iterations; golden ration; and a small
  !number which protects against trying to achieve fractional accuracy
  !for a minimum that happens to be exactly zero.
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  PARAMETER(ITMAX=100,CGOLD=.3819660,ZEPS=1.D-10)
  interface
     DOUBLE PRECISION PURE FUNCTION FUNC(ARG)
       DOUBLE PRECISION, intent(in) :: ARG(2)
     end function FUNC
  end interface


  A=MIN(AX,CX)
  B=MAX(AX,CX)
  V=BX
  W=V
  X=V
  E=0.
  FX=F1DIM(FUNC,X)
  FV=FX
  FW=FX
  DO 11 ITER=1,ITMAX	                                !main loop
     XM=0.5*(A+B)
     TOL1=TOL*ABS(X)+ZEPS
     TOL2=2.*TOL1
     IF (ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 3  !Test for done here
     IF (ABS(E).GT.TOL1) THEN     !Construct a trial parabolic fit
        R=(X-W)*(FX-FV)
        Q=(X-V)*(FX-FW)
        P=(X-V)*Q-(X-W)*R
        Q=.2*(Q-R)
        IF (Q.GT.0)  P=-P
        Q=ABS(Q)
        ETEMP=E
        E=D
        IF (ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(A-X).OR.  &
             P.GE.Q*(B-X))  GOTO 1
        !   The above conditions determine the acceptability of the 
        !   parabolic fit. Here it is o.k.:
        D=P/Q
        U=X+D
        IF (U-A.LT.TOL2.OR.B-U.LT.TOL2)  D=SIGN(TOL1,XM-X)
        GOTO 2
     ENDIF
1    IF (X.GE.XM) THEN
        E=A-X
     ELSE
        E=B-X
     ENDIF
     D=CGOLD*E
2    IF (ABS(D).GE.TOL1) THEN
        U=X+D
     ELSE
        U=X+SIGN(TOL1,D)
     ENDIF
     FU=F1DIM(FUNC,U)  !This the one function evaluation per iteration
     IF (FU.LE.FX) THEN
        IF (U.GE.X) THEN
           A=X
        ELSE
           B=X
        ENDIF
        V=W
        FV=FW
        W=X
        FW=FX
        X=U
        FX=FU
     ELSE
        IF (U.LT.X) THEN
           A=U
        ELSE
           B=U
        ENDIF
        IF (FU.LE.FW.OR.W.EQ.X) THEN
           V=W
           FV=FW
           W=U
           FW=FU
        ELSE IF (FU.LE.FV.OR.V.EQ.X.OR.V.EQ.W) THEN
           V=U
           FV=FU
        ENDIF
     ENDIF
11   CONTINUE
     !Pause ' Brent exceed maximum iterations.'
3    XMIN=X   !exit section
     BRENT=FX
     RETURN
  END
  
end module powell_routines
