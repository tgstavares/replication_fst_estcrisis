      SUBROUTINE DBSPDR(T,A,N,K,NDERIV,AD)
C***BEGIN PROLOGUE  DBSPDR
C***DATE WRITTEN   800901   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  E3,K6
C***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
C             SPLINE
C***AUTHOR  AMOS, D. E., (SNLA)
C***PURPOSE  Uses the B-representation to construct a divided difference
C            table ad preparatory to a (right) derivative calculation
C            in BSPEV.
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Reference
C         SIAM J. Numerical Analysis, 14, No. 3, June, 1977, pp.441-472.
C
C     Abstract     **** a double precision routine ****
C         DBSPDR is the BSPLDR routine of the reference.
C
C         DBSPDR uses the B-representation (T,A,N,K) to construct a
C         divided difference table ADIF preparatory to a (right)
C         derivative calculation in DBSPEV.  The lower triangular matrix
C         ADIF is stored in vector AD by columns.  The arrays are
C         related by
C
C           ADIF(I,J) = AD(I-J+1 + (2*N-J+2)*(J-1)/2)
C
C         I = J,N   ,   J=1,NDERIV.
C
C     Description of Arguments
C
C         Input      T,A are double precision
C          T       - knot vector of length N+K
C          A       - B-spline coefficient vector of length N
C          N       - number of B-spline coefficients
C                    N = sum of knot multiplicities-K
C          K       - order of the spline, K .GE. 1
C          NDERIV  - number of derivatives, 1 .LE. NDERIV .LE. K.
C                    NDERIV=1 gives the zero-th derivative =
C                    function value
C
C         Output     AD is double precision
C          AD      - table of differences in a vector of length
C                    (2*N-NDERIV+1)*NDERIV/2 for input to DBSPEV
C
C     Error Conditions
C         Improper input is a fatal error
C***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
C                 SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
C                 JUNE 1977, PP. 441-472.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  DBSPDR
C
C
      INTEGER I, ID, II, IPKMID, JJ, JM, K, KMID, N, NDERIV
      DOUBLE PRECISION A, AD, DIFF, FKMID, T
C     DIMENSION T(N+K), AD((2*N-NDERIV+1)*NDERIV/2)
      DIMENSION T(1), A(N), AD(1)
C***FIRST EXECUTABLE STATEMENT  DBSPDR
      IF(K.LT.1) GO TO 100
      IF(N.LT.K) GO TO 105
      IF(NDERIV.LT.1 .OR. NDERIV.GT.K) GO TO 110
      DO 10 I=1,N
        AD(I) = A(I)
   10 CONTINUE
      IF (NDERIV.EQ.1) RETURN
      KMID = K
      JJ = N
      JM = 0
      DO 30 ID=2,NDERIV
        KMID = KMID - 1
        FKMID = DBLE(FLOAT(KMID))
        II = 1
        DO 20 I=ID,N
          IPKMID = I + KMID
          DIFF = T(IPKMID) - T(I)
          IF (DIFF.NE.0.0D0) AD(II+JJ) = (AD(II+JM+1)-AD(II+JM))/
     1     DIFF*FKMID
          II = II + 1
   20   CONTINUE
        JM = JJ
        JJ = JJ + N - ID + 1
   30 CONTINUE
      RETURN
C
C
  100 CONTINUE
      CALL XERROR( ' DBSPDR,  K DOES NOT SATISFY K.GE.1', 35, 2, 1)
      RETURN
  105 CONTINUE
      CALL XERROR( ' DBSPDR,  N DOES NOT SATISFY N.GE.K', 35, 2, 1)
      RETURN
  110 CONTINUE
      CALL XERROR( ' DBSPDR,  NDERIV DOES NOT SATISFY 1.LE.NDERIV.LE.K',
     1 50, 2, 1)
      RETURN
      END
