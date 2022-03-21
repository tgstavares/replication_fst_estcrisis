      SUBROUTINE SD(X,N,IWRITE,XSD)
C
C     PURPOSE--THIS SUBROUTINE COMPUTES THE
C              STANDARD DEVIATION (WITH DENOMINATOR N)
C              OF THE DATA IN THE INPUT VECTOR X. 
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
C                                (UNSORTED OR SORTED) OBSERVATIONS.
C                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
C                                IN THE VECTOR X. 
C                     --IWRITE = AN INTEGER FLAG CODE WHICH 
C                                (IF SET TO 0) WILL SUPPRESS
C                                THE PRINTING OF THE
C                                STANDARD DEVIATION
C                                AS IT IS COMPUTED;
C                                OR (IF SET TO SOME INTEGER 
C                                VALUE NOT EQUAL TO 0),
C                                LIKE, SAY, 1) WILL CAUSE
C                                THE PRINTING OF THE
C                                STANDARD DEVIATION
C                                AT THE TIME IT IS COMPUTED.
C     OUTPUT ARGUMENTS--XSD    = THE SINGLE PRECISION VALUE OF THE
C                                COMPUTED STANDARD DEVIATION.
C     OUTPUT--THE COMPUTED SINGLE PRECISION VALUE OF THE
C             STANDARD DEVIATION
C     PRINTING--NONE, UNLESS IWRITE HAS BEEN SET TO A NON-ZERO
C               INTEGER, OR UNLESS AN INPUT ARGUMENT ERROR
C               CONDITION EXISTS.
C     RESTRICTIONS--THERE IS NO RESTRICTION ON THE MAXIMUM VALUE
C                   OF N FOR THIS SUBROUTINE.
C     OTHER DATAPAC   SUBROUTINES NEEDED--NONE.
C     FORTRAN LIBRARY SUBROUTINES NEEDED--SQRT.
C     MODE OF INTERNAL OPERATIONS--SINGLE PRECISION.
C     LANGUAGE--ANSI FORTRAN. 
C     REFERENCES--SNEDECOR AND COCHRAN, STATISTICAL METHODS,
C                 EDITION 6, 1967, PAGE 44.
C               --DIXON AND MASSEY, INTRODUCTION TO STATISTICAL
C                 ANALYSIS, EDITION 2, 1957, PAGES 19, 76.
C     WRITTEN BY--JAMES J. FILLIBEN
C                 STATISTICAL ENGINEERING LABORATORY (205.03)
C                 NATIONAL BUREAU OF STANDARDS
C                 WASHINGTON, D. C. 20234
C                 PHONE:  301-921-2315
C     ORIGINAL VERSION--JUNE      1972. 
C     UPDATED         --SEPTEMBER 1975. 
C     UPDATED         --NOVEMBER  1975. 
C
C---------------------------------------------------------------------
C
      DIMENSION X(1)
C
      IPR=6
C
C     CHECK THE INPUT ARGUMENTS FOR ERRORS
C
      AN=N
      IF(N.LT.1)GOTO50
      IF(N.EQ.1)GOTO55
      HOLD=X(1)
      DO60I=2,N
      IF(X(I).NE.HOLD)GOTO90
   60 CONTINUE
c$$$     WRITE(IPR, 9)HOLD
      XSD=0.0
      GOTO201
   50 WRITE(IPR,15) 
      WRITE(IPR,47)N
      RETURN
   55 WRITE(IPR,18) 
      XSD=0.0
      GOTO201
   90 CONTINUE
    9 FORMAT(1H ,109H***** NON-FATAL DIAGNOSTIC--THE FIRST  INPUT ARGUME
     1NT (A VECTOR) TO THE SD     SUBROUTINE HAS ALL ELEMENTS = ,E15.8,6
     1H *****)
   15 FORMAT(1H , 91H***** FATAL ERROR--THE SECOND INPUT ARGUMENT TO THE
     1 SD     SUBROUTINE IS NON-POSITIVE *****)
   18 FORMAT(1H ,100H***** NON-FATAL DIAGNOSTIC--THE SECOND INPUT ARGUME
     1NT TO THE SD     SUBROUTINE HAS THE VALUE 1 *****)
   47 FORMAT(1H , 35H***** THE VALUE OF THE ARGUMENT IS ,I8   ,6H *****)
C
C-----START POINT-----------------------------------------------------
C
      SUM=0.0
      DO100I=1,N
      SUM=SUM+X(I)
  100 CONTINUE
      XMEAN=SUM/AN
      SUM=0.0
      DO200I=1,N
      SUM=SUM+(X(I)-XMEAN)**2 
  200 CONTINUE
      VAR=SUM/AN
      XSD=SQRT(VAR) 
C
  201 IF(IWRITE.EQ.0)RETURN
      WRITE(IPR,999)
      WRITE(IPR,205)N,XSD
  205 FORMAT(1H ,37HTHE SAMPLE STANDARD DEVIATION OF THE ,I6,17H OBSERVA
     1TIONS IS ,E15.8)
  999 FORMAT(1H )
      RETURN
      END 