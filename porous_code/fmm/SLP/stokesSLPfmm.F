#include <fintrf.h>
C
C stokesSLPfmm.F 
C
C This subroutine is the main gateway to MATLAB.  When a MEX function
C  is executed MATLAB calls the MEXFUNCTION subroutine in the corresponding
C  MEX file.  
C
C
      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER PLHS(*), PRHS(*)
C-----------------------------------------------------------------------
C
c      CHARACTER*120 LINE 
c      INTEGER*4 MEXPRINTF,K
      INTEGER NLHS, NRHS
      INTEGER ERROR
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER MXCREATEDOUBLEMATRIX, MXGETPR

C-----------------------------------------------------------------------
C

      MWSIZE MXGETM, MXGETN
C
C KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS FOR USE
C IN ALL YOUR FORTRAN MEX FILES.
C---------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER F1,F2,XS,YS,XT,YT,IEQUAL,IPREC
      MWPOINTER U1,U2
C-----------------------------------------------------------------------
C

      MWSIZE MS, MT, NS
      REAL *8 RIEQUAL,RIPREC
      REAL *8, ALLOCATABLE :: RF1(:),RF2(:)
      REAL *8, ALLOCATABLE :: RXS(:),RYS(:)
      REAL *8, ALLOCATABLE :: RXT(:),RYT(:)
      REAL *8, ALLOCATABLE :: RU1(:),RU2(:)

#if defined MSWIND
C     For Windows only!
C     This resets the floating point exception to allow divide by zero,
C     overflow and invalid numbers. 
C     
      INTEGER(2) CONTROL
      CALL GETCONTROLFPQQ(CONTROL)
      CONTROL = CONTROL .OR. FPCW$ZERODIVIDE
      CONTROL = CONTROL .OR. FPCW$INVALID
      CONTROL = CONTROL .OR. FPCW$OVERFLOW
      CALL SETCONTROLFPQQ(CONTROL)
#endif

C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
c      IF (NRHS .NE. 8) THEN
c        CALL MEXERRMSGTXT('FMM requires four input arguments')
c      ELSEIF (NLHS .NE. 2) THEN
c        CALL MEXERRMSGTXT('FMM Requires two output argument')
c      ENDIF

c     GET THE NUMBER OF SOURCES.
      MS = MXGETM(PRHS(1))
      NS = MXGETN(PRHS(1));
      MS = MAX(MS,NS)
c     MS IS THE NUMBER OF SOURCES

c     GET THE NUMBER OF TARGETS.
      MT = MXGETM(PRHS(5))
      NT = MXGETN(PRHS(5));
      MT = MAX(MT,NT)
c     MT IS THE NUMBER OF TARGETS

      ALLOCATE(RF1(MS),STAT=ERROR)
      ALLOCATE(RF2(MS),STAT=ERROR)
      ALLOCATE(RXS(MS),STAT=ERROR)
      ALLOCATE(RYS(MS),STAT=ERROR)
      ALLOCATE(RXT(MT),STAT=ERROR)
      ALLOCATE(RYT(MT),STAT=ERROR)
      ALLOCATE(RU1(MT),STAT=ERROR)
      ALLOCATE(RU2(MT),STAT=ERROR)
C
C CREATE A MATRIX FOR RETURN ARGUMENT
      PLHS(1) = MXCREATEDOUBLEMATRIX(MT,1,0)
      PLHS(2) = MXCREATEDOUBLEMATRIX(MT,1,0)
C
C ASSIGN POINTERS TO THE VARIOUS PARAMETERS
C
      F1 = MXGETPR(PRHS(1))
      F2 = MXGETPR(PRHS(2))
      XS = MXGETPR(PRHS(3))
      YS = MXGETPR(PRHS(4))
      XT = MXGETPR(PRHS(5))
      YT = MXGETPR(PRHS(6))
      IEQUAL = MXGETPR(PRHS(7))
      IPREC = MXGETPR(PRHS(8))
C
      U1 = MXGETPR(PLHS(1))
      U2 = MXGETPR(PLHS(2))
C
C COPY RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES
      CALL MXCOPYPTRTOREAL8(F1, RF1, MS)
      CALL MXCOPYPTRTOREAL8(F2, RF2, MS)
      CALL MXCOPYPTRTOREAL8(XS, RXS, MS)
      CALL MXCOPYPTRTOREAL8(YS, RYS, MS)
      CALL MXCOPYPTRTOREAL8(XT, RXT, MT)
      CALL MXCOPYPTRTOREAL8(YT, RYT, MT)
      CALL MXCOPYPTRTOREAL8(IEQUAL, RIEQUAL, 1)
      CALL MXCOPYPTRTOREAL8(IPREC, RIPREC, 1)

C
C DO THE ACTUAL COMPUTATIONS IN A SUBROUTINE
C       CREATED ARRAYS.
C

c      WRITE(LINE,*) 'ARG = ',MT
c      K = MEXPRINTF(LINE//ACHAR(13))
      CALL stokesSLP(RF1,RF2,RXS,RYS,MS,RXT,RYT,MT,
     $      RIEQUAL,RIPREC,RU1,RU2)
c      WRITE(LINE,*) 'ARG = ',MT
c      K = MEXPRINTF(LINE//ACHAR(13))

C
C COPY OUTPUT WHICH IS STORED IN LOCAL ARRAY TO MATRIX OUTPUT
      CALL MXCOPYREAL8TOPTR(RU1, U1, MT)
      CALL MXCOPYREAL8TOPTR(RU2, U2, MT)
C

      DEALLOCATE(RF1)
      DEALLOCATE(RF2)
      DEALLOCATE(RXS)
      DEALLOCATE(RYS)
      DEALLOCATE(RXT)
      DEALLOCATE(RYT)
      DEALLOCATE(RU1)
      DEALLOCATE(RU2)

      RETURN
      END
