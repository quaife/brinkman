#include <fintrf.h>
C
C stokesDLPnewFMM.F
C
C
      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
C-------------------------------------------------------------------------
C   (pointer) Replace integer by integer*8 on 64-bit platformsi
C
      IMPLICIT NONE

      MWPOINTER PLHS(*), PRHS(*)

C------------------------------------------------------------------------
C

      INTEGER NLHS, NRHS
      INTEGER ERROR

C
C-----------------------------------------------------------------------
C (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER MXGETPR
      MWPOINTER MXGETPI
      MWPOINTER MXCREATEDOUBLEMATRIX

C-----------------------------------------------------------------------
C

      MWPOINTER MXGETM, MXGETN

C KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS FOR USE
C IN ALL YOUR FORTRAN MEX FILES.
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
C
      MWPOINTER REDIP1,IMDIP1,REDIP2,IMDIP2,XS,YS
      MWPOINTER REVEL,IMVEL

C-----------------------------------------------------------------------
C

      MWSIZE MS, NS, NEL
      REAL *8, ALLOCATABLE :: LXS(:),LYS(:)
      COMPLEX *16, ALLOCATABLE :: LDIP1(:),LDIP2(:)
      COMPLEX *16, ALLOCATABLE :: LVEL(:)
C----------------------------------------------------------------------
C
C     GET THE NUMBER OF SOURCES
      MS = MXGETM(PRHS(1))
      NS = MXGETN(PRHS(1))
      MS = MAX(MS,NS)
C
      ALLOCATE(LXS(MS),STAT=ERROR)
      ALLOCATE(LYS(MS),STAT=ERROR)
      ALLOCATE(LDIP1(MS),STAT=ERROR) 
      ALLOCATE(LDIP2(MS),STAT=ERROR) 
      ALLOCATE(LVEL(MS),STAT=ERROR)
C
C CREATE A MATRIX FOR RETURN ARGUMENT
      PLHS(1) = MXCREATEDOUBLEMATRIX(MS,1,1) 
C ASSIGN POINTERS TO THE VARIOUS PARAMETERS
C
      REDIP1 = MXGETPR(PRHS(1))
      IMDIP1 = MXGETPI(PRHS(1))
      REDIP2 = MXGETPR(PRHS(2))
      IMDIP2 = MXGETPI(PRHS(2))
      XS     = MXGETPR(PRHS(3))
      YS     = MXGETPR(PRHS(4))
C      
      REVEL = MXGETPR(PLHS(1))
      IMVEL = MXGETPI(PLHS(1))
C      
C
C COPY RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES
      CALL MXCOPYPTRTOREAL8(XS,LXS, MS)
      CALL MXCOPYPTRTOREAL8(YS,LYS, MS)
      CALL MXCOPYPTRTOCOMPLEX16(REDIP1,IMDIP1,LDIP1,MS)
      CALL MXCOPYPTRTOCOMPLEX16(REDIP2,IMDIP2,LDIP2,MS)
      
C APPLY FMM
      CALL STOKESDLPNEW(MS,LXS,LYS,LDIP1,LDIP2,LVEL)
C      
C COPY OUTPUT WHICH IS STORED IN LOCAL ARRAY TO MATRIX OUTPUT
C
      CALL MXCOPYCOMPLEX16TOPTR(LVEL,REVEL,IMVEL,MS)
C
      DEALLOCATE(LXS)
      DEALLOCATE(LYS) 
      DEALLOCATE(LDIP1)
      DEALLOCATE(LDIP2)
      DEALLOCATE(LVEL)      
     
      RETURN
      END SUBROUTINE


           

            

      


      

