*
      SUBROUTINE BI_F77_MPI_INITIALIZED(FLAG, IERR)
      INCLUDE 'mpif.h'
      INTEGER FLAG, IERR
      LOGICAL FLAG2
*
      CALL MPI_INITIALIZED(FLAG2, IERR)
      IF (FLAG2) THEN
         FLAG = 1
      ELSE
         FLAG = 0
      END IF
*
      RETURN
      END
