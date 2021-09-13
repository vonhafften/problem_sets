PROGRAM hello

  IMPLICIT NONE

  REAL :: temp
  REAL :: temp2

  READ(*,*)  temp
  READ(*,*)  temp2
  WRITE(*,*) 'the answer is'
  WRITE(*,*) temp * temp2

END PROGRAM hello
