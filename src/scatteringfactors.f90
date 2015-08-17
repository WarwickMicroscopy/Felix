SUBROUTINE ScatteringFactors 

  USE MyNumbers
  USE WriteToScreen

  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  REAL(RKIND) :: &
     Kirkland(103,13), Peng(103,8), DoyleAndTurner(103,9)
  
  DATA Kirkland(1,13)/
  DATA Kirkland(2,13)/
  
