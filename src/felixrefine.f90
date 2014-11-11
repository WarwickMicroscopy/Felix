!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixrefine
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all right reserved
!
! Version: :VERSION:
! Date:    :DATE:
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build:   :BUILD:
! Author:  :AUTHOR:
! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  This file is part of felixrefine.
!
!  felixrefine is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  felixrefine is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with felixrefine.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: Felixrefine.f90,v 1.89 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PROGRAM Felixrefine
 
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  !--------------------------------------------------------------------
  ! local variable definitions
  !--------------------------------------------------------------------
  
  IMPLICIT NONE


  INTEGER(IKIND) :: &
       IHours,IMinutes,ISeconds,IErr
  REAL(RKIND) :: &
       StartTime, CurrentTime, Duration, TotalDurationEstimate
  

  !-------------------------------------------------------------------
  ! constants
  !-------------------------------------------------------------------

  CALL Init_Numbers
  
  !-------------------------------------------------------------------
  ! set the error value to zero, will change upon error
  !-------------------------------------------------------------------

  IErr=0

  !--------------------------------------------------------------------
  ! MPI initialization
  !--------------------------------------------------------------------

  ! Initialise MPI  
  CALL MPI_Init(IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error in MPI_Init()"
     GOTO 9999
  ENDIF

  ! Get the rank of the current process
  CALL MPI_Comm_rank(MPI_COMM_WORLD,my_rank,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error in MPI_Comm_rank()"
     GOTO 9999
  ENDIF

  ! Get the size of the current communicator
  CALL MPI_Comm_size(MPI_COMM_WORLD,p,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error in MPI_Comm_size()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! protocal feature startup
  !--------------------------------------------------------------------
  
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"--------------------------------------------------------------"
     PRINT*,"Felixrefine: ", RStr
     PRINT*,"          ", DStr
     PRINT*,"          ", AStr
     PRINT*,"          on rank= ", my_rank, " of ", p, " in total."
     PRINT*,"--------------------------------------------------------------"
  END IF

  !--------------------------------------------------------------------
  ! timing startup
  !--------------------------------------------------------------------

  CALL cpu_time(StartTime)

  !--------------------------------------------------------------------
  ! INPUT section 
  !--------------------------------------------------------------------
  
  ISoftwareMode = 0 ! felixrefinemode
  
  !Read from input files
  CALL ReadInput (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error in ReadInput()"
     GOTO 9999
  ENDIF

  CALL FelixFunction(IErr)
  
  !--------------------------------------------------------------------
  ! finish off
  !--------------------------------------------------------------------
    
  CALL cpu_time(CurrentTime)
  Duration=(CurrentTime-StartTime)
  IHours = FLOOR(Duration/3600.0D0)
  IMinutes = FLOOR(MOD(Duration,3600.0D0)/60.0D0)
  ISeconds = MOD(Duration,3600.0D0)-IMinutes*60.0D0

  PRINT*, "Felixrefine(", my_rank, ") ", RStr, ", used time=", IHours, "hrs ",IMinutes,"mins ",ISeconds,"Seconds "

  !--------------------------------------------------------------------
  ! Shut down MPI
  !--------------------------------------------------------------------W
9999 &
  CALL MPI_Finalize(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error ", IErr, " in MPI_Finalize()"
     STOP
  ENDIF

  ! clean shutdown
  STOP
  
!!$800 PRINT*,"Felixrefine(", my_rank, "): ERR in CLOSE()"
!!$  IErr= 1
!!$  RETURN

END PROGRAM Felixrefine
