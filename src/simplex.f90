!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-17, all rights reserved
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
!  Felix is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  Felix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with Felix.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! A lot of these are currently out of use, but may be useful

! All procedures conatained in this file:
! NDimensionalDownhillSimplex()
! SimplexExtrapolate()
! OpenSimplexOutput()
! WriteOutSimplex()
! SaveSimplex() 


!>
!! Procedure-description: 
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE NDimensionalDownhillSimplex(RSimplexVariable,y,mp,np,ndim,ftol,iter,RStandardDeviation,RMean,IErr)

  USE MyNumbers
    
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE terminal_output
  USE MPI
  USE MyMPI

  USE felixfunction_mod

  IMPLICIT NONE

  INTEGER(IKIND) :: iter,mp,ndim,np,NMAX,ITMAX,IErr
  REAL(RKIND) :: ftol,RSimplexVariable(mp,np),y(mp),SimplexExtrapolate,RSendPacket(ndim+2),RExitFlag
  REAL(RKIND) :: rtol,Rsum,swap,ysave,Rytry,psum(ndim),amotry,RStandardDeviation,RMean,RStandardError,RStandardTolerance
  PARAMETER (NMAX=1000,ITMAX=50000)

  INTEGER(IKIND) :: i,ihi,ilo,inhi,j,m,n,IExitFlag
  
  Rytry=ZERO!initial value, has no significance

  IF(my_rank.EQ.0) THEN!why is this here, should be outside the subroutine?
1   DO n=1,ndim !enter here when starting or have just overall contracted
      Rsum=0 !recalculate psum
      DO m=1,ndim+1
        Rsum=Rsum+RSimplexVariable(m,n)
      ENDDO
      psum(n)=Rsum
    ENDDO
2   ilo=1 !enter here when have just changed a single point
    ysave=Rytry
    IF (y(1).GT.y(2)) THEN !Determine which point is highest (worst), next highest, and lowest (best)
      ihi=1
      inhi=2
    ELSE
      ihi=2
      inhi=1
    END IF
    DO i=1,ndim+1 !by looping over points in the simplex
      IF(y(i).LE.y(ilo)) THEN
        ilo=i!the best point
      END IF
      IF(y(i).GT.y(ihi)) THEN
        inhi=ihi
        ihi=i!the worst point
      ELSE IF(y(i).GT.y(inhi)) THEN
        IF(i.NE.ihi) inhi=i
      END IF
    ENDDO
	  !compute the range from highest to lowest
    rtol=2.*ABS(y(ihi)-y(ilo))/(ABS(y(ihi))+ABS(y(ilo)))

    IF(rtol.LT.ftol) THEN !returning, put the best point in slot 1
      swap=y(1)
      y(1)=y(ilo)
      y(ilo)=swap
      DO n=1,ndim
        swap=RSimplexVariable(1,n)
        RSimplexVariable(1,n)=RSimplexVariable(ilo,n)
        RSimplexVariable(ilo,n)=swap
      END DO
      psum = RESHAPE(RSimplexVariable(MAXLOC(y),:),SHAPE(psum)) ! psum = simplex point with highest correlation
      RSendPacket = [-10000.0_RKIND, psum, REAL(iter,RKIND)]
      CALL MPI_BCAST(RSendPacket,ndim+2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
      !Simulate-------------------------------------------------------------
      CALL SimulateAndFit(psum,iter,1,IErr)
      Rytry=RFigureofMerit
      RETURN
    END IF
     
    IF (iter.GE.ITMAX) THEN!We have reached the iteration limit, finish off
      psum = RESHAPE(RSimplexVariable(MAXLOC(y),:),SHAPE(psum)) ! psum = simplex point with highest correlation
      IErr = 1
      RSendPacket = [-10000.0_RKIND, psum, REAL(iter,RKIND)]
      CALL MPI_BCAST(RSendPacket,ndim+2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
      CALL message( LS, "Simplex halted after  iterations = ",iter )
      RETURN
    END IF
     
    !CALL SaveSimplex(RSimplexVariable,y,np,RStandardDeviation,RMean,iter,IErr)
    CALL message( LM, "------------------------------------------")

    CALL message( LM, "Iteration = ",iter)
    CALL message( LM, " current best fit = ",y(ilo) )

    CALL message( LM, "Simplex range ",rtol)
    CALL message( LM, "    will end at ",ftol)

    CALL message( LM, "------------------------------------------")
	
    iter=iter+2
    !begin a new iteration, reflect the simplex from the high point
    Rytry = SimplexExtrapolate(RSimplexVariable,y,psum,mp,np,ndim,ihi,-1.0D0,iter,IErr)
    CALL message( LM, "Simplex reflection:" )
    IF (Rytry.LE.y(ilo)) THEN !the reflected point is better than the best point, extrapolate by 2 times again
      Rytry = SimplexExtrapolate(RSimplexVariable,y,psum,mp,np,ndim,ihi,2.0D0,iter,IErr)
      CALL message( LM, "Extrapolation:" )
    ELSEIF (Rytry.GE.y(inhi)) THEN !the reflected point is worse than the second-highest, so look for an intermediate lower point 
      ysave=y(ihi)
      Rytry=SimplexExtrapolate(RSimplexVariable,y,psum,mp,np,ndim,ihi,0.5D0,iter,IErr)
      CALL message( LM, "Interpolation:" )
      IF(Rytry.GE.ysave) THEN !can't get rid of the highest point, contract about the best point
        CALL message( LM, "-----------------------------------------------------")
        CALL message( LM, "Entering Contraction Phase, Expect number Simulations = ",ndim+1 )
        CALL message( LM, "-----------------------------------------------------")
        DO i=1,ndim+1
          CALL message( LM, "Contraction Simulation",i )
          IF(i.NE.ilo) THEN
            DO j=1,ndim
              psum(j)=0.5*(RSimplexVariable(i,j)+RSimplexVariable(ilo,j))
              RSimplexVariable(i,j)=psum(j)
            ENDDO
            RSendPacket = [10000.0_RKIND, psum, REAL(iter,RKIND)]
            CALL MPI_BCAST(RSendPacket,ndim+2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
            !Simulate-------------------------------------------------------------
            CALL SimulateAndFit(psum,iter,0,IErr)
            y(i)=RFigureofMerit
          END IF
        ENDDO
        iter=iter+ndim
        GOTO 1 !go back to see if we are done and the next iteration
      END IF
    ELSE
      iter=iter-1
    END IF
    GOTO 2
  ELSE
    DO!Latch to loop cores other than zero waiting for MPI_BCAST (is it really necessary) 
      CALL MPI_BCAST(RSendPacket,ndim+2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
      RExitFlag = RSendPacket(1)                
      IF(RExitFlag.LT.ZERO) THEN
        IExitFLAG = 1
      ELSE
        IExitFLAG = 0
      END IF
      psum = RSendPacket(2:(ndim+1))
      iter = NINT(RSendPacket(ndim+2),KIND=IKIND)
      !Simulate-------------------------------------------------------------
      CALL SimulateAndFit(psum,iter,IExitFLAG,IErr)
      Rytry=RFigureofMerit
      IF(IExitFLAG.EQ.1) RETURN
    END DO

  END IF
  
END SUBROUTINE NDimensionalDownhillSimplex

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!>
!! Procedure-description: 
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
REAL(RKIND) FUNCTION SimplexExtrapolate(RSimplexVariable,y,psum,mp,np,ndim,ihi,fac,iter,IErr)

  USE MyNumbers
    
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  USE felixfunction_mod

  IMPLICIT NONE
  
  INTEGER(IKIND) :: ihi,mp,ndim,np,NMAX,IErr,iter,j
  REAL(RKIND) :: fac,RSimplexVariable(mp,np),psum(np),y(mp),RSendPacket(ndim+2)
  REAL(RKIND) :: fac1,fac2,Rytry,ptry(ndim)
  PARAMETER(NMAX=1000)
  
  !extrapolates by a factor fac through the face of the simplex across from the high point, tries it, and replaces the high point if the new point is better
  fac1=(1.0-fac)/ndim
  fac2=fac1-fac
  DO j=1,ndim
     ptry(j)=psum(j)*fac1-RSimplexVariable(ihi,j)*fac2
  ENDDO
  RSendPacket = [10000.0_RKIND, ptry, REAL(iter,RKIND)]
  CALL MPI_BCAST(RSendPacket,ndim+2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
  !Simulate-------------------------------------------------------------
  CALL SimulateAndFit(ptry,iter,0,IErr)
  Rytry=RFigureofMerit
      
  IF (Rytry.LT.y(ihi)) THEN
     y(ihi)=Rytry
     DO j=1,ndim
        psum(j)=psum(j)-RSimplexVariable(ihi,j)+ptry(j)
        RSimplexVariable(ihi,j)=ptry(j)
     ENDDO
  END IF

  SimplexExtrapolate=Rytry

  RETURN
END FUNCTION SimplexExtrapolate

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!>
!! Procedure-description: 
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE OpenSimplexOutput(IErr)
!ffs, 1 line of code
  USE MyNumbers

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr

  CHARACTER*200 :: filename

  WRITE(filename,*) "fr-Simplex.txt"

  OPEN( UNIT=IChOutSimplex,STATUS='UNKNOWN',FILE=TRIM( ADJUSTL(filename)) )

END SUBROUTINE OpenSimplexOutput

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!>
!! Procedure-description: 
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE WriteOutSimplex(RSimplexVariable,RSimplexFoM,IDimensions,RStandardDeviation,RMean,IIterations,IErr)

  USE MyNumbers

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,IDimensions,ind,IIterations
  REAL(RKIND),DIMENSION(IDimensions+1,IDimensions),INTENT(IN) :: RSimplexVariable
  REAL(RKIND),DIMENSION(IDimensions+1),INTENT(IN) :: RSimplexFoM
  REAL(RKIND),DIMENSION(IDimensions+1) :: RData
  REAL(RKIND) :: RStandardDeviation,RMean
  CHARACTER*200 :: CSizeofData,SFormatString
  
  WRITE(CSizeofData,*) IDimensions+1
  WRITE(SFormatString,*) "("//TRIM(ADJUSTL(CSizeofData))//"(1F6.3,1X),A1)"

  DO ind = 1,(IDimensions+1)
     RData = (/RSimplexVariable(ind,:), RSimplexFoM(ind)/)
     WRITE(IChOutSimplex,FMT=SFormatString) RData
  END DO

  WRITE(IChOutSimplex,FMT="(2(1F6.3,1X),I5.1,I5.1,A1)") RStandardDeviation,RMean,IStandardDeviationCalls,IIterations

  CLOSE(IChOutSimplex)

END SUBROUTINE WriteOutSimplex

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!>
!! Procedure-description: 
!!
!! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
!!
SUBROUTINE SaveSimplex(RSimplexVariable,RSimplexFoM,IDimensions,RStandardDeviation,RMean,IIterations,IErr)
  !what a useless subroutine, just calls two others
  USE MyNumbers

  USE IConst; USE RConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,IDimensions,IIterations
  REAL(RKIND),DIMENSION(IDimensions+1,IDimensions),INTENT(IN) :: RSimplexVariable
  REAL(RKIND),DIMENSION(IDimensions+1),INTENT(IN) :: RSimplexFoM
  REAL(RKIND) :: RStandardDeviation,RMean

  CALL OpenSimplexOutput(IErr)
 
  CALL WriteOutSimplex(RSimplexVariable,RSimplexFoM,IDimensions,RStandardDeviation,RMean,IIterations,IErr)

END SUBROUTINE SaveSimplex
