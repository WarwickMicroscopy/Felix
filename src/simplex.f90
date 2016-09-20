SUBROUTINE NDimensionalDownhillSimplex(RSimplexVariable,y,mp,np,ndim,ftol,iter,RStandardDeviation,RMean,IErr)

  USE MyNumbers
    
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: iter,mp,ndim,np,NMAX,ITMAX,IErr
  REAL(RKIND) :: ftol,RSimplexVariable(mp,np),y(mp),SimplexExtrapolate,RSendPacket(ndim+2),RExitFlag
  REAL(RKIND) :: rtol,Rsum,swap,ysave,Rytry,psum(ndim),amotry,RStandardDeviation,RMean,RStandardError,RStandardTolerance
  PARAMETER (NMAX=1000,ITMAX=50000)

  INTEGER(IKIND) :: i,ihi,ilo,inhi,j,m,n,IExitFlag
  CHARACTER*200 :: SPrintString
  
  Rytry=ZERO!initial value, has no significance

  IF(my_rank.EQ.0) THEN
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
      CALL SimulateAndFit(Rytry,psum,iter,1,IErr)
      RETURN
    END IF
     
    IF (iter.GE.ITMAX) THEN!We have reached the iteration limit, finish off
      psum = RESHAPE(RSimplexVariable(MAXLOC(y),:),SHAPE(psum)) ! psum = simplex point with highest correlation
      IErr = 1
      RSendPacket = [-10000.0_RKIND, psum, REAL(iter,RKIND)]
      CALL MPI_BCAST(RSendPacket,ndim+2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
      WRITE(SPrintString,FMT='(A22,I3,A11)') "Simplex halted after  ",iter," iterations"
      PRINT*,TRIM(ADJUSTL(SPrintString))
      RETURN
    END IF
     
    !CALL SaveSimplex(RSimplexVariable,y,np,RStandardDeviation,RMean,iter,IErr)
    PRINT*,"------------------------------------------"
    IF (iter.EQ.1) THEN    
      WRITE(SPrintString,FMT='(A15)') "First iteration"
	ELSE IF (iter.LT.10) THEN
      WRITE(SPrintString,FMT='(A10,I1,A18,F7.5)')&
	  "Iteration ",iter,", best fit so far ",y(ilo)
	ELSE IF (iter.LT.100) THEN
      WRITE(SPrintString,FMT='(A10,I2,A18,F7.5)')&
	  "Iteration ",iter,", best fit so far ",y(ilo)
	ELSE IF (iter.LT.1000) THEN
      WRITE(SPrintString,FMT='(A10,I3,A18,F7.5)')&
	  "Iteration ",iter,", best fit so far ",y(ilo)
	ELSE
      WRITE(SPrintString,FMT='(A10,I5,A18,F7.5)')&
	  "Iteration ",iter,", best fit so far ",y(ilo)
	END IF
    PRINT*,TRIM(ADJUSTL(SPrintString))
    WRITE(SPrintString,FMT='(A14,F7.5,A14,F7.5)') "Simplex range ",rtol,", will end at ",ftol
    PRINT*,TRIM(ADJUSTL(SPrintString))
    PRINT*,"------------------------------------------"
	
    iter=iter+2
    !begin a new iteration, reflect the simplex from the high point
    Rytry = SimplexExtrapolate(RSimplexVariable,y,psum,mp,np,ndim,ihi,-1.0D0,iter,IErr)
    PRINT*,"Simplex reflection of worst point:"
    IF (Rytry.LE.y(ilo)) THEN !the reflected point is better than the best point, extrapolate by 2 times again
      Rytry = SimplexExtrapolate(RSimplexVariable,y,psum,mp,np,ndim,ihi,2.0D0,iter,IErr)
      PRINT*,"Best fit so far, going further:"
    ELSEIF (Rytry.GE.y(inhi)) THEN !the reflected point is worse than the second-highest, so look for an intermediate lower point 
      ysave=y(ihi)
      Rytry=SimplexExtrapolate(RSimplexVariable,y,psum,mp,np,ndim,ihi,0.5D0,iter,IErr)
      PRINT*,"Poor result, interpolating for a better one:"
      IF(Rytry.GE.ysave) THEN !can't get rid of the highest point, contract about the best point
        PRINT*,"-----------------------------------------------------"
        PRINT*,"Entering Contraction Phase, Expect",ndim+1,"Simulations"
        PRINT*,"-----------------------------------------------------"
        DO i=1,ndim+1
          PRINT*,"Contraction Simulation",i
          IF(i.NE.ilo) THEN
            DO j=1,ndim
              psum(j)=0.5*(RSimplexVariable(i,j)+RSimplexVariable(ilo,j))
              RSimplexVariable(i,j)=psum(j)
            ENDDO
            RSendPacket = [10000.0_RKIND, psum, REAL(iter,RKIND)]
            CALL MPI_BCAST(RSendPacket,ndim+2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
            CALL SimulateAndFit(y(i),psum,iter,0,IErr)
          ENDIF
        ENDDO
        iter=iter+ndim
        GOTO 1 !go back to see if we are done and the next iteration
      ENDIF
    ELSE
      iter=iter-1
    ENDIF
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
      CALL SimulateAndFit(Rytry,psum,iter,IExitFLAG,IErr) ! Doesnt matter what this result is
      IF(IExitFLAG.EQ.1) RETURN
    END DO

  END IF
  
END SUBROUTINE NDimensionalDownhillSimplex

!!$----------------------------------------------------------------------------

REAL(RKIND) FUNCTION SimplexExtrapolate(RSimplexVariable,y,psum,mp,np,ndim,ihi,fac,iter,IErr)

  USE MyNumbers
    
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: ihi,mp,ndim,np,NMAX,IErr,iter,j
  REAL(RKIND) :: fac,RSimplexVariable(mp,np),psum(np),y(mp),RSendPacket(ndim+2)
  REAL(RKIND) :: fac1,fac2,Rytry,ptry(ndim)
  PARAMETER(NMAX=1000)
  CHARACTER*200 :: SPrintString
  
  !extrapolates by a factor fac through the face of the simplex across from the high point, tries it, and replaces the high point if the new point is better
  fac1=(1.0-fac)/ndim
  fac2=fac1-fac
  DO j=1,ndim
     ptry(j)=psum(j)*fac1-RSimplexVariable(ihi,j)*fac2
  ENDDO
  RSendPacket = [10000.0_RKIND, ptry, REAL(iter,RKIND)]
  CALL MPI_BCAST(RSendPacket,ndim+2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
  CALL SimulateAndFit(Rytry,ptry,iter,0,IErr)
  
  IF (Rytry.LT.y(ihi)) THEN
     y(ihi)=Rytry
     DO j=1,ndim
        psum(j)=psum(j)-RSimplexVariable(ihi,j)+ptry(j)
        RSimplexVariable(ihi,j)=ptry(j)
     ENDDO
  ENDIF

  SimplexExtrapolate=Rytry

  RETURN
END FUNCTION SimplexExtrapolate


!!$----------------------------------------------------------------------------

SUBROUTINE OpenSimplexOutput(IErr)

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


!!$----------------------------------------------------------------------------

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

!!$----------------------------------------------------------------------------

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


!!$----------------------------------------------------------------------------
