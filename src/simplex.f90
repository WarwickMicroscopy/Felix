SUBROUTINE NDimensionalDownhillSimplex(RSimplexVariable,y,mp,np,ndim,ftol,iter,RStandardDeviation,RMean,IErr)

  USE MyNumbers
    
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: iter,mp,ndim,np,NMAX,ITMAX,IErr
  REAL(RKIND) :: ftol,RSimplexVariable(mp,np),y(mp),SimplexFunction,SimplexExtrapolate,RSendPacket(ndim+2),RExitFlag
  PARAMETER (NMAX=1000,ITMAX=50000)

  INTEGER(IKIND) :: i,ihi,ilo,inhi,j,m,n,IExitFlag
  REAL(RKIND) :: rtol,sum,swap,ysave,ytry,psum(ndim),amotry,&
       RStandardDeviation,RMean,RStandardError,RStandardTolerance
  CHARACTER*200 :: SPrintString

  IF(my_rank.EQ.0) THEN
1    DO n = 1,ndim
        sum = 0
        DO m=1,ndim+1
           sum=sum+RSimplexVariable(m,n)
        ENDDO
        psum(n) = sum
     ENDDO
2    ilo = 1
     ysave = ytry
     IF (y(1).GT.y(2)) THEN
        ihi=1
        inhi=2
     ELSE
        ihi=2
        inhi=1
     END IF
     DO i=1,ndim+1
        IF(y(i).LE.y(ilo)) ilo=i
        IF(y(i).GT.y(ihi)) THEN
           inhi=ihi
           ihi=i
        ELSE IF(y(i).GT.y(inhi)) THEN
           IF(i.NE.ihi) inhi=i
        END IF
     ENDDO

     rtol=2.*ABS(y(ihi)-y(ilo))/(ABS(y(ihi))+ABS(y(ilo)))

     RStandardTolerance = RStandardError(RStandardDeviation,RMean,ytry,IErr)

     WRITE(SPrintString,FMT='(A14,F7.5,A14,F7.5)') "Simplex range ",rtol,", will end at ",ftol
     PRINT*,TRIM(ADJUSTL(SPrintString))
     IF(rtol.LT.ftol) THEN
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
        ytry = SimplexFunction(psum,iter,1,IErr)
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
     PRINT*,"--------------------------------"
     IF (iter.EQ.1) THEN    
       WRITE(SPrintString,FMT='(A15)') "First iteration"
	 ELSE IF (iter.LT.10) THEN
       WRITE(SPrintString,FMT='(A10,I1,A18,F7.5)') "Iteration ",iter,", figure of merit ",ytry
	 ELSE IF (iter.LT.100) THEN
       WRITE(SPrintString,FMT='(A10,I2,A18,F7.5)') "Iteration ",iter,", figure of merit ",ytry
	 ELSE IF (iter.LT.1000) THEN
       WRITE(SPrintString,FMT='(A10,I3,A18,F7.5)') "Iteration ",iter,", figure of merit ",ytry
	 ELSE
       WRITE(SPrintString,FMT='(A10,I5,A18,F7.5)') "Iteration ",iter,", figure of merit ",ytry
	 END IF
     PRINT*,TRIM(ADJUSTL(SPrintString))
     PRINT*,"--------------------------------"
     iter=iter+2
    
     ytry = SimplexExtrapolate(RSimplexVariable,y,psum,mp,np,ndim,ihi,-1.0D0,iter,IErr)

     
     IF (ytry.LE.y(ilo).OR.my_rank.NE.0) THEN
        ytry = SimplexExtrapolate(RSimplexVariable,y,psum,mp,np,ndim,ihi,2.0D0,iter,IErr)
     ELSEIF (ytry.GE.y(inhi)) THEN
        ysave=y(ihi)
        ytry=SimplexExtrapolate(RSimplexVariable,y,psum,mp,np,ndim,ihi,0.5D0,iter,IErr)
        IF(ytry.GE.ysave) THEN
           PRINT*,"-----------------------------------------------------"
           PRINT*,"Entering Expansion Phase Expect",ndim+1,"Simulations"
           PRINT*,"-----------------------------------------------------"
           DO i=1,ndim+1
              PRINT*,"Expansion Simulation",i
              IF(i.NE.ilo) THEN
                 DO j=1,ndim
                    psum(j)=0.5*(RSimplexVariable(i,j)+RSimplexVariable(ilo,j))
                    RSimplexVariable(i,j)=psum(j)
                 ENDDO
                 RSendPacket = [10000.0_RKIND, psum, REAL(iter,RKIND)]
                 CALL MPI_BCAST(RSendPacket,ndim+2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
                 y(i)=SimplexFunction(psum,iter,0,IErr)
              ENDIF
           ENDDO
           iter=iter+ndim
           GOTO 1
        ENDIF
     ELSE
        iter=iter-1
     ENDIF
     GOTO 2
  ELSE
     DO!We have reached the exit criteria, finish off
        RSendPacket = [-10000.0_RKIND, psum, REAL(iter,RKIND)]!RB added this line but no idea what I'm doing 
        CALL MPI_BCAST(RSendPacket,ndim+2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
        RExitFlag = RSendPacket(1)                
        IF(RExitFlag.LT.ZERO) THEN
           IExitFLAG = 1
        ELSE
           IExitFLAG = 0
        END IF
        psum = RSendPacket(2:(ndim+1))
        iter = NINT(RSendPacket(ndim+2),KIND=IKIND)
        ytry = SimplexFunction(psum,iter,IExitFLAG,IErr) ! Doesnt matter what this result is
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
  REAL(RKIND) :: fac,RSimplexVariable(mp,np),psum(np),y(mp),SimplexFunction,RSendPacket(ndim+2)
  REAL(RKIND) :: fac1,fac2,ytry,ptry(ndim)
  PARAMETER(NMAX=1000)
  CHARACTER*200 :: SPrintString
  
  fac1=(1.0-fac)/ndim
  fac2=fac1-fac
  DO j=1,ndim
     ptry(j)=psum(j)*fac1-RSimplexVariable(ihi,j)*fac2
  ENDDO
  RSendPacket = [10000.0_RKIND, ptry, REAL(iter,RKIND)]
  CALL MPI_BCAST(RSendPacket,ndim+2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
  
  ytry=SimplexFunction(ptry,iter,0,IErr)
  
  IF (ytry.LT.y(ihi)) THEN
     y(ihi)=ytry
     DO j=1,ndim
        psum(j)=psum(j)-RSimplexVariable(ihi,j)+ptry(j)
        RSimplexVariable(ihi,j)=ptry(j)
     ENDDO
  ENDIF

  SimplexExtrapolate=ytry

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
