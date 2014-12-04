SUBROUTINE NDimensionalDownhillSimplex(RSimplexVolume,y,mp,np,ndim,ftol,iter,IErr)

  USE MyNumbers
    
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       iter,mp,ndim,np,NMAX,ITMAX,IErr
  REAL(RKIND) :: &
       ftol,RSimplexVolume(mp,np),y(mp),SimplexFunction,SimplexExtrapolate
  PARAMETER (NMAX=100,ITMAX=5000)
!!$  EXTERNAL SimplexFunction

  INTEGER(IKIND) :: &
       i,ihi,ilo,inhi,j,m,n
  REAL(RKIND) :: &
       rtol,sum,swap,ysave,ytry,psum(NMAX),amotry
  
  IF(my_rank.EQ.0) THEN
     PRINT*,"Beginning Simplex",IErr
  END IF

  iter = 0
  
1 DO n = 1,ndim
     sum = 0
     DO m=1,ndim+1
        sum=sum+RSimplexVolume(m,n)
     ENDDO 
     psum(n) = sum
  ENDDO 
  
2 ilo = 1
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
     ENDIF
  ENDDO 
  
  rtol=2.*ABS(y(ihi)-y(ilo))/(ABS(y(ihi))+ABS(y(ilo)))

  PRINT*,"Current Tolerance",rtol,ftol
  IF(rtol.LT.ftol) THEN
     swap=y(1)
     y(1)=y(ilo)
     y(ilo)=swap
     DO n=1,ndim
        swap=RSimplexVolume(1,n)
        RSimplexVolume(1,n)=RSimplexVolume(ilo,n)
        RSimplexVolume(ilo,n)=swap
     END DO
     RETURN
  END IF

  IF (iter.GE.ITMAX) THEN
     IErr = 1
     RETURN
     PRINT*,"ITMAX exceeded in NDimensionalDownhillSimplex"
  END IF

  iter=iter+2

  ytry = SimplexExtrapolate(RSimplexVolume,y,psum,mp,np,ndim,ihi,-1.0D0,iter,IErr)
  IF (ytry.LE.y(ilo)) THEN
     ytry = SimplexExtrapolate(RSimplexVolume,y,psum,mp,np,ndim,ihi,2.0D0,iter,IErr)
  ELSEIF (ytry.GE.y(inhi)) THEN
     ysave=y(ihi)
     ytry=SimplexExtrapolate(RSimplexVolume,y,psum,mp,np,ndim,ihi,0.5D0,iter,IErr)
     IF(ytry.GE.ysave) THEN
        DO i=1,ndim+1
           IF(i.NE.ilo) THEN
              DO j=1,ndim
                 psum(j)=0.5*(RSimplexVolume(i,j)+RSimplexVolume(ilo,j))
                 RSimplexVolume(i,j)=psum(j)
              ENDDO
              y(i)=SimplexFunction(psum,iter,IErr)
           ENDIF
        ENDDO
        iter=iter+ndim
        GOTO 1
     ENDIF
  ELSE
     iter=iter-1
  ENDIF
  GOTO 2

END SUBROUTINE NDimensionalDownhillSimplex

REAL(RKIND) FUNCTION SimplexExtrapolate(RSimplexVolume,y,psum,mp,np,ndim,ihi,fac,iter,IErr)

  USE MyNumbers
    
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: &
       ihi,mp,ndim,np,NMAX,IErr,iter
  REAL(RKIND) :: &
       SimplexExtrapolate,fac,RSimplexVolume(mp,np),psum(np),y(mp),SimplexFunction 
  PARAMETER(NMAX=100)
!!$  EXTERNAL SimplexFunction

  INTEGER(IKIND) :: &
       j
  REAL(RKIND) :: &
       fac1,fac2,ytry,ptry(ndim)

  fac1=(1.0-fac)/ndim
  fac2=fac1-fac
  DO j=1,ndim
     ptry(j)=psum(j)*fac1-RSimplexVolume(ihi,j)*fac2
  ENDDO

  ytry=SimplexFunction(ptry,iter,IErr)
  
  IF (ytry.LT.y(ihi)) THEN
     y(ihi)=ytry
     DO j=1,ndim
        psum(j)=psum(j)-RSimplexVolume(ihi,j)+ptry(j)
        RSimplexVolume(ihi,j)=ptry(j)
     ENDDO
  ENDIF
  SimplexExtrapolate=ytry

  
  PRINT*,"-----------------------------------------------------"
  PRINT*,"Iteration",iter
  PRINT*,"Configuration",ptry
  PRINT*,"Figure of Merit",ytry
  PRINT*,"-----------------------------------------------------"


  RETURN
END FUNCTION SimplexExtrapolate