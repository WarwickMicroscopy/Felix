SUBROUTINE NDimensionalDownhillSimplex(p,y,mp,np,ndim,ftol,funk,iter)

  USE MyNumbers
    
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       iter,mp,ndim,np,NMAX,ITMAX
  REAL(KIND) :: &
       ftol,p(mp,np),y(mp),funk
  PARAMETER (NMAX=20,ITMAX=5000)
  EXTERNAL funk

  INTEGER(IKIND) :: &
       i,ihi,ilo,inhi,j,m,n
  REAL(RKIND) :: &
       rtol,sum,swap,ysave,ytry,psum(NMAX),amotry
  
  iter = 0
  
1 DO 12 n = 1,ndim
     sum = 0
     DO 11 m=1,ndim+1
        sum=sum+p(m,n)
     ENDDO 11
     psum(n) = sum
  ENDDO 12
  
2 ilo = 1
  IF (y(1).GT.y(2)) THEN
     ihi=1
     inhi=2
  ELSE
     ihi=2
     inhi=1
  END IF
  DO 13 i=1,ndim+1
     IF(y(i).LE.y(ilo)) ilo=i
     IF(y(i).GT.y(ihi)) THEN
        inhi=ihi
        ihi=i
     ELSE IF(y(i).GT.y(inhi)) THEN
        IF(i.NE.ihi) inhi=i
     ENDIF
  ENDDO 13
  
  rtol=2.*ABS(y(ihi)-y(ilo))/(ABS(y(ihi))+ABS(y(ilo)))
  IF(rtol.LT.ftol) THEN
     swap=y(1)
     y(1)=y(ilo)
     y(ilo)=swap
     DO 14 n=1,ndim
        swap=p(1,n)
        p(1,n)=p(ilo,n)
        p(ilo,n)=swap
     END DO 14
     RETURN
  END IF

  IF (iter.GE.ITMAX) pause 'ITMAX exceeded in NDimensionalDownhillSimplex'

  iter=iter+2

  ytry = SimplexExtrapolate(p,y,psum,mp,np,ndim,funk,ihi,-1.0D0)
  IF (ytry.LE.y(ilo)) THEN
     ytry = SimplexExtrapolate(p,y,psum,mp,np,ndim,funk,ihi,2.0D0)
  ELSEIF
     ysave=y(ihi)
     ytry=SimplexExtrapolate(p,y,psum,mp,np,ndim,funk,ihi,0.5D0)
     IF(ytry.GE.ysave) THEN
        DO 16 i=1,ndim+1
           IF(i.NE.ilo) THEN
              DO 15 j=1,ndim
                 psum(j)=0.5*(p(i,j)+p(ilo,j))
                 p(i,j)=psum(j)
              ENDDO 15
              y(i)=funk(psum)
           ENDIF
        ENDDO 16
        iter=iter+ndim
        GOTO 1
     ENDIF
  ELSE
     iter=iter-1
  ENDIF
  GOTO 2
END SUBROUTINE NDimensionalDownhillSimplex

REAL(RKIND) FUNCTION SimplexExtrapolate(p,y,psum,mp,np,ndim,funk,ihi,fac)

  USE MyNumbers
    
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: &
       ihi,mp,ndim,np,NMAX
  REAL(RKIND) :: &
       SimplexExtrapolate,fac,p(mp,np),psum(np),y(mp),funk 
  PARAMETER(NMAX=20)
  EXTERNAL funk

  INTEGER(IKIND) :: &
       j
  REAL(RKIND) :: &
       fac1,fac2,ytry,ptry(NMAX)
  fac1=(1.0-fac)/ndim
  fac2=fac1-fac
  DO 11 j=1,ndim
     ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
  ENDDO 11

  ytry=funk(ptry)
  
  IF (ytry.LT.y(ihi)) THEN
     y(ihi)=ytry
     DO 12 j=1,ndim
        psum(j)=psum(j)-p(ihi,j)+ptry(j)
        p(ihi,j)=ptry(j)
     ENDDO 12
  ENDIF
  SimplexExtrapolate=ytry
  RETURN
END FUNCTION SimplexExtrapolate
