

!>>> added simple hard-coded timing & profiling
  !certain Ccoeff loops take particularly long
  !GetUniqueSubset seems small, with odd long but still small realtive to whole Ccoeff

! These timings are used alongside a duplicate module to compare changes

!>>> added new method for r summation in Ccoeff
  ! multiples previous summation element instead of recalculating
  ! unfortunately negligable performance boost

!>>> peformance improvement, only call Ccoeff when zero diagonal element not included
  ! seems to be converging correctly and can view it approach

!>>> made new GetUniqueSubset2 assuming B values unique

!>>> Made Ccoeff more recursive

MODULE test_koch_mod

! equation 19
! A useful expansion of the exponential of the sum of two non-commuting matrices, one of which is diagonal
! By - Christoph T Koch and John C H Spence
! J. Phys. A: Math. Gen. 36 (2003) 803–816

! Dcoeff from equation 18

! bprime(0:u) are the unique B(l(0:q)), and there is a degeneracy d(k) for each bprime
! (bprime(k) unique -> d(k) = 0)
! bprime(k) used instead of bprime(l(k))

USE MyNumbers

IMPLICIT NONE

PRIVATE
PUBLIC :: CalculateElementS, GetCombinations

! used for profiling and timing
INTEGER(16) ::  time,time2,time3,time4,time5,time6,time7,time8,time9,time10,time11,time12,&
                timesum,timesum2,timesum3,timesum4,timesum5

CONTAINS

  ! S_n,m = e ** ( lambda ( A + B ) )
  ! where A, B matrices, B diagonal, lambda a scalar

  SUBROUTINE CalculateElementS ( lambda, A, Bmatrix, nnd, mnd, max_q, S )

    ! INPUTS
    INTEGER(4),INTENT(IN) :: nnd, mnd, max_q
    COMPLEX(CKIND),INTENT(IN) :: lambda, Bmatrix(:,:), A(:,:)
    ! OUTPUTS
    COMPLEX(CKIND),INTENT(OUT) :: S
    ! local variables
    COMPLEX(CKIND) :: Ccoeff, B(SIZE(Bmatrix,1)),&
                  sumproduct,& ! to store an iterative product
                  element
    INTEGER(4) :: ind,jnd, & ! a generic looping index
                  q, N, l(0:max_q)
    CHARACTER(100) :: formatting ! for writing terminal output

    ! check A, B square and same size
    IF(.NOT.(SIZE(A,1).EQ.SIZE(A,2).AND.SIZE(Bmatrix,1).EQ.SIZE(Bmatrix,2)&
          .AND.SIZE(A,1).EQ.SIZE(Bmatrix,1))) RETURN
    ! check diagonal elements of A equal zero
    IF(ANY([ (ABS(A(ind,ind)).GT.1e-40,ind=1,N) ])) RETURN
    ! check A matrix unique
    ! check B diagonals unique

    N = SIZE(A,1)

    !WRITE(*,'(a)')' -------------------------------------------------------------'  
!    !PRINT A, B, lambda
!    WRITE(*,'(a,i3,a,i3,a,(F12.7,SP,F12.7,"i"))') " n = ", nnd, "     m = ", mnd, "    lamda = ", lambda 
!    DO ind = 1,N
!!      IF(ind.EQ.1) THEN
!!        WRITE(formatting,'(a,i0,a,i0,a)') "(a,", N, "2(F12.7,1x),4x,a,", N, "2(F12.7,1x),4x,a)"
!!        WRITE(*,formatting) " A = [ ", A(:,ind), "]    B = [ ", Bmatrix(:,ind), "]"
!!      ELSE 
!      WRITE(formatting,'(a,i0,a,i0,a)') '(a,', N-1, '(F12.7,SP,F12.7,"i,"),(F12.7,SP,F12.7,"i"),a)'
!      WRITE(*,formatting) "     { ", A(:,ind), "},"
!    END DO
!    WRITE(*,*) '-------------------------------------------------------------'

!    WRITE(*,'(a,i3,a,i3,a,(F8.3,SP,F8.3,"i"))') " n = ", nnd, "     m = ", mnd, "    lamda = ", lambda 
!    DO ind = 1,10
!!      IF(ind.EQ.1) THEN
!!        WRITE(formatting,'(a,i0,a,i0,a)') "(a,", N, "2(F6.3,1x),4x,a,", N, "2(F6.3,1x),4x,a)"
!!        WRITE(*,formatting) " A = [ ", A(:,ind), "]    B = [ ", Bmatrix(:,ind), "]"
!!      ELSE 
!      WRITE(formatting,'(a,i0,a,i0,a)') "(a,", 5, '(F6.3,SP,F6.3,"i, "),a,', 5, '(F6.3,SP,F6.3,"i, "),a)'
!      WRITE(*,formatting) "     [ ", A(1:5,ind), "]        [ ", Bmatrix(1:5,ind), "]"
!    END DO
!    WRITE(*,*) '-------------------------------------------------------------'

    ! use single dim N array for B instead of N x N diagonal matrix
    DO ind = 1,N
      B(ind) = Bmatrix(ind,ind)
    END DO

    S = CMPLX(0,0,CKIND)

    ! + e ** ( lambda b_n ) * delta_n,m 
    IF(nnd.EQ.mnd) S = S + EXP(lambda*B(nnd))

    ! ------------------------------------
    ! summation over q
    ! ------------------------------------

    DO q = 1,max_q
      CALL system_clock(time)
      timesum=0;timesum2=0;timesum3=0;timesum4=0;
      !WRITE(*,'(a,i0)') 'q = ',q

      IF(q.EQ.1) THEN
        ! simply one term in summation a_n,m * Ccoeff
        l(0) = nnd
        l(1) = mnd
        IF(nnd.NE.mnd) THEN 
          CALL CalculateCcoeff ( B, lambda, l, q, Ccoeff )
          S = S + A(nnd,mnd) * Ccoeff
        END IF
      ELSEIF(q.EQ.2) THEN
        ! summation becomes a_n,l1 * a_l1,m * Ccoeff with l1 = 0 to l1 = N
        l(0) = nnd
        l(2) = mnd
        DO ind = 1,N
          l(1) = ind
          IF(ALL([ (l(jnd).NE.l(jnd+1),jnd=0,q-1) ])) THEN
            CALL CalculateCcoeff ( B, lambda, l, q, Ccoeff )
            S = S + A(nnd,l(1)) * A(l(1),mnd) * Ccoeff
          END IF
        END DO
      ELSE
        ! ------------------------------------
        ! q.GE.3
        ! ------------------------------------ 

        
        l = 1
        l(0) = nnd
        l(q) = mnd

        ! assume that can simply sum all of the final products together
        ! incrementing low indices l(1) first and l(q-1) last
  ! ---------------------------------------------------------------------------------
        DO WHILE (l(q-1).LT.N+1)

          !WRITE(*,*) l(0:q)
          !WRITE(*,'(a)') "--------------------------------------------------"
          !WRITE(formatting,'(a,i0,a)') "(a,", q-1, "(i4,1x))"
          !WRITE(*,formatting) 'l(1:q+1) = ', l(1:q+1)

          ! multiply a_n,l1 , a_l1,l2 , a_l2,l3, ... , a_lq-1,lq , a_lq,m
          ! A(l(0),l(1)) ... A(l(q-1),l(q))
          sumproduct = CMPLX(1,0,CKIND)
          ! assume the diagonals of A are 0
          IF(ALL([ (l(ind).NE.l(ind+1),ind=0,q-1) ])) THEN
            DO ind = 0,q-1
                sumproduct = sumproduct * A(l(ind),l(ind+1))
            END DO

            !WRITE(*,'(*(i2,1x))') l(0:q)
            !IF(ALL([ (l(ind).NE.l(ind+1),ind=0,q-1) ]).AND.ABS(sumproduct).LT.1e-40) WRITE(*,*) sumproduct, l(0:q), (/ ( A(l(jnd),l(jnd+1)),jnd=0,q-1 ) /)

            CALL SYSTEM_CLOCK(time7)
            CALL CalculateCcoeff ( B, lambda, l, q, Ccoeff ) 
            CALL SYSTEM_CLOCK(time8)
            timesum3 = timesum3 + time8 - time7
            S = S + sumproduct * Ccoeff
            !WRITE(*,*) sumproduct * Ccoeff
          END IF


          !WRITE(*,*) 'sumproduct'   
          !WRITE(*,'(F6.2)') Ccoeff
          !IF(ABS(Ccoeff).GT.1) WRITE(*,*) 'Ccoeff', Ccoeff, ' l(0:q) = ',l(0:q)
          !IF(q.EQ.7) WRITE(*,*) 'Ccoeff', Ccoeff,'||| q', q, ' l(0:q) = ',l(0:q)
          !WRITE(*,*) 'S',S

          ! iterate l_1 by 1 and check through each summation,
          ! if a summation has reached N, reset and increment summation above
          l(1) = l(1) + 1
          DO ind = 2,q-1
            IF(l(ind-1).EQ.N+1) THEN
              !IF(q.GE.5.AND.ind.EQ.(q-1)) WRITE(*,*) 'S',S
              !IF(q.GE.5.AND.ind.EQ.(q-1).AND.ALL([ (l(ind).NE.l(ind+1),ind=0,q-1) ])) WRITE(*,*) 'S',S, 'sumproduct * Ccoeff', sumproduct * Ccoeff 
              l(ind-1) = 1
              l(ind) = l(ind) + 1
            ELSE
              EXIT
            END IF
          END DO

        END DO
  ! ---------------------------------------------------------------------------------
      END IF
      !WRITE(*,'(a,(F11.5,SP,F11.5,"i"))') 'S = ',S
      CALL system_clock(time2)

      IF(q.GE.5) THEN
        !WRITE(*,*) 'total Time elapsed = ',time2 - time
        !WRITE(*,*) 'approx. time elapsed doing Ccoeff = ',timesum3
        !WRITE(*,*) 'approx. time elapsed doing inside Ccoeff = ',timesum5
        !WRITE(*,*) 'approx. time elapsed doing GetUniqueSubset = ',timesum2
        !WRITE(*,*) 'approx. time elapsed doing Dcoeff = ',timesum
        !WRITE(*,*) 'approx. time elapsed doing r sum in Ccoeff = ',timesum4
      END IF  
      !WRITE(*,'(a)')' -------------------------------------------------------------'  
    END DO
    !WRITE(*,'(a)')' -------------------------------------------------------------'  

  END SUBROUTINE


  ! print combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
  ! converted into fortran from python https://docs.python.org/2/library/itertools.html
  SUBROUTINE GetCombinations ()

    INTEGER(4),PARAMETER :: bigN = 4, q = 3
    INTEGER(4),PARAMETER :: r = q+1
    INTEGER(4) :: indices(0:r-1),i,n
    !INTEGER(4) :: pool(0:2)
    INTEGER(4) :: pool(0:bigN-1)
    pool = [( i,i=1,bigN )]

    !pool=[1,2,3]
    n = SIZE(pool)
    indices = 0
    CALL UseUniqueList( pool(indices(:)))

    DO WHILE (.TRUE.)
      DO i = r-1,0,-1
        IF(indices(i).NE.n-1) THEN
          EXIT
        ELSEIF(i.EQ.0) THEN
          RETURN
        END IF
      END DO
      
      !WRITE(*,*) i
      indices(i:r-1) = indices(i) + 1
      CALL UseUniqueList( pool(indices(:)) )
    END DO 

  END SUBROUTINE 


  ! https://rosettacode.org/wiki/Permutations#Fortran
  SUBROUTINE UseUniqueList( poolsubset )

    implicit none

    INTEGER(4),INTENT(IN) :: poolsubset(:)
 
    integer :: n, i
    integer :: a(SIZE(poolsubset))

    WRITE(*,*) poolsubset
    WRITE(*,*) '--------------------------------'

!    !read *, n
!    a = poolsubset
!    n = size(a)
!    call perm(1)

!    CONTAINS

!      recursive subroutine perm(i)
!          integer :: i, j, t
!          if (i == n) then
!              print *, a
!          else
!              do j = i, n
!                  t = a(i)
!                  a(i) = a(j)
!                  a(j) = t
!                  call perm(i + 1)
!                  t = a(i)
!                  a(i) = a(j)
!                  a(j) = t
!              end do
!          end if
!      end subroutine

  END SUBROUTINE


  SUBROUTINE CalculateCcoeff ( B, lambda, l, q, Ccoeff )

    COMPLEX(CKIND),INTENT(OUT) :: Ccoeff
    COMPLEX(CKIND),INTENT(IN) :: B(:), lambda
    INTEGER(4),INTENT(IN) :: l(0:), q
    COMPLEX(CKIND) :: bprime(0:q), Dcoeff, CIterationProduct, CSummationR, CPreviousValue, &
                      CLambdaBprimeValue, CLambdaJprimeFactorial, CExpLambdaBprimeValue
    INTEGER(4) :: d(0:q), u, jprime, k, r, ind
    CHARACTER(50) :: formatting ! for writing terminal output

    CALL SYSTEM_CLOCK(time11)
    CALL SYSTEM_CLOCK(time3)
    CALL GetUniqueSubset2 ( B, l, q, bprime, d, u )
    CALL SYSTEM_CLOCK(time4)
    timesum2 = timesum2 + time4 - time3

    !WRITE(*,*) 'Time elapsed GetUniqueSubset() = ',time4 - time3

    !WRITE(*,'(a)') "--------------------------------------------------"
    !WRITE(formatting,'(a,i0,a,i0,a)') "(a,", q+1, "(i4,1x)'|||'", q+1, "(F5.2,1x))"
    !WRITE(*,formatting) 'l(0:q), ( B(l(ind)), ind=1,q+1 ) ', l(0:q), ( B(l(ind)), ind=0,q )
    !WRITE(formatting,'(a,i0,a,i0,a)') '(a,', u+1, '(F5.2,1x)"|||",', u+1, '(i4,1x))'
    !WRITE(*,formatting) 'bprime(0:u), d(0:u)', bprime(0:u), d(0:u)
    !IF(q.EQ.3) THEN
    !  IF(ALL(l(0:q).EQ.[1,3,2,2])) THEN
    !    WRITE(formatting,'(a,i0,a,i0,a)') "(a,", q+1, "(i4,1x)'|||'", q+1, "(F5.2,1x))"
    !    WRITE(*,formatting) 'l(0:q), ( B(l(ind)), ind=1,q+1 ) ', l(0:q), ( B(l(ind)), ind=0,q )
    !    WRITE(formatting,'(a,i0,a,i0,a)') '(a,', u+1, '(F5.2,1x)"|||",', u+1, '(i4,1x))'
    !    WRITE(*,formatting) 'bprime(0:u), d(0:u)', bprime(0:u), d(0:u)
    !    WRITE(*,'(a)') '---------- SUM over k and jprime ----------------------------------- ' 
    !  END IF
    !END IF

    Ccoeff = CMPLX(0,0,CKIND)
    DO k = 0,u

      ! this remains constant independent of jprime value
      CLambdaBprimeValue = lambda * bprime(k)
      CExpLambdaBprimeValue = EXP( CLambdaBprimeValue )

      ! first recursive value
      CLambdaJprimeFactorial = 1/lambda
      DO jprime = 0,d(k)

        IF(q.EQ.3.AND.ALL(l(0:3).EQ.[1,1,1,1])) WRITE(*,*) '1,1,1,1 q',q,'jprime',jprime,'k',k,'d(k)',d(k),'brpime',bprime(0:u)
        CALL SYSTEM_CLOCK(time3)

        ! summation from r = 0 to r = q - jprime - 1
        CALL SYSTEM_CLOCK(time9)
        IF(q - jprime - 1.GE.0) THEN
          ! r = 0 value
          CPreviousValue = 1
          ! start product with r = 0 contribution
          CSummationR  = CMPLX(1,0,CKIND)
        ELSE
          WRITE(*,*) 'q - jprime - 1.LT.0'
          WRITE(*,*) 'q',q,'jprime',jprime,'k',k,'d(k)',d(k),'brpime',bprime(0:u)
          WRITE(*,*) 'l',l(0:q)
          CSummationR = CMPLX(0,0,CKIND)
        END IF
        DO r = 1, q - jprime - 1         
          ! each element of summation is simply a multiplication of previous
          CPreviousValue = CPreviousValue * CLambdaBprimeValue / r
          !CIterationProduct = CIterationProduct + CPreviousValue
          CSummationR = CSummationR + CPreviousValue
        END DO
        CALL SYSTEM_CLOCK(time10)
        timesum4 = timesum4 + time10 - time9
        !WRITE(*,*) 'Time elapsed doing r sum in Ccoeff = ',time10 - time9, 'q - jprime - 1', q - jprime - 1, 'CLambdaBprimeValue', CLambdaBprimeValue, 'CIterationProduct',CIterationProduct

        CALL SYSTEM_CLOCK(time5)
        CALL CalculateDcoeff ( B, l, q, bprime, k, d, u, jprime, Dcoeff )
        CALL SYSTEM_CLOCK(time6)
        timesum = timesum + time6 - time5
        !WRITE(*,*) 'Time elapsed doing Dcoeff = ',time6 - time5

        CLambdaJprimeFactorial = CLambdaJprimeFactorial * lambda / CMPLX(MAX(1,jprime),0,CKIND)
        Ccoeff = Ccoeff + ( -CSummationR + CExpLambdaBprimeValue ) * CLambdaJprimeFactorial * Dcoeff

        !IF(q.EQ.3) THEN
        !  IF(ALL(l(0:q).EQ.[1,3,2,2])) THEN
        !    WRITE(*,*) 'jprime = ',jprime, ' k = ',k, 'factorial & lambda', &
        !    lambda**real(jprime,KIND(8)) / real(factorial(jprime),kind(8))
        !    WRITE(*,*) 'l(0:q) ', l(0:q)        
        !    WRITE(*,*) 'Dcoeff', Dcoeff, 'k = ', k, ' jprime = ',jprime
        !    WRITE(*,'(a)') "--------------------------------------------------"
        !  END IF
        !END IF
        !WRITE(*,*) 'Dcoeff', Dcoeff, 'k = ', k, ' jprime = ',jprime
        !WRITE(*,*) ' exp( lambda * bprime(k) )', exp( lambda * bprime(k) )
        !WRITE(*,*) 'CIterationProduct', CIterationProduct
        !WRITE(*,*) 'k, jprime', k, jprime
        !WRITE(*,*) 'q - jprime - 1', q - jprime - 1
        !IF(ABS(Dcoeff).GT.1) WRITE(*,*) 'Dcoeff ', Dcoeff
        !IF(ABS(CIterationProduct).GT.1) WRITE(*,*) 'CIterationProduct ', CIterationProduct

        CALL SYSTEM_CLOCK(time4)
        !WRITE(*,*) 'Time elapsed doing a Ccoeff loop = ',time4 - time3
        !IF((time4 - time3).GT.100000.OR..TRUE.) WRITE(*,*) 'k', k, 'jprime', jprime
      END DO

    END DO

    !WRITE(*,*) 'Ccoeff', Ccoeff
    !IF(q.EQ.7) THEN
    !  IF(ALL(l(0:q).EQ.[1,2,2,2,2,1,1,1])) THEN
    !    WRITE(*,*) 'Ccoeff', Ccoeff
    !  END IF
    !END IF
    CALL SYSTEM_CLOCK(time12)
    timesum5 = timesum5 + time12 - time11
    
  END SUBROUTINE CalculateCcoeff



  SUBROUTINE GetUniqueSubset ( B, l, q, bprime, d, u )

    COMPLEX(CKIND),INTENT(IN) :: B(:)
    COMPLEX(CKIND),INTENT(INOUT) :: bprime(0:)
    INTEGER(4),INTENT(INOUT) :: d(0:)
    INTEGER(4),INTENT(IN) :: l(0:), q
    INTEGER(4) :: ind, k,&
                  u ! number of unique elements minus 1
    LOGICAL :: IsUnique

    u = -1
    d = 0 ! degenearacy, b(l(k)) unique -> d(k) = 0
    DO ind = 0, q
      IsUnique = .TRUE.
      DO k = 0,u
        IF(B(l(ind)).EQ.bprime(k)) THEN
          IsUnique = .FALSE.
          d(k) = d(k) + 1
          EXIT
        END IF
      END DO
      IF(IsUnique) THEN
        bprime(u + 1) = B(l(ind))
        u = u + 1
      END IF         
    END DO

  END SUBROUTINE

  SUBROUTINE GetUniqueSubset2 ( B, l, q, bprime, d, u )

    COMPLEX(CKIND),INTENT(IN) :: B(:)
    COMPLEX(CKIND),INTENT(INOUT) :: bprime(0:)
    INTEGER(4),INTENT(INOUT) :: d(0:)
    INTEGER(4),INTENT(IN) :: l(0:), q
    INTEGER(4) :: ind, jnd, k, &
                  u, & ! number of unique elements minus 1
                  UniqueListValuesL (0:q)
    LOGICAL :: IsUnique

    u = 0
    UniqueListValuesL(0) = l(0)
    d(0) = COUNT( [( l(0).EQ.l(jnd),jnd=0,q )] ) - 1
    DO ind = 1, q
      IF( ALL( [( l(ind).NE.UniqueListValuesL(jnd),jnd=0,u )] ) ) THEN ! unique l value
        u = u + 1
        UniqueListValuesL(u) = l(ind)
        ! degenearacy, b(l(k)) unique -> d(k) = 0
        d(u) = COUNT( [( l(ind).EQ.l(jnd),jnd=0,q )] ) - 1
      END IF
    END DO

    ! assume each b value unique
    bprime(0:u) = B(UniqueListValuesL(0:u))

  END SUBROUTINE



  SUBROUTINE CalculateDcoeff ( B, l, q, bprime, k, d, u, jprime, Dcoeff )

    COMPLEX(CKIND),INTENT(IN) :: B(:), bprime(0:)
    INTEGER(4),INTENT(IN) :: l(0:), q, k, d(0:), u, jprime
    COMPLEX(CKIND),INTENT(OUT) :: Dcoeff
    INTEGER(4) :: ind, INoOfPermittedr, r_index
    INTEGER(4) :: r(d(k) - jprime), permitted_r_values(q+1), r_permitted_referance(d(k) - jprime)
    COMPLEX(CKIND) :: rsum(d(k) - jprime)
    LOGICAL :: productdefined, NotFinished
    CHARACTER(50) :: formatting ! for writing terminal output

    ! -----------------------------------------------------------
    ! sign (-1) ** and capitial pie product from r = 0 to r = u
    ! -----------------------------------------------------------
    
    Dcoeff = CMPLX( (-1) ** ( d(k) - jprime ), 0, CKIND )

    DO r_index = 0, u

      !WRITE(*,*) 'r_index, k', r_index, k

      IF(r_index.NE.k) THEN
        productdefined = .true.
        Dcoeff = Dcoeff * ( bprime(k) - bprime(r_index) ) ** (-CMPLX( d(r_index) + 1, 0, CKIND) )

        !IF(q.EQ.7.AND.ALL(l(0:q).EQ.[1,2,2,2,2,1,1,1]).AND.k.EQ.0.AND.jprime.EQ.0) THEN
        !  WRITE(*,*) 'bprime(k) = ', bprime(k), ' bprime(r_index) = ',bprime(r_index)
        !  WRITE(*,*) 'fraction and d(r) power', ( bprime(k) - bprime(r_index) ) ** (-REAL( d(r_index) + 1, KIND(8) )) 
        !END iF

      END IF
    END DO

    IF(productdefined.EQV..FALSE.) THEN
      Dcoeff = CMPLX(1,0,CKIND)
      RETURN
    END IF
    ! NB assumed that both this and embedded sum should have value = 1 for null case

    ! -----------------------------------------------------------
    ! large loop for r_1, ..., r_(d_k - j') summation & product
    ! -----------------------------------------------------------
  
  ! ---------------------------------------------------------------------------------        
  ! ---------------------------------------------------------------------------------        
  ! ---------------------------------------------------------------------------------        
    IF((d(k) - jprime).GT.0) THEN

      ! NB this always seems to lead to INoOfPermittedr.GT.0

      ! find permitted r values (this list will not contain duplicate r values)
      permitted_r_values = -1
      INoOfPermittedr = 0
      DO ind = 0,q

        !IF(ALL(l(0:q).EQ.[1,3,2,2]).AND.k.EQ.2.AND.jprime.EQ.0) THEN
          !WRITE(*,*) 'B(l(ind))', B(l(ind)), '||| bprime(k)', bprime(k)
        !END iF
        !WRITE(*,*) 'B(l(ind)) , bprime(k)', B(l(ind)), '|||', bprime(k)

        IF(B(l(ind)).NE.bprime(k)) THEN
          INoOfPermittedr = INoOfPermittedr + 1
          permitted_r_values(INoOfPermittedr) = ind
        END IF
      END DO

      !IF(ALL(l(0:q).EQ.[1,3,2,2]).AND.k.EQ.2.AND.jprime.EQ.0) THEN
        !WRITE(*,*) 'permitted_r_values', permitted_r_values
        !WRITE(*,*) 'd(k) - jprime', d(k) - jprime
      !END iF


      IF(INoOfPermittedr.GT.0) THEN

        r = permitted_r_values( 1 )
        r_permitted_referance =  1

        !WRITE(*,*) 'Iteration with INoOfPermittedr.GT.0'
        !WRITE(*,*) 'q, permitted_r_values', q, '|||', permitted_r_values
        !WRITE(*,*) 'r', r(:)
        !WRITE(formatting,'(a,i0,a,i0,a)') "(a,", q+1, "(i4,1x)','", q+1, "(F5.2,1x))"
        !WRITE(*,formatting) 'l(0:q), ( B(l(ind)), ind=0,q ) ', l(0:q), ( B(l(ind)), ind=0,q )
        !WRITE(formatting,'(a,i0,a,i0,a)') '(a,', u+1, '(F5.2,1x)",",', u+1, '(i4,1x))'
        !WRITE(*,formatting) 'bprime(0:u), d(0:u)', bprime(0:u), d(0:u)

        rsum = CMPLX(0,0,CKIND)
        IF((d(k) - jprime).EQ.1) THEN
          DO ind = 1,INoOfPermittedr

            !IF(ALL(l(0:q).EQ.[1,3,2,2]).AND.k.EQ.2.AND.jprime.EQ.0) THEN
              !WRITE(*,*) 'r', r, '||| all permitted ', permitted_r_values(1:INoOfPermittedr), '||| r ref ', r_permitted_referance
            !END iF
            !WRITE(*,*) 'r stuff', r, '|||', permitted_r_values(INoOfPermittedr), '|||', r_permitted_referance

            rsum(1) = rsum(1) + CMPLX(1,0,CKIND)/(bprime(k) - B(l(r(1))))
            r(1) = permitted_r_values( ind+1 )
          END DO
        ELSE ! d(k) - jprime .GE. 2
  ! --------------------------------------------------------------------------------- 
          NotFinished = .TRUE.        
          DO WHILE (NotFinished)

            !IF(ALL(l(0:q).EQ.[1,3,2,2]).AND.k.EQ.2.AND.jprime.EQ.0) THEN
            !  WRITE(*,*) 'r', r, '||| all permitted ', permitted_r_values(1:INoOfPermittedr), '||| r ref ', r_permitted_referance
            !END iF
            !WRITE(*,*) 'NotFinished',NotFinished
            !WRITE(*,*) 'r stuff', r, '|||', permitted_r_values(INoOfPermittedr), '|||', r_permitted_referance

            rsum(d(k) - jprime) = rsum(d(k) - jprime) + CMPLX(1,0,CKIND)/(bprime(k) - B(l(r(d(k) - jprime))))

            ! moving outwards from innermost sum, check if max index of each sum has been reached
            DO ind = d(k) - jprime, 2, -1
              IF(r(ind).EQ.r(ind-1)) THEN
                ! contribute this summation to sum above
                rsum(ind-1) = rsum(ind-1) + CMPLX(1,0,CKIND)/(bprime(k) - B(l(r(ind-1)))) * rsum(ind)          

                IF(ind.EQ.2) THEN
                  IF(r(1).LT.permitted_r_values(INoOfPermittedr)) THEN ! increment r(1)
                    r(1) = permitted_r_values( r_permitted_referance(1) + 1 )
                    r_permitted_referance(1) = r_permitted_referance(1) + 1
                  ELSE
                    NotFinished = .FALSE.
                  END IF
                END IF

                ! set this r to first permitted value
                rsum(ind) = CMPLX(0,0,CKIND)
                r(ind) = permitted_r_values( 1 )
                r_permitted_referance(ind) = 1    
              ELSE ! r(ind).LT.r(ind-1)
                ! iterate sum to next permitted r value
                r(ind) = permitted_r_values( r_permitted_referance(ind) + 1 )
                r_permitted_referance(ind) = r_permitted_referance(ind) + 1      
                EXIT
              END IF
            END DO
          
          END DO
  ! ---------------------------------------------------------------------------------        
        END IF

        !IF(q.EQ.7.AND.ALL(l(0:q).EQ.[1,2,2,2,2,1,1,1]).AND.k.EQ.0.AND.jprime.EQ.0) THEN
        !  WRITE(*,*) 'rsum(1)', rsum(1)
        !END iF        
        !WRITE(*,*) rsum

        Dcoeff = Dcoeff * rsum(1)

      END IF

    END IF
  ! ---------------------------------------------------------------------------------        
  ! ---------------------------------------------------------------------------------        
  ! ---------------------------------------------------------------------------------         

  END SUBROUTINE



  RECURSIVE FUNCTION factorial ( n ) result ( f )
    INTEGER(4) :: f, n
    IF(n.EQ.1.OR.n.EQ.0) THEN
      f = 1
    ELSE
      f = n * factorial( n - 1_4 )
    END IF
  END FUNCTION

    

END MODULE
! to show in the meeting!
! test suite for range of values
! move inputs to function
! COMPLEX numbers or directly intensity
! impliment current debugging procedure
! impliment into felix
! HUGE PERFORMANCE IMPROVEMENTS
! max_q dynamical error handling
! impliment cleverer way to do permitted values using previous calculations, e.g. (q-1 - d(k)), u
! could impliment checking B diagonal, currently ignoring non-diagonal elements

! gfortran -o test_koch_mod test.f90 test_koch_mod.f90 -fbounds-check
! doing clever permutation stuff for non-unique multiplication
! google whether this has been programmed before
! affect of real byte size on performance
! internal functions have access to external scope but can overwrite
! variable scope for clarity?
! pure/elemental procedures
! use pointers to referance B elements?
! B diagonal so store as single dim array
