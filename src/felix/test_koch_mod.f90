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

MODULE test_koch_mod

! equation 19
! A useful expansion of the exponential of the sum of two non-commuting matrices, one of which is diagonal
! By - Christoph T Koch and John C H Spence
! J. Phys. A: Math. Gen. 36 (2003) 803–816

! Dcoeff from equation 18

! bprime(0:u) are the unique B(l(0:q)), and there is a degeneracy d(k) for each bprime
! (bprime(k) unique -> d(k) = 0)
! bprime(k) used instead of bprime(l(k))

! for testing, for example gfortran -o test_koch_mod test.f90 test_koch_mod.f90 -fbounds-check

USE MyNumbers

IMPLICIT NONE

PRIVATE
PUBLIC :: CalculateElementS, GetCombinations

! used for profiling and timing
INTEGER(8) ::  time,time2,time3,time4,time5,time6,time7,time8,time9,time10,time11,time12,&
                time13,time14,&
                timesum,timesum2,timesum3,timesum4,timesum5,timesum6

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

    N = SIZE(A,1)
    ! check A, B square and same size
    IF(.NOT.(SIZE(A,1).EQ.SIZE(A,2).AND.SIZE(Bmatrix,1).EQ.SIZE(Bmatrix,2)&
          .AND.SIZE(A,1).EQ.SIZE(Bmatrix,1))) RETURN
    ! check diagonal elements of A equal zero
    IF(ANY([ (ABS(A(ind,ind)).GT.1e-40,ind=1,N) ])) RETURN
    ! todo - check A matrix offdiagonal unique
    ! todo - check B diagonals unique

    N = SIZE(A,1)

    WRITE(*,'(a)')' -------------------------------------------------------------'  

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
      timesum=0;timesum2=0;timesum3=0;timesum4=0;timesum5=0;timesum6=0;
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
!        DO WHILE (l(q-1).LT.N+1)

!          ! multiply a_n,l1 , a_l1,l2 , a_l2,l3, ... , a_lq-1,lq , a_lq,m
!          ! A(l(0),l(1)) ... A(l(q-1),l(q))
!          sumproduct = CMPLX(1,0,CKIND)
!          ! assume the diagonals of A are 0
!          IF(ALL([ (l(ind).NE.l(ind+1),ind=0,q-1) ])) THEN
!            DO ind = 0,q-1
!                sumproduct = sumproduct * A(l(ind),l(ind+1))
!            END DO

!            CALL SYSTEM_CLOCK(time7)
!            CALL CalculateCcoeff ( B, lambda, l, q, Ccoeff ) 
!            CALL SYSTEM_CLOCK(time8)
!            timesum3 = timesum3 + time8 - time7
!            S = S + sumproduct * Ccoeff

!          END IF

!          ! iterate l_1 by 1 and check through each summation,
!          ! if a summation has reached N, reset and increment summation above
!          l(1) = l(1) + 1
!          DO ind = 2,q-1
!            IF(l(ind-1).EQ.N+1) THEN
!              !IF(q.GE.5.AND.ind.EQ.(q-1)) WRITE(*,*) 'S',S
!              l(ind-1) = 1
!              l(ind) = l(ind) + 1
!            ELSE
!              EXIT
!            END IF
!          END DO

!        END DO
  ! ---------------------------------------------------------------------------------
        CALL GetCombinations ( N, q, nnd, mnd, A, B, lambda, S )

      END IF
      WRITE(*,'(a,(F11.5,SP,F11.5,"i"))') 'S = ',S
      CALL system_clock(time2)

      IF(q.GE.5) THEN
        WRITE(*,'(a,i0,a,i0)') ' total Time elapsed = ',time2 - time,' for q = ',q
        WRITE(*,'(a,i0)') 'approx. time elapsed getting combintions and calculating = ',timesum6
        !WRITE(*,'(a,i0)') 'approx. time elapsed doing Ccoeff = ',timesum3
        WRITE(*,'(a,i0)') 'approx. time elapsed doing inside Ccoeff = ',timesum5
        WRITE(*,'(a,i0)') 'approx. time elapsed doing r sum in Ccoeff = ',timesum4
        WRITE(*,'(a,i0)') 'approx. time elapsed doing GetUniqueSubset = ',timesum2
        WRITE(*,'(a,i0)') 'approx. time elapsed doing Dcoeff = ',timesum
      END IF  
      !WRITE(*,'(a)')' -------------------------------------------------------------'  
    END DO
    !WRITE(*,'(a)')' -------------------------------------------------------------'  

  END SUBROUTINE


  ! This if for only calculating Ccoeff once and reusing for repeats, which should
  ! significantly reduce calculation times.
  ! print combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
  ! converted into fortran from python https://docs.python.org/2/library/itertools.html
  SUBROUTINE GetCombinations ( bigN, q, nnd, mnd, A, B, lambda, S )

    INTEGER(4),INTENT(IN) :: bigN, q, nnd, mnd
    COMPLEX(CKIND),INTENT(IN) :: A(:,:), B(:), lambda ! constant inputs !?? could be global
    COMPLEX(CKIND),INTENT(INOUT) :: S
    INTEGER(4) :: r
    INTEGER(4) :: indices(0:q-1-1),i,n
    !INTEGER(4) :: pool(0:2)
    INTEGER(4) :: pool(0:bigN-1), l(0:q)

    CALL SYSTEM_CLOCK(time13)

    pool = [( i,i=1,bigN )]
    l(0) = nnd
    l(q) = mnd

    r = q-1

    !pool=[1,2,3]
    n = SIZE(pool)
    indices = 0
    l(1:q-1) = pool(indices(:))
    CALL UseUniqueList( l, q, mnd, nnd, A, B, lambda, S )

    DO WHILE (.TRUE.)
      DO i = r-1,0,-1
        IF(indices(i).NE.n-1) THEN
          EXIT
        ELSEIF(i.EQ.0) THEN
          CALL SYSTEM_CLOCK(time14)
          timesum6 = timesum6 + time14 - time13
          RETURN
        END IF
      END DO
      
      !WRITE(*,*) i
      indices(i:r-1) = indices(i) + 1
      l(1:q-1) = pool(indices(:))
      CALL UseUniqueList( l, q, mnd, nnd, A, B, lambda, S )
      !WRITE(*,*) 'S', S
    END DO

  END SUBROUTINE 

  ! https://rosettacode.org/wiki/Permutations#Fortran
  SUBROUTINE UseUniqueList( l, q, mnd, nnd, A, B, lambda, S )

    implicit none

    INTEGER(4),INTENT(IN) :: l(0:) ! this acts as l
    INTEGER(4),INTENT(IN) :: q, mnd, nnd
    COMPLEX(CKIND),INTENT(IN) :: A(:,:), B(:), lambda ! constant inputs !?? could be global
    COMPLEX(CKIND),INTENT(INOUT) :: S
    COMPLEX(CKIND) :: Ccoeff, sumproduct
 
    integer :: n, i, INoUniquePerms, IUniquePerms( SIZE(l)**SIZE(l), SIZE(l) )
    integer :: littleA(SIZE(l))

!    WRITE(*,*) '--------------------------------'
!    WRITE(*,*) l
!    WRITE(*,*) '--------------------------------'

    CALL SYSTEM_CLOCK(time7)
    CALL CalculateCcoeff ( B, lambda, l, q, Ccoeff ) 
    CALL SYSTEM_CLOCK(time8)
    timesum3 = timesum3 + time8 - time7

    !read *, n
    littleA = l
    n = size(littleA)
    IUniquePerms = -1
    INoUniquePerms = 0
    call perm(1)

!    WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!    DO i = 1,INoUniquePerms
!      WRITE(*,*) IUniquePerms(i,:)
!    END DO
!    WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

    CONTAINS

      recursive subroutine perm(i)
          integer :: i, j, k, t
          if (i == n .AND. ALL([(ANY(littleA(:).NE.IUniquePerms(k,:)),k=1,INoUniquePerms)])) then

            INoUniquePerms = INoUniquePerms + 1
            IUniquePerms(INoUniquePerms,:) = littleA(:)
            !WRITE(*,*) littleA
              
            ! multiply a_n,l1 , a_l1,l2 , a_l2,l3, ... , a_lq-1,lq , a_lq,m
            ! A(l(0),l(1)) ... A(l(q-1),l(q))
            sumproduct = CMPLX(1,0,CKIND)
            ! assume the diagonals of A are 0
            IF(ALL([ (l(k).NE.l(k+1),k=0,q-1) ])) THEN
              DO k = 0,q-1
                  sumproduct = sumproduct * A(l(k),l(k+1))
              END DO

              S = S + sumproduct * Ccoeff
            END IF

          else
              do j = i, n
                  t = littleA(i)
                  littleA(i) = littleA(j)
                  littleA(j) = t
                  call perm(i + 1)
                  t = littleA(i)
                  littleA(i) = littleA(j)
                  littleA(j) = t
              end do
          end if
      end subroutine

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
    CALL GetUniqueSubset ( B, l, q, bprime, d, u )
    CALL SYSTEM_CLOCK(time4)
    timesum2 = timesum2 + time4 - time3

    Ccoeff = CMPLX(0,0,CKIND)
    DO k = 0,u

      ! this remains constant independent of jprime value
      CLambdaBprimeValue = lambda * bprime(k)
      CExpLambdaBprimeValue = EXP( CLambdaBprimeValue )

      ! first recursive value
      CLambdaJprimeFactorial = 1/lambda
      DO jprime = 0,d(k)

        !IF(q.EQ.3.AND.ALL(l(0:3).EQ.[1,1,1,1])) WRITE(*,*) '1,1,1,1 q',q,'jprime',jprime,'k',k,'d(k)',d(k),'brpime',bprime(0:u)
        CALL SYSTEM_CLOCK(time3)

        ! summation from r = 0 to r = q - jprime - 1
        CALL SYSTEM_CLOCK(time9)
        IF(q - jprime - 1.GE.0) THEN
          ! r = 0 value
          CPreviousValue = 1
          ! start product with r = 0 contribution
          CSummationR  = CMPLX(1,0,CKIND)
        ELSE ! This shouldn't happen anymore with the zero diagonals of A
          !WRITE(*,*) 'q - jprime - 1.LT.0'
          !WRITE(*,*) 'q',q,'jprime',jprime,'k',k,'d(k)',d(k),'brpime',bprime(0:u)
          !WRITE(*,*) 'l',l(0:q)
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

        CALL SYSTEM_CLOCK(time5)
        CALL CalculateDcoeff ( B, l, q, bprime, k, d, u, jprime, Dcoeff )
        CALL SYSTEM_CLOCK(time6)
        timesum = timesum + time6 - time5

        CLambdaJprimeFactorial = CLambdaJprimeFactorial * lambda / CMPLX(MAX(1,jprime),0,CKIND)
        Ccoeff = Ccoeff + ( -CSummationR + CExpLambdaBprimeValue ) * CLambdaJprimeFactorial * Dcoeff

        ! test values for single l list
        !IF(q.EQ.3) THEN
        !  IF(ALL(l(0:q).EQ.[1,3,2,2])) THEN
        !    WRITE(*,*) 'jprime = ',jprime, ' k = ',k, 'factorial & lambda', &
        !    lambda**real(jprime,KIND(8)) / real(factorial(jprime),kind(8))
        !    WRITE(*,*) 'l(0:q) ', l(0:q)        
        !    WRITE(*,*) 'Dcoeff', Dcoeff, 'k = ', k, ' jprime = ',jprime
        !    WRITE(*,'(a)') "--------------------------------------------------"
        !  END IF
        !END IF

        CALL SYSTEM_CLOCK(time4)

      END DO

    END DO

    CALL SYSTEM_CLOCK(time12)
    timesum5 = timesum5 + time12 - time11
    
  END SUBROUTINE CalculateCcoeff



  SUBROUTINE GetUniqueSubset ( B, l, q, bprime, d, u )

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

        rsum = CMPLX(0,0,CKIND)
        IF((d(k) - jprime).EQ.1) THEN
          DO ind = 1,INoOfPermittedr

            rsum(1) = rsum(1) + CMPLX(1,0,CKIND)/(bprime(k) - B(l(r(1))))
            r(1) = permitted_r_values( ind+1 )
          END DO
        ELSE ! d(k) - jprime .GE. 2
  ! --------------------------------------------------------------------------------- 
          NotFinished = .TRUE.        
          DO WHILE (NotFinished)

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

        Dcoeff = Dcoeff * rsum(1)

      END IF

    END IF
  ! ---------------------------------------------------------------------------------        
  ! ---------------------------------------------------------------------------------        
  ! ---------------------------------------------------------------------------------         

  END SUBROUTINE



  PURE RECURSIVE FUNCTION factorial ( n ) result ( f )
    INTEGER(4),INTENT(IN) :: n
    INTEGER(4) :: f
    IF(n.EQ.1.OR.n.EQ.0) THEN
      f = 1
    ELSE
      f = n * factorial( n - 1_4 )
    END IF
  END FUNCTION

    

END MODULE
