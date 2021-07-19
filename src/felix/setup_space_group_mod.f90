!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans, Jo Bithell & Rudolf A Roemer
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: smodules.f90,v 1.63 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>
!! Module-description: 
!!
MODULE setup_space_group_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: PreferredBasis, DetermineAllowedMovements, CountAllowedMovements, ConvertSpaceGroupToNumber

  CONTAINS
  
  
  !>
  !! Procedure-description: 
  !!
  !! Major-Authors: Richard Beanland (2018)
  !!
  SUBROUTINE PreferredBasis(IErr)
  !Called from felixrefine ~line 175
  !Changes the basis so that atomic coordinate refinement is possible
  !This has to be done for cases where there are several symmetrically equivalent
  !coordinates, and the atoms on the different sites move in different directions.
  !The dimension of movement is given in CountAllowedMovements.
  !The direction of movement is given in DetermineAllowedMovements.
  !Here we have to change the atomic coordinate to be the appropriate one
  !for the movements listed in those two subroutines.

    USE MyNumbers
    USE message_mod
    
    ! global inputs
    USE SPARA, ONLY : SSpaceGrp, SWyckoffSymbol
    USE RPARA, ONLY : RBasisAtomPosition

    IMPLICIT NONE

    CHARACTER(1) :: SWyckoff
    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND) ::  ind,jnd,knd,ISpaceGrp,IBasisAtoms

    IErr=0

    CALL ConvertSpaceGroupToNumber(ISpaceGrp,IErr)!Uses global variable SSpaceGrp
    IF(l_alert(IErr,"PreferredBasis","ConvertSpaceGroupToNumber")) RETURN
    
    IBasisAtoms=SIZE(RBasisAtomPosition,1)
    CALL ChangeOrigin(ISpaceGrp,IErr)!#only implemented for #142
    IF(l_alert(IErr,"PreferredBasis","ChangeOrigin")) RETURN
    
    DO ind=1,IBasisAtoms
      SWyckoff=SWyckoffSymbol(ind)
    
      SELECT CASE(ISpaceGrp)
      CASE(1)!P1
        SELECT CASE (SWyckoff)
        CASE('a')!point symmetry 1, coordinate [x,y,z], no reassignment
        
        CASE DEFAULT
          IErr = 1
          IF(l_alert(IErr,"PreferredBasis",&
              "Wyckoff Symbol for space group 1, P1, not recognised")) RETURN
        END SELECT
!!$  CASE(2)
!!$  CASE(3)
!!$  CASE(4)
!!$  CASE(5)
!!$  CASE(6)
!!$  CASE(7)
!!$  CASE(8)
!!$  CASE(9)
!!$  CASE(10)
!!$  CASE(11)
!!$  CASE(12)
!!$  CASE(13)
!!$  CASE(14)
!!$  CASE(15)
      CASE(15)!C2/c
        SELECT CASE (SWyckoff)
        CASE('a')!point symmetry -1, coordinate [0,0,0], no reassignment

        CASE('b')!point symmetry -1, coordinate [0,1/2,0], no reassignment

        CASE('c')!point symmetry -1, coordinate [1/4,1/4,0], no reassignment

        CASE('d')!point symmetry -1, coordinate [1/4,1/4,1/2], no reassignment

        CASE('e')!point symmetry 2, coordinate [0,y,1/4], no reassignment

        CASE('f')!point symmetry 1, coordinate [x,y,z], no reassignment

        CASE DEFAULT
          IErr = 1
          IF(l_alert(IErr,"PreferredBasis",&
              "Wyckoff Symbol for space group 15, C 2/c, not recognised")) RETURN
        END SELECT

!!$  CASE(16)
!!$  CASE(17)
!!$  CASE(18)
!!$  CASE(19)
!!$  CASE(20)
!!$  CASE(21)
!!$  CASE(22)
!!$  CASE(23)
!!$  CASE(24)
!!$  CASE(25)
!!$  CASE(26)
!!$  CASE(27)
!!$  CASE(28)
!!$  CASE(29)
!!$  CASE(30)
!!$  CASE(31)
!!$  CASE(32)
!!$  CASE(33)
!!$  CASE(34)
!!$  CASE(35)
      CASE(36)!C m c 21
        SELECT CASE (SWyckoff)
        CASE('a')!point symmetry m, coordinate [0,y,z], no reassignment
 
        CASE('b')!point symmetry 1, coordinate [x,y,z], no reassignment

        CASE DEFAULT
          IErr = 1
          IF(l_alert(IErr,"PreferredBasis",&
              "Wyckoff Symbol for space group 36, C m c 21, not recognised")) RETURN
        END SELECT
!!$  CASE(37)
!!$  CASE(38)
!!$  CASE(39)
!!$  CASE(40)
!!$  CASE(41)
!!$  CASE(42)
!!$  CASE(43)
!!$  CASE(44)
!!$  CASE(45)
!!$  CASE(46)
!!$  CASE(47)
!!$  CASE(48)
!!$  CASE(49)
!!$  CASE(50)
!!$  CASE(51)
!!$  CASE(52)
!!$  CASE(53)
!!$  CASE(54)
!!$  CASE(55)
!!$  CASE(56)
!!$  CASE(57)
!!$  CASE(58)
!!$  CASE(59)
!!$  CASE(60)
!!$  CASE(61)
!!$  CASE(62)
     CASE(63)!C m c m.  No equivalent positions swap axes so no reassignments
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m, coordinate [0,0,0] & eq, no reassignment

      CASE('b')!point symmetry 2/m, coordinate [0,1/2,0] & eq, no reassignment

      CASE('c')!point symmetry mm2, coordinate [0,y,1/4] & eq, no reassignment

      CASE('d')!point symmetry -1, coordinate [1/4,1/4,0] & eq, no reassignment

      CASE('e')!point symmetry 2(x), coordinate [x,0,0] & eq, no reassignment

      CASE('f')!point symmetry m(x), coordinate [0,y,z] & eq, no reassignment

      CASE('g')!point symmetry m(z), coordinate [x,y,1/4] & eq, no reassignment

      CASE('h')!point symmetry 1, coordinate [x,y,z] & eq, no reassignment

      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"PreferredBasis",&
              "Wyckoff Symbol for space group 63, C m c m, not recognised")) RETURN
      END SELECT

     CASE(64)!C m c a.  No equivalent positions swap axes so no reassignments
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m, coordinate [0,0,0] & eq, no reassignment

      CASE('b')!point symmetry 2/m, coordinate [1/2,0,0] & eq, no reassignment

      CASE('c')!point symmetry -1, coordinate [1/4,1/4,0] & eq, no reassignment

      CASE('d')!point symmetry 2(x), coordinate [x,0,0] & eq, no reassignment

      CASE('e')!point symmetry 2(y), coordinate [1/4,y,1/4] & eq, no reassignment

      CASE('f')!point symmetry m(x), coordinate [0,y,z] & eq, no reassignment

      CASE('g')!point symmetry 1, coordinate [x,y,z] & eq, no reassignment

      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"PreferredBasis",&
              "Wyckoff Symbol for space group 64, C m c a, not recognised")) RETURN    
      END SELECT

!!$  CASE(65)
!!$  CASE(66)
!!$  CASE(67)
     CASE(68)!C c c a.  No equivalent positions swap axes so no reassignments
      !N.B. multiple origin choices allowed, here origin at 222, -1 at[1/4,0,1/4]
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222, coordinate [0,0,0] & eq, no reassignment

      CASE('b')!point symmetry 222, coordinate [0,0,1/2] & eq, no reassignment

      CASE('c')!point symmetry -1, coordinate [1/4,0,1/4] & eq, no reassignment

      CASE('d')!point symmetry -1, coordinate [0,1/4,1/4] & eq, no reassignment

      CASE('e')!point symmetry 2(x), coordinate [x,0,0] & eq, no reassignment

      CASE('f')!point symmetry 2(y), coordinate [0,y,0] & eq, no reassignment

      CASE('g')!point symmetry 2(z), coordinate [0,0,z] & eq, no reassignment

      CASE('h')!point symmetry 2(z), coordinate [1/4,1/4,z] & eq, no reassignment

      CASE('i')!point symmetry 1, coordinate [x,y,z] & eq, no reassignment

      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"PreferredBasis",&
              "Wyckoff Symbol for space group 68, C c c a, not recognised")) RETURN
      END SELECT

!!$  CASE(69)
!!$  CASE(70)
!!$  CASE(71)
!!$  CASE(72)
!!$  CASE(73)
!!$  CASE(74)
!!$  CASE(75)
!!$  CASE(76)
!!$  CASE(77)
!!$  CASE(78)
!!$  CASE(79)
!!$  CASE(80)
!!$  CASE(81)
!!$  CASE(82)
!!$  CASE(83)
!!$  CASE(84)
!!$  CASE(85)
!!$  CASE(86)
!!$  CASE(87)
!!$  CASE(88)
!!$  CASE(89)
!!$  CASE(90)
!!$  CASE(91)
!!$  CASE(92)
!!$  CASE(93)
!!$  CASE(94)
!!$  CASE(95)
!!$  CASE(96)
!!$  CASE(97)
!!$  CASE(98)
     CASE(99)!P 4 m m
       SELECT CASE (SWyckoff)
       CASE('a')!point symmetry 4mm, coordinate [0,0,z], no reassignment

       CASE('b')!point symmetry 4mm, coordinate [1/2,1/2,z], no reassignment

       CASE('c')!point symmetry mm, coordinate [1/2,0,z] & eq, no reassignment

       CASE('d')!point symmetry m, coordinate [x,x,z] & eq
         !Change equivalent coordinate [x,-x,z] to [x,x,z]
         RBasisAtomPosition(ind,2)=RBasisAtomPosition(ind,1)

       CASE('e')!point symmetry m, coordinate [x,0,z] & eq
         !Change equivalent coordinate [0,x,z] to [x,0,z]
         IF (ABS(RBasisAtomPosition(ind,1)).LT.TINY) THEN
           RBasisAtomPosition(ind,1)=RBasisAtomPosition(ind,2)
           RBasisAtomPosition(ind,2)=ZERO
         END IF

       CASE('f')!point symmetry m, coordinate [x,1/2,z] & eq
         !Change equivalent coordinate [1/2,x,z] to [x,1/2,z]
         IF (ABS(RBasisAtomPosition(ind,1)-HALF).LT.TINY) THEN
           RBasisAtomPosition(ind,1)=RBasisAtomPosition(ind,2)
           RBasisAtomPosition(ind,2)=HALF
         END IF

       CASE('g')!point symmetry 1, coordinate [x,y,z] & eq, no reassignment

       CASE DEFAULT
         IErr = 1
         IF(l_alert(IErr,"PreferredBasis",&
              "Wyckoff Symbol for space group 99, P 4 m m, not recognised")) RETURN
              
       END SELECT
     
!!$  CASE(100)
!!$  CASE(101)
!!$  CASE(102)
!!$  CASE(103)
!!$  CASE(104)
!!$  CASE(105)
!!$  CASE(106)
!!$  CASE(107)
!!$  CASE(108)
!!$  CASE(109)
!!$  CASE(110)
!!$  CASE(111)
!!$  CASE(112)
!!$  CASE(113)
!!$  CASE(114)
!!$  CASE(115)
!!$  CASE(116)
!!$  CASE(117)
!!$  CASE(118)
!!$  CASE(119)
!!$  CASE(120)
!!$  CASE(121)
!!$  CASE(122)
!!$  CASE(123)
!!$  CASE(124)
!!$  CASE(125)
!!$  CASE(126)
!!$  CASE(127)
!!$  CASE(128)
!!$  CASE(129)
!!$  CASE(130)
!!$  CASE(131)
!!$  CASE(132)
!!$  CASE(133)
!!$  CASE(134)
!!$  CASE(135)
!!$  CASE(136)
!!$  CASE(137)
!!$  CASE(138)
     CASE(139)!I4/m m m
       SELECT CASE (SWyckoff)
       CASE('a')!point symmetry 4/mmm, coordinate [0,0,0], no reassignment
      
       CASE('b')!point symmetry 4/mmm, coordinate [0,0,1/2], no reassignment
      
       CASE('c')!point symmetry mmm, coordinate [0,1/2,0] or [1/2,0,0], no reassignment
      
       CASE('d')!point symmetry -4m2, coordinate [0,1/2,1/4] or [1/2,0,1/4], no reassignment
      
       CASE('e')!point symmetry 4mm, coordinate [0,0,z], no reassignment

       CASE('f')!point symmetry 2/m, coordinate [1/4,1/4,1/4] & eq, no reassignment
        
       CASE('g')!point symmetry mm, coordinate [0,1/2,z] & eq, no reassignment

       CASE('h')!point symmetry mm, coordinate [x,x,0]
         !Change to [x,x,0]
         RBasisAtomPosition(ind,2)=RBasisAtomPosition(ind,1)
         
       CASE('i')!point symmetry mm, coordinate [x,0,0] & eq
         !Change equivalent coordinate [0,y,0]
         IF (ABS(RBasisAtomPosition(ind,1)).LT.TINY) THEN
           RBasisAtomPosition(ind,1)=RBasisAtomPosition(ind,2)
           RBasisAtomPosition(ind,2)=ZERO
         END IF
         
       CASE('j')!point symmetry mm, coordinate [x,1/2,0] & eq
         !Change equivalent coordinate [1/2,y,0]
         IF (ABS(RBasisAtomPosition(ind,1)-HALF).LT.TINY) THEN
           RBasisAtomPosition(ind,1)=RBasisAtomPosition(ind,2)
           RBasisAtomPosition(ind,2)=HALF
         END IF
         
       CASE('k')!point symmetry 2, coordinate [x,1/2+x,1/4] & eq
         !Change to [x,1/2+x,1/4]
         RBasisAtomPosition(ind,2)=RBasisAtomPosition(ind,1)+HALF
       CASE('l')!point symmetry m, coordinate [x,y,0] & eq - no change

       CASE('m')!point symmetry m, coordinate [x,x,z] & eq
         !Change to [x,x,z]
         RBasisAtomPosition(ind,2)=RBasisAtomPosition(ind,1)
         
       CASE('n')!point symmetry m, coordinate [x,0,z] & eq
         !Change equivalent coordinate [0,y,z]
         IF (ABS(RBasisAtomPosition(ind,1)).LT.TINY) THEN
           RBasisAtomPosition(ind,1)=RBasisAtomPosition(ind,2)
           RBasisAtomPosition(ind,2)=ZERO
         END IF
         
       CASE('o')!point symmetry 1, coordinate [x,y,z] & eq, no reassignment

       CASE DEFAULT
         IErr = 1
         IF(l_alert(IErr,"PreferredBasis",&
             "Wyckoff Symbol for space group 139, I4/m m m, not recognised")) RETURN   
       END SELECT
!!$  CASE(140)
!!$  CASE(141)
     CASE(142)!I41/acd
       SELECT CASE (SWyckoff)

       CASE('a')!point symmetry -4, coordinate [0,1/4,3/8], no reassignment
     
       CASE('b')!point symmetry 222, coordinate [0,1/4,1/8], no reassignment
     
       CASE('c')!point symmetry -1, coordinate [0,0,0], no reassignment

       CASE('d')!point symmetry 2, coordinate [0,1/4,z], no reassignment

       CASE('e')!point symmetry 2, coordinate [x,0,1/4], & eq
         IF (ABS(RBasisAtomPosition(ind,3)).LT.TINY) THEN
           !change equivalent coordinate [1/4,1/4-x,0] & [1/4,3/4-x,0]
           IF (ABS(RBasisAtomPosition(ind,1)-0.25).LT.TINY) THEN
             RBasisAtomPosition(ind,1)=0.25-RBasisAtomPosition(ind,2)
             RBasisAtomPosition(ind,2)=ZERO
             RBasisAtomPosition(ind,2)=0.25
           END IF
           !change equivalent coordinate [3/4,3/4+x,0] & [3/4,1/4+x,0]
           IF (ABS(RBasisAtomPosition(ind,1)-0.75).LT.TINY) THEN
             RBasisAtomPosition(ind,1)=RBasisAtomPosition(ind,2)-0.25
             RBasisAtomPosition(ind,2)=ZERO
             RBasisAtomPosition(ind,2)=0.25
           END IF
         END IF
 
       CASE('f')!point symmetry 2, coordinate [x,x+1/4,1/8], & eq
         !atoms at height 3/8,5/8 have a 2-fold about [x,-x,.] and need to be changed
         IF (ABS(RBasisAtomPosition(ind,3)-0.375).LT.TINY) THEN
           RBasisAtomPosition(ind,2)=RBasisAtomPosition(ind,1)+0.25
           RBasisAtomPosition(ind,3)=0.125
         END IF
         IF (ABS(RBasisAtomPosition(ind,3)-0.625).LT.TINY) THEN
           RBasisAtomPosition(ind,2)=RBasisAtomPosition(ind,1)+0.25
           RBasisAtomPosition(ind,2)=0.125
         END IF

       CASE('g')!point symmetry 1, coordinate [x,y,z], no reassignment

       CASE DEFAULT
         IErr = 1
         IF(l_alert(IErr,"PreferredBasis",&
             "Wyckoff Symbol for space group 142, I41/a c d, not recognised")) RETURN
      END SELECT
!!$  CASE(143)
!!$  CASE(144)
!!$  CASE(145)
!!$  CASE(146)
!!$  CASE(147)
!!$  CASE(148)
!!$  CASE(149)
!!$  CASE(150)
!!$  CASE(151)
!!$  CASE(152)
!!$  CASE(153)
!!$  CASE(154)
!!$  CASE(155)
!!$  CASE(156)
!!$  CASE(157)
!!$  CASE(158)
!!$  CASE(159)
!!$  CASE(160)
!!$  CASE(161)
!!$  CASE(162)
!!$  CASE(163)
!!$  CASE(164)
!!$  CASE(165)
!!$  CASE(166)
!!$  CASE(167)
!!$  CASE(168)
!!$  CASE(169)
!!$  CASE(170)
!!$  CASE(171)
!!$  CASE(172)
!!$  CASE(173)
!!$  CASE(174)
!!$  CASE(175)
!!$  CASE(176)
!!$  CASE(177)
!!$  CASE(178)
!!$  CASE(179)
!!$  CASE(180)
!!$  CASE(181)
!!$  CASE(182)
!!$  CASE(183)
!!$  CASE(184)
!!$  CASE(185)
!!$  CASE(186)
!!$  CASE(187)
!!$  CASE(188)
!!$  CASE(189)
!!$  CASE(190)
!!$  CASE(191)
!!$  CASE(192)
!!$  CASE(193)
!!$  CASE(194)
!!$  CASE(195)
!!$  CASE(196)
!!$  CASE(197)
!!$  CASE(198)
!!$  CASE(199)
!!$  CASE(200)
!!$  CASE(201)
!!$  CASE(202)
!!$  CASE(203)
!!$  CASE(204)
!!$  CASE(205)
!!$  CASE(206)
!!$  CASE(207)
!!$  CASE(208)
!!$  CASE(209)
!!$  CASE(210)
!!$  CASE(211)
!!$  CASE(212)
!!$  CASE(213)
!!$  CASE(214)
!!$  CASE(215)
    CASE(216)!F-43m
        CALL Flattice(RBasisAtomPosition(ind,:),IErr)!puts the basis next to the origin
        SELECT CASE (SWyckoff)
        CASE('a')!point symmetry -43m, coordinate [0,0,0], no reassignment

        CASE('b')!point symmetry -43m, coordinate [1/2,1/2,1/2], no reassignment

        CASE('c')!point symmetry -43m, coordinate [1/4,1/4,1/4], no reassignment

        CASE('d')!point symmetry -43m, coordinate [3/4,3/4,3/4], no reassignment

        CASE('e')!point symmetry 3m, coordinate [x,x,x] & eq
          RBasisAtomPosition(ind,2)=RBasisAtomPosition(ind,1)
          RBasisAtomPosition(ind,3)=RBasisAtomPosition(ind,1)

        CASE('f')!point symmetry mm, coordinate [x,0,0] & eq
          !convert to [x,0,0]
          IF (ABS(RBasisAtomPosition(ind,1)).LE.TINY) THEN
            IF (ABS(RBasisAtomPosition(ind,2)).GT.TINY) &!it's[0,y,0]
              RBasisAtomPosition(ind,1)=RBasisAtomPosition(ind,2)
            IF (ABS(RBasisAtomPosition(ind,3)).GT.TINY) &!it's [0,0,z]
              RBasisAtomPosition(ind,1)=RBasisAtomPosition(ind,3)
            RBasisAtomPosition(ind,2)=ZERO
            RBasisAtomPosition(ind,3)=ZERO
          END IF

        CASE('g')!point symmetry mm, coordinate [x,1/4,1/4] & eq
          !move non-fixed coordinate to x by changing cyclically
          DO jnd=1,3!find the coord that isn't 1/4 or 3/4
            IF (ABS(MOD(4*RBasisAtomPosition(ind,jnd),1.0)).GT.TINY) knd=jnd
          END DO
          !cyclically move the coords
          RBasisAtomPosition(ind,1)=RBasisAtomPosition(ind,knd)
          RBasisAtomPosition(ind,2)=RBasisAtomPosition(ind,MOD(knd+1,3))
          RBasisAtomPosition(ind,3)=RBasisAtomPosition(ind,MOD(knd+2,3))

        CASE('h')!point symmetry m, coordinate [x,x,z] & eq
          IF(ABS(RBasisAtomPosition(ind,2))-ABS(RBasisAtomPosition(ind,3)).GT.TINY) THEN!it's zxx or -zx-x
            IF(ABS(RBasisAtomPosition(ind,2)-RBasisAtomPosition(ind,3)).GT.TINY) THEN!zxx
              RBasisAtomPosition(ind,3)=RBasisAtomPosition(ind,1)
              RBasisAtomPosition(ind,1)=RBasisAtomPosition(ind,2)!now it's xxz
            ELSE!-zx-x
              RBasisAtomPosition(ind,3)=-RBasisAtomPosition(ind,1)
              RBasisAtomPosition(ind,1)=RBasisAtomPosition(ind,2)!now it's xxz
            END IF
          END IF
          IF(ABS(RBasisAtomPosition(ind,1))-ABS(RBasisAtomPosition(ind,3)).GT.TINY) THEN!it's xzx or x-z-x
            IF(ABS(RBasisAtomPosition(ind,1)-RBasisAtomPosition(ind,3)).GT.TINY) THEN!xzx
              RBasisAtomPosition(ind,3)=RBasisAtomPosition(ind,2)
              RBasisAtomPosition(ind,2)=RBasisAtomPosition(ind,1)!now it's xxz
            ELSE!x-z-x
              RBasisAtomPosition(ind,3)=-RBasisAtomPosition(ind,2)
              RBasisAtomPosition(ind,2)=RBasisAtomPosition(ind,1)!now it's xxz
            END IF
          END IF

        CASE('i')!point symmetry 1, coordinate [x,y,z], no reassignment


        END SELECT
!!$  CASE(217)
!!$  CASE(218)
!!$  CASE(219)
!!$  CASE(220)
!!$  CASE(221)
!!$  CASE(222)
!!$  CASE(223)
!!$  CASE(224)
!!$  CASE(225)
!!$  CASE(226)
!!$  CASE(227)
!!$  CASE(228)
!!$  CASE(229)
!!$  CASE(230)
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"PreferredBasis","space group not recognised")) RETURN 
      END SELECT
    END DO
  END SUBROUTINE PreferredBasis

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>Put coordinates next to the origin in a face-centred lattice
  SUBROUTINE Flattice(Rcoord,IErr)

    USE MyNumbers
    USE message_mod

    IMPLICIT NONE

    INTEGER(IKIND), INTENT(OUT) :: IErr
    REAL(RKIND), INTENT(INOUT) :: Rcoord(3)

    IErr=0!I can't see how this subroutine can go wrong...
    !make sure all coords are positive
    Rcoord=MOD(Rcoord,1.0)
    ![1/2,1/2,0]
    IF(((Rcoord(1)-HALF).GT.TINY).AND.((Rcoord(2)-HALF).GT.TINY)) THEN
      Rcoord(1)=Rcoord(1)-HALF
      Rcoord(2)=Rcoord(2)-HALF
    END IF
    ![1/2,0,1/2]
    IF(((Rcoord(1)-HALF).GT.TINY).AND.((Rcoord(3)-HALF).GT.TINY)) THEN
      Rcoord(1)=Rcoord(1)-HALF
      Rcoord(3)=Rcoord(3)-HALF
    END IF
    ![0,1/2,1/2]
    IF(((Rcoord(2)-HALF).GT.TINY).AND.((Rcoord(3)-HALF).GT.TINY)) THEN
      Rcoord(2)=Rcoord(2)-HALF
      Rcoord(3)=Rcoord(3)-HALF
    END IF

  END SUBROUTINE Flattice

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: 
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE CountAllowedMovements(ISpaceGrp,SWyckoff,IVectors,IErr)

    USE MyNumbers
    USE message_mod

    IMPLICIT NONE

    INTEGER(IKIND),INTENT(IN) :: ISpaceGrp
    CHARACTER(1), INTENT(IN) :: SWyckoff
    INTEGER(IKIND), INTENT(OUT) :: IVectors, IErr
    
    IErr=0
    SELECT CASE(ISpaceGrp)
!!$  CASE(1)
    CASE(1)!P 1
      SELECT CASE (SWyckoff)
        CASE('a')
          IVectors = 3
        CASE DEFAULT
          IErr = 1
          IF(l_alert(IErr,"DetermineAllowedMovements",&
         "Wyckoff Symbol for space group 1, P1, not recognised")) RETURN
        END SELECT
!!$  CASE(2)
    CASE(2)!P -1
      SELECT CASE (SWyckoff)
        CASE('a')!point symmetry -1
          IVectors = 0
        CASE('b')!point symmetry -1
          IVectors = 0
        CASE('c')!point symmetry -1
          IVectors = 0
        CASE('d')!point symmetry -1
          IVectors = 0
        CASE('e')!point symmetry -1
          IVectors = 0
        CASE('f')!point symmetry -1
          IVectors = 0
        CASE('g')!point symmetry -1
          IVectors = 0
        CASE('h')!point symmetry -1
          IVectors = 0
        CASE('i')!point symmetry 1
          IVectors = 3
        CASE DEFAULT
          IErr = 1
          IF(l_alert(IErr,"PreferredBasis",&
                "Wyckoff Symbol for space group 2, P -1, not recognised")) RETURN
        END SELECT
!!$  CASE(3)
    CASE(3)!P 1 1 2
      SELECT CASE (SWyckoff)
        CASE('a')!point symmetry 2
          IVectors = 1
        CASE('b')!point symmetry 2
          IVectors = 1
        CASE('c')!point symmetry 2
          IVectors = 1
        CASE('d')!point symmetry 2
          IVectors = 1
        CASE('e')!point symmetry 1
          IVectors = 3
        CASE DEFAULT
          IErr = 1
          IF(l_alert(IErr,"PreferredBasis",&
                "Wyckoff Symbol for space group 3, P 1 1 2, not recognised")) RETURN
        END SELECT
!!$  CASE(4)
    CASE(4)!P 1 1 21
      SELECT CASE (SWyckoff)
        CASE('a')!point symmetry 1
          IVectors = 3
        CASE DEFAULT
          IErr = 1
          IF(l_alert(IErr,"PreferredBasis",&
                "Wyckoff Symbol for space group 4, P 1 1 21, not recognised")) RETURN
        END SELECT
!!$  CASE(5)
    CASE(5)!B 1 1 2
      SELECT CASE (SWyckoff)
        CASE('a')!point symmetry 2
          IVectors = 1
        CASE('b')!point symmetry 2
          IVectors = 1
        CASE('c')!point symmetry 1
          IVectors = 3
        CASE DEFAULT
          IErr = 1
          IF(l_alert(IErr,"PreferredBasis",&
                "Wyckoff Symbol for space group 5, B 1 1 2, not recognised")) RETURN
        END SELECT
!!$  CASE(6)
    CASE(6)!P 1 1 m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry m
        IVectors = 2
      CASE('b')!point symmetry m
        IVectors = 2
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"PreferredBasis",&
                "Wyckoff Symbol for space group 6, P 1 1 m, not recognised")) RETURN
      END SELECT
!!$  CASE(7)
    CASE(7)!P 1 1 b
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"PreferredBasis",&
              "Wyckoff Symbol for space group 7, P 1 1 b, not recognised")) RETURN
      END SELECT
!!$  CASE(8)
    CASE(8)!B 1 1 m
      SELECT CASE (SWyckoff)
        CASE('a')!point symmetry m
          IVectors = 2
        CASE('b')!point 1
          IVectors = 3
        CASE DEFAULT
          IErr = 1
          IF(l_alert(IErr,"DetermineAllowedMovements",&
                "Wyckoff Symbol for space group 8, B 1 1 m, not recognised")) RETURN  
        END SELECT 
!!$  CASE(9)
    CASE(9)!B 1 1 b 
      SELECT CASE (SWyckoff)
        CASE('a')!point symmetry 1
          IVectors = 3
        CASE DEFAULT
          IErr = 1
          IF(l_alert(IErr,"DetermineAllowedMovements",&
                "Wyckoff Symbol for space group 9, B 1 1 b, not recognised")) RETURN  
        END SELECT 
!!$  CASE(10)
    CASE(10)!P 1 1 2/m
      SELECT CASE (SWyckoff)
        CASE('a')!point symmetry 2/m
          IVectors = 0
        CASE('b')!point 2/m
          IVectors = 0
        CASE('c')!point symmetry 2/m
          IVectors = 0
        CASE('d')!point symmetry 2/m
          IVectors = 0
        CASE('e')!point symmetry 2/m
          IVectors = 0
        CASE('f')!point symmetry 2/m
          IVectors = 0
        CASE('g')!point symmetry 2/m
          IVectors = 0
        CASE('h')!point 2/m
          IVectors = 0
        CASE('i')!point symmetry 2
          IVectors = 1
        CASE('j')!point symmetry 2
          IVectors = 1
        CASE('k')!point symmetry 2
          IVectors = 1
        CASE('l')!point symmetry 2
          IVectors = 1
        CASE('m')!point symmetry m
          IVectors = 2
        CASE('n')!point m
          IVectors = 2
        CASE('o')!point symmetry 1
          IVectors = 3
        CASE DEFAULT
          IErr = 1
          IF(l_alert(IErr,"DetermineAllowedMovements",&
                "Wyckoff Symbol for space group 10, P 1 1 2/m, not recognised")) RETURN  
        END SELECT 
!!$  CASE(11)
    CASE(11)!P 1 1 21/m
      SELECT CASE (SWyckoff)
        CASE('a')!point symmetry -1
          IVectors = 0
        CASE('b')!point -1
          IVectors = 0
        CASE('c')!point symmetry -1
          IVectors = 0
        CASE('d')!point symmetry -1
          IVectors = 0
        CASE('e')!point symmetry m
          IVectors = 2
        CASE('f')!point symmetry 1
          IVectors = 3
        CASE DEFAULT
          IErr = 1
          IF(l_alert(IErr,"DetermineAllowedMovements",&
                "Wyckoff Symbol for space group 11, P 1 1 21/m, not recognised")) RETURN  
        END SELECT 
!!$  CASE(12)
    CASE(12)!B 1 1 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m
        IVectors = 0
      CASE('b')!point symmetry 2/m
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry -1
        IVectors = 0
      CASE('f')!point symmetry -1
        IVectors = 0
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry m
        IVectors = 2
      CASE('j')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 12, B 1 1 2/m, not recognised")) RETURN     
      END SELECT    
!!$  CASE(13)
    CASE(13)!P 1 1 2/b
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -1
        IVectors = 0
      CASE('b')!point symmetry -1
        IVectors = 0
      CASE('c')!point symmetry -1
        IVectors = 0
      CASE('d')!point symmetry -1
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 13, P 1 1 2/b, not recognised")) RETURN     
      END SELECT    
!!$  CASE(14)
    CASE(14)!P 1 1 21/b
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -1
        IVectors = 0
      CASE('b')!point symmetry -1
        IVectors = 0
      CASE('c')!point symmetry -1
        IVectors = 0
      CASE('d')!point symmetry -1
        IVectors = 0
      CASE('e')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 14, P 1 1 21/b, not recognised")) RETURN     
      END SELECT    
!!$  CASE(15)
    CASE(15)!C 1 2/c 1
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -1
        IVectors = 0
      CASE('b')!point symmetry -1
        IVectors = 0
      CASE('c')!point symmetry -1
        IVectors = 0
      CASE('d')!point symmetry -1
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 15, C 2/c, not recognised")) RETURN     
      END SELECT    
!!$  CASE(16)
    CASE(16)!P 2 2 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 222
        IVectors = 0
      CASE('f')!point symmetry 222
        IVectors = 0
      CASE('g')!point symmetry 222
        IVectors = 0
      CASE('h')!point symmetry 222
        IVectors = 0
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 2
        IVectors = 1
      CASE('m')!point symmetry 2
        IVectors = 1
      CASE('n')!point symmetry 2
        IVectors = 1
      CASE('o')!point symmetry 2
        IVectors = 1
      CASE('p')!point symmetry 2
        IVectors = 1
      CASE('q')!point symmetry 2
        IVectors = 1
      CASE('r')!point symmetry 2
        IVectors = 1
      CASE('s')!point symmetry 2
        IVectors = 1
      CASE('t')!point symmetry 2
        IVectors = 1
      CASE('u')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 16, P 2 2 2, not recognised")) RETURN     
      END SELECT    
!!$  CASE(17)
    CASE(17)!P 2 2 21
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 17, P 2 2 21, not recognised")) RETURN     
      END SELECT    
!!$  CASE(18)
    CASE(18)!P 21 21 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 18, P 21 21 2, not recognised")) RETURN     
      END SELECT    
!!$  CASE(19)
    CASE(19)!P 21 21 21
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 19, P 21 21 21, not recognised")) RETURN     
      END SELECT    
!!$  CASE(20)
    CASE(20)!C 2 2 21
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 20, C 2 2 21, not recognised")) RETURN     
      END SELECT    
!!$  CASE(21)
    CASE(21)!C 2 2 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 21, C 2 2 2, not recognised")) RETURN     
      END SELECT    
!!$  CASE(22)
    CASE(22)!F 2 2 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 22, F 2 2 2, not recognised")) RETURN     
      END SELECT    
!!$  CASE(23)
    CASE(23)!I 2 2 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 23, I 2 2 2, not recognised")) RETURN     
      END SELECT    
!!$  CASE(24)
    CASE(24)!I 21 21 21
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 24, I 21 21 21, not recognised")) RETURN     
      END SELECT    
!!$  CASE(25)
    CASE(25)!P m m 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mm
        IVectors = 1
      CASE('b')!point symmetry mm
        IVectors = 1
      CASE('c')!point symmetry mm
        IVectors = 1
      CASE('d')!point symmetry mm
        IVectors = 1
      CASE('e')!point symmetry m
        IVectors = 2
      CASE('f')!point symmetry m
        IVectors = 2
      CASE('g')!point symmetry m
        IVectors = 2
      CASE('h')!point symmetry m
        IVectors = 2
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 25, P m m 2, not recognised")) RETURN     
      END SELECT    
!!$  CASE(26)
    CASE(26)!P m c 21
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry m
        IVectors = 2
      CASE('b')!point symmetry m
        IVectors = 2
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 26, I4/m m m, not recognised")) RETURN     
      END SELECT    
!!$  CASE(27)
    CASE(27)!P c c 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 27, P c c 2, not recognised")) RETURN     
      END SELECT    
!!$  CASE(28)
    CASE(28)!P m a 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry m
        IVectors = 2
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 28, P m a 2, not recognised")) RETURN     
      END SELECT    
!!$  CASE(29)
    CASE(29)!P c a 21
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 29, P c a 21, not recognised")) RETURN     
      END SELECT    
!!$  CASE(30)
    CASE(30)!P n c 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 30, P n c 2, not recognised")) RETURN     
      END SELECT    
!!$  CASE(31)
    CASE(31)!P m n 21
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry m
        IVectors = 2
      CASE('b')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 31, P m n 21, not recognised")) RETURN     
      END SELECT    
!!$  CASE(32)
    CASE(32)!P b a 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 32, P b a 2, not recognised")) RETURN     
      END SELECT    
!!$  CASE(33)
    CASE(33)!P n a 21
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 33, P n a 21, not recognised")) RETURN     
      END SELECT    
!!$  CASE(34)
    CASE(34)!P n n 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 34, P n n 2, not recognised")) RETURN     
      END SELECT    
!!$  CASE(35)
    CASE(35)!C m m 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mm
        IVectors = 1
      CASE('b')!point symmetry mm
        IVectors = 1
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry m
        IVectors = 2
      CASE('e')!point symmetry m
        IVectors = 2
      CASE('f')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 35, C m m 2, not recognised")) RETURN     
      END SELECT    
!!$  CASE(36)   
    CASE(36)!C m c 21
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry m
        IVectors = 2
      CASE('b')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 36, C m c 21, not recognised")) RETURN     
      END SELECT
!!$  CASE(37)
    CASE(37)!C c c 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 37, C c c 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(38)
    CASE(38)!A m m 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mm
        IVectors = 1
      CASE('b')!point symmetry mm
        IVectors = 1
      CASE('c')!point symmetry m
        IVectors = 2
      CASE('d')!point symmetry m
        IVectors = 2
      CASE('e')!point symmetry m
        IVectors = 2
      CASE('f')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 38, A m m 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(39)
    CASE(39)!A b m 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry m
        IVectors = 2
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 39, A b m 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(40)
    CASE(40)!A m a 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry m
        IVectors = 2
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 40, A m a 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(41)
    CASE(41)!A b a 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 41, A b a 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(42)
    CASE(42)!F m m 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mm
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry m
        IVectors = 2
      CASE('d')!point symmetry m
        IVectors = 2
      CASE('e')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 42, F m m 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(43)
    CASE(43)!F d d 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 43, F d d 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(44)
    CASE(44)!I m m 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mm
        IVectors = 1
      CASE('b')!point symmetry mm
        IVectors = 1
      CASE('c')!point symmetry m
        IVectors = 2
      CASE('d')!point symmetry m
        IVectors = 2
      CASE('e')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 44, I m m 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(45)
    CASE(45)!I b a 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 45, I b a 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(46)
    CASE(46)!I m a 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry m
        IVectors = 2
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 46, I m a 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(47)
    CASE(47)!P 2/m 2/m 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mmm
        IVectors = 0
      CASE('b')!point symmetry mmm
        IVectors = 0
      CASE('c')!point symmetry mmm
        IVectors = 0
      CASE('d')!point symmetry mmm
        IVectors = 0
      CASE('e')!point symmetry mmm
        IVectors = 0
      CASE('f')!point symmetry mmm
        IVectors = 0
      CASE('g')!point symmetry mmm
        IVectors = 0
      CASE('h')!point symmetry mmm
        IVectors = 0
      CASE('i')!point symmetry mm
        IVectors = 1
      CASE('j')!point symmetry mm
        IVectors = 1
      CASE('k')!point symmetry mm
        IVectors = 1
      CASE('l')!point symmetry mm
        IVectors = 1
      CASE('m')!point symmetry mm
        IVectors = 1
      CASE('n')!point symmetry mm
        IVectors = 1
      CASE('o')!point symmetry mm
        IVectors = 1
      CASE('p')!point symmetry mm
        IVectors = 1
      CASE('q')!point symmetry mm
        IVectors = 1
      CASE('r')!point symmetry mm
        IVectors = 1
      CASE('s')!point symmetry mm
        IVectors = 1
      CASE('t')!point symmetry mm
        IVectors = 1
      CASE('u')!point symmetry m
        IVectors = 2
      CASE('v')!point symmetry m
        IVectors = 2
      CASE('w')!point symmetry m
        IVectors = 2
      CASE('x')!point symmetry m
        IVectors = 2
      CASE('y')!point symmetry m
        IVectors = 2
      CASE('z')!point symmetry m
        IVectors = 2
      CASE('alpha')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 47, P 2/m 2/m 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(48)
    CASE(48)!P 2/n 2/n 2/n
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry -1
        IVectors = 0
      CASE('f')!point symmetry -1
        IVectors = 0
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 2
        IVectors = 1
      CASE('m')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 48, P 2/n 2/n 2/n, not recognised")) RETURN     
      END SELECT
!!$  CASE(49)
    CASE(49)!P 2/c 2/c 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m
        IVectors = 0
      CASE('b')!point symmetry 2/m
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry 222
        IVectors = 0
      CASE('f')!point symmetry 222
        IVectors = 0
      CASE('g')!point symmetry 222
        IVectors = 0
      CASE('h')!point symmetry 222
        IVectors = 0
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 2
        IVectors = 1
      CASE('m')!point symmetry 2
        IVectors = 1
      CASE('n')!point symmetry 2
        IVectors = 1
      CASE('o')!point symmetry 2
        IVectors = 1
      CASE('p')!point symmetry 2
        IVectors = 1
      CASE('q')!point symmetry m
        IVectors = 2
      CASE('r')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 49, P 2/c 2/c 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(50)
    CASE(50)!P 2/b 2/a 2/n
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry -1
        IVectors = 0
      CASE('f')!point symmetry -1
        IVectors = 0
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 2
        IVectors = 1
      CASE('m')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 50, P 2/b 2/a 2/n, not recognised")) RETURN     
      END SELECT
!!$  CASE(51)
    CASE(51)!P 21/m 2/m 2/a
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m
        IVectors = 0
      CASE('b')!point symmetry 2/m
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry m
        IVectors = 2
      CASE('f')!point symmetry m
        IVectors = 2
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry m
        IVectors = 2
      CASE('j')!point symmetry m
        IVectors = 2
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 51, P 21/m 2/m 2/a, not recognised")) RETURN     
      END SELECT
!!$  CASE(52)
    CASE(52)!P 2/n 21/n 2/a
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -1
        IVectors = 0
      CASE('b')!point symmetry -1
        IVectors = 0
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 52, P 2/n 21/n 2/a, not recognised")) RETURN     
      END SELECT
!!$  CASE(53)
    CASE(53)!P 2/m 2/n 21/a
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m
        IVectors = 0
      CASE('b')!point symmetry 2/m
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry m
        IVectors = 2
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 53, P 2/m 2/n 21/a, not recognised")) RETURN     
      END SELECT
!!$  CASE(54)
    CASE(54)!P 21/c 2/c 2/a
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -1
        IVectors = 0
      CASE('b')!point symmetry -1
        IVectors = 0
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 54, P 21/c 2/c 2/a, not recognised")) RETURN     
      END SELECT
!!$  CASE(55)
    CASE(55)!P 21/b 21/a 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m
        IVectors = 0
      CASE('b')!point symmetry 2/m
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry m
        IVectors = 2
      CASE('h')!point symmetry m
        IVectors = 2
      CASE('i')!point symmetry m
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 55, P 21/b 21/a 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(56)
    CASE(56)!P 21/c 21/c 2/n
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -1
        IVectors = 0
      CASE('b')!point symmetry -1
        IVectors = 0
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 56, P 21/c 21/c 2/n, not recognised")) RETURN     
      END SELECT
!!$  CASE(57)
    CASE(57)!P 2/b 21/c 21/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -1
        IVectors = 0
      CASE('b')!point symmetry -1
        IVectors = 0
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry m
        IVectors = 2
      CASE('e')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 57, P 2/b 21/c 21/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(58)
    CASE(58)!P 21/n 21/n 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m
        IVectors = 0
      CASE('b')!point symmetry 2/m
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 0
      CASE('f')!point symmetry 2
        IVectors = 0
      CASE('g')!point symmetry m
        IVectors = 2
      CASE('h')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 58, P 21/n 21/n 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(59)
    CASE(59)!P 21/m 21/m 2/n
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mm
        IVectors = 1
      CASE('b')!point symmetry mm
        IVectors = 1
      CASE('c')!point symmetry -1
        IVectors = 0
      CASE('d')!point symmetry -1
        IVectors = 0
      CASE('e')!point symmetry m
        IVectors = 2
      CASE('f')!point symmetry m
        IVectors = 2
      CASE('g')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 59, P 21/m 21/m 2/n, not recognised")) RETURN     
      END SELECT
!!$  CASE(60)
    CASE(60)!P 21/b 2/c 21/n
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -1
        IVectors = 0
      CASE('b')!point symmetry -1
        IVectors = 0
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 60, P 21/b 2/c 21/n, not recognised")) RETURN     
      END SELECT
!!$  CASE(61)
    CASE(61)!P 21/b 21/c 21/a
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -1
        IVectors = 0
      CASE('b')!point symmetry -1
        IVectors = 0
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 61, P 21/b 21/c 21/a, not recognised")) RETURN     
      END SELECT
!!$  CASE(62)
    CASE(62)!P 21/n 21/m 21/a
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -1
        IVectors = 0
      CASE('b')!point symmetry -1
        IVectors = 0
      CASE('c')!point symmetry m
        IVectors = 2
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 62, P 21/n 21/m 21/a, not recognised")) RETURN     
      END SELECT
!!$  CASE(63)     
    CASE(63)!C 2/m 2/c 21/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m
        IVectors = 0
      CASE('b')!point symmetry 2/m
        IVectors = 0
      CASE('c')!point symmetry mm
        IVectors = 0
      CASE('d')!point symmetry -1
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry m
        IVectors = 2
      CASE('g')!point symmetry m
        IVectors = 2
      CASE('h')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 63, C 2/m 2/c 21/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(64)
    CASE(64)!C 2/m 2/c 21/a
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m
        IVectors = 0
      CASE('b')!point symmetry 2/m
        IVectors = 0
      CASE('c')!point symmetry -1
        IVectors = 0
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry m
        IVectors = 2
      CASE('g')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 64, C 2/m 2/c 21/a, not recognised")) RETURN     
      END SELECT
!!$  CASE(65)
    CASE(65)!C 2/m 2/m 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mmm
        IVectors = 0
      CASE('b')!point symmetry mmm
        IVectors = 0
      CASE('c')!point symmetry mmm
        IVectors = 0
      CASE('d')!point symmetry mmm
        IVectors = 0
      CASE('e')!point symmetry 2/m
        IVectors = 0
      CASE('f')!point symmetry 2/m
        IVectors = 0
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry mm
        IVectors = 1
      CASE('i')!point symmetry mm
        IVectors = 1
      CASE('j')!point symmetry mm
        IVectors = 1
      CASE('k')!point symmetry mm
        IVectors = 1
      CASE('l')!point symmetry mm
        IVectors = 1
      CASE('m')!point symmetry 2
        IVectors = 1
      CASE('n')!point symmetry m
        IVectors = 2
      CASE('o')!point symmetry m
        IVectors = 2
      CASE('p')!point symmetry m
        IVectors = 2
      CASE('q')!point symmetry m
        IVectors = 2
      CASE('r')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 65, C 2/m 2/m 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(66)
    CASE(66)!C 2/c 2/c 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry 2/m
        IVectors = 0
      CASE('f')!point symmetry 2/m
        IVectors = 0
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry m
        IVectors = 2
      CASE('m')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 66, C 2/c 2/c 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(67)
    CASE(67)!C 2/m 2/m 2/a
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry 2/m
        IVectors = 0
      CASE('f')!point symmetry 2/m
        IVectors = 0
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 2
        IVectors = 1
      CASE('m')!point symmetry m
        IVectors = 2
      CASE('n')!point symmetry m
        IVectors = 2
      CASE('o')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 67, C 2/m 2/m 2/a, not recognised")) RETURN     
      END SELECT
!!$  CASE(68)
    CASE(68)!C 2/c 2/c 2/a
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry -1
        IVectors = 0
      CASE('d')!point symmetry -1
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 68, C 2/c 2/c 2/a, not recognised")) RETURN     
      END SELECT
!!$  CASE(69)
    CASE(69)!F 2/m 2/m 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mmm
        IVectors = 0
      CASE('b')!point symmetry mmm
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry 2/m
        IVectors = 0
      CASE('f')!point symmetry 222
        IVectors = 0
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry mm
        IVectors = 1
      CASE('i')!point symmetry mm
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 2
        IVectors = 1
      CASE('m')!point symmetry m
        IVectors = 2
      CASE('n')!point symmetry m
        IVectors = 2
      CASE('o')!point symmetry m
        IVectors = 2    
      CASE('p')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 69, F 2/m 2/m 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(70)
    CASE(70)!F 2/d 2/d 2/d
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry -1
        IVectors = 0
      CASE('d')!point symmetry -1
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 70, F 2/d 2/d 2/d, not recognised")) RETURN     
      END SELECT
!!$  CASE(71)
    CASE(71)!I 2/m 2/m 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mmm
        IVectors = 0
      CASE('b')!point symmetry mmm
        IVectors = 0
      CASE('c')!point symmetry mmm
        IVectors = 0
      CASE('d')!point symmetry mmm
        IVectors = 0
      CASE('e')!point symmetry mm
        IVectors = 1
      CASE('f')!point symmetry mm
        IVectors = 1
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry mm
        IVectors = 1
      CASE('i')!point symmetry mm
        IVectors = 1
      CASE('j')!point symmetry mm
        IVectors = 1
      CASE('k')!point symmetry -1
        IVectors = 0
      CASE('l')!point symmetry m
        IVectors = 2
      CASE('m')!point symmetry m
        IVectors = 2
      CASE('n')!point symmetry m
        IVectors = 2
      CASE('o')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 71, I 2/m 2/m 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(72)
    CASE(72)!I 2/b 2/a 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry -1
        IVectors = 0
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry m
        IVectors = 2
      CASE('k')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 72, I 2/b 2/a 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(73)
    CASE(73)!I 2/b 2/c 2/a
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -1
        IVectors = 0
      CASE('b')!point symmetry -1
        IVectors = 0
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 73, !I 2/b 2/c 2/a, not recognised")) RETURN     
      END SELECT
!!$  CASE(74)
    CASE(74)!I 2/m 2/m 2/a
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m
        IVectors = 0
      CASE('b')!point symmetry 2/m
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry mm
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry m
        IVectors = 2
      CASE('i')!point symmetry m
        IVectors = 2
      CASE('j')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 74, I 2/m 2/m 2/a, not recognised")) RETURN     
      END SELECT
!!$  CASE(75)
    CASE(75)!P 4
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 4
        IVectors = 1
      CASE('b')!point symmetry 4
        IVectors = 1
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 75, P 4, not recognised")) RETURN     
      END SELECT
!!$  CASE(76)
    CASE(76)!P 41
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 76, P 41, not recognised")) RETURN     
      END SELECT
!!$  CASE(77)
    CASE(77)!P 42
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 77, P 42, not recognised")) RETURN     
      END SELECT
!!$  CASE(78)
    CASE(78)!P 43
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 78, P 43, not recognised")) RETURN     
      END SELECT
!!$  CASE(79)
    CASE(79)!I 4
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 4
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 179, I 4, not recognised")) RETURN     
      END SELECT
!!$  CASE(80)
    CASE(80)!I 41
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 80, I 41, not recognised")) RETURN     
      END SELECT
!!$  CASE(81)
    CASE(81)!P -4
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -4
        IVectors = 0
      CASE('b')!point symmetry -4
        IVectors = 0
      CASE('c')!point symmetry -4
        IVectors = 0
      CASE('d')!point symmetry -4
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 81, P -4, not recognised")) RETURN     
      END SELECT
!!$  CASE(82)
    CASE(82)!I -4
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -4
        IVectors = 0
      CASE('b')!point symmetry -4
        IVectors = 0
      CASE('c')!point symmetry -4
        IVectors = 0
      CASE('d')!point symmetry -4
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 82, I -4, not recognised")) RETURN     
      END SELECT
!!$  CASE(83)
    CASE(83)!P 4/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 4/m
        IVectors = 0
      CASE('b')!point symmetry 4/m
        IVectors = 0
      CASE('c')!point symmetry 4/m
        IVectors = 0
      CASE('d')!point symmetry 4/m
        IVectors = 0
      CASE('e')!point symmetry 2/m
        IVectors = 0
      CASE('f')!point symmetry 2/m
        IVectors = 0
      CASE('g')!point symmetry 4
        IVectors = 1
      CASE('h')!point symmetry 4
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry m
        IVectors = 2
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry 1
        IVectors = 3        
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 83, P 4/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(84)
    CASE(84)!P 42/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m
        IVectors = 0
      CASE('b')!point symmetry 2/m
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry -4
        IVectors = 0
      CASE('f')!point symmetry -4
        IVectors = 0
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry m
        IVectors = 2
      CASE('k')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 84, P 42/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(85)
    CASE(85)!P 4/n
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -4
        IVectors = 0
      CASE('b')!point symmetry -4
        IVectors = 1
      CASE('c')!point symmetry 4
        IVectors = 1
      CASE('d')!point symmetry -1
        IVectors = 0
      CASE('e')!point symmetry -1
        IVectors = 0
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 85, P 4/n, not recognised")) RETURN     
      END SELECT
!!$  CASE(86)
    CASE(86)!P 42/n
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -4
        IVectors = 0
      CASE('b')!point symmetry -4
        IVectors = 0
      CASE('c')!point symmetry -1
        IVectors = 0
      CASE('d')!point symmetry -1
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 86, P 42/n, not recognised")) RETURN     
      END SELECT
!!$  CASE(87)
    CASE(87)!I 4/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 4/m
        IVectors = 0
      CASE('b')!point symmetry 4/m
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry -4
        IVectors = 0
      CASE('e')!point symmetry 4
        IVectors = 1
      CASE('f')!point symmetry -1
        IVectors = 0
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry m
        IVectors = 2
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 87, I 4/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(88)
    CASE(88)!I 41/a
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -4
        IVectors = 0
      CASE('b')!point symmetry -4
        IVectors = 0
      CASE('c')!point symmetry -1
        IVectors = 0
      CASE('d')!point symmetry -1
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 88, I 41/a, not recognised")) RETURN     
      END SELECT
!!$  CASE(89)
    CASE(89)!P 4 2 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 42
        IVectors = 0
      CASE('b')!point symmetry 42
        IVectors = 0
      CASE('c')!point symmetry 42
        IVectors = 0
      CASE('d')!point symmetry 42
        IVectors = 0
      CASE('e')!point symmetry 222
        IVectors = 0
      CASE('f')!point symmetry 222
        IVectors = 0
      CASE('g')!point symmetry 4
        IVectors = 1
      CASE('h')!point symmetry 4
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 2
        IVectors = 1
      CASE('m')!point symmetry 2
        IVectors = 1
      CASE('n')!point symmetry 2
        IVectors = 1
      CASE('o')!point symmetry 2
        IVectors = 1
      CASE('p')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 89, P 4 2 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(90)
    CASE(90)!P 4 21 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 4
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 90, P 4 21 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(91)
    CASE(91)!P 41 2 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 91, P 41 2 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(92)
    CASE(92)!P 41 21 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 92, P 41 21 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(93)
    CASE(93)!P 42 2 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 222
        IVectors = 0
      CASE('f')!point symmetry 222
        IVectors = 0
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 2
        IVectors = 1
      CASE('m')!point symmetry 2
        IVectors = 1
      CASE('n')!point symmetry 2
        IVectors = 1
      CASE('o')!point symmetry 2
        IVectors = 1
      CASE('p')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 93, P 42 2 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(94)
    CASE(94)!P 42 21 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 94, P 42 21 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(95)
    CASE(95)!P 43 2 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 95, P 43 2 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(96)
    CASE(96)!P 43 21 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 96, P 43 21 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(97)
    CASE(97)!I 4 2 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 42
        IVectors = 0
      CASE('b')!point symmetry 42
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 4
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 97, I 4 2 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(98)
    CASE(98)!I 41 2 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 98, I 41 2 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(99)
    CASE(99)!P 4 m m 
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 4mm
        IVectors = 1
      CASE('b')!point symmetry 4mm
        IVectors = 1
      CASE('c')!point symmetry mm
        IVectors = 1
      CASE('d')!point symmetry m
        IVectors = 2
      CASE('e')!point symmetry m
        IVectors = 2
      CASE('f')!point symmetry m
        IVectors = 2
      CASE('g')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 99, P 4 m m , not recognised")) RETURN     
      END SELECT
!!$  CASE(100)
    CASE(100)!P 4 b m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 4
        IVectors = 1
      CASE('b')!point symmetry mm
        IVectors = 1
      CASE('c')!point symmetry m
        IVectors = 2
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 100, P 4 b m, not recognised")) RETURN     
      END SELECT
!!$  CASE(101)
    CASE(101)!P 42 c m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mm
        IVectors = 1
      CASE('b')!point symmetry mm
        IVectors = 1
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry m
        IVectors = 2
      CASE('e')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 101, P 42 c m, not recognised")) RETURN     
      END SELECT
!!$  CASE(102)
    CASE(102)!P 42 n m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mm
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry m
        IVectors = 2
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 102, P 42 n m, not recognised")) RETURN     
      END SELECT
!!$  CASE(103)
    CASE(103)!P 4 c c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 4
        IVectors = 1
      CASE('b')!point symmetry 4
        IVectors = 1
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 103, P 4 c c, not recognised")) RETURN     
      END SELECT
!!$  CASE(104)
    CASE(104)!P 4 n c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 4
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 104, P 4 n c, not recognised")) RETURN     
      END SELECT
!!$  CASE(105)
    CASE(105)!P 42 m c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mm
        IVectors = 1
      CASE('b')!point symmetry mm
        IVectors = 1
      CASE('c')!point symmetry mm
        IVectors = 1
      CASE('d')!point symmetry m
        IVectors = 2
      CASE('e')!point symmetry m
        IVectors = 2
      CASE('f')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 105, P 42 m c, not recognised")) RETURN     
      END SELECT
!!$  CASE(106)
    CASE(106)!P 42 b c 
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 106, P 42 b c , not recognised")) RETURN     
      END SELECT
!!$  CASE(107)
    CASE(107)!I 4 m m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 4mm
        IVectors = 1
      CASE('b')!point symmetry mm
        IVectors = 1
      CASE('c')!point symmetry m
        IVectors = 2
      CASE('d')!point symmetry m
        IVectors = 2
      CASE('e')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 107, I4 m m, not recognised")) RETURN     
      END SELECT
!!$  CASE(108)
    CASE(108)!I 4 c m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 4
        IVectors = 1
      CASE('b')!point symmetry mm
        IVectors = 1
      CASE('c')!point symmetry m
        IVectors = 2
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 108, I 4 c m, not recognised")) RETURN     
      END SELECT
!!$  CASE(109)
    CASE(109)!I 41 m d
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mm
        IVectors = 1
      CASE('b')!point symmetry m
        IVectors = 2
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 109, I 41 m d, not recognised")) RETURN     
      END SELECT
!!$  CASE(110)
    CASE(110)!I 41 c d
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 110, I 41 c d, not recognised")) RETURN     
      END SELECT
!!$  CASE(111)
    CASE(111)!P -4 2 m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -42m
        IVectors = 0
      CASE('b')!point symmetry -42m
        IVectors = 0
      CASE('c')!point symmetry -42m
        IVectors = 0
      CASE('d')!point symmetry -42m
        IVectors = 0
      CASE('e')!point symmetry 222
        IVectors = 0
      CASE('f')!point symmetry 222
        IVectors = 0
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry mm
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 2
        IVectors = 1
      CASE('m')!point symmetry 2
        IVectors = 1
      CASE('n')!point symmetry m
        IVectors = 2
      CASE('o')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 111, P -4 2 m, not recognised")) RETURN     
      END SELECT
!!$  CASE(112)
    CASE(112)!P -4 2 c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry -4
        IVectors = 0
      CASE('f')!point symmetry -4
        IVectors = 0
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 2
        IVectors = 1
      CASE('m')!point symmetry 2
        IVectors = 1
      CASE('n')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 112, P -4 2 c, not recognised")) RETURN     
      END SELECT
!!$  CASE(113)
    CASE(113)!P -4 21 m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -4
        IVectors = 0
      CASE('b')!point symmetry -4
        IVectors = 0
      CASE('c')!point symmetry mm
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry m
        IVectors = 2
      CASE('f')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 113, P -4 21 m, not recognised")) RETURN     
      END SELECT
!!$  CASE(114)
    CASE(114)!P -4 21 c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -4
        IVectors = 0
      CASE('b')!point symmetry -4
        IVectors = 0
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 114, P -4 21 c, not recognised")) RETURN     
      END SELECT
!!$  CASE(115)
    CASE(115)!P -4 m 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -42m
        IVectors = 0
      CASE('b')!point symmetry -42m
        IVectors = 0
      CASE('c')!point symmetry -42m
        IVectors = 0
      CASE('d')!point symmetry -42m
        IVectors = 0
      CASE('e')!point symmetry mm
        IVectors = 1
      CASE('f')!point symmetry mm
        IVectors = 1
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry m
        IVectors = 2
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 115, P -4 m 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(116)
    CASE(116)!P -4 c 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry -4
        IVectors = 0
      CASE('d')!point symmetry -4
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 116, P -4 c 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(117)
    CASE(117)!P -4 b 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -4
        IVectors = 0
      CASE('b')!point symmetry -4
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 117, P -4 b 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(118)
    CASE(118)!P -4 n 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -4
        IVectors = 0
      CASE('b')!point symmetry -4
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 118, P -4 n 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(119)
    CASE(119)!I -4 m 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -42m
        IVectors = 0
      CASE('b')!point symmetry -42m
        IVectors = 0
      CASE('c')!point symmetry -42m
        IVectors = 0
      CASE('d')!point symmetry -42m
        IVectors = 0
      CASE('e')!point symmetry mm
        IVectors = 1
      CASE('f')!point symmetry mm
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry m
        IVectors = 2
      CASE('j')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 119, I -4 m 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(120)
    CASE(120)!I -4 c 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry -4
        IVectors = 0
      CASE('c')!point symmetry -4
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 120, I -4 c 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(121)
    CASE(121)!I -4 2 m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -42m
        IVectors = 0
      CASE('b')!point symmetry -42m
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry -4
        IVectors = 0
      CASE('e')!point symmetry mm
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry m
        IVectors = 2
      CASE('j')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 121, I -4 2 m, not recognised")) RETURN     
      END SELECT
!!$  CASE(122)
    CASE(122)!I -4 2 d
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -4
        IVectors = 0
      CASE('b')!point symmetry -4
        IVectors = 0
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 122, I -4 2 d, not recognised")) RETURN     
      END SELECT
!!$  CASE(123)
    CASE(123)!P 4/m 2/m 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 4/mmm
        IVectors = 0
      CASE('b')!point symmetry 4/mmm
        IVectors = 0
      CASE('c')!point symmetry 4/mmm
        IVectors = 0
      CASE('d')!point symmetry 4/mmm
        IVectors = 0
      CASE('e')!point symmetry mmm
        IVectors = 0
      CASE('f')!point symmetry mmm
        IVectors = 0
      CASE('g')!point symmetry 4mm
        IVectors = 1
      CASE('h')!point symmetry 4mm
        IVectors = 1
      CASE('i')!point symmetry mm
        IVectors = 1
      CASE('j')!point symmetry mm
        IVectors = 1
      CASE('k')!point symmetry mm
        IVectors = 1
      CASE('l')!point symmetry mm
        IVectors = 1
      CASE('m')!point symmetry mm
        IVectors = 1
      CASE('n')!point symmetry mm
        IVectors = 1
      CASE('o')!point symmetry mm
        IVectors = 1
      CASE('p')!point symmetry m
        IVectors = 2
      CASE('q')!point symmetry m
        IVectors = 2
      CASE('r')!point symmetry m
        IVectors = 2
      CASE('s')!point symmetry m
        IVectors = 2
      CASE('t')!point symmetry m
        IVectors = 2
      CASE('u')!point symmetry 1
        IVectors = 3        
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 123, P 4/m 2/m 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(124)
    CASE(124)!P 4/m 2/c 2/c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 42
        IVectors = 0
      CASE('b')!point symmetry 4/m
        IVectors = 0
      CASE('c')!point symmetry 42
        IVectors = 0
      CASE('d')!point symmetry 4/m
        IVectors = 0
      CASE('e')!point symmetry 2/m
        IVectors = 0
      CASE('f')!point symmetry 222
        IVectors = 0
      CASE('g')!point symmetry 4
        IVectors = 1
      CASE('h')!point symmetry 4
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 2
        IVectors = 1
      CASE('m')!point symmetry m
        IVectors = 2
      CASE('n')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 124, P 4/m 2/c 2/c, not recognised")) RETURN     
      END SELECT
!!$  CASE(125)
    CASE(125)!P 4/n 2/b 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 42
        IVectors = 0
      CASE('b')!point symmetry 42
        IVectors = 0
      CASE('c')!point symmetry -42m
        IVectors = 0
      CASE('d')!point symmetry -42m
        IVectors = 0
      CASE('e')!point symmetry 2/m
        IVectors = 0
      CASE('f')!point symmetry 2/m
        IVectors = 0
      CASE('g')!point symmetry 4
        IVectors = 1
      CASE('h')!point symmetry mm
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 2
        IVectors = 1
      CASE('m')!point symmetry m
        IVectors = 2
      CASE('n')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 125, P 4/n 2/b 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(126)
    CASE(126)!P 4/n 2/n 2/c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 42
        IVectors = 0
      CASE('b')!point symmetry 42
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry -4
        IVectors = 0
      CASE('e')!point symmetry 4
        IVectors = 1
      CASE('f')!point symmetry -1
        IVectors = 0
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 126, P 4/n 2/n 2/c, not recognised")) RETURN     
      END SELECT
!!$  CASE(127)
    CASE(127)!P 4/m 21/b 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 4/m
        IVectors = 0
      CASE('b')!point symmetry 4/m
        IVectors = 0
      CASE('c')!point symmetry mmm
        IVectors = 0
      CASE('d')!point symmetry mmm
        IVectors = 0
      CASE('e')!point symmetry 4
        IVectors = 1
      CASE('f')!point symmetry mm
        IVectors = 1
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry mm
        IVectors = 1
      CASE('i')!point symmetry m
        IVectors = 2
      CASE('j')!point symmetry m
        IVectors = 2
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 127, P 4/m 21/b 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(128)
    CASE(128)!P 4/m 21/n 2/c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 4/m
        IVectors = 0
      CASE('b')!point symmetry 4/m
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 4
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry m
        IVectors = 2
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 128, P 4/m 21/n 2/c, not recognised")) RETURN     
      END SELECT
!!$  CASE(129)
    CASE(129)!P 4/n 21/m 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -42m
        IVectors = 0
      CASE('b')!point symmetry -42m
        IVectors = 0
      CASE('c')!point symmetry 4mm
        IVectors = 1
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry 2/m
        IVectors = 0
      CASE('f')!point symmetry mm
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry m
        IVectors = 2
      CASE('j')!point symmetry m
        IVectors = 2
      CASE('k')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 129, P 4/n 21/m 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(130)
    CASE(130)!P 4/n 21/c 2/c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry -4
        IVectors = 0
      CASE('c')!point symmetry 4
        IVectors = 1
      CASE('d')!point symmetry -1
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 130, P 4/n 21/c 2/c, not recognised")) RETURN     
      END SELECT
!!$  CASE(131)
    CASE(131)!P 42/m 2/m 2/c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mmm
        IVectors = 0
      CASE('b')!point symmetry mmm
        IVectors = 0
      CASE('c')!point symmetry mmm
        IVectors = 0
      CASE('d')!point symmetry mmm
        IVectors = 0
      CASE('e')!point symmetry -42m
        IVectors = 1
      CASE('f')!point symmetry -42m
        IVectors = 0
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry mm
        IVectors = 1
      CASE('i')!point symmetry mm
        IVectors = 1
      CASE('j')!point symmetry mm
        IVectors = 1
      CASE('k')!point symmetry mm
        IVectors = 1
      CASE('l')!point symmetry mm
        IVectors = 1
      CASE('m')!point symmetry mm
        IVectors = 1
      CASE('n')!point symmetry 2
        IVectors = 1
      CASE('o')!point symmetry m
        IVectors = 2
      CASE('p')!point symmetry m
        IVectors = 2
      CASE('q')!point symmetry m
        IVectors = 2
      CASE('r')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 131, P 42/m 2/m 2/c, not recognised")) RETURN     
      END SELECT
!!$  CASE(132)
    CASE(132)!P 42/m 2/c 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mmm
        IVectors = 0
      CASE('b')!point symmetry -42m
        IVectors = 0
      CASE('c')!point symmetry mmm
        IVectors = 0
      CASE('d')!point symmetry -42m
        IVectors = 0
      CASE('e')!point symmetry 222
        IVectors = 0
      CASE('f')!point symmetry 2/m
        IVectors = 0
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry mm
        IVectors = 1
      CASE('i')!point symmetry mm
        IVectors = 1
      CASE('j')!point symmetry mm
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 2
        IVectors = 1
      CASE('m')!point symmetry 2
        IVectors = 1
      CASE('n')!point symmetry m
        IVectors = 2
      CASE('o')!point symmetry m
        IVectors = 2
      CASE('p')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 132, P 42/m 2/c 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(133)
    CASE(133)!P 42/n 2/b 2/c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry -4
        IVectors = 0
      CASE('e')!point symmetry -1
        IVectors = 0
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 133, P 42/n 2/b 2/c, not recognised")) RETURN     
      END SELECT
!!$  CASE(134)
    CASE(134)!P 42/n 2/n 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -42m
        IVectors = 0
      CASE('b')!point symmetry -42m
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 2/m
        IVectors = 0
      CASE('f')!point symmetry 2/m
        IVectors = 0
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 2
        IVectors = 1
      CASE('m')!point symmetry m
        IVectors = 2
      CASE('n')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 134, P 42/n 2/n 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(135)
    CASE(135)!P 42/m 21/b 2/c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m
        IVectors = 0
      CASE('b')!point symmetry -4
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry m
        IVectors = 2
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 135, P 42/m 21/b 2/c, not recognised")) RETURN     
      END SELECT
!!$  CASE(136)
    CASE(136)!P 42/m 21/n 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry mmm
        IVectors = 0
      CASE('b')!point symmetry mmm
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry -4
        IVectors = 0
      CASE('e')!point symmetry mm
        IVectors = 1
      CASE('f')!point symmetry mm
        IVectors = 1
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry m
        IVectors = 2
      CASE('j')!point symmetry m
        IVectors = 2
      CASE('k')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 136, P 42/m 21/n 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(137)
    CASE(137)!P 42/n 21/m 2/c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -42m
        IVectors = 0
      CASE('b')!point symmetry -42m
        IVectors = 0
      CASE('c')!point symmetry mm
        IVectors = 1
      CASE('d')!point symmetry mm
        IVectors = 1
      CASE('e')!point symmetry -1
        IVectors = 0
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry m
        IVectors = 2
      CASE('h')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 137, P 42/n 21/m 2/c, not recognised")) RETURN     
      END SELECT
!!$  CASE(138)
    CASE(138)!P 42/n 21/c 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry -4
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry mm
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry m
        IVectors = 2
      CASE('j')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 138, P 42/n 21/c 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(139)
    CASE(139)!I 4/m 2/m 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 4/mmm
        IVectors = 0
      CASE('b')!point symmetry 4/mmm
        IVectors = 0
      CASE('c')!point symmetry mmm
        IVectors = 0
      CASE('d')!point symmetry -4m2
        IVectors = 0
      CASE('e')!point symmetry 4mm
        IVectors = 0
      CASE('f')!point symmetry 2/m
        IVectors = 0
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry mm
        IVectors = 1
      CASE('i')!point symmetry mm
        IVectors = 1
      CASE('j')!point symmetry mm
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry m
        IVectors = 2
      CASE('m')!point symmetry m
        IVectors = 2
      CASE('n')!point symmetry m
        IVectors = 2
      CASE('o')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 139, I4/m m m, not recognised")) RETURN     
      END SELECT    
!!$  CASE(140)
    CASE(140)!I 4/m 2/c 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 42
        IVectors = 0
      CASE('b')!point symmetry -42m
        IVectors = 0
      CASE('c')!point symmetry 4/m
        IVectors = 1
      CASE('d')!point symmetry mmm
        IVectors = 0
      CASE('e')!point symmetry 2/m
        IVectors = 0
      CASE('f')!point symmetry 4
        IVectors = 1
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry mm
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry m
        IVectors = 2
      CASE('m')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 140, I 4/m 2/c 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(141)
    CASE(141)!I 41/a 2/m 2/d
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -42m
        IVectors = 0
      CASE('b')!point symmetry -42m
        IVectors = 0
      CASE('c')!point symmetry 2/m
        IVectors = 0
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry mm
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry m
        IVectors = 2
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 141, I 41/a 2/m 2/d, not recognised")) RETURN     
      END SELECT
!!$  CASE(142)
    CASE(142)!I 41/a 2/c 2/d
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -4
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry -1
        IVectors = 0
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 142, I 41/a 2/c 2/d, not recognised")) RETURN     
      END SELECT
!!$  CASE(143)
    CASE(143)!P 3
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 3
        IVectors = 1
      CASE('b')!point symmetry 3
        IVectors = 1
      CASE('c')!point symmetry 3
        IVectors = 1
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 143, P 3, not recognised")) RETURN     
      END SELECT
!!$  CASE(144)
    CASE(144)!P 31
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 144, P 31, not recognised")) RETURN     
      END SELECT
!!$  CASE(145)
    CASE(145)!P 32
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 145, P 32, not recognised")) RETURN     
      END SELECT
!!$  CASE(146)
    CASE(146)!R 3
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 3
        IVectors = 1
      CASE('b')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 146, R 3, not recognised")) RETURN     
      END SELECT
!!$  CASE(147)
    CASE(147)!P -3
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -3
        IVectors = 0
      CASE('b')!point symmetry -3
        IVectors = 0
      CASE('c')!point symmetry 3
        IVectors = 1
      CASE('d')!point symmetry 3
        IVectors = 1
      CASE('e')!point symmetry -1
        IVectors = 0
      CASE('f')!point symmetry -1
        IVectors = 0
      CASE('g')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 147, P -3, not recognised")) RETURN     
      END SELECT
!!$  CASE(148)
    CASE(148)!R -3
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -3
        IVectors = 0
      CASE('b')!point symmetry -3
        IVectors = 0
      CASE('c')!point symmetry 3
        IVectors = 1
      CASE('d')!point symmetry -1
        IVectors = 0
      CASE('e')!point symmetry -1
        IVectors = 0
      CASE('f')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 148, R -3, not recognised")) RETURN     
      END SELECT
!!$  CASE(149)
    CASE(149)!P 3 1 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 32
        IVectors = 0
      CASE('b')!point symmetry 32
        IVectors = 0
      CASE('c')!point symmetry 32
        IVectors = 0
      CASE('d')!point symmetry 32
        IVectors = 0
      CASE('e')!point symmetry 32
        IVectors = 0
      CASE('f')!point symmetry 32
        IVectors = 0
      CASE('g')!point symmetry 3
        IVectors = 1
      CASE('h')!point symmetry 3
        IVectors = 1
      CASE('i')!point symmetry 3
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 149, P 3 1 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(150)
    CASE(150)!P 3 2 1
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 32
        IVectors = 0
      CASE('b')!point symmetry 32
        IVectors = 0
      CASE('c')!point symmetry 3
        IVectors = 1
      CASE('d')!point symmetry 3
        IVectors = 1
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 150, P 3 2 1, not recognised")) RETURN     
      END SELECT
!!$  CASE(151)
    CASE(151)!P 31 1 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 151, P 31 1 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(152)
    CASE(152)!P 31 2 1
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 152, P 31 2 1, not recognised")) RETURN     
      END SELECT
!!$  CASE(153)
    CASE(153)!P 32 1 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 153, P 32 1 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(154)
    CASE(154)!P 32 2 1
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 154, P 32 2 1, not recognised")) RETURN     
      END SELECT
!!$  CASE(155)
    CASE(155)!R 3 2 
      SELECT CASE (SWyckoff)
        CASE('a')!point symmetry 32
          IVectors = 0
        CASE('b')!point symmetry 32
          IVectors = 0
        CASE('c')!point symmetry 3
          IVectors = 1
        CASE('d')!point symmetry 2
          IVectors = 1
        CASE('e')!point symmetry 2
          IVectors = 1
        CASE('f')!point symmetry 1
          IVectors = 3
        CASE DEFAULT
          IErr = 1
          IF(l_alert(IErr,"DetermineAllowedMovements",&
                "Wyckoff Symbol for space group 155, R 3 2, not recognised")) RETURN  
        END SELECT 
!!$  CASE(156)
    CASE(156)!P 3 m 1
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 3m
        IVectors = 1
      CASE('b')!point symmetry 3m
        IVectors = 1
      CASE('c')!point symmetry 3m
        IVectors = 1
      CASE('d')!point symmetry m
        IVectors = 2
      CASE('e')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 156, P 3 m 1, not recognised")) RETURN     
      END SELECT
!!$  CASE(157)
    CASE(157)!P 3 1 m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 3m
        IVectors = 1
      CASE('b')!point symmetry 3
        IVectors = 1
      CASE('c')!point symmetry m
        IVectors = 2
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 157, P 3 1 m, not recognised")) RETURN     
      END SELECT
!!$  CASE(158)
    CASE(158)!P 3 c 1
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 3
        IVectors = 1
      CASE('b')!point symmetry 3
        IVectors = 1
      CASE('c')!point symmetry 3
        IVectors = 1
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 158, P 3 c 1, not recognised")) RETURN     
      END SELECT
!!$  CASE(159)
    CASE(159)!P 3 1 c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 3
        IVectors = 1
      CASE('b')!point symmetry 3
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 159, P 3 1 c, not recognised")) RETURN     
      END SELECT
!!$  CASE(160)
    CASE(160)!R 3 m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 3m
        IVectors = 1
      CASE('b')!point symmetry m
        IVectors = 2
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 160, R 3 m, not recognised")) RETURN     
      END SELECT
!!$  CASE(161)
    CASE(161)!R 3 c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 3
        IVectors = 1
      CASE('b')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 161, R 3 c, not recognised")) RETURN     
      END SELECT
!!$  CASE(162)
    CASE(162)!P -3 1 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -3m
        IVectors = 0
      CASE('b')!point symmetry -3m
        IVectors = 0
      CASE('c')!point symmetry 32
        IVectors = 0
      CASE('d')!point symmetry 32
        IVectors = 0
      CASE('e')!point symmetry 3m
        IVectors = 1
      CASE('f')!point symmetry 2/m
        IVectors = 0
      CASE('g')!point symmetry 2/m
        IVectors = 0
      CASE('h')!point symmetry 3
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 162, P -3 1 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(163)
    CASE(163)!P -3 1 2/c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 32
        IVectors = 0
      CASE('b')!point symmetry -3
        IVectors = 0
      CASE('c')!point symmetry 32
        IVectors = 0
      CASE('d')!point symmetry 32
        IVectors = 0
      CASE('e')!point symmetry 3
        IVectors = 1
      CASE('f')!point symmetry 3
        IVectors = 1
      CASE('g')!point symmetry -1
        IVectors = 0
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 163, P -3 1 2/c, not recognised")) RETURN     
      END SELECT
!!$  CASE(164)
    CASE(164)!P -3 2/m 1
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -3m
        IVectors = 0
      CASE('b')!point symmetry -3m
        IVectors = 0
      CASE('c')!point symmetry 3m
        IVectors = 1
      CASE('d')!point symmetry 3m
        IVectors = 1
      CASE('e')!point symmetry 2/m
        IVectors = 0
      CASE('f')!point symmetry 2/m
        IVectors = 0
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry m
        IVectors = 2
      CASE('j')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 164, P -3 2/m 1, not recognised")) RETURN     
      END SELECT
!!$  CASE(165)
    CASE(165)!P -3 2/c 1
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 32
        IVectors = 0
      CASE('b')!point symmetry -3
        IVectors = 0
      CASE('c')!point symmetry 3
        IVectors = 1
      CASE('d')!point symmetry 3
        IVectors = 1
      CASE('e')!point symmetry -1
        IVectors = 0
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 165, P -3 2/c 1, not recognised")) RETURN     
      END SELECT
!!$  CASE(166)
    CASE(166)!R -3 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -3m
        IVectors = 0
      CASE('b')!point symmetry -3m
        IVectors = 0
      CASE('c')!point symmetry 3m
        IVectors = 1
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry 2/m
        IVectors = 0
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry m
        IVectors = 2
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 166, R -3 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(167)
    CASE(167)!R -3 2/c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 32
        IVectors = 0
      CASE('b')!point symmetry -3
        IVectors = 0
      CASE('c')!point symmetry 3
        IVectors = 1
      CASE('d')!point symmetry -1
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 167, R -3 2/c, not recognised")) RETURN     
      END SELECT
!!$  CASE(168)
    CASE(168)!P 6
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 6
        IVectors = 1
      CASE('b')!point symmetry 3
        IVectors = 1
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 168, P 6 not recognised")) RETURN     
      END SELECT
!!$  CASE(169)
    CASE(169)!P 61
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 169, P 61, not recognised")) RETURN     
      END SELECT
!!$  CASE(170)
    CASE(170)!P 65
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 170, P 65, not recognised")) RETURN     
      END SELECT
!!$  CASE(171)
    CASE(171)!P 62
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 171, P 62, not recognised")) RETURN     
      END SELECT
!!$  CASE(172)
    CASE(172)!P 64
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 172, P 64, not recognised")) RETURN     
      END SELECT
!!$  CASE(173)
    CASE(173)!P 63
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 3
        IVectors = 1
      CASE('b')!point symmetry 3
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 173, P 63, not recognised")) RETURN     
      END SELECT
!!$  CASE(174)
    CASE(174)!P -6
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -6
        IVectors = 0
      CASE('b')!point symmetry -6
        IVectors = 0
      CASE('c')!point symmetry -6
        IVectors = 0
      CASE('d')!point symmetry -6
        IVectors = 0
      CASE('e')!point symmetry -6
        IVectors = 0
      CASE('f')!point symmetry -6
        IVectors = 0
      CASE('g')!point symmetry 3
        IVectors = 1
      CASE('h')!point symmetry 3
        IVectors = 1
      CASE('i')!point symmetry 3
        IVectors = 1
      CASE('j')!point symmetry m
        IVectors = 2
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 174, P -6, not recognised")) RETURN     
      END SELECT
!!$  CASE(175)
    CASE(175)!P 6/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 6/m
        IVectors = 0
      CASE('b')!point symmetry 6/m
        IVectors = 0
      CASE('c')!point symmetry -6
        IVectors = 0
      CASE('d')!point symmetry -6
        IVectors = 0
      CASE('e')!point symmetry 6
        IVectors = 1
      CASE('f')!point symmetry 2/m
        IVectors = 0
      CASE('g')!point symmetry 2/m
        IVectors = 0
      CASE('h')!point symmetry 3
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry m
        IVectors = 2
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 175, P 6/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(176)
    CASE(176)!P 63/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -6
        IVectors = 0
      CASE('b')!point symmetry -3
        IVectors = 0
      CASE('c')!point symmetry -6
        IVectors = 0
      CASE('d')!point symmetry -6
        IVectors = 0
      CASE('e')!point symmetry 3
        IVectors = 1
      CASE('f')!point symmetry 3
        IVectors = 1
      CASE('g')!point symmetry -1
        IVectors = 0
      CASE('h')!point symmetry m
        IVectors = 2
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 176, P 63/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(177)
    CASE(177)!P 6 2 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 622
        IVectors = 0
      CASE('b')!point symmetry 622
        IVectors = 0
      CASE('c')!point symmetry 32
        IVectors = 0
      CASE('d')!point symmetry 32
        IVectors = 0
      CASE('e')!point symmetry 6
        IVectors = 1
      CASE('f')!point symmetry 222
        IVectors = 0
      CASE('g')!point symmetry 222
        IVectors = 0
      CASE('h')!point symmetry 3
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 2
        IVectors = 1
      CASE('m')!point symmetry 2
        IVectors = 1
      CASE('n')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 177, P 6 2 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(178)
    CASE(178)!P 61 2 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 178, P 61 2 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(179)
    CASE(179)!P 65 2 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 179, P 65 2 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(180)
    CASE(180)!P 62 2 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 180, P 62 2 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(181)
    CASE(181)!P 64 2 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 181, P 64 2 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(182)
    CASE(182)!P 63 2 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 32
        IVectors = 0
      CASE('b')!point symmetry 32
        IVectors = 0
      CASE('c')!point symmetry 32
        IVectors = 0
      CASE('d')!point symmetry 32
        IVectors = 0
      CASE('e')!point symmetry 3
        IVectors = 1
      CASE('f')!point symmetry 3
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 182, P 63 2 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(183)
    CASE(183)!P 6 m m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 6mm
        IVectors = 1
      CASE('b')!point symmetry 3m
        IVectors = 1
      CASE('c')!point symmetry mm
        IVectors = 1
      CASE('d')!point symmetry m
        IVectors = 2
      CASE('e')!point symmetry m
        IVectors = 2
      CASE('f')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 183, P 6 m m, not recognised")) RETURN     
      END SELECT
!!$  CASE(184)
    CASE(184)!P 6 c c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 6
        IVectors = 1
      CASE('b')!point symmetry 3
        IVectors = 1
      CASE('c')!point symmetry 2
        IVectors = 1
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 184, P 6 c c, not recognised")) RETURN     
      END SELECT
!!$  CASE(185)
    CASE(185)!P 63 c m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 3m
        IVectors = 1
      CASE('b')!point symmetry 3
        IVectors = 1
      CASE('c')!point symmetry m
        IVectors = 2
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 185, P 63 c m, not recognised")) RETURN     
      END SELECT
!!$  CASE(186)
    CASE(186)!P 63 m c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 3m
        IVectors = 1
      CASE('b')!point symmetry 3m
        IVectors = 1
      CASE('c')!point symmetry m
        IVectors = 2
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 186, 63 m c, not recognised")) RETURN     
      END SELECT
!!$  CASE(187)
    CASE(187)!P -6 m 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -6m2
        IVectors = 0
      CASE('b')!point symmetry -6m2
        IVectors = 0
      CASE('c')!point symmetry -6m2
        IVectors = 0
      CASE('d')!point symmetry -6m2
        IVectors = 0
      CASE('e')!point symmetry -6m2
        IVectors = 0
      CASE('f')!point symmetry -6m2
        IVectors = 0
      CASE('g')!point symmetry 3m
        IVectors = 1
      CASE('h')!point symmetry 3m
        IVectors = 1
      CASE('i')!point symmetry 3m
        IVectors = 1
      CASE('j')!point symmetry mm
        IVectors = 1
      CASE('k')!point symmetry mm
        IVectors = 1
      CASE('l')!point symmetry m
        IVectors = 2
      CASE('m')!point symmetry m
        IVectors = 2
      CASE('n')!point symmetry m
        IVectors = 2
      CASE('o')!point symmetry 1
        IVectors = 3        
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 187, P -6 m 2, not recognised")) RETURN     
      END SELECT
!!$  CASE(188)
    CASE(188)!P -6 c 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 32
        IVectors = 0
      CASE('b')!point symmetry -6
        IVectors = 0
      CASE('c')!point symmetry 32
        IVectors = 0
      CASE('d')!point symmetry -6
        IVectors = 0
      CASE('e')!point symmetry 32
        IVectors = 0
      CASE('f')!point symmetry -6
        IVectors = 0
      CASE('g')!point symmetry 3
        IVectors = 1
      CASE('h')!point symmetry 3
        IVectors = 1
      CASE('i')!point symmetry 3
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 188, P -6 c 2, not recognised")) RETURN    
      END SELECT
!!$  CASE(189)
    CASE(189)!P -6 2 m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -6m2
        IVectors = 0
      CASE('b')!point symmetry -6m2
        IVectors = 0
      CASE('c')!point symmetry -6
        IVectors = 0
      CASE('d')!point symmetry -6
        IVectors = 0
      CASE('e')!point symmetry 3m
        IVectors = 1
      CASE('f')!point symmetry mm
        IVectors = 1
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry 3
        IVectors = 1
      CASE('i')!point symmetry m
        IVectors = 2
      CASE('j')!point symmetry m
        IVectors = 2
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry 1
        IVectors = 3        
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 189, P -6 2 m, not recognised")) RETURN     
      END SELECT
!!$  CASE(190)
    CASE(190)!P -6 2 c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 32
        IVectors = 0
      CASE('b')!point symmetry -6
        IVectors = 0
      CASE('c')!point symmetry -6
        IVectors = 0
      CASE('d')!point symmetry -6
        IVectors = 0
      CASE('e')!point symmetry 3
        IVectors = 1
      CASE('f')!point symmetry 3
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry m
        IVectors = 2
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 190, P -6 2 c, not recognised")) RETURN     
      END SELECT
!!$  CASE(191)
    CASE(191)!P 6/m 2/m 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 6/mmm
        IVectors = 0
      CASE('b')!point symmetry 6/mmm
        IVectors = 0
      CASE('c')!point symmetry -6m2
        IVectors = 0
      CASE('d')!point symmetry -6m2
        IVectors = 0
      CASE('e')!point symmetry 6mm
        IVectors = 1
      CASE('f')!point symmetry mmm
        IVectors = 0
      CASE('g')!point symmetry mmm
        IVectors = 0
      CASE('h')!point symmetry 3m
        IVectors = 1
      CASE('i')!point symmetry mm
        IVectors = 1
      CASE('j')!point symmetry mm
        IVectors = 1
      CASE('k')!point symmetry mm
        IVectors = 1
      CASE('l')!point symmetry mm
        IVectors = 1
      CASE('m')!point symmetry mm
        IVectors = 1
      CASE('n')!point symmetry m
        IVectors = 2
      CASE('o')!point symmetry m
        IVectors = 2
      CASE('p')!point symmetry m
        IVectors = 2
      CASE('q')!point symmetry m
        IVectors = 2
      CASE('r')!point symmetry 1
        IVectors = 3            
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 191, P 6/m 2/m 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(192)
    CASE(192)!P 6/m 2/c 2/c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 62
        IVectors = 0
      CASE('b')!point symmetry 6/m
        IVectors = 0
      CASE('c')!point symmetry 32
        IVectors = 0
      CASE('d')!point symmetry -6
        IVectors = 0
      CASE('e')!point symmetry 6
        IVectors = 1
      CASE('f')!point symmetry 222
        IVectors = 0
      CASE('g')!point symmetry 2/m
        IVectors = 0
      CASE('h')!point symmetry 3
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry m
        IVectors = 2
      CASE('m')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 192, P 6/m 2/c 2/c, not recognised")) RETURN     
      END SELECT
!!$  CASE(193)
    CASE(193)!P 63/m 2/c 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -6m2
        IVectors = 0
      CASE('b')!point symmetry -3m
        IVectors = 0
      CASE('c')!point symmetry -6
        IVectors = 0
      CASE('d')!point symmetry 32
        IVectors = 0
      CASE('e')!point symmetry 6
        IVectors = 0
      CASE('f')!point symmetry 2/m
        IVectors = 0
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry 3
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry m
        IVectors = 2
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 193, P 63/m 2/c 2/m, not recognised")) RETURN     
      END SELECT
!!$  CASE(194)
    CASE(194)!P 63/m 2/m 2/c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -3m
        IVectors = 0
      CASE('b')!point symmetry -6m2
        IVectors = 0
      CASE('c')!point symmetry -6m2
        IVectors = 0
      CASE('d')!point symmetry -6m2
        IVectors = 0
      CASE('e')!point symmetry 3m
        IVectors = 1
      CASE('f')!point symmetry 3m
        IVectors = 1
      CASE('g')!point symmetry 2/m
        IVectors = 0
      CASE('h')!point symmetry mm
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry m
        IVectors = 2
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry 1
        IVectors = 3        
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 194, P 63/m 2/m 2/c, not recognised")) RETURN     
      END SELECT
!!$  CASE(195)
    CASE(195)!P 2 3
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 23
        IVectors = 0
      CASE('b')!point symmetry 23
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 3
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 1
        IVectors = 3            
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 195, P 2 3, not recognised")) RETURN     
      END SELECT
!!$  CASE(196)
    CASE(196)!F 2 3
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 23
        IVectors = 0
      CASE('b')!point symmetry 23
        IVectors = 0
      CASE('c')!point symmetry 23
        IVectors = 0
      CASE('d')!point symmetry 23
        IVectors = 0
      CASE('e')!point symmetry 3
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 1
        IVectors = 3        
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 196, F 2 3, not recognised")) RETURN     
      END SELECT
!!$  CASE(197)
    CASE(197)!I 2 3
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 23
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry 3
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 2
        IVectors = 1
      CASE('f')!point symmetry 1
        IVectors = 3    
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 197, I 2 3, not recognised")) RETURN     
      END SELECT
!!$  CASE(198)
    CASE(198)!P 21 3
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 3
        IVectors = 1
      CASE('b')!point symmetry 1
        IVectors = 3            
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 198, P 21 3, not recognised")) RETURN     
      END SELECT
!!$  CASE(199)
    CASE(199)!I 21 3
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 3
        IVectors = 1
      CASE('b')!point symmetry 2
        IVectors = 1
      CASE('c')!point symmetry 1
        IVectors = 3            
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 199, I 21 3, not recognised")) RETURN     
      END SELECT
!!$  CASE(200)
    CASE(200)!P 2/m -3
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry m3
        IVectors = 0
      CASE('b')!point symmetry m3
        IVectors = 0
      CASE('c')!point symmetry mmm
        IVectors = 0
      CASE('d')!point symmetry mmm
        IVectors = 0
      CASE('e')!point symmetry mm
        IVectors = 1
      CASE('f')!point symmetry mm
        IVectors = 1
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry mm
        IVectors = 1
      CASE('i')!point symmetry 3
        IVectors = 1
      CASE('j')!point symmetry m
        IVectors = 2
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry 1
        IVectors = 3            
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 200, P 2/m -3, not recognised")) RETURN     
      END SELECT
!!$  CASE(201)
    CASE(201)!P 2/n -3
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 23
        IVectors = 0
      CASE('b')!point symmetry -3
        IVectors = 0
      CASE('c')!point symmetry -3
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 3
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 1
        IVectors = 3            
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 201, P 2/n -3, not recognised")) RETURN     
      END SELECT
!!$  CASE(202)
    CASE(202)!F 2/m -3
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry m3
        IVectors = 0
      CASE('b')!point symmetry m3
        IVectors = 0
      CASE('c')!point symmetry 23
        IVectors = 0
      CASE('d')!point symmetry 2/m
        IVectors = 0
      CASE('e')!point symmetry mm
        IVectors = 1
      CASE('f')!point symmetry 3
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry m
        IVectors = 2
      CASE('i')!point symmetry 1
        IVectors = 3            
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 202, F 2/m -3, not recognised")) RETURN     
      END SELECT
!!$  CASE(203)
    CASE(203)!F 2/d -3
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 23
        IVectors = 0
      CASE('b')!point symmetry 23
        IVectors = 0
      CASE('c')!point symmetry -3
        IVectors = 0
      CASE('d')!point symmetry -3
        IVectors = 0
      CASE('e')!point symmetry 3
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 1
        IVectors = 3            
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 139, I4/m m m, not recognised")) RETURN     
      END SELECT
!!$  CASE(204)
    CASE(204)!I 2/m -3
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry m3
        IVectors = 0
      CASE('b')!point symmetry mmm
        IVectors = 0
      CASE('c')!point symmetry -3
        IVectors = 0
      CASE('d')!point symmetry mm
        IVectors = 1
      CASE('e')!point symmetry mm
        IVectors = 1
      CASE('f')!point symmetry 3
        IVectors = 1
      CASE('g')!point symmetry m
        IVectors = 2
      CASE('h')!point symmetry 1
        IVectors = 3        
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 139, I4/m m m, not recognised")) RETURN     
      END SELECT
!!$  CASE(205)
    CASE(205)!P 21/a -3
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -3
        IVectors = 0
      CASE('b')!point symmetry -3
        IVectors = 0
      CASE('c')!point symmetry 3
        IVectors = 1
      CASE('d')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 139, I4/m m m, not recognised")) RETURN     
      END SELECT
!!$  CASE(206)
    CASE(206)!I 21/a -3
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -3
        IVectors = 0
      CASE('b')!point symmetry -3
        IVectors = 0
      CASE('c')!point symmetry 3
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 1
        IVectors = 3            
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 139, I4/m m m, not recognised")) RETURN     
      END SELECT
!!$  CASE(207)
    CASE(207)!P 4 3 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 43
        IVectors = 0
      CASE('b')!point symmetry 43
        IVectors = 0
      CASE('c')!point symmetry 42
        IVectors = 0
      CASE('d')!point symmetry 42
        IVectors = 0
      CASE('e')!point symmetry 4
        IVectors = 1
      CASE('f')!point symmetry 4
        IVectors = 1
      CASE('g')!point symmetry 3
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 139, I4/m m m, not recognised")) RETURN     
      END SELECT
!!$  CASE(208)
    CASE(208)!P 42 3 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 23
        IVectors = 0
      CASE('b')!point symmetry 32
        IVectors = 0
      CASE('c')!point symmetry 32
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 222
        IVectors = 0
      CASE('f')!point symmetry 222
        IVectors = 0
      CASE('g')!point symmetry 3
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry 2
        IVectors = 1
      CASE('l')!point symmetry 2
        IVectors = 1
      CASE('m')!point symmetry 1
        IVectors = 3        
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 139, I4/m m m, not recognised")) RETURN     
      END SELECT
!!$  CASE(209)
    CASE(209)!F 4 3 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 43
        IVectors = 0
      CASE('b')!point symmetry 43
        IVectors = 0
      CASE('c')!point symmetry 23
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 4
        IVectors = 1
      CASE('f')!point symmetry 3
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 1
        IVectors = 3            
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 139, I4/m m m, not recognised")) RETURN     
      END SELECT
!!$  CASE(210)
    CASE(210)!F 41 3 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 23
        IVectors = 0
      CASE('b')!point symmetry 23
        IVectors = 0
      CASE('c')!point symmetry 32
        IVectors = 0
      CASE('d')!point symmetry 32
        IVectors = 0
      CASE('e')!point symmetry 3
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 1
        IVectors = 3        
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 139, I4/m m m, not recognised")) RETURN     
      END SELECT
!!$  CASE(211)
    CASE(211)!I 4 3 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 43
        IVectors = 0
      CASE('b')!point symmetry 42
        IVectors = 0
      CASE('c')!point symmetry 32
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 4
        IVectors = 1
      CASE('f')!point symmetry 3
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 1
        IVectors = 3        
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 139, I4/m m m, not recognised")) RETURN     
      END SELECT
!!$  CASE(212)
    CASE(212)!P 43 3 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 32
        IVectors = 0
      CASE('b')!point symmetry 32
        IVectors = 0
      CASE('c')!point symmetry 3
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 139, I4/m m m, not recognised")) RETURN     
      END SELECT
!!$  CASE(213)
    CASE(213)!P 41 3 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 32
        IVectors = 0
      CASE('b')!point symmetry 32
        IVectors = 0
      CASE('c')!point symmetry 3
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 1
        IVectors = 3    
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 139, I4/m m m, not recognised")) RETURN     
      END SELECT
!!$  CASE(214)
    CASE(214)!I 41 3 2
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 32
        IVectors = 0
      CASE('b')!point symmetry 32
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry 222
        IVectors = 0
      CASE('e')!point symmetry 3
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 1
        IVectors = 3            
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 139, I4/m m m, not recognised")) RETURN     
      END SELECT
!!$  CASE(215)
    CASE(215)!P -4 3 m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -43m
        IVectors = 0
      CASE('b')!point symmetry -43m
        IVectors = 0
      CASE('c')!point symmetry -42m
        IVectors = 0
      CASE('d')!point symmetry -42m
        IVectors = 0
      CASE('e')!point symmetry 3m
        IVectors = 1
      CASE('f')!point symmetry mm
        IVectors = 1
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry m
        IVectors = 2
      CASE('j')!point symmetry 1
        IVectors = 3            
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 215, P -4 3 m, not recognised")) RETURN     
      END SELECT
!!$  CASE(216)
    CASE(216)!F -4 3 m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -43m
        IVectors = 0
      CASE('b')!point symmetry -43m
        IVectors = 0
      CASE('c')!point symmetry -43m
        IVectors = 0
      CASE('d')!point symmetry -43m
        IVectors = 0
      CASE('e')!point symmetry 3m
        IVectors = 1
      CASE('f')!point symmetry mm
        IVectors = 1
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry m
        IVectors = 2
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements","space group not recognised")) RETURN
      END SELECT
!!$  CASE(217)
    CASE(217)!I -4 3 m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -43m
        IVectors = 0
      CASE('b')!point symmetry -43m
        IVectors = 0
      CASE('c')!point symmetry 3m
        IVectors = 1
      CASE('d')!point symmetry -4
        IVectors = 0
      CASE('e')!point symmetry mm
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry m
        IVectors = 2
      CASE('h')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements","space group not recognised")) RETURN
      END SELECT
!!$  CASE(218)
    CASE(218)!P -4 3 n
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 23
        IVectors = 0
      CASE('b')!point symmetry 222
        IVectors = 0
      CASE('c')!point symmetry -4
        IVectors = 0
      CASE('d')!point symmetry -4
        IVectors = 0
      CASE('e')!point symmetry 3
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements","space group not recognised")) RETURN
      END SELECT
!!$  CASE(219)
    CASE(219)!F -4 3 c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 23
        IVectors = 0
      CASE('b')!point symmetry 23
        IVectors = 0
      CASE('c')!point symmetry -4
        IVectors = 0
      CASE('d')!point symmetry -4
        IVectors = 0
      CASE('e')!point symmetry 3
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements","space group not recognised")) RETURN
      END SELECT
!!$  CASE(220)
    CASE(220)!I -4 3 d
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -4
        IVectors = 0
      CASE('b')!point symmetry -4
        IVectors = 0
      CASE('c')!point symmetry 3
        IVectors = 1
      CASE('d')!point symmetry 2
        IVectors = 1
      CASE('e')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements","space group not recognised")) RETURN
      END SELECT
!!$  CASE(221)
    CASE(221)!P 4/m -3 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry m3m
        IVectors = 0
      CASE('b')!point symmetry m3m
        IVectors = 0
      CASE('c')!point symmetry 4/mmm
        IVectors = 0
      CASE('d')!point symmetry 4/mmm
        IVectors = 0
      CASE('e')!point symmetry 4mm
        IVectors = 1
      CASE('f')!point symmetry 4mm
        IVectors = 1
      CASE('g')!point symmetry 3m
        IVectors = 1
      CASE('h')!point symmetry mm
        IVectors = 1
      CASE('i')!point symmetry mm
        IVectors = 1
      CASE('j')!point symmetry mm
        IVectors = 1
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry m
        IVectors = 2
      CASE('m')!point symmetry m
        IVectors = 2
      CASE('n')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements","space group not recognised")) RETURN
      END SELECT
!!$  CASE(222)
    CASE(222)!P 4/n -3 2/n
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 43
        IVectors = 0
      CASE('b')!point symmetry 42
        IVectors = 0
      CASE('c')!point symmetry -3
        IVectors = 0
      CASE('d')!point symmetry -4
        IVectors = 0
      CASE('e')!point symmetry 4
        IVectors = 1
      CASE('f')!point symmetry 3
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements","space group not recognised")) RETURN
      END SELECT
!!$  CASE(223)
    CASE(223)!P 42/m -3 2/n
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry m3
        IVectors = 0
      CASE('b')!point symmetry mmm
        IVectors = 0
      CASE('c')!point symmetry -42m
        IVectors = 0
      CASE('d')!point symmetry -42m
        IVectors = 0
      CASE('e')!point symmetry 32
        IVectors = 0
      CASE('f')!point symmetry mm
        IVectors = 1
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry mm
        IVectors = 1
      CASE('i')!point symmetry 3
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements","space group not recognised")) RETURN
      END SELECT
!!$  CASE(224)
    CASE(224)!P 42/n -3 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -43m
        IVectors = 0
      CASE('b')!point symmetry -3m
        IVectors = 0
      CASE('c')!point symmetry -3m
        IVectors = 0
      CASE('d')!point symmetry -42m
        IVectors = 0
      CASE('e')!point symmetry 3m
        IVectors = 1
      CASE('f')!point symmetry 222
        IVectors = 1
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry 2
        IVectors = 1
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements","space group not recognised")) RETURN
      END SELECT
!!$  CASE(225)
    CASE(225)!F 4/m -3 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry m3m
        IVectors = 0
      CASE('b')!point symmetry m3n
        IVectors = 0
      CASE('c')!point symmetry -43m
        IVectors = 0
      CASE('d')!point symmetry mmm
        IVectors = 0
      CASE('e')!point symmetry 4mm
        IVectors = 1
      CASE('f')!point symmetry 3m
        IVectors = 1
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry mm
        IVectors = 1
      CASE('i')!point symmetry mm
        IVectors = 1
      CASE('j')!point symmetry m
        IVectors = 2
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements","space group not recognised")) RETURN
      END SELECT
!!$  CASE(226)
    CASE(226)!F 4/m -3 2/c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 43
        IVectors = 0
      CASE('b')!point symmetry m3
        IVectors = 0
      CASE('c')!point symmetry -42m
        IVectors = 0
      CASE('d')!point symmetry 4/m
        IVectors = 0
      CASE('e')!point symmetry mm
        IVectors = 1
      CASE('f')!point symmetry 4
        IVectors = 1
      CASE('g')!point symmetry 3
        IVectors = 1
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry m
        IVectors = 2
      CASE('j')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements","space group not recognised")) RETURN
      END SELECT
!!$  CASE(227)
    CASE(227)!F 41/d -3 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -43m
        IVectors = 0
      CASE('b')!point symmetry -43m
        IVectors = 0
      CASE('c')!point symmetry -3m
        IVectors = 0
      CASE('d')!point symmetry -3m
        IVectors = 0
      CASE('e')!point symmetry 3m
        IVectors = 1
      CASE('f')!point symmetry mm
        IVectors = 1
      CASE('g')!point symmetry m
        IVectors = 2
      CASE('h')!point symmetry 2
        IVectors = 1
      CASE('i')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements","space group not recognised")) RETURN
      END SELECT
!!$  CASE(228)
    CASE(228)!F 41/d -3 2/c
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 23
        IVectors = 0
      CASE('b')!point symmetry 32
        IVectors = 0
      CASE('c')!point symmetry -3
        IVectors = 0
      CASE('d')!point symmetry -4
        IVectors = 0
      CASE('e')!point symmetry 3
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 1
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements","space group not recognised")) RETURN
      END SELECT
!!$  CASE(229)
    CASE(229)!I 4/m -3 2/m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry m3m
        IVectors = 0
      CASE('b')!point symmetry 4/mmm
        IVectors = 0
      CASE('c')!point symmetry -3m
        IVectors = 0
      CASE('d')!point symmetry -42m
        IVectors = 0
      CASE('e')!point symmetry 4mm
        IVectors = 0
      CASE('f')!point symmetry 3m
        IVectors = 1
      CASE('g')!point symmetry mm
        IVectors = 1
      CASE('h')!point symmetry mm
        IVectors = 3
      CASE('i')!point symmetry 2
        IVectors = 1
      CASE('j')!point symmetry m
        IVectors = 2
      CASE('k')!point symmetry m
        IVectors = 2
      CASE('l')!point symmetry 1
        IVectors = 3        
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements","space group not recognised")) RETURN
      END SELECT
!!$  CASE(230)
    CASE(230)!I 41/a -3 2/d
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -3
        IVectors = 0
      CASE('b')!point symmetry 32
        IVectors = 0
      CASE('c')!point symmetry 222
        IVectors = 0
      CASE('d')!point symmetry -4
        IVectors = 0
      CASE('e')!point symmetry 3
        IVectors = 1
      CASE('f')!point symmetry 2
        IVectors = 1
      CASE('g')!point symmetry 2
        IVectors = 1
      CASE('h')!point symmetry 1
        IVectors = 3
      CASE DEFAULT     
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements","space group not recognised")) RETURN   
      END SELECT
    END SELECT
  END SUBROUTINE CountAllowedMovements

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  !>
  !! Procedure-description: 
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE DetermineAllowedMovements(ISpaceGrp,SWyckoff,RMoveMatrix,IErr)
  !Called from felixrefine ~line 1600
    USE MyNumbers
    USE message_mod

    IMPLICIT NONE

    CHARACTER(1), INTENT(IN) :: SWyckoff
    INTEGER(IKIND), INTENT(IN) :: ISpaceGrp
    INTEGER(IKIND), INTENT(OUT) :: IErr
    REAL(RKIND), INTENT(OUT) :: RMoveMatrix(ITHREE,ITHREE)
    INTEGER(IKIND) ::  ind

    !Initialisation of RMoveMatrix should not be needed, this is just to catch
    !any sloppy coding where the move count has not been done properly
    RMoveMatrix=ZERO

    SELECT CASE(ISpaceGrp)
    CASE(1)!P1
      SELECT CASE (SWyckoff)
      CASE('x')
        RMoveMatrix(1,:) = (/ZERO, ONE, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ZERO, ONE/)
      CASE('b')
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ONE, ZERO/)
        RMoveMatrix(3,:) = (/ZERO, ZERO, ONE/)
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for this space group not recognised")) RETURN     
      END SELECT
!!$  CASE(2)
!!$  CASE(3)
!!$  CASE(4)
!!$  CASE(5)
!!$  CASE(6)
!!$  CASE(7)
!!$  CASE(8)
!!$  CASE(9)
!!$  CASE(10)
!!$  CASE(11)
!!$  CASE(12)
!!$  CASE(13)
!!$  CASE(14)
!!$  CASE(15)
    CASE(15)!C 1 2/c 1
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -1

      CASE('b')!point symmetry -1

      CASE('c')!point symmetry -1

      CASE('d')!point symmetry -1

      CASE('e')!point symmetry 2 along y
        RMoveMatrix(1,:) = (/ZERO, ONE, ZERO/)
      CASE('f')!point symmetry 1
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ONE, ZERO/)
        RMoveMatrix(3,:) = (/ZERO, ZERO, ONE/)
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 15, C 2/c, not recognised")) RETURN
      END SELECT
!!$  CASE(16)
!!$  CASE(17)
!!$  CASE(18)
!!$  CASE(19)
!!$  CASE(20)
!!$  CASE(21)
!!$  CASE(22)
!!$  CASE(23)
!!$  CASE(24)
!!$  CASE(25)
!!$  CASE(26)
!!$  CASE(27)
!!$  CASE(28)
!!$  CASE(29)
!!$  CASE(30)
!!$  CASE(31)
!!$  CASE(32)
!!$  CASE(33)
!!$  CASE(34)
!!$  CASE(35)
    CASE(36)!C m c 21
      !NEED TO CODE ALTERNATIVE ORIENTATIONS Ccm21,Bb21m,Bm21b,A21ma,A21am
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 1, coordinate [x,y,z],
        RMoveMatrix(1,:) = (/ZERO, ONE, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ZERO, ONE/)
      CASE('b')!point symmetry 1, coordinate [x,y,z],
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ONE, ZERO/)
        RMoveMatrix(3,:) = (/ZERO, ZERO, ONE/)
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for this space group not recognised")) RETURN     
      END SELECT
!!$  CASE(37)
!!$  CASE(38)
!!$  CASE(39)
!!$  CASE(40)
!!$  CASE(41)
!!$  CASE(42)
!!$  CASE(43)
!!$  CASE(44)
!!$  CASE(45)
!!$  CASE(46)
!!$  CASE(47)
!!$  CASE(48)
!!$  CASE(49)
!!$  CASE(50)
!!$  CASE(51)
!!$  CASE(52)
!!$  CASE(53)
!!$  CASE(54)
!!$  CASE(55)
!!$  CASE(56)
!!$  CASE(57)
!!$  CASE(58)
!!$  CASE(59)
!!$  CASE(60)
!!$  CASE(61)
!!$  CASE(62)
    CASE(63)!Cmcm
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m, coordinate [0,0,0] & eq, no movement

      CASE('b')!point symmetry 2/m, coordinate [0,1/2,0] & eq, no movement

      CASE('c')!point symmetry mm2, coordinate [0,y,1/4] & eq
        RMoveMatrix(1,:) = (/ZERO, ONE, ZERO/)
      CASE('d')!point symmetry -1, coordinate [1/4,1/4,0] & eq, no movement

      CASE('e')!point symmetry 2(x), coordinate [x,0,0] & eq
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
      CASE('f')!point symmetry m(x), coordinate [0,y,z] & eq
        RMoveMatrix(1,:) = (/ZERO, ONE, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ZERO, ONE/)
      CASE('g')!point symmetry m(z), coordinate [x,y,1/4] & eq
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ONE, ZERO/)
      CASE('h')!point symmetry 1, coordinate [x,y,z] & eq
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ONE, ZERO/)
        RMoveMatrix(3,:) = (/ZERO, ZERO, ONE/)
      CASE DEFAULT
          IErr = 1
          IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 63, C m c m, not recognised")) RETURN
      END SELECT

    CASE(64)!Cmca
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m, coordinate [0,0,0] & eq, no movement

      CASE('b')!point symmetry 2/m, coordinate [1/2,0,0] & eq, no movement

      CASE('c')!point symmetry -1, coordinate [1/4,1/4,0] & eq, no movement

      CASE('d')!point symmetry 2(x), coordinate [x,0,0] & eq
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
      CASE('e')!point symmetry 2(y), coordinate [1/4,y,1/4] & eq
        RMoveMatrix(1,:) = (/ZERO, ONE, ZERO/)
      CASE('f')!point symmetry m(x), coordinate [0,y,z] & eq
        RMoveMatrix(1,:) = (/ZERO, ONE, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ZERO, ONE/)
      CASE('g')!point symmetry 1, coordinate [x,y,z] & eq
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ONE, ZERO/)
        RMoveMatrix(3,:) = (/ZERO, ZERO, ONE/)
      CASE DEFAULT
          IErr = 1
          IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 64, C m c a, not recognised")) RETURN
      END SELECT

!!$  CASE(65)
!!$  CASE(66)
!!$  CASE(67)
     CASE(68)!C c c a
      !N.B. multiple origin choices allowed, here origin at 222, -1 at
      ![1/4,0,1/4]
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 222, coordinate [0,0,0] & eq, no movement

      CASE('b')!point symmetry 222, coordinate [0,0,1/2] & eq, no movement

      CASE('c')!point symmetry -1, coordinate [1/4,0,1/4] & eq, no movement

      CASE('d')!point symmetry -1, coordinate [0,1/4,1/4] & eq, no movement

      CASE('e')!point symmetry 2(x), coordinate [x,0,0] & eq
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
      CASE('f')!point symmetry 2(y), coordinate [0,y,0] & eq
        RMoveMatrix(1,:) = (/ZERO, ONE, ZERO/)
      CASE('g')!point symmetry 2(z), coordinate [0,0,z] & eq
        RMoveMatrix(1,:) = (/ZERO, ZERO, ONE/)
      CASE('h')!point symmetry 2(z), coordinate [1/4,1/4,z] & eq
        RMoveMatrix(1,:) = (/ZERO, ZERO, ONE/)
      CASE('i')!point symmetry 1, coordinate [x,y,z] & eq
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ONE, ZERO/)
        RMoveMatrix(3,:) = (/ZERO, ZERO, ONE/)
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"PreferredBasis",&
          "Wyckoff Symbol for space group 68, C c c a, not recognised")) RETURN
      END SELECT

!!$  CASE(69)
!!$  CASE(70)
!!$  CASE(71)
!!$  CASE(72)
!!$  CASE(73)
!!$  CASE(74)
!!$  CASE(75)
!!$  CASE(76)
!!$  CASE(77)
!!$  CASE(78)
!!$  CASE(79)
!!$  CASE(80)
!!$  CASE(81)
!!$  CASE(82)
!!$  CASE(83)
!!$  CASE(84)
!!$  CASE(85)
!!$  CASE(86)
!!$  CASE(87)
!!$  CASE(88)
!!$  CASE(89)
!!$  CASE(90)
!!$  CASE(91)
!!$  CASE(92)
!!$  CASE(93)
!!$  CASE(94)
!!$  CASE(95)
!!$  CASE(96)
!!$  CASE(97)
!!$  CASE(98)
     CASE(99)!P 4 m m 
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 4mm, coordinate [0,0,z], allowed movement along [001]
        RMoveMatrix(1,:) = (/ZERO, ZERO, ONE/)
      CASE('b')!point symmetry 4mm, coordinate [1/2,1/2,z], allowed movement along [001]
        RMoveMatrix(1,:) = (/ZERO, ZERO, ONE/)
      CASE('c')!point symmetry mm, coordinate [1/2,0,z] & eq, allowed movement along [001]
        RMoveMatrix(1,:) = (/ZERO, ZERO, ONE/)
      CASE('d')!point symmetry m, coordinate [x,x,z] & eq, allowed movement along [110] & [001]
        RMoveMatrix(1,:) = (/ONE, ONE, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ZERO, ONE/)
      CASE('e')!point symmetry m, coordinate [x,0,z] & eq, allowed movement along [100] & [001]
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ZERO, ONE/)
      CASE('f')!point symmetry m, coordinate [x,1/2,z] & eq, allowed movement along [100] & [001]
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ZERO, ONE/)
      CASE('g')!point symmetry 1, coordinate [x,y,z] & eq
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ONE, ZERO/)
        RMoveMatrix(3,:) = (/ZERO, ZERO, ONE/)
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"PreferredBasis",&
          "Wyckoff Symbol for space group 68, C c c a, not recognised")) RETURN
      END SELECT

!!$  CASE(100)
!!$  CASE(101)
!!$  CASE(102)
!!$  CASE(103)
!!$  CASE(104)
!!$  CASE(105)
!!$  CASE(106)
!!$  CASE(107)
!!$  CASE(108)
!!$  CASE(109)
!!$  CASE(110)
!!$  CASE(111)
!!$  CASE(112)
!!$  CASE(113)
!!$  CASE(114)
!!$  CASE(115)
!!$  CASE(116)
!!$  CASE(117)
!!$  CASE(118)
!!$  CASE(119)
!!$  CASE(120)
!!$  CASE(121)
!!$  CASE(122)
!!$  CASE(123)
!!$  CASE(124)
!!$  CASE(125)
!!$  CASE(126)
!!$  CASE(127)
!!$  CASE(128)
!!$  CASE(129)
!!$  CASE(130)
!!$  CASE(131)
!!$  CASE(132)
!!$  CASE(133)
!!$  CASE(134)
!!$  CASE(135)
!!$  CASE(136)
!!$  CASE(137)
!!$  CASE(138)
    CASE(139)!I4/m m m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 4/mmm, coordinate [0,0,0], no allowed movements
      
      CASE('b')!point symmetry 4/mmm, coordinate [0,0,1/2], no allowed movements
      
      CASE('c')!point symmetry mmm, coordinate [0,1/2,0] or [1/2,0,0], no allowed movements
      
      CASE('d')!point symmetry -4m2, coordinate [0,1/2,1/4] or [1/2,0,1/4], no allowed movements
      
      CASE('e')!point symmetry 4mm, coordinate [0,0,z] allowed movement along z
        RMoveMatrix(1,:) = (/ZERO, ZERO, ONE/)
      CASE('f')!point symmetry 2/m, coordinate [1/4,1/4,1/4] & eq, no allowed movements
        
      CASE('g')!point symmetry mm, coordinate [0,1/2,z] & eq, allowed movement along z
        RMoveMatrix(1,:) = (/ZERO, ZERO, ONE/)
      CASE('h')!point symmetry mm, coordinate [x,x,0] & eq, allowed movement along [110]
        RMoveMatrix(1,:) = (/ONE, ONE, ZERO/)
      CASE('i')!point symmetry mm, coordinate [x,0,0] & eq, allowed movement along [100]
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
      CASE('j')!point symmetry mm, coordinate [x,1/2,0] & eq, allowed movement along [100]
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
      CASE('k')!point symmetry 2, coordinate [x,1/2+x,1/4] & eq, allowed movement along [110]
        RMoveMatrix(1,:) = (/ONE, ONE, ZERO/)
      CASE('l')!point symmetry m, coordinate [x,y,0] & eq, allowed movement along [100] & [010]
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ONE, ZERO/)
      CASE('m')!point symmetry m, coordinate [x,x,z] & eq, allowed movement along [110] & [001]
        RMoveMatrix(1,:) = (/ONE, ONE, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ZERO, ONE/)
      CASE('n')!point symmetry m, coordinate [x,0,z] & eq, allowed movement along [100] & [001]
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ZERO, ONE/)
      CASE('o')!point symmetry 1, coordinate [x,y,z] & eq, allowed movement along [100], [010] & [001]
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ONE, ZERO/)
        RMoveMatrix(3,:) = (/ZERO, ZERO, ONE/)
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
          "Wyckoff Symbol for space group I4/m m m not recognised")) RETURN
      END SELECT
!!$  CASE(140)
!!$  CASE(141)
     CASE(142)!I41/acd
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -4, no allowed movements

      CASE('b')!point symmetry 222, no allowed movements

      CASE('c')!point symmetry -1, no allowed movements

      CASE('d')!point symmetry 2, allowed movement along z
        RMoveMatrix(1,:) = (/ZERO, ZERO, ONE/)
      CASE('e')!point symmetry 2, allowed movement along x
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
      CASE('f')!point symmetry 2, allowed movement along [x,x,0]
        RMoveMatrix(1,:) = (/ONE, ONE, ZERO/)
      CASE('g')!point symmetry 1
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ONE, ZERO/)
        RMoveMatrix(3,:) = (/ZERO, ZERO, ONE/)
      CASE DEFAULT
         IErr = 1
         IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 142, I 41/a 2/c 2/d, not recognised")) RETURN
      END SELECT
!!$  CASE(143)
!!$  CASE(144)
!!$  CASE(145)
!!$  CASE(146)
!!$  CASE(147)
!!$  CASE(148)
!!$  CASE(149)
!!$  CASE(150)
!!$  CASE(151)
!!$  CASE(152)
!!$  CASE(153)
!!$  CASE(154)
!!$  CASE(155)
!!$  CASE(156)
!!$  CASE(157)
!!$  CASE(158)
!!$  CASE(159)
!!$  CASE(160)
!!$  CASE(161)
!!$  CASE(162)
!!$  CASE(163)
!!$  CASE(164)
!!$  CASE(165)
!!$  CASE(166)
!!$  CASE(167)
!!$  CASE(168)
!!$  CASE(169)
!!$  CASE(170)
!!$  CASE(171)
!!$  CASE(172)
!!$  CASE(173)
!!$  CASE(174)
!!$  CASE(175)
!!$  CASE(176)
!!$  CASE(177)
!!$  CASE(178)
!!$  CASE(179)
!!$  CASE(180)
!!$  CASE(181)
!!$  CASE(182)
!!$  CASE(183)
!!$  CASE(184)
!!$  CASE(185)
!!$  CASE(186)
!!$  CASE(187)
!!$  CASE(188)
!!$  CASE(189)
!!$  CASE(190)
!!$  CASE(191)
!!$  CASE(192)
!!$  CASE(193)
!!$  CASE(194)
!!$  CASE(195)
!!$  CASE(196)
!!$  CASE(197)
!!$  CASE(198)
!!$  CASE(199)
!!$  CASE(200)
!!$  CASE(201)
!!$  CASE(202)
!!$  CASE(203)
!!$  CASE(204)
!!$  CASE(205)
!!$  CASE(206)
!!$  CASE(207)
!!$  CASE(208)
!!$  CASE(209)
!!$  CASE(210)
!!$  CASE(211)
!!$  CASE(212)
!!$  CASE(213)
!!$  CASE(214)
!!$  CASE(215)
    CASE(216)!F-43m
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry -43m, coordinate [0,0,0], no allowed movements

      CASE('b')!point symmetry -43m, coordinate [1/2,1/2,1/2], no allowed movements

      CASE('c')!point symmetry -43m, coordinate [1/4,1/4,1/4], no allowed movements

      CASE('d')!point symmetry -43m, coordinate [3/4,3/4,3/4], no allowed movements

      CASE('e')!point symmetry 3m, coordinate [x,x,x] allowed movement along [111]
        RMoveMatrix(1,:) = (/ONE, ONE, ONE/)
      CASE('f')!point symmetry mm, coordinate [x,0,0] allowed movements along x
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
      CASE('g')!point symmetry mm, coordinate [x,1/4,1/4] allowed movement along x
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
      CASE('h')!point symmetry m, coordinate [x,x,z], allowed movement along [110] and [001]
        RMoveMatrix(1,:) = (/ONE, ONE, ZERO/)
        RMoveMatrix(1,:) = (/ZERO, ZERO, ONE/)
      CASE('i')!point symmetry 1, coordinate [x,y,z], allowed movement along x,y,z
        RMoveMatrix(1,:) = (/ONE, ZERO, ZERO/)
        RMoveMatrix(2,:) = (/ZERO, ONE, ZERO/)
        RMoveMatrix(3,:) = (/ZERO, ZERO, ONE/)
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
          "Wyckoff Symbol for space group I4/m m m not recognised")) RETURN
      END SELECT

!!$  CASE(217)
!!$  CASE(218)
!!$  CASE(219)
!!$  CASE(220)
!!$  CASE(221)
!!$  CASE(222)
!!$  CASE(223)
!!$  CASE(224)
!!$  CASE(225)
!!$  CASE(226)
!!$  CASE(227)
!!$  CASE(228)
!!$  CASE(229)
!!$  CASE(230)
    CASE DEFAULT
      IErr = 1
      IF(l_alert(IErr,"DetermineAllowedMovements","space group not recognised")) RETURN 
    END SELECT
  END SUBROUTINE DetermineAllowedMovements

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Convert SSpacegrp to lower case and Compare
  !! SSpaceGrpNoSpaces with every space group
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE ChangeOrigin(ISpaceGrp,IErr)

    ! called from PreferredBasis to change to the preferred origin
    USE MyNumbers
    USE message_mod

    ! global inputs
    USE RPARA, ONLY : RBasisAtomPosition
    USE SPARA, ONLY : SWyckoffSymbol

    IMPLICIT NONE

    INTEGER(IKIND),INTENT(OUT) :: IErr
    INTEGER(IKIND) :: ind,ISpaceGrp,IBasisAtoms,IChangeFLAG

    IBasisAtoms=SIZE(RBasisAtomPosition,1)
    SELECT CASE(ISpaceGrp)
      CASE(142)!I41/acd
        !Change from choice 1 (origin -4 at [0,0,0]) to choice 2 (origin -1 at [0,0,0])
        IChangeFLAG = 0
        IErr = 1!Set the error flag in the case we don't get an 'a' site
        
        !Go through the atoms and look for an 'a' site incompatible with choice 2
        DO ind = 1,IBasisAtoms
          !Wyckoff 'a' is -4 at [000],[0,1/2,1/2],[0,1/2,1/4],[1/2,0,1/4] in 1
          ! and [0,1/4,3/8], [0,3/4,5/8], [1/2,1/4,5/8], [1/2,3/4,5/8] in 2
          IF (SWyckoffSymbol(ind).EQ.'a') THEN
            IErr = 0!We have an 'a' site
            !check for origin 1 by multiplying by 4 and checking it is an integer
            IF(MODULO(FOUR*RBasisAtomPosition(ind,3),1.0).LT.TINY) IChangeFLAG = 1
          END IF
        END DO
        IF (IChangeFLAG.EQ.1) THEN!add[[0,1/4,3/8]
          DO ind = 1,IBasisAtoms
            RBasisAtomPosition(ind,2)=MODULO((RBasisAtomPosition(ind,2)+0.25),1.0)
            RBasisAtomPosition(ind,3)=MODULO((RBasisAtomPosition(ind,2)+0.375),1.0)
          END DO
        END IF
    END SELECT
    
  END SUBROUTINE ChangeOrigin

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Convert SSpacegrp to lower case and Compare
  !! SSpaceGrpNoSpaces with every space group
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE ConvertSpaceGroupToNumber(ISpaceGrp,IErr)

    ! called from 1) PreferredBasis and 2) felixrefine to setup coordinate refinement
    USE MyNumbers
    USE message_mod

    ! global inputs
    USE SPARA, ONLY : SSpaceGrp
    USE SConst, ONLY : SAllSpaceGrp

    IMPLICIT NONE

    INTEGER(IKIND),INTENT(OUT) :: ISpaceGrp, IErr
    INTEGER(IKIND) :: jnd, IIndex, ind
    CHARACTER(LEN(SSpaceGrp)) :: SSpaceGrpNoSpaces
    CHARACTER(20) :: SSpaceGrpToCompare

    ! Push Spaces In SSpaceGrp to the end of the String
    
    jnd = 0
    ISpaceGrp = 0
    SSpaceGrpNoSpaces = ' '
    SSpaceGrpToCompare = ' '
    
    DO ind = 1,LEN(SSpaceGrp)
       IF(INDEX(STRING='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ01234567879-/',&
            SUBSTRING=SSpaceGrp(ind:ind)).NE.0) THEN
          jnd = jnd + 1
          SSpaceGrpNoSpaces(jnd:jnd) = SSpaceGrp(ind:ind)
       END IF
    END DO
    
    ! Convert SSpacegrp to lower case 

    CALL StrLowCase( SSpaceGrpNoSpaces,SSpaceGrpNoSpaces,IErr )

    ! Compare SSpaceGrpNoSpaces with every space group 

    DO ind = 1,SIZE(SAllSpaceGrp)

       CALL StrLowCase( SAllSpaceGrp(ind),SSpaceGrpToCompare,IErr )
       IIndex = INDEX(TRIM(ADJUSTL(SSpaceGrpToCompare)),TRIM(ADJUSTL(SSpaceGrpNoSpaces)))
       IF (IIndex.NE.0) THEN
          ISpaceGrp = ind
          EXIT
       END IF
    END DO

!DBG IF (my_rank.EQ.0) PRINT*, ISpaceGrp
    IF(ISpaceGrp.EQ.0) THEN
      IErr = 1
      IF(l_alert(IErr,"ConvertSpaceGroupToNumber",&
            "Space Group was not found. Check .cif file")) RETURN
    END IF
    
  END SUBROUTINE ConvertSpaceGroupToNumber

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: 
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE StrLowCase( Input_String,Output_String,IErr )
    ! used twice in ConvertSpaceGroupToNumber

    USE MyNumbers

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: Input_String
    CHARACTER(LEN(Input_String)), INTENT(OUT) :: Output_String
    INTEGER(IKIND), INTENT(OUT) :: IErr
    CHARACTER(*), PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz',&
         UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
    INTEGER(IKIND) :: ind, n

    IErr=0
    ! Copy input string
    Output_String = Input_String
    
    ! Convert case character by character
    DO ind = 1, LEN(Output_String,KIND=IKIND)
       n = INDEX(UPPER_CASE, Output_String(ind:ind))
       IF ( n.NE.0 ) Output_String(ind:ind) = LOWER_CASE(n:n)
    END DO
  END SUBROUTINE  StrLowCase

END MODULE setup_space_group_mod
