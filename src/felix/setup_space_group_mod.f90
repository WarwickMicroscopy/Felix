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
    USE MyNumbers
    USE message_mod
    
    ! global inputs
    USE SPARA, ONLY : SSpaceGrp, SWyckoffSymbol
    USE RPARA, ONLY : RBasisAtomPosition

    IMPLICIT NONE

    CHARACTER*1 :: SWyckoff
    INTEGER(IKIND), INTENT(OUT) :: IErr
    INTEGER(IKIND) ::  ind,ISpaceGrp,IBasisAtoms

    IErr=0

    CALL ConvertSpaceGroupToNumber(ISpaceGrp,IErr)!Uses global variable SSpaceGrp
    IF(l_alert(IErr,"PreferredBasis","ConvertSpaceGroupToNumber")) RETURN
    
    IBasisAtoms=SIZE(RBasisAtomPosition,1)
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
     CASE(63)
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

     CASE(64)!Cmca.  No equivalent positios swap axes so no reassignments in fact
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
!!$  CASE(68)
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
!!$  CASE(99)
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
!!$  CASE(142)
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
!!$  CASE(216)
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
    CHARACTER*1, INTENT(IN) :: SWyckoff
    INTEGER(IKIND), INTENT(OUT) :: IVectors, IErr
    
    IErr=0
    SELECT CASE(ISpaceGrp)
    CASE(1)!P1
      SELECT CASE (SWyckoff)
      CASE('a')
        IVectors = 3_IKIND
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
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
      CASE('a')
        IVectors = 2
      CASE('b')
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
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
    CASE(63)!Cmcm
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m, coordinate [0,0,0] & eq, no movement

      CASE('b')!point symmetry 2/m, coordinate [0,1/2,0] & eq, no movement

      CASE('c')!point symmetry mm2, coordinate [0,y,1/4] & eq
        IVectors = 1
      CASE('d')!point symmetry -1, coordinate [1/4,1/4,0] & eq, no movement

      CASE('e')!point symmetry 2(x), coordinate [x,0,0] & eq
        IVectors = 1
      CASE('f')!point symmetry m(x), coordinate [0,y,z] & eq
        IVectors = 2
      CASE('g')!point symmetry m(z), coordinate [x,y,1/4] & eq
        IVectors = 2
      CASE('h')!point symmetry 1, coordinate [x,y,z] & eq
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"CountAllowedMovements",&
              "Wyckoff Symbol for space group 63, C m c m, not recognised")) RETURN
      END SELECT

    CASE(64)!Cmca
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry 2/m, coordinate [0,0,0] & eq, no movement

      CASE('b')!point symmetry 2/m, coordinate [1/2,0,0] & eq, no movement

      CASE('c')!point symmetry -1, coordinate [1/4,1/4,0] & eq, no movement

      CASE('d')!point symmetry 2(x), coordinate [x,0,0] & eq
        IVectors = 1
      CASE('e')!point symmetry 2(y), coordinate [1/4,y,1/4] & eq
        IVectors = 1
      CASE('f')!point symmetry m(x), coordinate [0,y,z] & eq
        IVectors = 2
      CASE('g')!point symmetry 1, coordinate [x,y,z] & eq
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"CountAllowedMovements",&
              "Wyckoff Symbol for space group 64, C m c a, not recognised")) RETURN
      END SELECT
!!$  CASE(65)
!!$  CASE(66)
!!$  CASE(67)
!!$  CASE(68)
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
!!$  CASE(99)
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
        IVectors = 0
      CASE('b')!point symmetry 4/mmm, coordinate [0,0,1/2], no allowed movements
        IVectors = 0
      CASE('c')!point symmetry mmm, coordinate [0,1/2,0] or [1/2,0,0], no allowed movements
        IVectors = 0
      CASE('d')!point symmetry -4m2, coordinate [0,1/2,1/4] or [1/2,0,1/4], no allowed movements
        IVectors = 0
      CASE('e')!point symmetry 4mm, coordinate [0,0,z] allowed movement along z
        IVectors = 1
      CASE('f')!point symmetry 2/m, coordinate [1/4,1/4,1/4] & eq, no allowed movements
        IVectors = 0
      CASE('g')!point symmetry mm, coordinate [0,1/2,z] & eq, allowed movement along z
        IVectors = 1
      CASE('h')!point symmetry mm, coordinate [x,x,0] & eq, allowed movement along [110]
        IVectors = 1
      CASE('i')!point symmetry mm, coordinate [x,0,0] & eq, allowed movement along [100]
        IVectors = 1
      CASE('j')!point symmetry mm, coordinate [x,1/2,0] & eq, allowed movement along [100]
        IVectors = 1
      CASE('k')!point symmetry 2, coordinate [x,1/2+x,1/4] & eq, allowed movement along [100]
        IVectors = 1
      CASE('l')!point symmetry m, coordinate [x,y,0] & eq, allowed movement along [100] & [010]
        IVectors = 2
      CASE('m')!point symmetry m, coordinate [x,x,z] & eq, allowed movement along [110] & [001]
        IVectors = 2
      CASE('n')!point symmetry m, coordinate [x,0,z] & eq, allowed movement along [100] & [001]
        IVectors = 2
      CASE('o')!point symmetry 1, coordinate [x,y,z] & eq, allowed movement along [100], [010] & [001]
        IVectors = 3
      CASE DEFAULT
        IErr = 1
        IF(l_alert(IErr,"DetermineAllowedMovements",&
              "Wyckoff Symbol for space group 139, I4/m m m, not recognised")) RETURN     
      END SELECT    
!!$  CASE(140)
!!$  CASE(141)
!!$  CASE(142)
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
!!$  CASE(216)
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

    CHARACTER*1, INTENT(IN) :: SWyckoff
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
    CASE(36)
      !NEED TO CODE ALTERNATIVE ORIENTATIONS Ccm21,Bb21m,Bm21b,A21ma,A21am
      SELECT CASE (SWyckoff)
      CASE('a')!point symmetry m, coordinate [0,y,z],
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
!!$  CASE(68)
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
!!$  CASE(99)
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
        !Redefine basis set since [x,-x,0] is an equivalent coordinate
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
!!$  CASE(142)
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
!!$  CASE(216)
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
    CHARACTER*20 :: SSpaceGrpToCompare

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
