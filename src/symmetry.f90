!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixsim
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
!  This file is part of felixsim.
!
!  felixsim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  felixsim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with felixsim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: smodules.f90,v 1.63 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE DetermineAllowedMovements(ISpaceGrp,SWyckoffSymbol,RVector,IVector,IErr)

  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: ISpaceGrp,IErr,IVector,ind
  CHARACTER*1 :: SWyckoffSymbol
  REAL(RKIND),DIMENSION(IVector,ITHREE),INTENT(OUT) :: RVector

  !Tells user, entering DetermineAllowedMovements
  CALL Message("DetermineAllowedMovements",IMust,IErr)

  SELECT CASE(ISpaceGrp)
  CASE(1)
     SELECT CASE (SWyckoffSymbol)
     CASE('x')
       RVector(1,:) = (/ZERO, ONE, ZERO/)
       RVector(2,:) = (/ZERO, ZERO, ONE/)
     CASE('b')
       RVector(1,:) = (/ONE, ZERO, ZERO/)
       RVector(2,:) = (/ZERO, ONE, ZERO/)
       RVector(3,:) = (/ZERO, ZERO, ONE/)
     CASE DEFAULT
        PRINT*,"---------------------------------------"
        PRINT*,"This Wyckoff Symbol for this space group is not yet implemented"
        PRINT*,"You can correct this in symmetry.f90"
        PRINT*,"---------------------------------------"        
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
     SELECT CASE (SWyckoffSymbol)
     CASE('a')
       RVector(1,:) = (/ZERO, ONE, ZERO/)
       RVector(2,:) = (/ZERO, ZERO, ONE/)
     CASE('b')
       RVector(1,:) = (/ONE, ZERO, ZERO/)
       RVector(2,:) = (/ZERO, ONE, ZERO/)
       RVector(3,:) = (/ZERO, ZERO, ONE/)
     CASE DEFAULT
        PRINT*,"---------------------------------------"
        PRINT*,"This Wyckoff Symbol for this space group is not yet implemented"
        PRINT*,"You can correct this in symmetry.f90"
        PRINT*,"---------------------------------------"     
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
!!$  CASE(63)
!!$  CASE(64)
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
!!$  CASE(139)
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
     PRINT*,"---------------------------------------"
     PRINT*,"This space group is not yet implemented"
     PRINT*,"You can correct this in symmetry.f90"
     PRINT*,"---------------------------------------"
  END SELECT
END SUBROUTINE DetermineAllowedMovements

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE CountAllowedMovements(ISpaceGrp,SWyckoffSymbol,IVectors,IErr)

  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND),INTENT(IN) :: ISpaceGrp
  INTEGER(IKIND) :: IVectors,IErr
  CHARACTER*1 :: SWyckoffSymbol

  !Tells user, entering CountAllowedMovements
  CALL Message("CountAllowedMovements",IMust,IErr)
  
  SELECT CASE(ISpaceGrp)
  CASE(1)
     SELECT CASE (SWyckoffSymbol)
     CASE('x')
        IVectors = 2_IKIND
     CASE('b')
        IVectors = 3_IKIND
     CASE DEFAULT
        PRINT*,"---------------------------------------"
        PRINT*,"This Wyckoff Symbol for this space group is not yet implemented"
        PRINT*,"You can correct this in symmetry.f90"
        PRINT*,"---------------------------------------"     
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
     SELECT CASE (SWyckoffSymbol)
     CASE('a')
        IVectors = 2_IKIND
     CASE('b')
        IVectors = 3_IKIND
     CASE DEFAULT
        PRINT*,"---------------------------------------"
        PRINT*,"This Wyckoff Symbol for this space group is not yet implemented"
        PRINT*,"You can correct this in symmetry.f90"
        PRINT*,"---------------------------------------"     
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
!!$  CASE(63)
!!$  CASE(64)
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
!!$  CASE(139)
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
     PRINT*,"---------------------------------------"
     PRINT*,"This space group is not yet implemented"
     PRINT*,"You can correct this in symmetry.f90"
     PRINT*,"---------------------------------------"
  END SELECT
END SUBROUTINE CountAllowedMovements

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE ConvertSpaceGroupToNumber(ISpaceGrp,IErr)

  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,jnd,IIndex,ind
  INTEGER(IKIND),INTENT(OUT) :: ISpaceGrp
  CHARACTER(LEN(SSpaceGrp)) :: SSpaceGrpNoSpaces
  CHARACTER*20 :: SSpaceGrpToCompare
       
  
  !Tells user, entering ConvertSpaceGroupToNumber
  CALL Message("ConvertSpaceGroupToNumber",IMust,IErr)

!!$  Push Spaces In SSpaceGrp to the end of the String
  
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
  
!!$  Convert SSpacegrp to lower case 

  CALL StrLowCase( SSpaceGrpNoSpaces,SSpaceGrpNoSpaces,IErr )

!!$  Compare SSpaceGrpNoSpaces with every space group 

  DO ind = 1,SIZE(CSpaceGrp)

     CALL StrLowCase( CSpaceGrp(ind),SSpaceGrpToCompare,IErr )
     IIndex = INDEX(TRIM(ADJUSTL(SSpaceGrpToCompare)),TRIM(ADJUSTL(SSpaceGrpNoSpaces)))
     IF (IIndex.NE.0) THEN
        ISpaceGrp = ind
        EXIT
     END IF
  END DO

  IF(ISpaceGrp.EQ.0) THEN
     IErr = 1
     IF(my_rank.EQ.0) THEN
        PRINT*,"Your Space Group",SSpaceGrpNoSpaces,"Was not found, Check your .cif file"
     END IF
  END IF
  
END SUBROUTINE ConvertSpaceGroupToNumber

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE StrLowCase( Input_String,Output_String,IErr ) 

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  CHARACTER(*), INTENT(IN)     :: &
       Input_String
  CHARACTER(LEN(Input_String)) :: &
       Output_String
  CHARACTER(*),PARAMETER :: &
       LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz',&
       UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
  INTEGER(IKIND) :: &
       ind, n,IErr

  ! Copy input string
  Output_String = Input_String
  
  ! Convert case character by character

  DO ind = 1, LEN(Output_String,KIND=IKIND)
     n = INDEX(UPPER_CASE, Output_String(ind:ind))
     IF ( n.NE.0 ) Output_String(ind:ind) = LOWER_CASE(n:n)
  END DO
END SUBROUTINE  StrLowCase

SUBROUTINE MakeUnitVector(RVector,IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind
  REAL(RKIND),DIMENSION(:,:),INTENT(INOUT) :: &
       RVector
  
  DO ind = 1,SIZE(RVector,DIM=1)
     RVector(ind,:) = RVector(ind,:)/SQRT(SUM(RVector(ind,:)**2))
  END DO
END SUBROUTINE MakeUnitVector
