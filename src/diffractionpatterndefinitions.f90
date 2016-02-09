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
! $Id: diffractionpatterndefinitions.f90,v 1.11 2014/03/25 15:37:30 phsht Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE ReflectionDetermination( IErr )

  USE MyNumbers
  USE WriteToScreen
  USE ErrorCodes

  USE CConst; USE IConst
  USE IPara; USE RPara
  USE SPara
  USE CPara
  USE IChannels
  USE BlochPara

  USE MyMPI

  IMPLICIT NONE

  REAL(RKIND), DIMENSION(:,:), ALLOCATABLE :: RgDummyVecMat,RgVecMagLaue
  REAL(RKIND) :: dummy,RMaxAcceptanceGVecMag,RMinLaueZoneValue,RMaxLaueZoneValue,RGzUnitVec, &
       RLaueZoneGz,RLaueZoneElectronWaveVectorMag
  INTEGER(IKIND), DIMENSION(:),ALLOCATABLE :: IOriginGVecIdentifier
  INTEGER(IKIND) :: IErr,ind,jnd,icheck,ihklrun,IFind,IFound,knd,IMaxLaueZoneLevel, &
       ICutOff,ITotalLaueZoneLevel,ICounter,IHOLZGVecMagSize,INumInitReflections, &
       INumFinalReflections, IBSMaxLocGVecAmp, IZerothLaueZoneLevel,INumTotalReflections, &
       ILaueLevel
  CHARACTER*20 :: Sind,Sjnd

  CALL Message("ReflectionDetermination",IMust,IErr)

  icheck = 0

  IHKLMAXValue = 15!RB starting value, increments if necessary
  ihklrun = 0

  DO WHILE (icheck.EQ.0)
     ihklrun = ihklrun+1     

     CALL Message("ReflectionDetermination",IInfo,IErr,MessageVariable = "IHklrun", &
          IVariable = IHklrun)

     CALL NewHKLMake(IHKLMAXValue,RZDirC,TWODEG2RADIAN,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"Reflectiondetermination(", my_rank, ") error in NewHKLMake()"
        RETURN
     END IF

     IF(SIZE(Rhkl,DIM=1).LT.IMinReflectionPool) THEN
        IHKLMAXValue = IHKLMAXValue*2
        Deallocate(Rhkl,STAT=ierr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"ReflectionDetermination(",my_rank,") error",IErr,"deallocating Rhkl"
           RETURN
        END IF

        CYCLE

     ELSE
        icheck = 1
     END IF
  END DO

  CALL SortHKL(Rhkl,SIZE(Rhkl,1),IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReflectionDetermination(): error in SortHKL"
     RETURN
  END IF

  IF(IXDirectionFLAG.EQ.0) THEN
     IDiffractionFLAG = 1
     RXDirC(1) = Rhkl(2,1)
     RXDirC(2) = Rhkl(2,2)
     RXDirC(3) = Rhkl(2,3)
     CALL CrystallographyInitialisation( IErr )
     IF( IErr.NE.0 ) THEN
        PRINT*,"DiffractionPatternDefinitions(", my_rank, ") error",IErr, &
             "in CrystallographyInitialisation()"
        RETURN
     END IF

  END IF

  ALLOCATE(RgVecMatT(SIZE(Rhkl,DIM=1),THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternDefinitions(",my_rank,")error allocating RgVecMatT"
     RETURN
  END IF

  ALLOCATE(RgDummyVecMat(SIZE(Rhkl,DIM=1),THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternDefinitions(",my_rank,")error allocating RgDummyVecMat"
     RETURN
  END IF

  ALLOCATE(RgVecMag(SIZE(Rhkl,DIM=1)), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternDefinitions(",my_rank,")error allocating RgVecMag"
     RETURN
  END IF
  
  ALLOCATE(RgVecVec(SIZE(Rhkl,DIM=1)), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternDefinitions(",my_rank,")error allocating RgVecVec"
     RETURN
  END IF

!!$  For IZOLZFLAG=0, we want to identify the number of HOLZ in the Bloch Problem

!!$  Loop through the z-component of the g-vector matrix in the microscope frame
!!$  (Along the K-vector in the Bragg condition) 
!!$  Find the minimum g-z vector - RGzUnitVector (this is the quantisation condition)
!!$  ICutoff stops the conditional statement determining the minimum g-vector
!!$  Also populate a Dummy g-vector matrix  
  ICutOff = 1
  DO ind=1,SIZE(Rhkl,DIM=1)
     WRITE(Sind,'(I10.1)')ind
     DO jnd=1,THREEDIM
        RgVecMatT(ind,jnd)= &!RB what is RgVecMatT???
             Rhkl(ind,1)*RarVecM(jnd) + &
             Rhkl(ind,2)*RbrVecM(jnd) + &
             Rhkl(ind,3)*RcrVecM(jnd)
        RgDummyVecMat(ind,jnd)=RgVecMatT(ind,jnd)
     ENDDO
     CALL Message("ReflectionDetermination",IMoreInfo,IErr, &
          MessageVariable = "Rhkl(h,k,l)", &
          RVector = Rhkl(ind,:))
     IF((RgVecMatT(ind,3).GT.TINY.OR.RgVecMatT(ind,3).LT.-TINY).AND.ICutOff.NE.0) THEN
        RGzUnitVec=ABS(RgVecMatT(ind,3))
        ICutOff=0
     END IF
  ENDDO

  !No higher order Laue Zones with ZOLZFlag switched on
  IF(ICutOff.EQ.1) THEN
     RGzUnitVec=ZERO
  END IF
  CALL Message("ReflectionDetermination", IMoreInfo,IErr, &
       MessageVariable="RGzUnitVec", RVariable=RGzUnitVec)

!!$  Divide the non-zero Gz positions in the matrix by the mimimum Gz distance
!!$  This quantises the dummy matrix identifying the various Laue Zones
  WHERE(RgDummyVecMat(:,3).GT.TINY.OR.RgDummyVecMat(:,3).LT.-TINY)
     RgDummyVecMat(:,3)=RgDummyVecMat(:,3)/RGzUnitVec
  END WHERE

!!$     Write out the z component Dummy g-vector Matrix in DebugMODE
  DO ind=1,SIZE(Rhkl,DIM=1)
     WRITE(Sind,'(I10.1)')ind
     CALL Message("ReflectionDetermination", IAllInfo+IDEBUG,IErr, &
          MessageVariable="Reciprocal Vector" &
          //",RgDummyVecMat"//"("//TRIM(ADJUSTL(Sind))//",3)", &
          RVariable=RgDummyVecMat(ind,3))
  END DO

!!$     Find the maximum and minimum Laue Zone for the system in question
!!$     and print to screen
  RMaxLaueZoneValue=MAXVAL(RgDummyVecMat(:,3),DIM=1)
  RMinLaueZoneValue=MINVAL(RgDummyVecMat(:,3),DIM=1)
  CALL Message("ReflectionDetermination", IInfo+IDEBUG,IErr, &
       MessageVariable="Minimum Laue Zone,RMinLaueZoneValue", &
       RVariable=RMinLaueZoneValue)
  CALL Message("ReflectionDetermination", IInfo+IDEBUG,IErr, &
       MessageVariable="Maximum Laue Zone,RMaxLaueZoneValue", &
       RVariable=RMaxLaueZoneValue)
  CALL Message("ReflectionDetermination", IInfo,IErr, &
       MessageVariable="Minimum Laue Zone,RMinLaueZoneValue", &
       IVariable=NINT(RMinLaueZoneValue,IKIND))
  CALL Message("ReflectionDetermination", IInfo,IErr, &
       MessageVariable="Maximum Laue Zone,RMaxLaueZoneValue", &
       IVariable=NINT(RMaxLaueZoneValue,IKIND))

!!$     Identify the total Laue Zones in the system and print to screen
  ITotalLaueZoneLevel=NINT(RMaxLaueZoneValue+ABS(RMinLaueZoneValue)+1,IKIND)
  CALL Message("ReflectionDetermination", IInfo,IErr, &
       MessageVariable="Total Laue Zones in system, ITotalLaueZoneLevel", &
       IVariable=ITotalLaueZoneLevel)
  IZerothLaueZoneLevel=NINT(ABS(RMinLaueZoneValue)+1,IKIND)

!!$     HOLZ Acceptance Angle
!!$     Allocate enough space to determine Acceptance angle (in reciprocal Angstroms) 
!!$     For each Laue Zone
  ALLOCATE(RgVecMagLaue(SIZE(Rhkl,DIM=1),ITotalLaueZoneLevel),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"DiffractionPatternDefinitions(",my_rank,")error allocating RgVecMagLaue"
     RETURN
  END IF

!!$     Loop through the Laue Zones (first negative: -1,-2,...,minLaueZone,
!!$     Then positive values,0,1,..,MaxLauezone) 
!!$     Identify which values in RgVecMat are quantised in Gz, ie all -1 Laue zones
!!$     have the same gz vector component in the microscope frame
!!$     ZerothlaueZoneLevel indicates the index which corresponds to the zeroth order
!!$     LaueZone

!!$     Find the magnitude from the central K-Vector in the ideal Bragg case, this is just the
!!$     square root of x and y components squared of RgVecMat (the x and y component)

  IF(RAcceptanceAngle.NE.ZERO.AND.IZOLZFLAG.EQ.0) THEN
     INumtotalReflections=0

     DO ind=1,ITotalLaueZoneLevel
        WRITE(Sind,'(I10.1)')ind
        ILaueLevel=ind-IZerothLaueZoneLevel
        RLaueZoneGz=RGzUnitVec*ILaueLevel

        DO jnd=1,SIZE(Rhkl,DIM=1)
           WRITE(Sjnd,'(I10.1)')jnd
           IF(RgVecMatT(jnd,3).GE.(RLaueZoneGz-TINY).AND. &
                RgVecMatT(jnd,3).LE.(RLaueZoneGz+TINY)) THEN
              RgVecMagLaue(jnd,ind)=SQRT((RgVecMatT(jnd,1)**2)+(RgVecMatT(jnd,2)**2))              
              IF(ind.LT.IZerothLaueZoneLevel) THEN
                 CALL Message("ReflectionDetermination", IMoreInfo+IDEBUG,IErr, &
                      MessageVariable="Negative Laue Zone Reciprocal Vector" &
                      //",RgVecMagLaue"//"("//TRIM(ADJUSTL(Sjnd))//","//TRIM(ADJUSTL(Sind))//")", &
                      RVariable=RgVecMagLaue(jnd,ind))
              ELSE
                 CALL Message("ReflectionDetermination", IMoreInfo+IDEBUG,IErr, &
                      MessageVariable="Positive Laue Zone Reciprocal Vector" &
                      //",RgVecMagLaue"//"("//TRIM(ADJUSTL(Sjnd))//","//TRIM(ADJUSTL(Sind))//")", &
                      RVariable=RgVecMagLaue(jnd,ind))
              END IF
           ELSE
              RgVecMagLaue(jnd,ind)=NEGHUGE
           END IF

           CALL Message("ReflectionDetermination", IMoreInfo+IDEBUG,IErr, &
                MessageVariable="Reciprocal Vector" &
                //",RgVecMatT"//"("//TRIM(ADJUSTL(Sjnd))//",3)", &
                RVariable=RgVecMatT(jnd,3))
        END DO

        !At each Laue Zone, determine the magnitude from Kz-Gz, ie. from the K Vector 
        !incident on the central spot on each Laue Zone Plane (under the bragg condition)
        !I am currently using the Total Laue Level, dependent on the negative Laue zones           
        INumInitReflections=COUNT(RgVecMagLaue(:,ind).NE.NEGHUGE)
        RLaueZoneElectronWaveVectorMag=RElectronWaveVectorMagnitude-ABS(RLaueZoneGz)           
        CALL Message("ReflectionDetermination", IInfo,IErr, &
             MessageVariable="For Laue Zone", IVariable=ILaueLevel)
        CALL Message("ReflectionDetermination", IInfo,IErr, &
             MessageVariable="Reduced Laue Zone K-Vector,RLaueZoneElectronWaveVectorMag", &
             RVariable=RLaueZoneElectronWaveVectorMag)
        CALL Message("ReflectionDetermination", IInfo,IErr, &
             MessageVariable="Acceptance Angle has reduced Initial Reflection Pool from", &
             IVariable=INumInitReflections)

        !Make any reflection gvec magnitude greater than the Acceptance angle gVec equal 1D9 for each level,
        !the rest are restricted by the Acceptance angle - need to change to degrees 1E-3 is 
        !for milliradians conversion- for now, we take the absolute g-vector 
        !this will probably change according to different beam selection rules
        RMaxAcceptanceGVecMag=(RLaueZoneElectronWaveVectorMag*TAN(RAcceptanceAngle*DEG2RADIAN))

        WHERE(ABS(RgVecMagLaue(:,ind)).GT.RMaxAcceptanceGVecMag)
           RgVecMagLaue(:,ind)=NEGHUGE
        END WHERE

        INumFinalReflections=COUNT(RgVecMagLaue(:,ind).NE.NEGHUGE)
        CALL Message("ReflectionDetermination", IInfo,IErr, &
             MessageVariable="To", &
             IVariable=INumFinalReflections)
        INumTotalReflections=INumTotalReflections+INumInitReflections
     END DO

     CALL Message("ReflectionDetermination", IInfo,IErr, &
          MessageVariable="INumTotalReflections", &
          IVariable=INumTotalReflections)
     CALL Message("ReflectionDetermination", IInfo,IErr, &
          MessageVariable="Size of Rhkl", &
          IVariable=SIZE(Rhkl,DIM=1))

!!$        Error check here to ensure that quantised Laue Zones are the same as the total
!!$        number of Reflections
     IF(INumTotalReflections.NE.SIZE(Rhkl,DIM=1)) THEN
        CALL ErrorChecks("Reflection Determination","Reflection Determination", &
             IPotError,IReflectionMismatch)
     END IF

!!$        Each value in each Laue Zone is unique in Size(Rhkl)
!!$        The rest are filled with NEGHUGE (-1D9), we can then find all the indexes where the 
!!$        magnitudes are still inside the acceptance angle
!!$        Firstly find how many instances there are
     ICounter=0
     DO ind=1,SIZE(Rhkl,DIM=1)
        IF(SUM(RgVecMagLaue(ind,:))/REAL(ITotalLaueZoneLevel,RKIND).GT.NEGHUGE) THEN
           ICounter=ICounter+1
        END IF
     END DO

     IHOLZGVecMagSize=ICounter

     ALLOCATE(IOriginGVecIdentifier(IHOLZGVecMagSize),STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"DiffractionPatternDefinitions(",my_rank,")error allocating IOriginGVecIdentifier"
        RETURN
     END IF

     IOriginGVecIdentifier=0
     !Store all the indexes for the magnitudes of each Laue Level in OriginGVecIdentifier    
     ICounter=1

     DO ind=1,SIZE(Rhkl,DIM=1)
        IF((SUM(RgVecMagLaue(ind,:))/REAL(ITotalLaueZoneLevel,RKIND)).GT.NEGHUGE) THEN
           IOriginGVecIdentifier(ICounter)=ind
           ICounter=ICounter+1
        END IF
     END DO
  END IF


!!$  Calculate all gvectors magnitudes (x,y,z componenents) from the origin

  DO ind=1,SIZE(Rhkl,DIM=1)
     RgVecMag(ind)= SQRT(DOT_PRODUCT(RgVecMatT(ind,:),RgVecMatT(ind,:)))
  ENDDO


!!$  If the user specifies an acceptance angle, felix uses it to restrict the number
!!$  of reflections in IMinReflectionPool & nReflections
!!$  check to ensure everything is only in magnitudes (there are negative g-vectors) 

  IF(RAcceptanceAngle.NE.ZERO.AND.IZOLZFLAG.EQ.1) THEN
     CALL Message("ReflectionDetermination",IInfo,IErr, &
          MessageString="Determining restriction effect of Acceptance Angle" &
          // " (in k-space for ZOLZ only). If unintended, please cancel the simulation and" &
          // " set the RAcceptance angle to 0.0, and/or switch IZOLZFLAG to 0")
     !Acceptance Angle 
     RMaxAcceptanceGVecMag=(RElectronWaveVectorMagnitude*TAN(RAcceptanceAngle*DEG2RADIAN))

     CALL Message("ReflectionDetermination",IInfo,IErr, &
          MessageVariable="RMaxGVecMAG",RVariable=RMaxAcceptanceGVecMag)   
     CALL Message("ReflectionDetermination",IInfo,IErr, &
          MessageVariable="RElectronWaveVector",RVariable=RElectronWaveVectorMagnitude)

     IF(RgVecMag(IMinReflectionPool).GT.RMaxAcceptanceGVecMag) THEN
        CALL Message("ReflectionDetermination",IInfo,IErr, &
             MessageString="Number of Reflections (IMinReflectionPool) Exceeds cut-off from" &
             //" Acceptance angle, calculating new cut-off value (reciprocal angstroms)")

        !New max Gvector amplitude is the G-vector specified by the acceptance angle
        RBSMaxGVecAmp = RMaxAcceptanceGVecMag 
     ELSE
        CALL Message("ReflectionDetermination",IInfo,IErr, &
             MessageString="Warning: Number of Reflections in Reflection pool does not" &
             //" exceed the Acceptance angle (reciprocal angstroms)" &
             //" continuing in normal mode, for full range increase IMinReflectionPool")

!!$     Normal Reflection Pool value
        RBSMaxGVecAmp = RgVecMag(IMinReflectionPool)
     END IF

  ELSEIF(RAcceptanceAngle.NE.ZERO.AND.IZOLZFLAG.EQ.0) THEN
     CALL Message("ReflectionDetermination",IInfo,IErr, &
          MessageString="Determining restriction effect of Acceptance Angle" &
          // " (in k-space with HOLZ). If unintended, please cancel the simulation and" &
          // " set the RAcceptance angle to 0.0 and/or switch IZOLZFLAG to 1")
     IBSMaxLocGVecAmp=MAXVAL(IOriginGVecIdentifier)
!!$     DO ind=1,ICounter
!!$        RgVecHOLZMag(ind)=RgVecMag(IOriginGVecIdentifier(ind))
!!$     END DO
     RBSMaxGVecAmp=RgVecMag(IBSMaxLocGVecAmp)
     IF(RgVecMag(IBSMaxLocGVecAmp).LT.RgVecMag(IMinreflectionPool)) THEN
        CALL Message("ReflectionDetermination",IInfo,IErr, &
             MessageString="Number of Reflections (IMinReflectionPool) Exceeds cut-off from" &
             //" Acceptance angle,")
        CALL Message("ReflectionDetermination",IInfo,IErr, &
             MessageString="using Acceptance angle cut-off value (reciprocal Angstroms)")
     ELSE
        CALL Message("ReflectionDetermination",IInfo,IErr, &
             MessageString="Warning: Number of Reflections in Reflection pool does not")
        CALL Message("ReflectionDetermination",IInfo,IErr, &
             MessageString=" exceed the Acceptance angle (reciprocal angstroms)")
        CALL Message("ReflectionDetermination",IInfo,IErr, &
             MessageString=" continuing in normal mode, for full range, increase IMinReflectionPool")
        !Normal Reflection Pool value
        RBSMaxGVecAmp = RgVecMag(IMinReflectionPool)
     END IF
  ELSE
     RBSMaxGVecAmp = RgVecMag(IMinReflectionPool)
  END IF

  nReflections = 0
  nStrongBeams = 0
  nWeakBeams = 0

  DO ind=1,SIZE(Rhkl,DIM=1)
     IF (ABS(RgVecMag(ind)).LE.RBSMaxGVecAmp) THEN
        WRITE(Sind,'(I10.1)')ind
        CALL Message("ReflectionDetermination", IAllInfo+IDEBUG,IErr, &
             MessageVariable="Allowed Reciprocal g-vector magnitude" &
             //",RgVecMag" // "("//ADJUSTL(TRIM(Sind))//")", &
             RVariable=RgVecMag(ind))
        nReflections = nReflections + 1
     END IF
  ENDDO

  CALL Message("ReflectionDetermination", IInfo,IErr, &
       MessageVariable="Reflection Pool changed from, IMinReflectionpool", &
       IVariable=IMinReflectionPool)
  CALL Message("ReflectionDetermination", IInfo,IErr, &
       MessageVariable="to nReflections",IVariable=nReflections)

!RB Deallocate local variables	   
  DEALLOCATE(RgDummyVecMat,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReflectionDetermination(",my_rank,") error deallocating RgDummyVecMat"
     RETURN
  ENDIF
END SUBROUTINE ReflectionDetermination

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
SUBROUTINE SpecificReflectionDetermination (IErr)
    
    USE MyNumbers
    USE WriteToScreen
    USE IConst

    USE IPara; USE RPara

    USE MyMPI

    IMPLICIT NONE

    INTEGER(IKIND) :: IFind,IFound,ind,jnd,knd,IErr

    CALL Message("SpecificReflectionDetermination",IMust,IErr)
    
  IFind = 0
  
  IF(IHKLSelectFLAG.EQ.1) THEN
     DO ind = 1,SIZE(RInputHKLs,DIM=1)
        DO jnd = 1,SIZE(Rhkl,DIM=1)
           IF(ABS(Rhkl(jnd,1)-RInputHKLs(ind,1)).LE.RTolerance.AND.&
                ABS(Rhkl(jnd,2)-RInputHKLs(ind,2)).LE.RTolerance.AND.&
                ABS(Rhkl(jnd,3)-RInputHKLs(ind,3)).LE.RTolerance) THEN
              IFound = 0
              DO knd = 1,IFind
                 IF(ABS(Rhkl(IOutputReflections(knd),1)-RInputHKLs(ind,1)).LE.RTolerance.AND.&
                      ABS(Rhkl(IOutputReflections(knd),2)-RInputHKLs(ind,2)).LE.RTolerance.AND.&
                      ABS(Rhkl(IOutputReflections(knd),3)-RInputHKLs(ind,3)).LE.RTolerance) THEN
                    IFound = 1
                    EXIT
                 END IF
              END DO

              IF(IFound.EQ.0) THEN
                 IFind = IFind +1
                 IOutputReflections(IFind) = jnd
                 CALL Message("SpecificReflectionDetermination",IMoreInfo,IErr, &
                      MessageVariable = "At",IVariable = jnd)
              ELSE 
                 CALL Message("SpecificReflectionDetermination",IMoreInfo,IErr, &
                      MessageVariable = "At",IVariable = jnd)
                 CALL Message("SpecificReflectionDetermination",IMoreInfo,IErr, &
                      MessageString = "However it is not unique")
                 
              END IF
              EXIT
           ELSE
              IF((jnd.EQ.SIZE(Rhkl,DIM=1).AND.IWriteFLAG.GE.3.AND.my_rank.EQ.0).or.&
                   (jnd.EQ.SIZE(Rhkl,DIM=1).AND.IWriteFLAG.GE.10)) THEN
                 PRINT*,"DiffractionPatternDefinitions(",my_rank,&
                      ") Could Not Find Requested HKL ",&
                      RInputHKLs(ind,:)," Will Ignore and Continue"
              END IF
              CYCLE
           END IF
        END DO
     END DO
     
     IF(IFind.LE.0) THEN
        IErr = 1
        IF( IErr.NE.0 ) THEN
           PRINT*,"DiffractionPatternDefinitions(", my_rank, ") error ", IErr, &
                " No requested HKLs are allowed using the purposed geometry"
           RETURN
        ENDIF
     END IF
     IF(INoOfLacbedPatterns.NE.IFind) THEN
        INoOfLacbedPatterns = IFind
     END IF
  END IF
  
END SUBROUTINE SpecificReflectionDetermination

!!$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE DiffractionPatternCalculation (IErr)
  
  USE MyNumbers
  USE WriteToScreen
  USE IConst
  
  USE IPara; USE RPara;
  
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: ind,IErr
  CHARACTER*20 :: Sind
  
  CALL Message("DiffractionPatternCalculation",IMust,IErr)
  
  RNormDirM = RNormDirM/sqrt(DOT_PRODUCT(RNormDirM,RNormDirM))
  
  DO ind =1,SIZE(Rhkl,DIM=1)
     RgVecVec(ind) = DOT_PRODUCT(RgVecMatT(ind,:),RNormDirM)
  END DO
  
  CALL Message("DiffractionPatternCalculation",IInfo,IErr, &
       MessageVariable = "GMax", RVariable = RBSMaxGVecAmp)     
  
  ! smallest g is gmag(2) IF 000 beam is included !!!add error catch here
  
  RMinimumGMag = RgVecMag(2)
  
  CALL Message("DiffractionPatternCalculation",IInfo,IErr, &
       MessageVariable = "MinimumGMag", RVariable = RMinimumGMag)
  
  IF (nReflections.LT.INoOfLacbedPatterns) THEN
     INoOfLacbedPatterns = nReflections
  END IF
  
  CALL Message("DiffractionPatternCalculation",IInfo,IErr, &
       MessageVariable = "No. of Reflections", IVariable = nReflections)
  
  ! resolution in k space
  RDeltaK = RMinimumGMag*RConvergenceAngle/REAL(IPixelCount,RKIND)

  CALL Message("DiffractionPatternCalculation",IInfo,IErr, &
       MessageVariable = "RDeltaK", RVariable = RDeltaK)

  RETURN

END SUBROUTINE DiffractionPatternCalculation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE NewHKLMake(Ihklmax,Rhkl0Vec,RHOLZAcceptanceAngle,IErr)
  
  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE IChannels
  
  USE MyMPI
  USE MPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr, Ihklmax,ind,jnd,knd,INhkl
  REAL(RKIND) :: RHOLZAcceptanceAngle
  REAL(RKIND), DIMENSION(THREEDIM) :: Rhkl0Vec,RhklDummyUnitVec,RhklDummyVec,Rhkl0UnitVec

  CALL Message("NewHKLMake",IMust,IErr)

  INhkl = 0
  
  Rhkl0UnitVec= Rhkl0Vec/SQRT(DOT_PRODUCT(REAL(Rhkl0Vec,RKIND),REAL(Rhkl0Vec,RKIND)))

!RB first count the number of reflections in the acceptance angle
  DO ind=-Ihklmax,Ihklmax,1
     DO jnd=-Ihklmax,Ihklmax,1
        DO knd=-Ihklmax,Ihklmax,1
           
           RhklDummyVec= REAL((/ ind,jnd,knd /),RKIND)
           
           IF(ind.NE.0.AND.jnd.NE.0.AND.knd.NE.0) THEN
              RhklDummyUnitVec= RhklDummyVec / &
                   SQRT(DOT_PRODUCT(REAL(RhklDummyVec,RKIND),REAL(RhklDummyVec,RKIND)))
           ELSE
              RhklDummyUnitVec = RhklDummyVec
           END IF
           
           SELECT CASE(SSpaceGroupName)
              
           CASE("F") !Face Centred
              IF(((ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),TWO)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),TWO))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),TWO))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),TWO))).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),TWO))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),TWO))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),TWO))).GT.TINY))) THEN
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
			  
           CASE("I")! Body Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
			  
           CASE("A")! A-Face Centred
              IF(ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
			  
           CASE("B")! B-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1
                 ENDIF
              END IF
			  
           CASE("C")! C-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),TWO)).LE.TINY) THEN
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
			  
           CASE("R")! Rhombohedral Reverse
              IF(ABS(MOD(RhklDummyVec(1)-RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
			  
           CASE("V")! Rhombohedral Obverse
              IF(ABS(MOD(-RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                      .LE.SIN(RHOLZAcceptanceAngle)) THEN
                    INhkl = INhkl +1       
                 ENDIF
              END IF
			  
           CASE("P")! Primitive
              IF(IZolzFLAG.EQ.1) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    INhkl=INhkl+1
                 END IF
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) &
                   .LE.SIN(RHOLZAcceptanceAngle)) THEN
                 INhkl = INhkl +1       
              ENDIF
           CASE DEFAULT
              PRINT*,"HKLMake(): unknown space group", SSpaceGroupName, "--- aborting"
              IErr=1
              RETURN
           END SELECT

        END DO
     END DO
  END DO

!RB now allocate the hkl list...  
  ALLOCATE(Rhkl((INhkl),THREEDIM),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"hklMake(",my_rank,")error allocating Rhkl"
     RETURN
  ENDIF

!RB ...and calculate it all again, filling Rhkl  
  INhkl = 0

  DO ind=-Ihklmax,Ihklmax,1
     DO jnd=-Ihklmax,Ihklmax,1
        DO knd=-Ihklmax,Ihklmax,1

           RhklDummyVec= REAL((/ ind,jnd,knd /),RKIND)

           IF(ind.NE.0.AND.jnd.NE.0.AND.knd.NE.0) THEN
              RhklDummyUnitVec= RhklDummyVec / &
                   SQRT(DOT_PRODUCT(REAL(RhklDummyVec,RKIND),REAL(RhklDummyVec,RKIND)))
           ELSE
              RhklDummyUnitVec = RhklDummyVec
           END IF
           
           SELECT CASE(SSpaceGroupName)
		   
           CASE("F") !Face Centred
              IF(((ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),TWO)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY).AND.&
                   (ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),TWO))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),TWO))).LE.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),TWO))).LE.TINY)).OR.&
                   (((ABS(MOD(RhklDummyVec(1),TWO))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(2),TWO))).GT.TINY).AND.&
                   ((ABS(MOD(RhklDummyVec(3),TWO))).GT.TINY))) THEN
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("I")! Body Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("A")! A-Face Centred
              IF(ABS(MOD(RhklDummyVec(2)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("B")! B-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(3),TWO)).LE.TINY) THEN
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("C")! C-Face Centred
              IF(ABS(MOD(RhklDummyVec(1)+RhklDummyVec(2),TWO)).LE.TINY) THEN
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("R")! Rhombohedral Reverse
              IF(ABS(MOD(RhklDummyVec(1)-RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                 IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("V")! Rhombohedral Obverse
              IF(ABS(MOD(-RhklDummyVec(1)+RhklDummyVec(2)+RhklDummyVec(3),THREE)).LE.TINY) THEN
                IF(IZolzFLAG.EQ.1) THEN
                    IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                       INhkl=INhkl+1
                       Rhkl(INhkl,:)= RhklDummyVec
                    END IF
                 ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                    INhkl =  INhkl + 1
                    Rhkl(INhkl,:) = RhklDummyVec                 
                 END IF
              END IF
			  
           CASE("P")! Primitive

		   IF(IZolzFLAG.EQ.1) THEN
                 IF( ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)) .LE. TINY ) THEN
                    INhkl=INhkl+1
                    Rhkl(INhkl,:)= RhklDummyVec
                 END IF
              ELSEIF (ABS(DOT_PRODUCT(RhklDummyUnitVec,Rhkl0UnitVec)).LE.sin(RHOLZAcceptanceAngle)) THEN
                 INhkl =  INhkl + 1
                 Rhkl(INhkl,:) = RhklDummyVec                 
              END IF
           CASE DEFAULT
              PRINT*,"HKLMake(): unknown space group", SSpaceGroupName, "--- aborting"
              IErr=1
              RETURN
           END SELECT

        END DO
     END DO
  END DO

END SUBROUTINE NewHKLmake
