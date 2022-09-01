!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-19, all rights reserved
!
! Version: 2.0
! Date: 31-08-2022 
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build: cRED
! Author:  r.beanland@warwick.ac.uk
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

!>
!! Module-description: Holds main top level simulating SUBROUTINEs considering
!! different thicknesses
!!
MODULE refinementcontrol_mod
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Simulate

  CONTAINS
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Simulates and produces images for each thickness
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!  
  SUBROUTINE Simulate(IErr)

    USE MyNumbers
    USE IConst, ONLY : ITHREE
    USE MyMPI
    USE message_mod

    USE bloch_mod

    !global outputs
    USE RPara, ONLY : RImageSimi, RSimulatedPatterns
    ! RImageSimi(x_coordinate, y_coordinate y, LACBED_pattern_ID , thickness_ID )
    ! RSimulatedPatterns( Pixel_ID, LACBED_pattern_ID , thickness_ID )
    ! RSimulatedPatterns has a long list of pixel instead of a 2D image matrix
    USE RPARA, ONLY : RIndividualReflections
    USE IPara, ONLY : IInitialSimulationFLAG, IPixelComputed
    
    !global inputs
    USE RPARA, ONLY : RBlurRadius
    USE IPARA, ONLY : ICount,IDisplacements,ILocalPixelCountMax,INoOfLacbedPatterns,&
          ILocalPixelCountMin,IPixelLocations,ISizeX,ISizeY,IThicknessCount, nBeams

    IMPLICIT NONE

    INTEGER(IKIND) :: IErr, ind,jnd,knd,pnd,IIterationFLAG,IYPixelIndex,IXPixelIndex
    REAL(RKIND) :: RKn,RThickness

    ! Reset simuation   
    RIndividualReflections = ZERO

    ! Simulation (different local pixels for each core)
    CALL message(LS,"Bloch wave calculation...")
    DO knd = ILocalPixelCountMin,ILocalPixelCountMax,1
      IYPixelIndex = IPixelLocations(knd,1)
      IXPixelIndex = IPixelLocations(knd,2)
      ! fills array for each pixel number knd (x & y coordinates are IXPixelIndex & IYPixelIndex)
      CALL BlochCoefficientCalculation(IYPixelIndex,IXPixelIndex,knd, &
              ILocalPixelCountMin, nBeams, RThickness,RKn, IErr)
      IF(l_alert(IErr,"Simulate","BlochCoefficientCalculation")) RETURN
    END DO

    !===================================== ! MPI gatherv into RSimulatedPatterns
    CALL MPI_GATHERV(RIndividualReflections,SIZE(RIndividualReflections),MPI_DOUBLE_PRECISION,&
         RSimulatedPatterns,ICount,IDisplacements,MPI_DOUBLE_PRECISION,&
         root,MPI_COMM_WORLD,IErr)
    !=====================================
    IF(l_alert(IErr,"SimulateAndFit","MPI_GATHERV")) RETURN

    ! put 1D array RSimulatedPatterns into 2D image RImageSimi
    ! remember dimensions of RSimulatedPatterns(INoOfLacbedPatterns,IThicknessCount,IPixelTotal)
    ! and RImageSimi(height, width, INoOfLacbedPatterns,IThicknessCount )
    RImageSimi = ZERO
    ind = 0
    DO IYPixelIndex = 1,ISizeY
      DO IXPixelIndex = 1,ISizeX
        ind = ind+1
        RImageSimi(IYPixelIndex,IXPixelIndex,:,:) = RSimulatedPatterns(:,:,ind)
      END DO
    END DO
    ! Gaussian blur to match experiment using global variable RBlurRadius
    IF (RBlurRadius.GT.TINY) THEN
      DO ind=1,INoOfLacbedPatterns
        DO jnd=1,IThicknessCount
          CALL BlurG(RImageSimi(:,:,ind,jnd),ISizeX,ISizeY,RBlurRadius,IErr)
        END DO
      END DO
    END IF

  END SUBROUTINE Simulate

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !>
  !! Procedure-description: Performs a 2D Gaussian blur on the input image using 
  !! global variable RBlurRadius and renormalises the output image to have the
  !! same min and max as the input image
  !!
  !! Closed procedure, no access to global variables
  !!
  !! Major-Authors: Richard Beanland (2016)
  !! 
  SUBROUTINE BlurG(RImageToBlur,IPixX,IPixY,RBlurringRadius,IErr)

    USE MyNumbers
    USE MPI
    USE message_mod

    IMPLICIT NONE

    REAL(RKIND),DIMENSION(IPixX,IPixY),INTENT(INOUT) :: RImageToBlur
    INTEGER(IKIND),INTENT(IN) :: IPixX,IPixY
    REAL(RKIND),INTENT(IN) :: RBlurringRadius
    INTEGER(IKIND),INTENT(OUT) :: IErr

    REAL(RKIND),DIMENSION(IPixX,IPixY) :: RTempImage,RShiftImage
    INTEGER(IKIND) :: ind,jnd,IKernelRadius,IKernelSize
    REAL(RKIND),DIMENSION(:), ALLOCATABLE :: RGauss1D
    REAL(RKIND) :: Rind,Rsum,Rmin,Rmax

    ! get min and max of input image
    Rmin=MINVAL(RImageToBlur)
    Rmax=MAXVAL(RImageToBlur)

    ! set up a 1D kernel of appropriate size  
    IKernelRadius=NINT(3*RBlurringRadius)
    ALLOCATE(RGauss1D(2*IKernelRadius+1),STAT=IErr)!ffs
    Rsum=0
    DO ind=-IKernelRadius,IKernelRadius
      Rind=REAL(ind)
      RGauss1D(ind+IKernelRadius+1)=EXP(-(Rind**2)/(2*(RBlurringRadius**2)))
      Rsum=Rsum+RGauss1D(ind+IKernelRadius+1)
      IF(ind==0) IErr=78 
    END DO
    RGauss1D=RGauss1D/Rsum!normalise
    RTempImage=RImageToBlur*0_RKIND !reset the temp image 

    ! apply the kernel in direction 1
    DO ind = -IKernelRadius,IKernelRadius
       IF (ind.LT.0) THEN
          RShiftImage(1:IPixX+ind,:)=RImageToBlur(1-ind:IPixX,:)
          DO jnd = 1,1-ind!edge fill on right
             RShiftImage(IPixX-jnd+1,:)=RImageToBlur(IPixX,:)
          END DO
       ELSE
          RShiftImage(1+ind:IPixX,:)=RImageToBlur(1:IPixX-ind,:)
          DO jnd = 1,1+ind!edge fill on left
             RShiftImage(jnd,:)=RImageToBlur(1,:)
          END DO
       END IF
       RTempImage=RTempImage+RShiftImage*RGauss1D(ind+IKernelRadius+1)
    END DO

    ! make the 1D blurred image the input for the next direction
    RImageToBlur=RTempImage
    RTempImage=RImageToBlur*0_RKIND ! reset the temp image

    ! apply the kernel in direction 2  
    DO ind = -IKernelRadius,IKernelRadius
       IF (ind.LT.0) THEN
          RShiftImage(:,1:IPixY+ind)=RImageToBlur(:,1-ind:IPixY)
          DO jnd = 1,1-ind!edge fill on bottom
             RShiftImage(:,IPixY-jnd+1)=RImageToBlur(:,IPixY)
          END DO
       ELSE
          RShiftImage(:,1+ind:IPixY)=RImageToBlur(:,1:IPixY-ind)
          DO jnd = 1,1+ind!edge fill on top
             RShiftImage(:,jnd)=RImageToBlur(:,1)
          END DO
       END IF
       RTempImage=RTempImage+RShiftImage*RGauss1D(ind+IKernelRadius+1)
    END DO
    DEALLOCATE(RGauss1D,STAT=IErr)

    ! set intensity range of output image to match that of the input image
    RTempImage=RTempImage-MINVAL(RTempImage)
    RTempImage=RTempImage*(Rmax-Rmin)/MAXVAL(RTempImage)+Rmin
    ! return the blurred image
    RImageToBlur=RTempImage;

  END SUBROUTINE BlurG        

END MODULE refinementcontrol_mod
