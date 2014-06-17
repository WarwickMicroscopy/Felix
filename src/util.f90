!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! FelixSim
!
! Richard Beanland, Keith Evans and Rudolf A Roemer
!
! (C) 2013/14, all rights reserved
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  This file is part of FelixSim.
!
!  FelixSim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  FelixSim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with FelixSim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: util.f90,v 1.95 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Log: util.f90,v $
! Revision 1.95  2014/04/28 12:26:19  phslaz
! Fixed minor bugs with the new reflection pool request
!
! Revision 1.93  2014/04/24 10:48:50  phslaz
! Resorthkl corrected to sort by gmax and not dot(hkl,hkl), language used in resorthkl needs attention, continuous changing of integers is horrible. Codes makes and runs
!
! Revision 1.90  2014/04/09 13:45:39  phslaz
! cleaned up the write flags also added in some of the amplitude/phase imaging
!
! Revision 1.89  2014/03/27 21:39:14  phslaz
! Added two new flags IImageFlag and IOutputFlag
! IImageFLAG = 0 (Montage) 1 (Montage and Reflections)
! IOutputFLAG = 0 (nothing) 1 (EigenSpectra) 2 (UgMat) 3 (Wavefunctions)
! Have also put many Print statments into IWriteflag control
! code compiles and runs
!
! Revision 1.88  2014/03/27 18:13:43  phslaz
! MPI_Reduce attempt, compiles but fails to fun
!
! Revision 1.87  2014/03/26 17:04:52  phslaz
! Felix now creates images
!
! Revision 1.84  2014/03/25 15:45:31  phslaz
! conflict resolution
!
! Revision 1.83  2014/03/25 15:35:34  phsht
! included consistent start of each source file and GPL statements
!
! Revision 1.82  2014/03/25 14:03:39  phslaz
! fixed allocation issue with FelixCoefficients subroutine, allocated memory wasnt being deallocated within loop
!
! Revision 1.81  2014/03/24 12:50:31  phslaz
! fixed write flag 2 and rewrote lacbed to read file once instead of every thickness
!
! Revision 1.80  2014/03/21 15:55:36  phslaz
! New Lacbed code Working
!
! Revision 1.79  2014/03/13 18:10:32  phslaz
! Seg Fault due to 32bit suspected constraint
!
! Revision 1.78  2014/03/07 10:49:45  phslaz
! Corrected issues with inpcif, should now work with badly structured cifs
!
! Revision 1.76  2014/02/21 15:26:42  phslaz
! collapsed a few parts of main.f90 into subroutines in util.f90
!
! Revision 1.75  2014/02/07 14:33:05  phslaz
! LACBED code now reads eigen spectra output
!
! Revision 1.74  2014/01/23 18:50:55  phslaz
! Absorption and 2Pi convention installed and checked with matlab
!
! Revision 1.73  2014/01/20 18:33:59  phslaz
! Isotropic Debye Waller Factor and Atomic Site Partial Occupancy done
!
! Revision 1.71  2014/01/17 15:43:10  phslaz
! Bug Fix : Periodic Table now compiles
!
! Revision 1.69  2014/01/16 16:12:42  phsht
! work on scattering factors
!
! Revision 1.68  2014/01/13 13:46:20  phslaz
! *** empty log message ***
!
! Revision 1.65  2014/01/07 17:11:52  phslaz
! Bug Fix : Hklmake was being called before the basevector assignment and as such, being called with lots of zeros...
!
! Revision 1.64  2014/01/07 11:49:17  phsht
! transformation to microscope reference included, needs testing
!
! Revision 1.63  2013/12/19 16:30:28  phsht
! new version of HKLMake(), using .cif information
!
! Revision 1.62  2013/12/19 14:58:57  phsht
! symmetry operations now correctly interpreted from .cif;
! structure correctly read in
!
! Revision 1.61  2013/12/17 17:40:53  phsht
! make inpcif.f90 which now seems to work
!
! Revision 1.60  2013/11/27 12:30:21  phsht
! more consistency in MPI output routines and their error handling
!
! Revision 1.59  2013/10/21 15:56:31  phslaz
! Changed Variable names
!
! Revision 1.58  2013/10/03 15:48:57  phsht
! added routine HKLMake() in util.f90 to make BWM compatible;
! checked for 64 pixels with input file of this version

! Revision 1.57  2013/10/03 12:52:25  phsht
! two errors corrected for beam selection parameter
!
! Revision 1.56  2013/10/03 11:15:16  phsht
! new input file structure
!
! Revision 1.55  2013/10/02 20:43:28  phsht
! replaced old input names with new names;
! structure of input file still needs changing
!
! Revision 1.54  2013/09/24 16:29:10  phslaz
! Ugmateffective Fixed Now, images now match matlab :)
!
! Revision 1.53  2013/09/24 15:10:07  phslaz
! Added Near Zero Checks
!
! Revision 1.52  2013/09/19 11:19:33  phsht
! eigenvector calculation now working
!
! Revision 1.51  2013/09/18 16:08:42  phsht
! work with Keiths to include changes done with Richard yesterday and
! also to continue until intensities; reached these, but don't quite agree
! yet. Nevertheless, up to the diagonalization, things do work! And eigenvalues
! are also correct, just eigenvectors not yet.
!
! Revision 1.50  2013/09/10 16:51:54  phsht
! before implementation of UgMatEffective
!
! Revision 1.49  2013/09/09 14:09:46  phsht
! works up to the potential
!
! Revision 1.48  2013/09/09 10:58:20  phsht
! subroutines up to this stage seems to produce same output as matlab code
!
! Revision 1.47  2013/09/05 15:20:55  phslaz
! Bug Fix : Didnt declare kstep
!
! Revision 1.46  2013/09/05 15:19:15  phslaz
! Added calculation of kstep
!
! Revision 1.45  2013/09/05 15:08:23  phslaz
! hklpositions is causing issues and isnt massive important right now so ive commented it out for the time being
!
! Revision 1.44  2013/09/05 14:54:46  phslaz
! Bug Fix : Was getting an Unclassifed Function Error from hklpositions, hopefully this will fix it
!
! Revision 1.43  2013/09/05 14:43:39  phslaz
! Added Calculation of positions of reflections in the final image
!
! Revision 1.42  2013/09/05 14:22:48  phslaz
! Beam Selection Works Perfectly
!
! Revision 1.41  2013/09/05 14:21:24  phslaz
! MIssed a comment
!
! Revision 1.40  2013/09/05 14:16:19  phslaz
! Bug Fix : Forgot to uncomment the error
!
! Revision 1.39  2013/09/05 14:04:20  phslaz
! Debugging : Seems to the Sg Calculation Causing the issue, have commented
!
! Revision 1.38  2013/09/05 14:02:56  phslaz
! Debugging : Added Print Lines
!
! Revision 1.37  2013/09/05 11:10:55  phslaz
! Bug Fix : Getting a Segementatiion Fault from beam selection, trying to track it
!
! Revision 1.36  2013/09/05 11:07:40  phslaz
! Bug Fix : forgot an =
!
! Revision 1.35  2013/09/05 11:05:31  phslaz
! Bug Fix i forgot a THEN, i thknk thats the issue
!
! Revision 1.34  2013/09/05 11:01:38  phslaz
! Added Beam Selection Criteria
!
! Revision 1.33  2013/09/05 09:07:33  phslaz
! Resorthkl being called with size(hkl) = 7101 which is the size of the array not the number of reflections (which is 1/3 that)
!
! Revision 1.32  2013/09/05 08:50:38  phslaz
! There was a second call of ReSortHKL but was calling with iErr not Size(HKL) changed it, we'll see what happens
!
! Revision 1.31  2013/09/05 08:46:21  phslaz
! ResortHKL still doing something funny, have added some print lines
!
! Revision 1.30  2013/09/05 08:32:37  phslaz
! Loop wasnt looping over the whole hkl matrix
!
! Revision 1.29  2013/09/05 08:06:10  phsht
! ReSortHKL() now sorts smaller entries first
!
! Revision 1.28  2013/09/05 00:05:03  phslaz
! Resort is backwards, need to invert
!
! Revision 1.27  2013/09/04 23:50:05  phslaz
! Messing around with variables in ILength and ind
!
! Revision 1.26  2013/09/04 23:44:55  phslaz
! Sorting seems to work backwards :)
!
! Revision 1.25  2013/09/04 19:21:46  phsht
! corrected calling of ReSortHKL(); but not sure that result is ok
!
! Revision 1.24  2013/09/04 16:16:06  phslaz
! Changed HKLarray from (N,THREEDIMENSION) to HKLarray(ILength,THREEDIMENSION)
!
! Revision 1.23  2013/09/04 16:12:55  phslaz
! I think ive broke resortHKl
!
! Revision 1.22  2013/09/04 16:09:50  phslaz
! Yeah my previous update broke everything, ive commented out the printstatement
!
! Revision 1.21  2013/09/04 16:07:49  phslaz
! Added Some Prints in, still trying to figure out how the resortHKL works, should it be called with Size(HKL)?
!
! Revision 1.20  2013/09/04 15:38:47  phslaz
! Bug Fix : Trying to diagnose the ResortHKL subroutine, see if theres anything wrong
!
! Revision 1.19  2013/09/04 15:22:30  phslaz
!
! Added a PRINT line in resorthkl trying to bugfix
!
! Revision 1.18  2013/09/04 15:16:57  phslaz
! Removed the definition of gVecMag here
!
! Revision 1.17  2013/09/04 14:59:35  phslaz
! Attempting to Print out HKL after sort
!
! Revision 1.16  2013/09/04 14:52:13  phslaz
! Added A print for minimumG
!
! Revision 1.15  2013/09/04 14:27:16  phslaz
! i wrote int instead of ind, because im an idiot
!
! Revision 1.14  2013/09/04 14:24:18  phslaz
! Modifed a line to output 8 g vector magnitudes not all of them (for debuging)
!
! Revision 1.13  2013/09/03 16:31:33  phslaz
! Bug Fix : Values looked good but i think its running on the unsorted HKL so ive called it prior to the function
!
! Revision 1.12  2013/09/03 16:23:34  phslaz
! Added gVecMag output (to check functionality)
!
! Revision 1.11  2013/09/03 16:20:03  phslaz
! Bug Fix : HopeFully fixed gVecMag calculation
!
! Revision 1.10  2013/09/03 16:10:49  phslaz
! Uncommented gVecMag calculation and added real variable gvecmag(Size(HKL,DIM=1)) (not sure if this will work)
!
! Revision 1.9  2013/09/03 16:05:07  phslaz
! Bug Fix : End - Endif
!
! Revision 1.8  2013/09/03 16:03:05  phslaz
! Added "minimum G Mag"
!
! Revision 1.7  2013/09/02 15:42:42  phsht
! code checked with matlab version up to and including BraggCentral
!
! Revision 1.6  2013/08/29 11:25:15  phsht
! small change
!
! Revision 1.5  2013/06/11 14:53:08  phsht
! more work in the Diffraction defs part
!
! Revision 1.4  2013/06/10 15:08:02  phsht
! work this morning
!
! Revision 1.3  2013/06/10 08:20:28  phsht
! more work before realizing that MATLAB file is not quite the latest version
!
! Revision 1.2  2013/06/07 07:11:28  phsht
! more ground work
!
! Revision 1.1  2013/04/03 19:35:45  phsht
! first installation of basic Fortran routines/structure
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!--------------------------------------------------------------------
!	ReSort:
!
!	sort the Lyapunov eigenvalues s.t. the largest comes first. RESORT()
!	is based on ShellSort from "Numerical Recipes", routine SHELL().
!---------------------------------------------------------------------

SUBROUTINE ReSortHKL( RHKLarray, N )

  USE MyNumbers
    
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER N
  REAL(RKIND) RHKLarray(N,THREEDIM)
  REAL(RKIND) RhklarraySearch(THREEDIM), RhklarrayCompare(THREEDIM)
  
  REAL(KIND=RKIND) ALN2I, LocalTINY
  PARAMETER (ALN2I=1.4426950D0, LocalTINY=1.D-5)
  
  INTEGER NN,M,L,K,J,I,LOGNB2, index
  REAL(KIND=RKIND) dummy

  IF((IWriteFLAG.EQ.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"ReSort()"
  END IF
  
  LOGNB2=INT(LOG(REAL(N))*ALN2I+LocalTINY)
  M=N
  DO 12 NN=1,LOGNB2
     M=M/2
     K=N-M
     DO 11 J=1,K
        I=J
3       CONTINUE
        L=I+M
        RhklarraySearch = Rhklarray(L,1)*RarVecO + &
             Rhklarray(L,2)*RbrVecO + &
             Rhklarray(L,3)*RcrVecO    
        RhklarrayCompare = Rhklarray(I,1)*RarVecO + &
             Rhklarray(I,2)*RbrVecO + &
             Rhklarray(I,3)*RcrVecO
        IF( &
             DOT_PRODUCT(RHKLarraySearch(:),RHKLarraySearch(:)) .LT. &
             DOT_PRODUCT(RHKLarrayCompare(:),RHKLarrayCompare(:))) THEN
!!$           
!!$             DOT_PRODUCT(RHKLarray(L,:),RHKLarray(L,:)) .LT. &
!!$             DOT_PRODUCT(RHKLarray(I,:),RHKLarray(I,:))) THEN
!!$           
           DO 100 index=1,THREEDIM
              dummy        = RHKLarray(I,index)
              RHKLarray(I,index)= RHKLarray(L,index)
              RHKLarray(L,index)= dummy
100        ENDDO
           
           I=I-M
           IF(I.GE.1) GOTO 3
        ENDIF
11   ENDDO
12 ENDDO
  
  !PRINT*,"Finishing ResortHKL"

  !	PRINT*,"array0(1),array0(N)",array0(1),array0(N)
  RETURN

END SUBROUTINE ReSortHKL

!---------------------------------------------------------------------
SUBROUTINE CONVERTAtomName2Number(name, number, IErr)

  USE IPara
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER IErr, ind, number
  CHARACTER*2 name

  INTEGER, PARAMETER :: NElements=103

  CHARACTER*2 A(NElements)

  DATA A/" H", "He", "Li", "Be", " B", " C", " N", "O", "F", "Ne", &
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", &
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", &
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", &
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", &
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", &
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", &
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",& 
        "Md","No","Lr"/

  DO ind=1,NElements
     IF(TRIM(name)==TRIM(A(ind))) THEN
        number= ind
        IF((IWriteFLAG.EQ.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
           PRINT*,"DBG: name, number ", name, number
        END IF
           RETURN
     ENDIF
  ENDDO

  PRINT*,"CONVERTAtomName2Number(): could not find index for atom of name ", name
  IErr=1
  RETURN

  PRINT*,"DBG: name, number ", name, number
END SUBROUTINE CONVERTAtomName2Number

!---------------------------------------------------------------------
SUBROUTINE ImageMask (IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER ind,jnd, ierr
  REAL(RKIND) :: Rradius, RImageRadius
  
  IF((IWriteFLAG.EQ.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"DBG: ImageMask()"
  END IF

  IPixelTotal =0
  SELECT CASE (IMaskFLAG)
  CASE(0) ! circle
     DO ind=1,2*IPixelCount
        DO jnd=1,2*IPixelCount
           Rradius= (ind-(REAL(IPixelCount,RKIND)+0.5))**2 + &
                (jnd-(REAL(IPixelCount,RKIND)+0.5))**2
           Rradius=SQRT(DBLE(Rradius))
           RImageRadius = IPixelCount+0.5
           IF(Rradius.LE.RImageRadius) THEN
              RMask(jnd,ind) = 1
              IPixelTotal= IPixelTotal + 1
              
           ELSE
              RMask(jnd,ind) = 0
           ENDIF
        ENDDO
     ENDDO
  CASE(1) ! square
     RMask = 1
     IPixelTotal = (2*IPixelCount)**2
  END SELECT
  
  ALLOCATE( &
       IPixelLocations(IPixelTotal,2), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Imagemask(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  IPixelTotal = 0
 
  SELECT CASE (IMaskFLAG)
  CASE(0) ! circle
     DO ind=1,2*IPixelCount
        DO jnd=1,2*IPixelCount
           Rradius= (ind-(REAL(IPixelCount,RKIND)+0.5))**2 + &
                (jnd-(REAL(IPixelCount,RKIND)+0.5))**2
           Rradius=SQRT(DBLE(Rradius))
           RImageRadius = IPixelCount+0.5
           IF(Rradius.LE.RImageRadius) THEN
              !RMask(jnd,ind) = 1
              IPixelTotal= IPixelTotal + 1
              IPixelLocations(IPixelTotal,1) = ind
              IPixelLocations(IPixelTotal,2) = jnd
           ELSE
              !RMask(jnd,ind) = 0
           ENDIF
        ENDDO
     ENDDO
  CASE(1) ! square
     IPixelTotal = 0
     DO ind = 1,2*IPixelCount
        DO jnd = 1,2*IPixelCount
           IPixelTotal = IPixelTotal+1
           IPixelLocations(IPixelTotal,1) = ind
           IPixelLocations(IPixelTotal,2) = jnd
        END DO
     END DO
  END SELECT
  
END SUBROUTINE ImageMask

!---------------------------------------------------------------------
SUBROUTINE CountTotalAtoms(IErr)

  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER ind,jnd,knd,hnd,ierr, ifullind, iuniind
  LOGICAL Lunique

  
  IF((IWriteFLAG.EQ.6.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"CountTotalAtoms()"
  END IF

  ALLOCATE( &
       RFullAtomicFracCoordVec( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1),&
       THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  ALLOCATE( &
       SFullAtomicNameVec( &
       SIZE(RSymVec,1)*SIZE(RAtomSiteFracCoordVec,1)), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  
  RFullAtomicFracCoordVec = ZERO
  
  IF((IWriteFLAG.GE.10.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"SIZE OF RFULLATOMICFRACCOORDVEC = ",SIZE(RFullAtomicFracCoordVec,1)
  END IF

  DO ind=1, SIZE(RSymVec,DIM=1)
     
     DO jnd=1, SIZE(RAtomSiteFracCoordVec,DIM=1)
       
        Ifullind= SIZE(RSymVec,1)*(jnd-1) + ind
        
        RFullAtomicFracCoordVec(Ifullind,:)= &
             MATMUL(RSymMat(ind,:,:),RAtomSiteFracCoordVec(jnd,:)) &
             + RSymVec(ind,:)
        SFullAtomicNameVec(Ifullind) = SAtomName(jnd)
        
        ! renormalize such that all values are non-negative
        DO knd=1,THREEDIM
           IF( RFullAtomicFracCoordVec(Ifullind,knd) .LT. ZERO) THEN
              RFullAtomicFracCoordVec(Ifullind,knd)= &
                   RFullAtomicFracCoordVec(Ifullind,knd)+1.D0
           ENDIF
        ENDDO
        
     ENDDO
     
  ENDDO

  DO ind = 1,SIZE(RFullAtomicFracCoordVec,DIM=1)
     DO jnd = 1,SIZE(RFullAtomicFracCoordVec,DIM=2)
        IF (RFullAtomicFracCoordVec(ind,jnd).LT.ZERO) THEN
           RFullAtomicFracCoordVec(ind,jnd) = &
                RFullAtomicFracCoordVec(ind,jnd) + ONE
        END IF
        IF (RFullAtomicFracCoordVec(ind,jnd).GE.ONE) THEN
           RFullAtomicFracCoordVec(ind,jnd) = &
                RFullAtomicFracCoordVec(ind,jnd) - ONE
        END IF
     END DO
  END DO

  
  ! Calculate the set of unique fractional atomic positions
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"CountTotalAtoms(",my_rank,") ITotalAtoms = ",ITotalAtoms
  END IF

  ALLOCATE( &
       MNP(1000,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  ALLOCATE( &
       SMNP(1000), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  
  MNP = ZERO
  
  MNP(1,:)= RFullAtomicFracCoordVec(1,:)
  SMNP(1)= SFullAtomicNameVec(1)
  
  Iuniind = 1
  
  IF(ITotalAtoms.EQ.0)THEN
     
     DO ind=2,SIZE(RFullAtomicFracCoordVec,1)
        
        DO jnd=1,Iuniind
           
           IF ( RFullAtomicFracCoordVec(ind,1) .EQ. MNP(jnd,1) .AND. &
                RFullAtomicFracCoordVec(ind,2) .EQ. MNP(jnd,2) .AND. &
                RFullAtomicFracCoordVec(ind,3) .EQ. MNP(jnd,3) .AND. &
                SFullAtomicNameVec(ind) .EQ. SMNP(jnd) ) THEN
              !this seems NOT a unique coordinate
              Lunique=.FALSE.
              EXIT
           ENDIF
           Lunique=.TRUE.
        ENDDO
        
        IF(Lunique .EQV. .TRUE.) THEN
           Iuniind=Iuniind+1
           MNP(Iuniind,:)= RFullAtomicFracCoordVec(ind,:)
           SMNP(Iuniind)= SFullAtomicNameVec(ind)
        ENDIF
        
     ENDDO
     
     ITotalAtoms = Iuniind
     
  END IF
  
  IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"CountTotalAtoms(",my_rank,") ITotalAtoms = ",ITotalAtoms
  END IF

  DEALLOCATE( &
       MNP,SMNP, &
       RFullAtomicFracCoordVec, &
       SFullAtomicNameVec,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"CountTotalAtoms(", my_rank, ") error ", IErr, &
          " in Deallocation"
     RETURN
  ENDIF
  
END SUBROUTINE CountTotalAtoms

SUBROUTINE AngularOffset (RImageSim,RImageExpi,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER ind,jnd, ierr
  REAL(RKIND),DIMENSION(IImageSizeXY(1),IImageSizeXY(2)) :: &
       RImageSim,RImageExpi  

END SUBROUTINE AngularOffset

SUBROUTINE CorrectCentreOffset (RUncorrImage,RCorrImage,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER ind,jnd, ierr
  REAL(RKIND),DIMENSION(IImageSizeXY(1),IImageSizeXY(2)) :: &
       RUnCorrImage
  REAL(RKIND),DIMENSION(IImageSizeXY(1)-IOffset(1),IImageSizeXY(1)-IOffset(1)),INTENT(OUT) :: &
       RCorrImage

  RCorrImage = RUncorrImage((IOffset(1)+1):,(IOffset(2)+1):) 
  
END SUBROUTINE CorrectCentreOffset

SUBROUTINE PhaseCorrelate(RImageSim,RImageExpi,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara
  USE IChannels
  
  USE MPI
  USE MyMPI
  USE MyFFTW

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr

  REAL(RKIND),DIMENSION(IImageSizeXY(1),IImageSizeXY(2)) :: &
       RImageExpi,RImageSim
  
  type(C_PTR) :: &
       Iplan

  real(C_DOUBLE), pointer :: &
       RImageSimDummy(:,:)
  type(C_PTR) :: & 
       p1,p2,p3,p4
  complex(C_DOUBLE_COMPLEX), pointer :: &
       CDummy1(:,:),CDummy2(:,:),CCorrelatedImage(:,:)
  integer(C_INT) :: IX,IY
  
  IX = IImageSizeXY(1)
  IY = IImageSizeXY(2)

  !PRINT*,"IX, IY =",IX,IY

  p1 = fftw_alloc_real(INT(IImageSizeXY(1)*IImageSizeXY(2), C_SIZE_T))
  p2 = fftw_alloc_complex(INT(IImageSizeXY(1)*IImageSizeXY(2), C_SIZE_T))
  p3 = fftw_alloc_complex(INT(IImageSizeXY(1)*IImageSizeXY(2), C_SIZE_T))
  p4 = fftw_alloc_complex(INT(IImageSizeXY(1)*IImageSizeXY(2), C_SIZE_T))
  call c_f_pointer(p1, RImageSimDummy, [IImageSizeXY(1),IImageSizeXY(2)])
  call c_f_pointer(p2, CDummy1, [IImageSizeXY(1),IImageSizeXY(2)])
  call c_f_pointer(p3, CDummy2, [IImageSizeXY(1),IImageSizeXY(2)])
  call c_f_pointer(p4, CCorrelatedImage, [IImageSizeXY(1),IImageSizeXY(2)])
  !...use arr and arr(i,j) as usual...
  
  ! Set the dummy array to the input simulated data

  RImageSimDummy = RImageSim


  !PRINT*,RImageSimDummy(:2,:2)

  ! Plan and Execute the fft of the Simulated Data 
  
  Iplan = FFTW_PLAN_DFT_r2c_2D(IX,IY,RImageSimDummy,CDummy1,FFTW_ESTIMATE)
  CALL FFTW_EXECUTE_DFT_R2C(Iplan,RImageSimDummy,CDummy1)
  CALL FFTW_DESTROY_PLAN(Iplan)

  ! Set the dummy array to the input experimental data

  RImageSimDummy = RImageExpi


  !PRINT*,RImageSimDummy(:2,:2)
 
  ! Plan and Execute the fft of the Experimental Data 

  Iplan = FFTW_PLAN_DFT_R2C_2D(IX,IY,RImageSimDummy,CDummy2,FFTW_ESTIMATE)

  CALL FFTW_EXECUTE_DFT_R2C(Iplan,RImageSimDummy,CDummy2)
  CALL FFTW_DESTROY_PLAN(Iplan)

  !Calculate the Phase Correlation

  CCorrelatedImage = (CDummy1*CONJG(CDummy2))/&
       (&
       ABS(CDummy1*CONJG(CDummy2)+&
       CMPLX(TINY,TINY,C_DOUBLE_COMPLEX)))

  
  ! Plan and Execute the inverse fft of the phase correlation

  Iplan = FFTW_PLAN_DFT_C2R_2D(IX,IY,CCorrelatedImage,RImageSimDummy,FFTW_ESTIMATE)

  CALL FFTW_EXECUTE_DFT_C2R(Iplan,CCorrelatedImage,RImageSimDummy)

  CALL FFTW_DESTROY_PLAN(Iplan)

  
  RCrossCorrelation = MAXVAL(RImageSimDummy)
  IOffset = MAXLOC(RImageSimDummy)
  
  !PRINT*,RImageSimDummy(:2,:2)

  call fftw_free(p1)
  call fftw_free(p2)
  call fftw_free(p3)
  call fftw_free(p4)
  
END SUBROUTINE PhaseCorrelate

!!$SUBROUTINE BiCubicResampling(RImin,IImSize,RImout,RScaleFactor,IErr)
!!$  
!!$USE MyNumbers
!!$  
!!$  USE CConst; USE IConst
!!$  USE IPara; USE RPara
!!$  USE IChannels
!!$  
!!$  USE MPI
!!$  USE MyMPI
!!$
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER(IKIND) :: &
!!$       ind,jnd, ierr
!!$  INTEGER(IKIND),DIMENSION(2),Intent(IN) :: &
!!$       IImSize
!!$  INTEGER(IKIND),DIMENSION(2):: &
!!$       IImNewSize
!!$  REAL(RKIND),DIMENSION(IImSize(2),IImSize(2)),INTENT(IN) :: &
!!$       RImin
!!$  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT) :: &
!!$       RImout
!!$  REAL(RKIND) :: &
!!$       RScaleFactor,RrowFunction,RTotalInterpolation
!!$  REAL(RKIND),DIMENSION(:) :: &
!!$       Ralist,Rblist
!!$
!!$  IF(RScaleFactor.LT.ONE) THEN
!!$     IImNewSize(1) = FLOOR(IImSize(1)*RScaleFactor)
!!$     IImNewSize(2) = FLOOR(IImSize(2)*RScaleFactor)
!!$  ELSE
!!$     IImNewSize(1) = CEILING(IImSize(1)*RScaleFactor)
!!$     IImNewSize(2) = CEILING(IImSize(2)*RScaleFactor)
!!$  END IF
!!$
!!$  ALLOCATE(&
!!$       Ralist(IImNewSize(1)),&
!!$       STAT=IErr)
!!$  IF( IErr.NE.0 ) THEN
!!$     PRINT*,"BiCubicResampling(", my_rank, ") error ", IErr, &
!!$          " in allocation"
!!$     RETURN
!!$  ENDIF
!!$  ALLOCATE(&
!!$       RImout(IImNewSize(1),IImNewSize(2)),&
!!$       STAT=IErr)
!!$  IF( IErr.NE.0 ) THEN
!!$     PRINT*,"BiCubicResampling(", my_rank, ") error ", IErr, &
!!$          " in allocation"
!!$     RETURN
!!$  ENDIF
!!$
!!$  ALLOCATE(&
!!$       Rblist(IImNewSize(2)),&
!!$       STAT=IErr)
!!$  IF( IErr.NE.0 ) THEN
!!$     PRINT*,"BiCubicResampling(", my_rank, ") error ", IErr, &
!!$          " in allocation"
!!$     RETURN
!!$  ENDIF
!!$
!!$  Ralist = (1:IImNewSize(1))*(1/RScaleFactor)
!!$  Rblist = (1:IImNewSize(2))*(1/RScaleFactor)
!!$
!!$  DO ind=1,IImNewSize(1)
!!$     DO jnd=1,IImNewSize(2)
!!$     RrowFunction = -Rblist(ind)*(1-Rblist(ind))**2*RImin() + &
!!$          (1-2*Rblist(ind)**2+Rblist(ind)**3)*RImin() + &
!!$          Rblist(ind)*(1+Rblist(ind)-Rblist(ind)**3)*RImin() - &
!!$          (Rblist(ind)**2)*(1-Rblist(ind))*RImin()
!!$     RImout(ind,jnd) = 
!!$     END DO
!!$  END DO
!!$  
!!$ 
!!$END SUBROUTINE BiCubicResampling
