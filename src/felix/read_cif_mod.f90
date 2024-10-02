!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-19, all rights reserved
!
! Version: 1.3
! Date: 13-05-2024
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build: g-vector limit 
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
!! Module-description: 
!!
MODULE read_cif_mod

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_cif 

  CONTAINS

  !>
  !! Procedure-description: Reads from felix.cif input file using toolbox
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE read_cif(IErr)

    ! -----------------------------------------------------------------------
    ! ReadCif: Read the input file
    !
    !	IErr	error code
    ! ----------------------------------------------------------------------
    !
    ! based on 
    ! 
    !                     CIF Tool Box Application 'tbx_ex.f'
    !                     ------------------------
    !                     Version June 1998
    !
    ! ----------------------------------------------------------------------
    
    !
    ! The tool box common variable file 'ciftbx.cmn' must be present 
    ! at the start of EACH routine using the CIFtbx functions.

    USE MyNumbers
    USE message_mod

    ! global outputs (or inout)
    USE SPARA, ONLY : SChemicalFormula, SSpaceGroupName, SBasisAtomLabel, &
          SBasisAtomName, SWyckoffSymbol, SSpaceGrp, SSymString
    USE RPARA, ONLY : RLengthX, RLengthY, RLengthZ, RAlpha, RBeta, RGamma, RVolume, &
          RAnisotropicDebyeWallerFactorTensor, RBasisAtomPosition, RBasisIsoDW, &
          RBasisOccupancy, RSymMat, RSymVec, RBasisAtomDelta
    USE IPARA, ONLY : IVolumeFLAG, IBasisAtomicNumber, IMaxPossibleNAtomsUnitCell, &
          ISymCount, IBasisAnisoDW,ISpaceGrp

    ! global inputs
    USE SConst, ONLY : SAllSpaceGrp,SElementSymbolMatrix
    USE RPARA, ONLY : RDebyeWallerConstant
    USE IPARA, ONLY : IAtomsToRefine,ISimFLAG,ILN
    USE SPARA, ONLY : SPrintString
    USE IConst
    
    IMPLICIT NONE

    INCLUDE       'ciftbx-f90.cmn'

    LOGICAL       f1,f2,f3
    CHARACTER(32)  name
    CHARACTER(32)  SChemForm
    CHARACTER(80)  line
    CHARACTER(4)   label(6)
    CHARACTER(1)   SAlphabetarray(52)
    CHARACTER(52)  alphabet
    CHARACTER(62)  alphabetnum
    CHARACTER(2)   rs
    CHARACTER(1)   slash
    CHARACTER(1)   SAtomChar2,SAtomChar3,SAtomChar4
    CHARACTER(40)  Stext
    REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RPrint
    REAL          cela,celb,celc,siga,sigb,sigc
    REAL          x,y,z,u,su,sx,sy,sz,B,sB,sOcc,Uso,suso,Occ,OccSum
    REAL          numb,sdev,dum
    REAL          xf(6),yf(6),zf(6),uij(6,6)
    INTEGER       i,j,nsite, iset, imark
    DATA SAlphabetarray /"A","B","C","D","E","F","G","H","I","J","K","L","M","N", &
         "O","P","Q","R","S","T","U","V","W","X","Y","Z","a","b","c","d","e", &
         "f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v", &
         "w","x","y","z"/
    DATA alphabet /"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"/
    DATA alphabetnum /"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890"/
    DATA          cela,celb,celc,siga,sigb,sigc/6*0./
    DATA          x,y,z,u,sx,sy,sz,su/8*0./
    DATA          xf,yf,zf,uij/54*0./
    DATA          rs/'\\'/

    INTEGER IErr, IAtomCount, ICommaPosLeft, ICommaPosRight, &
         Ipos,Idpos, IoneI,IFRACminus, Inum,Idenom,IAtomID,ind,jnd,knd
    INTEGER(IKIND),DIMENSION(:),ALLOCATABLE :: ILinkedOccupancies
    CHARACTER(32) Csym(ITHREE)

    ! fudge to deal with gfortran vs. g77
    slash = rs(1:1)

    ! Assign the CIFtbx files 
    f1 = init_( 1, 2, 3, 6 )

    ! Request dictionary validation check
    IF(.NOT.dict_('cif_core.dic','valid dtype')) THEN
      ! Restore the old clipping action for text fields
      clipt_ = .TRUE.
      pclipt_ = .TRUE.
    END IF
    
    ! Open the CIF to be accessed
    name='felix.cif'
    IF(.NOT.ocif_(name)) THEN
      IErr=1; IF(l_alert(IErr,"ReadCif","Cannot find .cif")) RETURN
    END IF
    ! Assign the data block to be accessed
    IF(.NOT.data_(' ')) THEN
      IErr=1; IF(l_alert(IErr,"ReadCif","No cif data_ statement found")) RETURN
    END IF
    ! Extracts crystal formula
    f1 = char_('_chemical_formula_structural', name)
    IF(.NOT.f1) THEN
      f1 = char_('_chemical_formula_iupac', name)
      IF(.NOT.f1) THEN
        f1 = char_('_chemical_formula_sum', name)
        IF(.NOT.f1) THEN
          IErr=1; IF(l_alert(IErr,"ReadCif","Chemical formula missing")) RETURN
        END IF
      END IF
    END IF
    ! strip spaces/brackets and set global variable SChemicalFormula
    ILN=0!Global variable with length of chemical formula string
    DO jnd = 1,LEN(TRIM(name))
      IF ( SCAN(alphabetnum,name(jnd:jnd)).NE.0) THEN!it is an allowed character, put it in 
        ILN=ILN+1
        SChemicalFormula(ILN:ILN) = name(jnd:jnd)
      END IF
    END DO
    
    ! Extract some cell dimensions; test all is OK
    ! NEED TO PUT IN A CHECK FOR LENGTH UNITS
    siga = 0.
    sigb = 0.
    sigc = 0.
    f1 = numb_('_cell_length_a', cela, siga)
    f2 = numb_('_cell_length_b', celb, sigb)
    f3 = numb_('_cell_length_c', celc, sigc)
    ! error check
    IF(.NOT.(f1.AND.f2.AND.f3)) THEN
        IErr=1; IF(l_alert(IErr,"ReadCif","Cell dimension(s) missing")) RETURN
    END IF
    RLengthX=cela; RLengthY=celb; RLengthZ=celc !global variables
 
    siga = 0.
    sigb = 0.
    sigc = 0.
    f1 = numb_('_cell_angle_alpha', cela, siga)
    f2 = numb_('_cell_angle_beta', celb, sigb)
    f3 = numb_('_cell_angle_gamma', celc, sigc)
    IF(.NOT.(f1.AND.f2.AND.f3)) THEN
        IErr=1; IF(l_alert(IErr,"ReadCif","Cell angles(s) missing")) RETURN
    END IF
    
    ! convert angles from degrees to radians
    IF (cela.GT.TWOPI) THEN!assume this angle is expressed in degrees
      RAlpha=cela*DEG2RADIAN;
    END IF
    IF (celb.GT.TWOPI) THEN!assume this angle is expressed in degrees
      RBeta=celb*DEG2RADIAN;
    END IF
    IF (celc.GT.TWOPI) THEN!assume this angle is expressed in degrees
      RGamma=celc*DEG2RADIAN;
    END IF
    CALL message( LXL, dbg14, "Unit cell angles alpha, beta, gamma", (/ RAlpha*RADIAN2DEG,RBeta*RADIAN2DEG,RGamma*RADIAN2DEG /) )
    f1 = numb_('_cell_volume', cela, siga)
    !Cell volume
    IF((f1) .EQV. .FALSE.) THEN
      IVolumeFLAG= 0
      RVolume= RLengthX*RLengthY*RLengthZ* &
            SQRT(ONE-COS(RAlpha)*COS(RAlpha)-COS(RBeta)*COS(RBeta)-COS(RGamma)*COS(RGamma) + &
            TWO*COS(RAlpha)*COS(RBeta)*COS(RGamma))
    ELSE 
      RVolume= cela
      IVolumeFLAG= 1
    END IF
    CALL message ( LXL, dbg14, "Unit cell volume", RVolume )

    ! Extract space group
    !try the number first
    f1 = numb_('_symmetry_Int_tables_number',numb,sx)
    IF (numb.LT.TINY) f1 = numb_('_space_group_IT_number',numb,sx)
    !if no number, look for a string
    IF (numb.LT.TINY) THEN
      f1 = char_('_symmetry_space_group_name_H-M', name)
      IF (SCAN(name,alphabet).EQ.0) f1 = char_('_symmetry_space_group_name_Hall',name)
      !If we still have nothing, we have run out of options
      IF (SCAN(name,alphabet).EQ.0) THEN
        IErr=1; IF(l_alert(IErr,"ReadCif","No Space Group")) RETURN
      ELSE!we have a string
        SSpaceGrp = TRIM(ADJUSTL(name))
        !Get ISpaceGrp (global variable)
        CALL ConvertSpaceGroupToNumber(IErr)
        IF (l_alert(IErr,"ReadCif","error converting sapce group")) RETURN
      END IF
    ELSE!we have a number
      ISpaceGrp = numb
      SSpaceGrp = SAllSpaceGrp(ISpaceGrp)
    END IF
    !SpaceGroupName is the first letter of the spacegroup, used to determine reflection rules
    SSpaceGroupName=TRIM(SSpaceGrp(1:1))
    IF (SCAN(alphabet,SSpaceGroupName).GT.26) SSpaceGroupName=SAlphabetarray(SCAN(alphabet,SSpaceGroupName)-26)

    WRITE(SPrintString,FMT='(A10,A,A2,A)') "Material: ",SChemicalFormula(1:ILN),", ",SSpaceGrp
    CALL message( LS, dbg3, SPrintString)
    
    ! ----------------------------------------------------------
    ! Extract atom site data
    ! count how many atoms
    IAtomCount=0
    DO 
      f1 = char_('_atom_site_label', name)
      IAtomCount= IAtomCount+1
      IF(loop_ .NEQV. .TRUE.) EXIT
    END DO

    WRITE(SPrintString,FMT='(A11,I3)') "IAtomCount ", IAtomCount
    CALL message( LM, dbg3, SPrintString)

    IF (ISimFLAG.EQ.0.AND.SIZE(IAtomsToRefine,DIM=1).GT.IAtomCount) THEN
      IErr=1; IF(l_alert(IErr,"ReadCif",&
            "Number of atomic sites to refine is larger than the number of atoms. "//&
            "Please correct in felix.inp")) RETURN
    END IF
    !allocate variables
    !coordinates of the basis
    ALLOCATE(RBasisAtomPosition(IAtomCount,ITHREE),STAT=IErr)
    ALLOCATE(RBasisAtomDelta(IAtomCount,ITHREE),STAT=IErr)
    ALLOCATE(SBasisAtomLabel(IAtomCount),STAT=IErr)
    ALLOCATE(SBasisAtomName(IAtomCount),STAT=IErr)
    ALLOCATE(IBasisAtomicNumber(IAtomCount),STAT=IErr)
    ALLOCATE(SWyckoffSymbol(IAtomCount),STAT=IErr)
    ALLOCATE(RBasisIsoDW(IAtomCount),STAT=IErr)
    ALLOCATE(RBasisOccupancy(IAtomCount),STAT=IErr)
    ALLOCATE(ILinkedOccupancies(IAtomCount),STAT=IErr)
    ALLOCATE(RAnisotropicDebyeWallerFactorTensor(IAtomCount,ITHREE,ITHREE),STAT=IErr)
    ALLOCATE(IBasisAnisoDW(IAtomCount),STAT=IErr)
    IF(l_alert(IErr,"ReadCif","Allocation error")) RETURN
    !initialise variables
    IBasisAtomicNumber = 0
    RAnisotropicDebyeWallerFactorTensor = ZERO
    IMaxPossibleNAtomsUnitCell = 0
    RBasisAtomDelta = ZERO
    ! input data loop
    DO ind=1,IAtomCount
      B = ZERO
      Uso = ZERO
      f1 = char_('_atom_site_label', name)
      SBasisAtomLabel(ind)=name
      f1 = char_('_atom_site_type_symbol', name)
      ! accommodate cifs without atom symbols by using the label
      IF(name.EQ."") THEN
        SBasisAtomName(ind)=SBasisAtomLabel(ind)
      ELSE
        SBasisAtomName(ind)=name(1:2)
      END IF
      ! checks on second letter of name
      SAtomChar2=TRIM(SBasisAtomName(ind)(2:2))
      ! remove numbers from single-letter elements (O,F etc.)
      IF (SCAN(SAtomChar2,"1234567890+-()").GT.0) THEN
         SAtomChar2= SBasisAtomName(ind)(1:1)
         WRITE(SBasisAtomName(ind),'(A1,A1)') SAtomChar2," "
      ENDIF
      SAtomChar2=TRIM(SBasisAtomName(ind)(2:2))
      IF (SAtomChar2.NE." ") THEN
        ! convert second letter to lower case
        IF (SCAN(alphabet,SAtomChar2).LT.26) THEN
          SAtomChar2=SAlphabetarray(SCAN(alphabet,SAtomChar2)+26)
          SAtomChar3=SBasisAtomName(ind)(1:1)
          WRITE(SBasisAtomName(ind),'(A1,A1)') SAtomChar3,SAtomChar2
        END IF
      END IF
      !checks on first letter of name
      SAtomChar2=TRIM(SBasisAtomName(ind)(1:1))
      ! convert first letter to upper case
      IF (SCAN(alphabet,SAtomChar2).GT.26) THEN
        SAtomChar2=SAlphabetarray(SCAN(alphabet,SAtomChar2)-26)
        SAtomChar3=SBasisAtomName(ind)(2:2)
        WRITE(SBasisAtomName(ind),'(A1,A1)') SAtomChar2,SAtomChar3
      END IF
      !get atomic number
      IBasisAtomicNumber(ind)=0
      DO jnd=1,INElements!NB must match SElementSymbolMatrix defined in smodules line 73
        IF(TRIM(SBasisAtomName(ind)).EQ.TRIM(SElementSymbolMatrix(jnd))) THEN
          IBasisAtomicNumber(ind)=jnd
        END IF
      END DO
      !Deuterium - replace with hydrogen
      IF(IBasisAtomicNumber(ind).EQ.104) IBasisAtomicNumber(ind)=1
      IF (IBasisAtomicNumber(ind).EQ.0) THEN
        WRITE(SPrintString,'(A,I0,A,A)') &
             "Could not find Z for atom ", ind ," with symbol ", &
             SBasisAtomName(ind)
        IErr=1
        IF(l_alert(IErr,"ReadCif",SPrintString)) RETURN
      END IF
      !Wyckoff symbol
      f1 = char_('_atom_site_Wyckoff_symbol',name)
      !If there is no Wyckoff symbol use 'x'.  To be picked up later if coord refinement is attempted!
      IF (name.NE."") THEN
        SWyckoffSymbol(ind) = name
      ELSE
        SWyckoffSymbol(ind) = "x"
      END IF
      !coordinates
      f2 = numb_('_atom_site_fract_x', x, sx)
      RBasisAtomPosition(ind,1)= x
      f2 = numb_('_atom_site_fract_y', y, sy)
      RBasisAtomPosition(ind,2)= y
      f2 = numb_('_atom_site_fract_z', z, sz)
      RBasisAtomPosition(ind,3)= z
      !Isotropic D-W factor
      f2 = numb_('_atom_site_B_iso_or_equiv',B,sB)
      f2 = numb_('_atom_site_U_iso_or_equiv',Uso,suso)
      IF(ABS(B).GT.TINY) THEN
        RBasisIsoDW(ind) = B
      ELSE
        IF(ABS(Uso).GT.TINY) THEN
          RBasisIsoDW(ind) = Uso*(8*PI**2)
        ELSE
          RBasisIsoDW(ind) = RDebyeWallerConstant
        END IF
      END IF
      ! if not given, occupancy is initially assumed to be 1.0 
      RBasisOccupancy(ind) = ONE
      f2 = numb_('_atom_site_occupancy',Occ, sOcc)
      IF(Occ.GT.TINY) RBasisOccupancy(ind) = Occ
      CALL message( LM, dbg7, "For Atom ",ind)
      CALL message( LM, dbg7, SBasisAtomLabel(ind)//SBasisAtomName(ind)//&
            " Z=",IBasisAtomicNumber(ind) )
      CALL message( LM, dbg7, "RBasisAtomPosition", RBasisAtomPosition(ind,:) )
      CALL message( LM, dbg7, "(DWF, occupancy) respectively = ",&
            (/ RBasisIsoDW(ind), RBasisOccupancy(ind) /) )
      
      IF(loop_ .NEQV. .TRUE.) EXIT
    END DO

    ! Occupancy check
    ! LinkedOccupancies lists the unique basis atom sites
    ILinkedOccupancies = 0
    knd = 0
    DO ind=1,IAtomCount
    ! if ILinkedOccupancies = 0 it's a new atom
      IF (ILinkedOccupancies(ind).EQ.0) THEN
        knd = knd + 1
        ILinkedOccupancies(ind) = knd
        ! look for atoms on the same site
        DO jnd = ind+1,IAtomCount
          IF (ABS(SUM(RBasisAtomPosition(ind,:)-RBasisAtomPosition(jnd,:))).LT.TINY) THEN
            ILinkedOccupancies(jnd) = knd
          END IF
        END DO
      END IF
    END DO
    WRITE(SPrintString,'(I,A22,I3,A13)') IAtomCount," atoms in the basis on ",knd," unique sites"
    CALL message( LS, TRIM(ADJUSTL(SPrintString)))
    ! now go through this list and check site occupancies are <=1
    DO ind=1,knd
      OccSum = ZERO
      DO jnd=1,IAtomCount
        IF (ILinkedOccupancies(jnd).EQ.ind) THEN  
          OccSum = OccSum + RBasisOccupancy(jnd)
        END IF
      END DO
      IF (OccSum.GT.ONE) THEN  ! oops
        IErr = 1
      END IF
    END DO
    IF(l_alert(IErr,"ReadCif","site occupancy >1, aborting: please correct the .cif")) RETURN

!    ! Anisotropic D-W factor
!    DO ind=1,IAtomCount
!      f2 = numb_('_atom_site_aniso_U_11',u,su) 
!      RAnisotropicDebyeWallerFactorTensor(ind,1,1) = u
!      f2 = numb_('_atom_site_aniso_U_22',u,su) 
!      RAnisotropicDebyeWallerFactorTensor(ind,2,2) = u
!      f2 = numb_('_atom_site_aniso_U_33',u,su) 
!      RAnisotropicDebyeWallerFactorTensor(ind,3,3) = u
!      f2 = numb_('_atom_site_aniso_U_23',u,su) 
!      RAnisotropicDebyeWallerFactorTensor(ind,2,3) = u
!      f2 = numb_('_atom_site_aniso_U_12',u,su) 
!      RAnisotropicDebyeWallerFactorTensor(ind,1,2) = u
!      f2 = numb_('_atom_site_aniso_U_13',u,su) 
!      RAnisotropicDebyeWallerFactorTensor(ind,1,3) = u
!      IBasisAnisoDW(ind) = ind
!    END DO

    ! count how many symmetry elements
    ISymCount=0 ! we assume that ONE of the below will work
    Stext = '_symmetry_equiv_pos_as_xyz'
    DO
      f1 = char_(Stext, name)
      IF(name.EQ."") EXIT
      DO
        f2 = char_(name, line)
        ISymCount=ISymCount+1
        IF(text_ .NEQV. .TRUE.) EXIT
       END DO
      IF(loop_ .NEQV. .TRUE.) EXIT
    END DO
    IF (ISymCount.EQ.0) THEN
      Stext = '_space_group_symop_operation_xyz'
      DO 
        f1 = char_(Stext, name)
        DO 
          f2 = char_(name, line)
          ISymCount=ISymCount+1
          IF(text_ .NEQV. .TRUE.) EXIT
         END DO
        IF(loop_ .NEQV. .TRUE.) EXIT
      END DO
    END IF
    ALLOCATE(SSymString(ISymCount),STAT=IErr)
    ALLOCATE(RSymMat(ISymCount,ITHREE,ITHREE),STAT=IErr)
    ALLOCATE(RSymVec(ISymcount,ITHREE),STAT=IErr)
    IF(l_alert(IErr,"ReadCif","allocate RSymMat")) RETURN
    
    RSymVec=ZERO
    RSymMat=ZERO
    
    ! Fill the symmetry matrix
    ind=0
    DO 
      f1 = char_(Stext, name)
      DO
        ind=ind+1
        f2 = char_(name, line)
        SSymString(ind)=TRIM(ADJUSTL(name))
        ICommaPosLeft = SCAN(name, ",")
        ICommaPosRight= SCAN(name, ",",.TRUE.)
        Csym(1)= name(1:ICommaPosLeft-1)
        Csym(2)= name(ICommaPosLeft+1:ICommaPosRight-1)
        Csym(3)= name(ICommaPosRight+1:LEN_TRIM(name))
        DO jnd=1,ITHREE
          IFRACminus=1
          name= Csym(jnd)
          Ipos= SCAN(name, "xX")
          IF(Ipos.GT.0) THEN ! there is an X
            IoneI=1
            IF(Ipos.GT.1) THEN
              IF(name(Ipos-1:Ipos-1).EQ."-") IoneI=-1
            END IF
            RSymMat(ind,jnd,1)=IoneI
          END IF
             
          Ipos= SCAN(name, "yY")
          IF(Ipos.GT.0) THEN ! there is a Y
            IoneI=1
            IF(Ipos.GT.1)THEN
              IF(name(Ipos-1:Ipos-1).EQ."-") IoneI=-1
            END IF
            RSymMat(ind,jnd,2)=IoneI
          END IF
             
          Ipos= SCAN(name, "zZ")
          IF(Ipos.GT.0) THEN ! there is a Z
            IoneI=1
            IF(Ipos.GT.1) THEN
              IF(name(Ipos-1:Ipos-1).EQ."-") IoneI=-1
            END IF
            RSymMat(ind,jnd,3)=IoneI
          END IF
          Ipos= SCAN(name, "/")
          IF(Ipos.GT.1) THEN
            IF(Ipos < LEN_TRIM(NAME) ) THEN ! there is an /
              Inum  = IACHAR(name(Ipos-1:Ipos-1))-48
              Idenom= IACHAR(name(Ipos+1:Ipos+1))-48
              IF(Ipos.GT.2) THEN
                IF(name(Ipos-2:Ipos-2)=="-") IFRACminus=-1
              END IF
            END IF
            RSymVec(ind,jnd)=IFRACminus*REAL(Inum)/REAL(Idenom)
!WRITE(SPrintString,'(I4, A11, F5.2,1X,F5.2,1X,F5.2)') ind,": RSymVec: ", RSymVec(ind,:)
!IF(my_rank.EQ.0) PRINT*, TRIM(ADJUSTL(SPrintString))
          END IF
!WRITE(SPrintString,'(A8,1X,A8,1X,A8)') Csym(:)
!IF(my_rank.EQ.0) PRINT*, TRIM(ADJUSTL(SPrintString))
!WRITE(SPrintString,'(I4,A9,F3.0,1X,F3.0,1X,F3.0,2X,F3.0,1X,F3.0,1X,F3.0,2X,F3.0,1X,F3.0,1X,F3.0)') ind,"RSymMat: ",RSymMat(ind,1,:),RSymMat(ind,2,:),RSymMat(ind,3,:)
!IF(my_rank.EQ.0) PRINT*, TRIM(ADJUSTL(SPrintString))
!WRITE(SPrintString,'(I4, A11, F3.0,1X,F3.0,1X,F3.0)') ind,": RSymVec: ", RSymVec(ind,:)
!IF(my_rank.EQ.0) PRINT*, TRIM(ADJUSTL(SPrintString))
        END DO
        IF(text_ .NEQV. .TRUE.) EXIT
      END DO

      IF(loop_ .NEQV. .TRUE.) EXIT
    END DO

    ! closes the cif file
    CALL close_

  END SUBROUTINE read_cif

  !>
  !! Procedure-description: Strips any character from string which is not simply
  !! numeric or alphabetic
  !!
  !! Major-Authors: Jacob Richardson (2017)
  !! 
  SUBROUTINE strip_chars(SString,SStripped)
  !now redundant
    !?? only used once but it is a nice pure utility function
    CHARACTER(*), INTENT(IN) :: SString
    CHARACTER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: SStripped
    CHARACTER(LEN_TRIM(SString)) :: SPaddedStripped
    INTEGER :: i, n
    n = 0
    DO i = 1,LEN_TRIM(SString) 
      IF ( VERIFY(SString(i:i), &
            "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890") == 0) THEN
        n = n + 1
        SPaddedStripped(n:n) = SString(i:i)
      END IF
    END DO
    SStripped = SPaddedStripped(1:n)
  END SUBROUTINE

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !>
  !! Procedure-description: Convert SSpacegrp to lower case and Compare
  !! SSpaceGrpNoSpaces with every space group
  !!
  !! Major-Authors: Keith Evans (2014), Richard Beanland (2016)
  !!
  SUBROUTINE ConvertSpaceGroupToNumber(IErr)

    ! called from read_cif
    USE MyNumbers
    USE message_mod

    ! global inputs
    USE SPARA, ONLY : SSpaceGrp
    USE SConst, ONLY : SAllSpaceGrp
    USE IPARA, ONLY : ISpaceGrp

    IMPLICIT NONE

    INTEGER(IKIND),INTENT(OUT) :: IErr
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

END MODULE read_cif_mod
