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
! $Id: ReadCifFile.f90,v 1.23 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE ReadCifFile(IErr)

  ! -----------------------------------------------------------------------
  ! ReadCifFile: Read the input file
  !
  !	IErr	error code
  ! ----------------------------------------------------------------------
  !
  ! based on 
  ! 
  !                     CIF Tool Box Application 'tbx_ex.f'
  !                     ------------------------
  !                     Version: June 1998
  !
  ! ----------------------------------------------------------------------
  
  !
  ! The tool box common variable file 'ciftbx.cmn' must be present 
  ! at the start of EACH routine using the CIFtbx functions.

  USE MyNumbers
  USE WriteToScreen
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE SPara
  USE CPara
  USE IChannels

  USE MyMPI
  
  IMPLICIT NONE

  INCLUDE       'ciftbx-f90.cmn'

  LOGICAL       f1,f2,f3
  CHARACTER*32  name,Sind
  CHARACTER*80  line
  CHARACTER*4   label(6)
  CHARACTER*1   SAlphabetarray(52)
  CHARACTER*52  alphabet
  CHARACTER*2   rs
  CHARACTER*1   slash
  CHARACTER string*(30)
  REAL          cela,celb,celc,siga,sigb,sigc
  REAL          x,y,z,u,su,sx,sy,sz,B,sB,sOcc,Uso,suso,Occ
  REAL          numb,sdev,dum
  REAL          xf(6),yf(6),zf(6),uij(6,6)
  INTEGER       i,j,nsite, iset, imark
  DATA SAlphabetarray /"A","B","C","D","E","F","G","H","I","J","K","L","M","N", &
       "O","P","Q","R","S","T","U","V","W","X","Y","Z","a","b","c","d","e", &
       "f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v", &
       "w","x","y","z"/
  DATA alphabet /"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"/
  DATA          cela,celb,celc,siga,sigb,sigc/6*0./
  DATA          x,y,z,u,sx,sy,sz,su/8*0./
  DATA          xf,yf,zf,uij/54*0./
  DATA          rs/'\\'/

  INTEGER IAtomCount, ICommaPosLeft, ICommaPosRight, &
       Ipos,Idpos, IXYZminus,IFRACminus, Inum,Idenom,IAtomID
  
  CHARACTER*32 Csym(THREEDIM)

  INTEGER IErr, ind

  ! fudge to deal with gfortran vs. g77

  slash = rs(1:1)

  ! Assign the CIFtbx files 
  f1 = init_( 1, 2, 3, 6 )
  
  !Tells user, entering ReadCIFFile
  CALL Message("ReadCIFFile",IMust,IErr)

  ! Request dictionary validation check
  IF(.NOT.dict_('cif_core.dic','valid dtype')) THEN
     CALL Message("ReadCIFFile",IInfo,IErr, &
          MessageString = "Requested core dictionary not present")

     ! Restore the old clipping action for text fields
     clipt_ = .TRUE.
     pclipt_ = .TRUE.
  END IF
  ! Open the CIF to be accessed

100 name='felix.cif'
  CALL Message("ReadCIFFile",IInfo,IErr, &
       MessageVariable = "Read data from CIF ", MessageString = name)
  IF(.NOT.ocif_(name)) THEN
     !maybe error message here?
     CALL Message("ReadCIFFile",IInfo,IErr, &
             MessageString = "CIF cannot be opened")
     IErr=1
     RETURN
  END IF
  
  ! Assign the data block to be accessed
120 IF(.NOT.data_(' ')) THEN
     IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"ReadCifFile(", my_rank, ") No data_ statement found"
     END IF
     IErr=1
  END IF
  
130 IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     CALL Message("ReadCIFFile",IInfo,IErr, &
          MessageVariable = "Access items in data block ", MessageString = bloc_)
  END IF
  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageString = "Cell length origin") 
  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageVariable = "RLengthX",RVariable = RLengthX)
  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageVariable = "RLengthY",RVariable = RLengthY)
  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageVariable = "RLengthZ",RVariable = RLengthZ)
  
  
  ! Extract some cell dimensions; test all is OK
  ! NEED TO PUT IN A CHECK FOR LENGTH UNITS
  siga = 0.
  sigb = 0.
  sigc = 0.
  
  f1 = numb_('_cell_length_a', cela, siga)
  f2 = numb_('_cell_length_b', celb, sigb)
  f3 = numb_('_cell_length_c', celc, sigc)
  !error call
  IF(.NOT.(f1.AND.f2.AND.f3)) THEN
     IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"ReadCifFile(", my_rank, ") Cell dimension(s) missing!"
     END IF
     IErr=1
  END IF

  RLengthX=cela; RLengthY=celb; RLengthZ=celc
  CALL Message("ReadCIFFile",IInfo,IErr,MessageString = "Cell length") 
  CALL Message("ReadCIFFile",IInfo,IErr,MessageVariable = "RLengthX",RVariable = RLengthX)
  CALL Message("ReadCIFFile",IInfo,IErr,MessageVariable = "RLengthY",RVariable = RLengthY)
  CALL Message("ReadCIFFile",IInfo,IErr,MessageVariable = "RLengthZ",RVariable = RLengthZ)

  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageString = "Standard deviation of length") 
  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageVariable = "siga",RVariable = REAL(siga,RKIND))
  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageVariable = "sigb",RVariable = REAL(sigb,RKIND))
  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageVariable = "sigc",RVariable = REAL(sigc,RKIND))


  siga = 0.
  sigb = 0.
  sigc = 0.
  f1 = numb_('_cell_angle_alpha', cela, siga)
  f2 = numb_('_cell_angle_beta', celb, sigb)
  f3 = numb_('_cell_angle_gamma', celc, sigc)
  
  IF(.NOT.(f1.AND.f2.AND.f3)) THEN
     IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
        PRINT*,"ReadCifFile(", my_rank, ") Cell angle(s) missing!"
     END IF
     IErr=1
  ENDIF

  ! convert angles from degrees to radians
  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageString = "Angle (input)")
  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageVariable = "cela",RVariable = REAL(cela,RKIND))
  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageVariable = "celb",RVariable = REAL(celb,RKIND))
  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageVariable = "celc",RVariable = REAL(celc,RKIND))

  IF (cela.GT.TWOPI) THEN

        CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageVariable = "Angle alpha", &
             RVariable = REAL(cela,RKIND))
        CALL Message("ReadCIFFile",IMoreInfo,IErr, &
             MessageString = "Which is greater than Two Pi, Program will assume this angle is expressed in degrees")

     RAlpha=cela*TWOPI/360.D0;
  END IF
 
  IF (celb.GT.TWOPI) THEN
     CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageVariable = "Angle beta", &
          RVariable = REAL(celb,RKIND))
     CALL Message("ReadCIFFile",IMoreInfo,IErr,& 
          MessageString = "Which is greater than Two Pi, Program will assume this angle is expressed in degrees")
 
     RBeta=celb*TWOPI/360.D0;
  END IF
  IF (celc.GT.TWOPI) THEN
     CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageVariable = "Angle gamma", &
          RVariable = REAL(celc,RKIND))  
     CALL Message("ReadCIFFile",IMoreInfo,IErr, &
          MessageString = "Which is greater than Two Pi, Program will assume this angle is expressed in degrees")
  
     RGamma=celc*TWOPI/360.D0;
  END IF

  CALL Message("ReadCIFFile",IInfo,IErr,MessageString = "Angle (radians)")
  CALL Message("ReadCIFFile",IInfo,IErr,MessageVariable = "RAlpha",RVariable = RAlpha)
  CALL Message("ReadCIFFile",IInfo,IErr,MessageVariable = "RBeta",RVariable = RBeta)
  CALL Message("ReadCIFFile",IInfo,IErr,MessageVariable = "RGamma",RVariable = RGamma)

  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageString = "Standard deviation of angle") 
  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageVariable = "siga",RVariable = REAL(siga,RKIND))
  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageVariable = "sigb",RVariable = REAL(sigb,RKIND))
  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageVariable = "sigc",RVariable = REAL(sigc,RKIND))


 ! IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
 !    PRINT*,"ReadCifFile(", my_rank, ")       ",siga,sigb,sigc
 ! END IF
  
  f1 = numb_('_cell_volume', cela, siga)
  
  !Error message
  IF((f1) .EQV. .FALSE.) THEN
     CALL Message("ReadCIFFile",IInfo,IErr, &
          MessageString = "Volume missing from felix.cif (_cell_volume), calculating from cell parameters (_cell_length etc.)")
!!$     IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
!!$        PRINT*,"ReadCifFile(", my_rank, ") Volume missing!"
!!$     END IF
     IVolumeFLAG= 0
     RVolume= RLengthX*RLengthY*RLengthZ* &
          SQRT(1.0D0 - &
          COS(RAlpha)*COS(RAlpha)-COS(RBeta)*COS(RBeta)-COS(RGamma)*COS(RGamma) + &
          2.0D0*COS(RAlpha)*COS(RBeta)*COS(RGamma))
  ELSE 
     RVolume= cela
     IVolumeFLAG= 1
  ENDIF

  CALL Message("ReadCIFFile",IInfo,IErr,MessageString = "Cell volume")
  CALL Message("ReadCIFFile",IInfo,IErr,MessageVariable = "RVolume", &
          RVariable = RVolume)  

  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageString = "Standard deviation of cell volume")
  CALL Message("ReadCIFFile",IMoreInfo,IErr,MessageVariable = "siga", &
       RVariable = REAL(siga,RKIND))
  
    
  ! Extract atom type symbol data in a loop
  CALL Message("ReadCIFFile",IInfo,IErr,MessageString = "Atom type", &
       RVariable = REAL(siga,RKIND))


  DO      
     f1 = char_('_atom_type_symbol', name)
     CALL Message("ReadCIFFile",IInfo,IErr,MessageString = name)
  
     IF(loop_ .NEQV. .TRUE.) EXIT
  ENDDO

  ! Extract space group notation (expected char string)
  
  f1 = char_('_symmetry_space_group_name_H-M', name)
  CALL Message("ReadCIFFile",IInfo,IErr,MessageVariable = "Space Group",MessageString = name)


  !different types of space groups as well as different phrasing of Hall space groups

  IF (SCAN(name,'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz').EQ.0) THEN
     f1 = char_('_symmetry_space_group_name_Hall',name)
     IF (SCAN(name,'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz').EQ.0) THEN
        f1 = numb_('_symmetry_Int_tables_number',numb,sx)
        IF (numb.LT.TINY) THEN
           f1 = numb_('_space_group_IT_number',numb,sx)
           IF (numb.LT.TINY) THEN
              !Error message
              IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
                 PRINT*,"No Space Group"
              END IF
              IErr = 1
           ELSE
              name = CSpaceGrp(NINT(numb))
           END IF
        ELSE 
           name = CSpaceGrp(NINT(numb))
        END IF
     END IF
  END IF

  CALL Message("ReadCIFFile",IInfo,IErr,MessageVariable = "Space Group", &
       MessageString =name(1:long_))

  SSpaceGroupName=TRIM(name(1:1))
  SSpaceGrp = TRIM(ADJUSTL(name))
  
  !sometimes space group is input in lowercase letters - below changes the first letter to
  !uppercase
  IF (SCAN(alphabet,SSpaceGroupName).GT.26) THEN
     SSpaceGroupName=SAlphabetarray(SCAN(alphabet,SSpaceGroupName)-26)
  END IF
 
  ! ----------------------------------------------------------
  ! Extract atom site data in a loop
  !-----------------------------------------------------------
  CALL Message("ReadCIFFile",IInfo,IErr,MessageString = "Atom sites")
  
  ! counting loop
  IAtomCount=0
  DO 
     f1 = char_('_atom_site_label', name)
     
     IAtomCount= IAtomCount+1

     IF(loop_ .NEQV. .TRUE.) EXIT
  ENDDO
  
  CALL Message("ReadCIFFile",IInfo,IErr,MessageVariable = "IAtomCount", IVariable = IAtomCount) 
  
  ALLOCATE( &
       RAtomSiteFracCoordVec(IAtomCount,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadCifFile(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  ALLOCATE( &
       SAtomName(IAtomCount), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadCifFile(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  ALLOCATE( &
       IAtomNumber(IAtomCount), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadCifFile(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  ALLOCATE( &
       RIsotropicDebyeWallerFactors(IAtomCount), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadCifFile(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  ALLOCATE( &
       RAtomicSitePartialOccupancy(IAtomCount), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadCifFile(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  ALLOCATE( &
       RAnisotropicDebyeWallerFactorTensor(IAtomCount,THREEDIM,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadCifFile(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  ALLOCATE( &
       IAnisotropicDWFTensor(IAtomCount), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadCifFile(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  RAnisotropicDebyeWallerFactorTensor = ZERO

  ! actual data loop

  ITotalAtoms=0
  DO ind=1,IAtomCount
   
     f1 = char_('_atom_site_label', name)
     SAtomName(ind)=name(1:2)
     ! remove the oxcidation state numbers
     Ipos=SCAN(SAtomName(ind),"1234567890")

     IF(Ipos>0) THEN
        WRITE(SAtomName(ind),'(A1,A1)') name(1:1)," "
     ENDIF

     WRITE(Sind,*) ind
     CALL CONVERTAtomName2Number(SAtomName(ind),IAtomNumber(ind), IErr)
     CALL Message("ReadCIFFile",IInfo,IErr,MessageVariable = "IAtomNumber("//TRIM(ADJUSTL(Sind))//")", IVariable = IAtomNumber(ind))

     IF(loop_ .NEQV. .TRUE.) EXIT
     
  ENDDO
  
  !----------------------------------------------------
  ! RESET
  !---------------------------------------------------
  
  CALL CifReset(IErr)
  
  DO ind=1,IAtomCount
     f2 = numb_('_atom_site_fract_x', x, sx)
     RAtomSiteFracCoordVec(ind,1)= x
  ENDDO
  
  !----------------------------------------------------
  ! RESET
  !---------------------------------------------------
  
  CALL CifReset(IErr)

  DO ind=1,IAtomCount
     
     f2 = numb_('_atom_site_fract_y', y, sy)
     RAtomSiteFracCoordVec(ind,2)= y

  ENDDO

  !----------------------------------------------------
  ! RESET
  !---------------------------------------------------

  CALL CifReset(IErr)

  DO ind=1,IAtomCount
  
     f2 = numb_('_atom_site_fract_z', z, sz)
     RAtomSiteFracCoordVec(ind,3)= z
     
  ENDDO 

  
     CALL CifReset(IErr)

  DO ind=1,IAtomCount
     
     B = 0.D0
     Uso = 0.D0
     f2 = numb_('_atom_site_B_iso_or_equiv',B,sB)
     f2 = numb_('_atom_site_U_iso_or_equiv', Uso, suso)
     IF(ABS(B).GT.TINY) THEN
        RIsotropicDebyeWallerFactors(ind) = B
     ELSEIF(ABS(Uso).GT.TINY) THEN
        RIsotropicDebyeWallerFactors(ind) = Uso*(8*PI**2)
     END IF
     
     IF(ABS(B).LT.TINY.AND.ABS(Uso).LT.TINY) THEN
        B = RDebyeWallerConstant
        CALL Message("ReadCIFFile",IInfo,IErr,MessageVariable = "For Atom", IVariable = ind, &
             MessageString =  "Debye Waller Factors missing in the CIF File ")
        IF (IMinReflectionPool.EQ.666) THEN
        CALL Message("ReadCIFFile",IMust,IErr, MessageVariable = "For ye magic vessel ", IVariable = ind, &
              MessageString = "Thar be no Debye Waller Factor in Yar Cif File matey")
        END IF
     END IF
     
     IF(loop_ .NEQV. .TRUE.) EXIT
     
  ENDDO

  !----------------------------------------------------
  ! RESET
  !---------------------------------------------------

  CALL CifReset(IErr)

  !----------------------------------------------------
  ! RESET
  !---------------------------------------------------

  CALL CifReset(IErr)
  
  DO ind=1,IAtomCount
     f2 = numb_('_atom_site_occupancy',Occ, sOcc)
     RAtomicSitePartialOccupancy(ind) = Occ
  ENDDO

  !----------------------------------------------------
  ! RESET
  !----------------------------------------------------

  CALL CifReset
  
  ALLOCATE( &
       SWyckoffSymbols(SIZE(IAtomicSitesToRefine)),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadCifFile(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF

  DO ind=1,SIZE(IAtomicSitesToRefine)
!!$     IAtomID = IAtomsToRefine(ind)
     f2 = char_('_atom_site_Wyckoff_symbol',name)
     SWyckoffSymbols(ind) = name
  ENDDO

  !----------------------------------------------------
  ! RESET
  !---------------------------------------------------

  CALL CifReset(IErr)
  
  DO ind=1,IAtomCount

     f2 = numb_('_atom_site_aniso_U_11',u,su) 
     RAnisotropicDebyeWallerFactorTensor(ind,1,1) = u
     f2 = numb_('_atom_site_aniso_U_22',u,su) 
     RAnisotropicDebyeWallerFactorTensor(ind,2,2) = u
     f2 = numb_('_atom_site_aniso_U_33',u,su) 
     RAnisotropicDebyeWallerFactorTensor(ind,3,3) = u
     f2 = numb_('_atom_site_aniso_U_23',u,su) 
     RAnisotropicDebyeWallerFactorTensor(ind,2,3) = u
     f2 = numb_('_atom_site_aniso_U_12',u,su) 
     RAnisotropicDebyeWallerFactorTensor(ind,1,2) = u
     f2 = numb_('_atom_site_aniso_U_13',u,su) 
     RAnisotropicDebyeWallerFactorTensor(ind,1,3) = u
     
     IAnisotropicDWFTensor(ind) = ind
  ENDDO

  IF(((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10).AND.IAnisoDebyeWallerFactorFlag.EQ.1) THEN
     PRINT*,"RAnisotropicDebyeWallerFactorTensor",RAnisotropicDebyeWallerFactorTensor
  END IF

  ! ----------------------------------------------------------
  ! Extract atom site data in a loop
  !-----------------------------------------------------------

  CALL Message("ReadCIFFile",IInfo,IErr, MessageString = "Symmetries")      
  !IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
  !   PRINT*,"ReadCifFile(", my_rank, ") Symmetries"
  !END IF

  ! counting loop
  ISymCount=0
  DO 
     f1 = char_('_symmetry_equiv_pos_as_xyz', name)

     DO 
        f2 = char_(name, line)
        ISymCount=ISymCount+1

        IF(text_ .NEQV. .TRUE.) EXIT
     ENDDO

     IF(loop_ .NEQV. .TRUE.) EXIT
  ENDDO

  CALL Message("ReadCIFFile",IInfo,IErr,MessageVariable = "found", &
       IVariable = ISymCount, MessageString = "symmetries")
  !IF((IWriteFLAG.GE.1.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
  !   PRINT*,"ReadCifFile(", my_rank, ") found", ISymCount, "symmetries"
  !END IF
  
  ALLOCATE( &
       RSymVec(ISymCount,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadCifFile(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  ALLOCATE( &
       RSymMat(ISymCount,THREEDIM,THREEDIM), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"ReadCifFile(", my_rank, ") error ", IErr, " in ALLOCATE()"
     RETURN
  ENDIF
  
  RSymVec=ZERO
  RSymMat=ZERO
  
  ! actual data loop
  ISymCount=0
  DO 
     f1 = char_('_symmetry_equiv_pos_as_xyz', name)

     DO 
        f2 = char_(name, line)
        ISymCount=ISymCount+1

        ICommaPosLeft = SCAN(name, ",")
        ICommaPosRight= SCAN(name, ",",.TRUE.)

        Csym(1)= name(1:ICommaPosLeft-1)
        Csym(2)= name(ICommaPosLeft+1:ICommaPosRight-1)
        Csym(3)= name(ICommaPosRight+1:LEN_TRIM(name))

        DO ind=1,THREEDIM

           IXYZminus=1
           IFRACminus=1
           
           name= Csym(ind)

           Ipos= SCAN(name, "xX")
           IF(Ipos > 0) THEN ! there is an X
              IF(Ipos>1) THEN
                 IF(name(Ipos-1:Ipos-1)=="-") IXYZminus=-1
              ENDIF
              RSymMat(ISymCount, ind,1)=IXYZminus
           ENDIF
           
           Ipos= SCAN(name, "yY")
           IF(Ipos > 0) THEN ! there is an Y
              IF(Ipos>1) THEN
                 IF(name(Ipos-1:Ipos-1)=="-") IXYZminus=-1
              END IF
              RSymMat(ISymCount, ind,2)=IXYZminus
           ENDIF
           
           Ipos= SCAN(name, "zZ")
           IF(Ipos > 0) THEN ! there is an Z
              IF(Ipos>1) THEN
                 IF(name(Ipos-1:Ipos-1)=="-") IXYZminus=-1
              END IF
              RSymMat(ISymCount, ind,3)=IXYZminus
           ENDIF
           
           Ipos= SCAN(name, "/")
           IF(Ipos > 1) THEN
              IF(Ipos < LEN_TRIM(NAME) ) THEN ! there is an /
                 Inum  = IACHAR(name(Ipos-1:Ipos-1))-48
                 Idenom= IACHAR(name(Ipos+1:Ipos+1))-48
                 
                 IF(Ipos>2) THEN
                    IF(name(Ipos-2:Ipos-2)=="-") IFRACminus=-1
                 ENDIF
                 
              ENDIF
              RSymVec(ISymCount,ind)=IFRACminus*REAL(Inum)/REAL(Idenom)
           ENDIF
        ENDDO
        IF(text_ .NEQV. .TRUE.) EXIT
     ENDDO

     IF(loop_ .NEQV. .TRUE.) EXIT
  ENDDO

  !closes the cif file
  CALL close_
  
  !moved return until after the Call counttotalatoms
  IF (ITotalAtoms.EQ.0) THEN
     CALL CountTotalAtoms(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"ReadCifFile(", my_rank, ") error in CountTotalAtoms()"
       ! GOTO 9999
     ENDIF
  END IF

!!$Reset Message Counter
IMessageCounter =0  
  
  RETURN

END SUBROUTINE ReadCifFile

SUBROUTINE CifReset(IErr)
  
  USE WriteToScreen
  USE IConst

  USE IPara

  IMPLICIT NONE
  INTEGER(IKIND):: IErr
  
  INCLUDE       'ciftbx-f90.cmn'
  
  CHARACTER*30 string
  
!!$  Only print message once
  DO WHILE (IMessageCounter.LT.2)
     CALL Message("CifReset",IMust,IErr)
     CALL Message("CifReset",IMust+IDebug,IErr,MessageString= "CifReset is called several times")
     IMessageCounter = IMessageCounter +1
  END DO

  IF (find_('_cell_angle_alpha','name',string)) THEN
  END IF
END SUBROUTINE CifReset
