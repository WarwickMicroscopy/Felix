PROGRAM main
  USE read_dm3_mod
  IMPLICIT NONE
  REAL(KIND(1.0D0)),ALLOCATABLE :: RImageArray(:,:,:) 
  INTEGER(4) :: IErr

  CALL ReadInDM3directory('../../samples/SrTiO3_short_dm3/dm3_images_672x668/',672,668,IErr,RImageArray)
  IF(IErr.NE.0) THEN
    WRITE(*,'(A)') 'Error found in ReadInDM3directory'
    RETURN
  ENDIF

  CALL WriteOutImageArrayDirectory('../../samples/SrTiO3_short_dm3/dm3_converted_672x672_binary/',RImageArray)

END PROGRAM
