!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix

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
!! Accounts for the possibility that the incident electron beam is not parallel to the surface normal of the specimen.
!! Proceeds to calculate the resulting wavefunction and the diffracted intensities from the beam

!! Finds the structure matrix that accounts for absorption and the angle of the beam
!! Finds the eigenvalues and eigenvectors of structure matrix
!! Forms the scattering matrix
!! Finds diffracted intensities from scattering matrix.
!! Based off "Structure refinement using precession electron diffraction tomography
!! and dynamical diffraction: theory and implementation" By Lukas Palatinus(2015)

MODULE angle_correction_mod
	
	IMPLICIT NONE
	
	PRIVATE
	PUBLIC :: components, structure_matrix_absorption_orientation, eigen_structure_matrix, &
				M_matrix_formation, weighting_coefficients, scattering_matrix, diffracted_intensities
	
	CONTAINS 
		
		! Calculates components of g vectors and K vector along surface normal vector
		SUBROUTINE components
		
			! Import modules and global variables
		
			USE MyNumbers
			USE MyMPI
			USE BlochPara ONLY : RTiltedK
			USE RPARA, ONLY : RNormDirM, RgDotNorm			
			USE IPara ONLY : nBeams ! nBeams after strong and weak beam determination
		
			IMPLICIT NONE
		
			! Define local variables
			REAL :: RKn
			
			! Now dot g vectors and K vector with surface normal vector
			
			! First do K vector
			RKn = DOT_PRODUCT(RTiltedK, RNormDirM)
			IF(l_alert(IErr,"angle_correction","components, K component")) CALL abort
			
			! gVectors dotted with normal vector are contained in RgDotNorm
			! So no need to calculate as below:
			
			! Then g vectors
			DO ind = 1, nBeams
				RGn(ind) = DOT_PRODUCT(RgMatrix(3, ind, :), RNormDirM)
			END DO
			IF(l_alert(IErr,"angle_correction","components, G component")) CALL abort
			
			RETURN RGn, RKn
		
		END SUBROUTINE components
		
		
		! Calculates the new structure matrix from Ug matrix, accounting for absorption
		SUBROUTINE structure_matrix_absorption_orientation(RKn)
			
			! Import modules and global variables
			
			USE MyNumbers
			USE MyMPI
			USE RPARA ONLY RMeanInnerPotential, RgMatrix, RgDotNorm
			USE CPara, ONLY : CUgMat ! This gives Ug matrix with absorption
			USE BlochPara ONLY : RTiltedK
			USE IPara ONLY : nBeams
			EXTERNAL NORM2
			
			IMPLICIT NONE
			
			! Define local variables
			
			COMPLEX :: CStructureMatrix(nBeams, nBeams), CDiagonal_element(nBeams) &
					   CElement_off(nBeams)
			REAL :: RKn
			INTEGER :: ind, jnd
			
			! Finding off and on diagonal elements of structure matrix
			! First loop iterates over row
			! Second loop over columns
			K = ((NORM2(RTiltedK))**2 + RMeanInnerPotential)**0.5 ! This is mod of K, in sample, RMeanInnerPotential is U0
			DO ind = 1, nBeams
				Kg = (NORM2(K + RgMatrix(3, ind, :)))**0.5 ! This is mod of K + g
				DO jnd = 1, nBeams
					IF ind = jnd THEN
						! On diagonal elements
						CDiagonal_element = (K**2 - Kg**2)/(1 + RgDotNorm(ind)/RKn)**0.5
						CStructureMatrix(ind, ind) = CDiagonal_element
					ELSE
						! Off diagonal elements
						Ugij = CUgMat(ind)(jnd) - CUgMat(jnd)(ind) ! Definitely need to check
						CElement_off = (Ugij)/(((1 + (RgDotNorm(ind))/(RKn))**0.5) * ((1 + (RgDotNorm(jnd))/(RKn))**0.5))
						CStructureMatrix(ind, jnd) = CElement_off
				END DO
			END Do
			
			IF(l_alert(IErr,"angle_correction","structure_matrix")) CALL abort
			
			RETURN CStructureMatrix
			
		END SUBROUTINE structure_matrix_absorption_orientation
	
	
		! Finds the eigenvalues and eigenvectors
		SUBROUTINE eigen_structure_matrix(CStructureMatrix)
			! Finds eigenvalues and eigenvectors of new structure matrix
			! Then forms matrix of eigenvectors(C in Kirkland, B in Palatinus)
			! Also forms matrix of eigenvalues
			
			! Import modules and global variables
			
			USE MyNumbers
			USE MyMPI
			USE bloch_mod
			USE IPara ONLY : nBeams
			
			IMPLICIT NONE
			
			! Define local variables
			
			REAL :: REigenvalues(nBeams), RThickness, &
					REigenvalue_matrix(nBeams, nBeams), RKn
			COMPLEX :: CEigenvector_matrix(nBeams, nBeams) &
					   CStructureMatrix(nBeams, nBeams)
			INTEGER :: ind, jnd
			
			! Solve for eigenvalues and eigenvectors of new structure matrix
			! Matrices are hermitian, so eigenvalues are real
			
			CALL EigenSpectrum(nBeams, CStructureMatrix)
			CEigenvector_matrix = VR ! VR is output from ZGEEV(in EigenSpectrum) containing eigenvectors along columns of matrix
			IF(l_alert(IErr,"angle_correction","eigen_structure_matrix, eigenvectors")) CALL abort
			REigenvalues = W ! W is output from ZGEEV(in EigenSpectrum) containing eigenvalues as an array
			IF(l_alert(IErr,"angle_correction","eigen_structure_matrix, eigenvalues")) CALL abort
			
			! Form matrix of eigenvalues, with exponent, thickness, K component, along diagonal
			DO ind = 1, nBeams
				REigenvalue_matrix(ind, ind) = EXP((((TWOPI * CIMAGONE * RThickness)/ (2 * RKn))) * REigenvalues(ind))
			END DO
			
			RETURN REigenvalue_matrix, CEigenvector_matrix
			
		END SUBROUTINE eigen_structure_matrix
		
		
		! Form the M matrix for calculating the scattering matrix
		SUBROUTINE M_matrix_formation(RKn)
			! This diagonal matrix takes into account the orientation of incident beam compared to surface normal
			! See equation (5) in Palatinus
			! Becomes identity matrix if g vectors are parallel to surface of crystal(dot product goes to 0)
			! Scattering matrix then becomes same as Kirkland
			
			! Import modules and global variables
			
			USE MyNumbers
			USE MyMPI
			USE IPara ONLY : nBeams
			USE RPARA ONLY : RgDotNorm
			
			IMPLICIT NONE
			
			! Define local variables
			
			REAL :: RKn
			REAL, ALLOCATABLE :: RElement(nBeams)
			REAL, ALLOCATABLE :: RMmatrix(nBeams, nBeams)
			REAL :: ind

			
			! Calculates each element on diagonal and inputs to matrix
			DO ind = 1, nBeams
				RElement = 1/((1 + (RgDotNorm(ind)/RKn))**0.5)
				RMmatrix(ind, ind) = RElement
			END DO
			
			IF(l_alert(IErr,"angle_correction","M_matrix_formation, form matrix")) CALL abort
			
			RETURN RMmatrix
			
		END SUBROUTINE M_matrix_formation
		
		
		
		! Takes eigenvalues and eigenvectors from structure matrix to get weighting coefficients for bloch waves
		SUBROUTINE weighting_coefficients(CEigenvector_matrix)
			! Takes eigenvectors matrix
			! Inverts it
			! Weighting coefficients can then be found from
			! Product of matrix and initial conditions on wavefunction
			
			! Import modules and global variables
			
			USE MyNumbers
			USE MyMPI
			USE bloch_mod
			USE IPara ONLY : nBeams
			
			!Define local variables
			
			IMPLICIT NONE
			
			COMPLEX :: CEigenvector_matrix(nBeams, nBeams), CInitial_wavefunction(nBeams) &
						CInverted_vectors(nBeams, nBeams), CWeighting_Coefficients(nBeams)
			
			! Boundary conditions on surface of sample
			! These particular conditions are only valid for singular incident plane wave
			CInitial_wavefunction = CZERO ! All diffracted beams are zero
			CInitial_wavefunction(1) = CONE ! 000 beam has unit amplitude
			
			! Invert matrix of eignvectors
			CInverted_vectors = CALL INVERT(nBeams, CEigenvector_matrix)
			IF(l_alert(IErr,"angle_correction","weighting_coefficients, inverted_vectors")) CALL abort
			
			! Finds weighting coefficients
			CWeighting_Coefficients = MATMUL(CInverted_vectors, CInitial_wavefunction)
			IF(l_alert(IErr,"angle_correction","weighting_coefficients, alpha")) CALL abort
			
			RETURN CWeighting_Coefficients		
			
		END SUBROUTINE weighting_coefficients
		
	
		! Computes the scattering matrix from the eigenvalue and eigenvector matrices
		SUBROUTINE scattering_matrix (CEigenvector_matrix, REigenvalue_matrix, RMmatrix)
			
			! Import modules and global variables
			
			USE MyNumbers
			USE MyMPI
			USE bloch_mod
			USE IPara ONLY : nBeams
			
			IMPLICIT NONE
			
			! Defining local variables
			COMPLEX :: CEigenvector_matrix(nBeams, nBeams)
			REAL :: RMmatrix(nBeams, nBeams), REigenvalue_matrix(nBeams, nBeams)
			COMPLEX :: CScatter_matrix(nBeams, nBeams), CAlt_scatter_matrix(nBeams, nBeams)
		
			! Eigenvector and M matrix need to be inverted for scattering matrix
			
			! First do eigenvector matrix
			CInverted_eigen = CALL INVERT(nBeams, CEigenvector_matrix)
			IF(l_alert(IErr,"angle_correction","scattering_matrix, inverted_eigenvectors")) CALL abort
			
			! Now do the M matrix
			RInverted_M = CALL INVERT(nBeams, RMmatrix)
			IF(l_alert(IErr,"angle_correction","scattering_matrix, inverted_M")) CALL abort
			
			!Calculating scattering matrix
			CScatter_matrix = MATMUL(MATMUL(MATMUL(RMmatrix, CEigenvector_matrix), &
							  REigenvalue_matrix), MATMUL(CInverted_eigen, RInverted_M))
			IF(l_alert(IErr,"angle_correction","scattering_matrix, scattering")) CALL abort
			
			!!!! Alternative Scattering matrix that operates on weighting coefficients instead of boundary conditions
			
			! Forming alternative scatter matrix
			CAlt_scatter_matrix = MATMUL(RMmatrix, MATMUL(CEigenvector_matrix, &
											MATMUL(REigenvalue_matrix, &
											MATMUL(MATMUL(CInverted_eigen, RInverted_M), &
											CEigenvector_matrix)))
			IF(l_alert(IErr,"angle_correction","scattering_matrix, alt_scattering")) CALL abort
			
			RETURN CScatter_matrix, CAlt_scatter_matrix
			
		END SUBROUTINE scattering_matrix
		
		
		! Finds the diffracted intensities from using scattering matrix
		SUBROUTINE diffracted_intensities(CScatter_matrix, CAlt_scatter_matrix, CWeighting_Coefficients)
			! 2 Methods can be used here
			! Calculates wavefunction column vector from
			! Product of scattering matrix and wavefunction on surface(boundary condition)
			! Take column vector of wavefunction
			! Inverse fourier transform of wavefunction to get in terms of all space
			! Get intensity from modulus squared of each element
			
			! Import modules and global variables
			
			USE MyNumbers
			USE MyMPI
			USE IPara ONLY : nBeams
			
			IMPLICIT NONE
			
			! Define local variables
			
			INTEGER :: ind
			REAL :: RIntensity_vector(nBeams), RMmatrix(nBeams, nBeams)
			COMPLEX :: CScatter_matrix(nBeams, nBeams), &
						CAlt_scatter_matrix(nBeams, nBeams), &
						CInitial_wavefunction(nBeams), &
						CFinal_wavefunction_z(nBeams), &
						CFinal_wavefunction(nBeams), &
						CWeighting_Coefficients(nBeams), &
			EXTERNAL IFFT2
			
			
			! Boundary conditions on surface of sample
			! These particular conditions are only valid for single incident plane wave
			CInitial_wavefunction = CZERO ! All diffracted beams are zero
			CInitial_wavefunction(1) = CONE ! 000 beam has amplitude of unity
			
			! First find wavefunction at any thickness in specimen, psi(z).
			CFinal_wavefunction_z = MATMUL(CScatter_matrix, CInitial_wavefunction)
			IF(l_alert(IErr,"angle_correction","diffracted_intensities, final_wavefunction_z")) CALL abort
			
			! Now need to find total wavefunction as a function of all space
			! Do 2D inverse fourier transform on CFinal_wavefunction_z
			CFinal_wavefunction = CALL IFFT2(CFinal_wavefunction_z)	! Need to check this
			IF(l_alert(IErr,"angle_correction","diffracted_intensities, final_wavefunction")) CALL abort
			
			! Modulus squared of each element of final wavefunction for intensity
			! Could maybe do this without a loop(?)
			DO ind = 1, nBeams
				RIntensity_vector(ind) = CFinal_wavefunction(ind) * CONJG(CFinal_wavefunction(ind))
			END DO
			
			IF(l_alert(IErr,"angle_correction","diffracted_intensities, RIntensity_vector")) CALL abort
			
			!!!! Alternative way to find wavefunction as function of depth
			
			! Wavefunction at any thickness
			CFinal_wavefunction_z = MATMUL(CAlt_scatter_matrix, CWeighting_Coefficients)
			
			! Total wavefunction is then found from inverse fourier transform
			CFinal_wavefunction = CALL IFFT2(CFinal_wavefunction_z)	! Need to check this
			IF(l_alert(IErr,"angle_correction","diffracted_intensities, alt_final_wavefunction")) CALL abort
			
			! Modulus squared of each element of final wavefunction for intensity
			! Again, could maybe do this without a loop(?)
			DO ind = 1, nBeams
				RIntensity_vector(ind) = CFinal_wavefunction(ind) * CONJG(CFinal_wavefunction(ind))
			END DO
			
			RETURN RIntensity_vector

		END SUBROUTINE diffracted_intensities
	
END MODULE angle_correction_mod
