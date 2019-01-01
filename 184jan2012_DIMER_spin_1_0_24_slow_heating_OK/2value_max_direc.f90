!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!     HOANG Danh Tai - Laboratoire de Physique Théorique et Modélisation 
!     UMR 8089 CNRS-Université de Cergy-Pontoise
!     2, Avenue Adolphe Chauvin, 95032 Cergy-Pontoise Cedex, France
!----------------------------------------------------------------------------------------------------!
!     PROGRAMME: MONTE CARLO TRANSPORT OF SPINS
!     !24.11.2010: Tim gia tri Mz_Tc, Cvmax, Ksimax tu cac gia tri tinh duoc boi simulation direc
!       
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!            

      PROGRAM average_thermal_direc_max

      IMPLICIT NONE

!!! BEGIN -----------------------------------------------------------------------------------------------
      CHARACTER (256) Ligne31

      INTEGER (KIND=4),PARAMETER :: n_T = 21            ! n_T = 21 x 1
      REAL (KIND=8),PARAMETER :: L = 80.              ! L = nx = ny

      INTEGER (KIND=4)           :: i_T

      REAL    (KIND=8) :: T,L1
      REAL    (KIND=8) :: Mz_Tc,Cvmax,Ksimax,Mz_Tc1,Cvmax1,Ksimax1

      REAL (KIND=8),DIMENSION(n_T) :: EN_direc,Mz_direc,Cv_direc,Ksi_direc
                 
!!! END Decalation -----------------------------------------------------------------------------------------


      Mz_Tc=0.
      Cvmax=0.
      Ksimax=0.

      OPEN(unit=11,file='average_thermal_direc.dat')

      DO i_T=1,n_T
            READ(11,*)T,EN_direc(i_T),Mz_direc(i_T),Cv_direc(i_T),Ksi_direc(i_T)

            IF (Cv_direc(i_T)>Cvmax) THEN
                  Cvmax=Cv_direc(i_T)      
            END IF

            IF (Ksi_direc(i_T)>Ksimax) THEN
                  Ksimax=Ksi_direc(i_T)
                  Mz_Tc=Mz_direc(i_T)      
            END IF
            
      END DO 

      CLOSE(11)

!!! Write log cua cac gia tri cuc dai
      L1=log(L)
      Mz_Tc1=log(Mz_Tc)
      Cvmax1=log(Cvmax)
      Ksimax1=log(Ksimax)

      OPEN (unit=31,file='average_thermal_direc_max_log.dat')
            WRITE(Ligne31,*)L1,Mz_Tc1,Cvmax1,Ksimax1
            WRITE(31,'(a)') trim(Ligne31)
      CLOSE (31)

      END PROGRAM average_thermal_direc_max
