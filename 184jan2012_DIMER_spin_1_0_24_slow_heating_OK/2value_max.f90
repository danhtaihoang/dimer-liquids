!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!     HOANG Danh Tai - Laboratoire de Physique Théorique et Modélisation 
!     UMR 8089 CNRS-Université de Cergy-Pontoise
!     2, Avenue Adolphe Chauvin, 95032 Cergy-Pontoise Cedex, France
!----------------------------------------------------------------------------------------------------!
!     PROGRAMME: MONTE CARLO TRANSPORT OF SPINS
!     !16/11/2010: Tinh gia tri TB cua E, M, EM, Cv, Ksi, V theo pp histogram
!       
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!            

      PROGRAM average_thermal_histo_max

      IMPLICIT NONE

!!! BEGIN -----------------------------------------------------------------------------------------------
      CHARACTER (256) Ligne21,Ligne31

      INTEGER (KIND=4) :: i_T,n_T,natx,naty,natz,n_Ti_histo,n_To_histo

      REAL    (KIND=8) :: T,L1,L
      REAL    (KIND=8) :: Mz_Tc,Cvmax,Ksimax,V1max,V2max,Mz_Tc1,Cvmax1,Ksimax1,V1max1,V2max1

      REAL    (KIND=8),DIMENSION(:),ALLOCATABLE :: EN_histo,Mz_histo,Cv_histo,Ksi_histo,V1_histo,V2_histo

                 
!!! END Decalation -----------------------------------------------------------------------------------------

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! Doc gia tri parameter tu file parameter.in
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      CHARACTER (LEN=150) :: tamp
      OPEN(11,file='1parameter.in')
      
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I5))')    tamp, natx
      READ(11, '(A30,(I5))')    tamp, naty
      READ(11, '(A30,(I5))')    tamp, natz
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp,n_Ti_histo
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp      
      READ(11, '(A30,(F7.4))')  tamp,T
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp,n_To_histo
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      CLOSE(11) 

      n_T=n_To_histo*n_Ti_histo
      L=natx

      ALLOCATE(EN_histo(n_T))
      ALLOCATE(Mz_histo(n_T))
      ALLOCATE(Cv_histo(n_T))
      ALLOCATE(Ksi_histo(n_T))
      ALLOCATE(V1_histo(n_T))
      ALLOCATE(V2_histo(n_T))

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! Tim gia tri lon nhat
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      Mz_Tc=0.
      Cvmax=0.
      Ksimax=0.
      V1max=0.
      V2max=0.

      OPEN(unit=12,file='average_thermal_histo.dat')

      DO i_T=1,n_T

            READ(12,*)T,EN_histo(i_T),Mz_histo(i_T),Cv_histo(i_T),Ksi_histo(i_T),V1_histo(i_T),V2_histo(i_T)

            IF (Cv_histo(i_T)>Cvmax) THEN
                  Cvmax=Cv_histo(i_T)      
            END IF

            IF (Ksi_histo(i_T)>Ksimax) THEN
                  Ksimax=Ksi_histo(i_T)
                  Mz_Tc=Mz_histo(i_T)      
            END IF
            
            IF (V1_histo(i_T)>V1max) THEN
                  V1max=V1_histo(i_T)      
            END IF

            IF (V2_histo(i_T)>V2max) THEN
                  V2max=V2_histo(i_T)      
            END IF

      END DO 

      CLOSE(12)


!!! Write cac gia tri cuc dai --------------------------------------------      
      OPEN (unit=21,file='average_thermal_histo_max.dat')
            WRITE(Ligne21,*)L,Mz_Tc,Cvmax,Ksimax,V1max,V2max
            WRITE(21,'(a)') trim(Ligne21)
      CLOSE (21)

!!! Write log cua cac gia tri cuc dai
      L1=log(L)
      Mz_Tc1=log(Mz_Tc)
      Cvmax1=log(Cvmax)
      Ksimax1=log(Ksimax)
      V1max1=abs(log(V1max))
      V2max1=abs(log(V2max))

      OPEN (unit=31,file='average_thermal_histo_max_log.dat')
            WRITE(Ligne31,*)L1,Mz_Tc1,Cvmax1,Ksimax1,V1max1,V2max1
            WRITE(31,'(a)') trim(Ligne31)
      CLOSE (31)

      END PROGRAM average_thermal_histo_max

