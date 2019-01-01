!!! ================================================================================================
!!! ================================================================================================  
!!! Chuong trinh tong hop ket qua tu cach tinh MP
!!! 24.11.2011: Su dung cho DIMER
!!! Chi su dung de tinh P_MP, con lai muon su dung phai sua lai cho chinh xac
!!! ================================================================================================
!!! ================================================================================================       
      PROGRAM test_open
      IMPLICIT NONE
      
      INTEGER (KIND=4), PARAMETER :: nF=20
      
      CHARACTER (LEN=150) :: tamp

      INTEGER (KIND=4) :: natx,naty,natz,i,j,nE
      REAL    (KIND=8) :: T

      REAL    (KIND=8) :: Cv,Ksi,ABC3,ABC4,ABC5,ABC6,ABC7
      REAL    (KIND=8) :: T2,E_MP,M_MP,Cv_MP,Ksi_MP,E_2_MP,M_2_MP
          
      REAL    (KIND=8),DIMENSION(nF) :: E_moy,M_moy,E_2_moy,M_2_moy
      
      REAL    (KIND=8),DIMENSION(:,:),ALLOCATABLE  :: E,P
      REAL    (KIND=8),DIMENSION(:),ALLOCATABLE :: P_MP
     
!!!=================================================================================================
!!!=================================================================================================
!!!=================================================================================================
      
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
      READ(11, '(A30,(I8))')    tamp,nE
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp

      CLOSE(11) 


      ALLOCATE(E(nF,nE))
      ALLOCATE(P(nF,nE))
      ALLOCATE(P_MP(nE))
      


!!!=================================================================================================
!!!=================================================================================================




!!! Doc gia tri P tu cac file       
      DO i=1,nF

            IF (i==1) THEN
                  OPEN(unit=12,file='1/value_P_at_To.dat')
                  OPEN(unit=13,file='1/average_thermal.dat')
            END IF

            IF (i==2) THEN
                  OPEN(unit=12,file='2/value_P_at_To.dat')
                  OPEN(unit=13,file='2/average_thermal.dat')
            END IF
                             
            IF (i==3) THEN
                  OPEN(unit=12,file='3/value_P_at_To.dat')
                  OPEN(unit=13,file='3/average_thermal.dat')
            END IF 
            
             IF (i==4) THEN
                  OPEN(unit=12,file='4/value_P_at_To.dat')
                  OPEN(unit=13,file='4/average_thermal.dat')
            END IF

            IF (i==5) THEN
                  OPEN(unit=12,file='5/value_P_at_To.dat')
                  OPEN(unit=13,file='5/average_thermal.dat')
            END IF
                             
            IF (i==6) THEN
                  OPEN(unit=12,file='6/value_P_at_To.dat')
                  OPEN(unit=13,file='6/average_thermal.dat')
            END IF 
            
             IF (i==7) THEN
                  OPEN(unit=12,file='7/value_P_at_To.dat')
                  OPEN(unit=13,file='7/average_thermal.dat')
            END IF

            IF (i==8) THEN
                  OPEN(unit=12,file='8/value_P_at_To.dat')
                  OPEN(unit=13,file='8/average_thermal.dat')
            END IF
                             
            IF (i==9) THEN
                  OPEN(unit=12,file='9/value_P_at_To.dat')
                  OPEN(unit=13,file='9/average_thermal.dat')
            END IF 
            
            IF (i==10) THEN
                  OPEN(unit=12,file='10/value_P_at_To.dat')
                  OPEN(unit=13,file='10/average_thermal.dat')
            END IF 


            IF (i==11) THEN
                  OPEN(unit=12,file='11/value_P_at_To.dat')
                  OPEN(unit=13,file='11/average_thermal.dat')
            END IF

            IF (i==12) THEN
                  OPEN(unit=12,file='12/value_P_at_To.dat')
                  OPEN(unit=13,file='12/average_thermal.dat')
            END IF
                             
            IF (i==13) THEN
                  OPEN(unit=12,file='13/value_P_at_To.dat')
                  OPEN(unit=13,file='13/average_thermal.dat')
            END IF 
            
             IF (i==14) THEN
                  OPEN(unit=12,file='14/value_P_at_To.dat')
                  OPEN(unit=13,file='14/average_thermal.dat')
            END IF

            IF (i==15) THEN
                  OPEN(unit=12,file='15/value_P_at_To.dat')
                  OPEN(unit=13,file='15/average_thermal.dat')
            END IF
                             
            IF (i==16) THEN
                  OPEN(unit=12,file='16/value_P_at_To.dat')
                  OPEN(unit=13,file='16/average_thermal.dat')
            END IF 
            
             IF (i==17) THEN
                  OPEN(unit=12,file='17/value_P_at_To.dat')
                  OPEN(unit=13,file='17/average_thermal.dat')
            END IF

            IF (i==18) THEN
                  OPEN(unit=12,file='18/value_P_at_To.dat')
                  OPEN(unit=13,file='18/average_thermal.dat')
            END IF
                             
            IF (i==19) THEN
                  OPEN(unit=12,file='19/value_P_at_To.dat')
                  OPEN(unit=13,file='19/average_thermal.dat')
            END IF 
            
            IF (i==20) THEN
                  OPEN(unit=12,file='20/value_P_at_To.dat')
                  OPEN(unit=13,file='20/average_thermal.dat')
            END IF 


            IF (i==21) THEN
                  OPEN(unit=12,file='21/value_P_at_To.dat')
                  OPEN(unit=13,file='21/average_thermal.dat')
            END IF

            IF (i==22) THEN
                  OPEN(unit=12,file='22/value_P_at_To.dat')
                  OPEN(unit=13,file='22/average_thermal.dat')
            END IF
                             
            IF (i==23) THEN
                  OPEN(unit=12,file='23/value_P_at_To.dat')
                  OPEN(unit=13,file='23/average_thermal.dat')
            END IF 
            
             IF (i==24) THEN
                  OPEN(unit=12,file='24/value_P_at_To.dat')
                  OPEN(unit=13,file='24/average_thermal.dat')
            END IF

            IF (i==25) THEN
                  OPEN(unit=12,file='25/value_P_at_To.dat')
                  OPEN(unit=13,file='25/average_thermal.dat')
            END IF
                             
            IF (i==26) THEN
                  OPEN(unit=12,file='26/value_P_at_To.dat')
                  OPEN(unit=13,file='26/average_thermal.dat')
            END IF 
            
             IF (i==27) THEN
                  OPEN(unit=12,file='27/value_P_at_To.dat')
                  OPEN(unit=13,file='27/average_thermal.dat')
            END IF

            IF (i==28) THEN
                  OPEN(unit=12,file='28/value_P_at_To.dat')
                  OPEN(unit=13,file='28/average_thermal.dat')
            END IF
                             
            IF (i==29) THEN
                  OPEN(unit=12,file='29/value_P_at_To.dat')
                  OPEN(unit=13,file='29/average_thermal.dat')
            END IF 
            
            IF (i==30) THEN
                  OPEN(unit=12,file='30/value_P_at_To.dat')
                  OPEN(unit=13,file='30/average_thermal.dat')
            END IF 

            IF (i==31) THEN
                  OPEN(unit=12,file='31/value_P_at_To.dat')
                  OPEN(unit=13,file='31/average_thermal.dat')
            END IF

            IF (i==32) THEN
                  OPEN(unit=12,file='32/value_P_at_To.dat')
                  OPEN(unit=13,file='32/average_thermal.dat')
            END IF
                             
            IF (i==33) THEN
                  OPEN(unit=12,file='33/value_P_at_To.dat')
                  OPEN(unit=13,file='33/average_thermal.dat')
            END IF 
            
             IF (i==34) THEN
                  OPEN(unit=12,file='34/value_P_at_To.dat')
                  OPEN(unit=13,file='34/average_thermal.dat')
            END IF

            IF (i==35) THEN
                  OPEN(unit=12,file='35/value_P_at_To.dat')
                  OPEN(unit=13,file='35/average_thermal.dat')
            END IF
                             
            IF (i==36) THEN
                  OPEN(unit=12,file='36/value_P_at_To.dat')
                  OPEN(unit=13,file='36/average_thermal.dat')
            END IF 
            
             IF (i==37) THEN
                  OPEN(unit=12,file='37/value_P_at_To.dat')
                  OPEN(unit=13,file='37/average_thermal.dat')
            END IF

            IF (i==38) THEN
                  OPEN(unit=12,file='38/value_P_at_To.dat')
                  OPEN(unit=13,file='38/average_thermal.dat')
            END IF
                             
            IF (i==39) THEN
                  OPEN(unit=12,file='39/value_P_at_To.dat')
                  OPEN(unit=13,file='39/average_thermal.dat')
            END IF 
            
            IF (i==40) THEN
                  OPEN(unit=12,file='40/value_P_at_To.dat')
                  OPEN(unit=13,file='40/average_thermal.dat')
            END IF 
            

            IF (i==41) THEN
                  OPEN(unit=12,file='41/value_P_at_To.dat')
                  OPEN(unit=13,file='41/average_thermal.dat')
            END IF

            IF (i==42) THEN
                  OPEN(unit=12,file='42/value_P_at_To.dat')
                  OPEN(unit=13,file='42/average_thermal.dat')
            END IF
                             
            IF (i==43) THEN
                  OPEN(unit=12,file='43/value_P_at_To.dat')
                  OPEN(unit=13,file='43/average_thermal.dat')
            END IF 
            
             IF (i==44) THEN
                  OPEN(unit=12,file='44/value_P_at_To.dat')
                  OPEN(unit=13,file='44/average_thermal.dat')
            END IF

            IF (i==45) THEN
                  OPEN(unit=12,file='45/value_P_at_To.dat')
                  OPEN(unit=13,file='45/average_thermal.dat')
            END IF
                             
            IF (i==46) THEN
                  OPEN(unit=12,file='46/value_P_at_To.dat')
                  OPEN(unit=13,file='46/average_thermal.dat')
            END IF 
            
             IF (i==47) THEN
                  OPEN(unit=12,file='47/value_P_at_To.dat')
                  OPEN(unit=13,file='47/average_thermal.dat')
            END IF

            IF (i==48) THEN
                  OPEN(unit=12,file='48/value_P_at_To.dat')
                  OPEN(unit=13,file='48/average_thermal.dat')
            END IF
                             
            IF (i==49) THEN
                  OPEN(unit=12,file='49/value_P_at_To.dat')
                  OPEN(unit=13,file='49/average_thermal.dat')
            END IF 
            
            IF (i==50) THEN
                  OPEN(unit=12,file='50/value_P_at_To.dat')
                  OPEN(unit=13,file='50/average_thermal.dat')
            END IF 

            IF (i==51) THEN
                  OPEN(unit=12,file='51/value_P_at_To.dat')
                  OPEN(unit=13,file='51/average_thermal.dat')
            END IF

            IF (i==52) THEN
                  OPEN(unit=12,file='52/value_P_at_To.dat')
                  OPEN(unit=13,file='52/average_thermal.dat')
            END IF
                             
            IF (i==53) THEN
                  OPEN(unit=12,file='53/value_P_at_To.dat')
                  OPEN(unit=13,file='53/average_thermal.dat')
            END IF 
            
             IF (i==54) THEN
                  OPEN(unit=12,file='54/value_P_at_To.dat')
                  OPEN(unit=13,file='54/average_thermal.dat')
            END IF

            IF (i==55) THEN
                  OPEN(unit=12,file='55/value_P_at_To.dat')
                  OPEN(unit=13,file='55/average_thermal.dat')
            END IF
                             
            IF (i==56) THEN
                  OPEN(unit=12,file='56/value_P_at_To.dat')
                  OPEN(unit=13,file='56/average_thermal.dat')
            END IF 
            
             IF (i==57) THEN
                  OPEN(unit=12,file='57/value_P_at_To.dat')
                  OPEN(unit=13,file='57/average_thermal.dat')
            END IF

            IF (i==58) THEN
                  OPEN(unit=12,file='58/value_P_at_To.dat')
                  OPEN(unit=13,file='58/average_thermal.dat')
            END IF
                             
            IF (i==59) THEN
                  OPEN(unit=12,file='59/value_P_at_To.dat')
                  OPEN(unit=13,file='59/average_thermal.dat')
            END IF 
            
            IF (i==60) THEN
                  OPEN(unit=12,file='60/value_P_at_To.dat')
                  OPEN(unit=13,file='60/average_thermal.dat')
            END IF 

            IF (i==61) THEN
                  OPEN(unit=12,file='61/value_P_at_To.dat')
                  OPEN(unit=13,file='61/average_thermal.dat')
            END IF

            IF (i==62) THEN
                  OPEN(unit=12,file='62/value_P_at_To.dat')
                  OPEN(unit=13,file='62/average_thermal.dat')
            END IF
                             
            IF (i==63) THEN
                  OPEN(unit=12,file='63/value_P_at_To.dat')
                  OPEN(unit=13,file='63/average_thermal.dat')
            END IF 
            
             IF (i==64) THEN
                  OPEN(unit=12,file='64/value_P_at_To.dat')
                  OPEN(unit=13,file='64/average_thermal.dat')
            END IF

            IF (i==65) THEN
                  OPEN(unit=12,file='65/value_P_at_To.dat')
                  OPEN(unit=13,file='65/average_thermal.dat')
            END IF
                             
            IF (i==66) THEN
                  OPEN(unit=12,file='66/value_P_at_To.dat')
                  OPEN(unit=13,file='66/average_thermal.dat')
            END IF 
            
             IF (i==67) THEN
                  OPEN(unit=12,file='67/value_P_at_To.dat')
                  OPEN(unit=13,file='67/average_thermal.dat')
            END IF

            IF (i==68) THEN
                  OPEN(unit=12,file='68/value_P_at_To.dat')
                  OPEN(unit=13,file='68/average_thermal.dat')
            END IF
                             
            IF (i==69) THEN
                  OPEN(unit=12,file='69/value_P_at_To.dat')
                  OPEN(unit=13,file='69/average_thermal.dat')
            END IF 
            
            IF (i==70) THEN
                  OPEN(unit=12,file='70/value_P_at_To.dat')
                  OPEN(unit=13,file='70/average_thermal.dat')
            END IF 

            IF (i==71) THEN
                  OPEN(unit=12,file='71/value_P_at_To.dat')
                  OPEN(unit=13,file='71/average_thermal.dat')
            END IF

            IF (i==72) THEN
                  OPEN(unit=12,file='72/value_P_at_To.dat')
                  OPEN(unit=13,file='72/average_thermal.dat')
            END IF
                             
            IF (i==73) THEN
                  OPEN(unit=12,file='73/value_P_at_To.dat')
                  OPEN(unit=13,file='73/average_thermal.dat')
            END IF 
            
             IF (i==74) THEN
                  OPEN(unit=12,file='74/value_P_at_To.dat')
                  OPEN(unit=13,file='74/average_thermal.dat')
            END IF

            IF (i==75) THEN
                  OPEN(unit=12,file='75/value_P_at_To.dat')
                  OPEN(unit=13,file='75/average_thermal.dat')
            END IF
                             
            IF (i==76) THEN
                  OPEN(unit=12,file='76/value_P_at_To.dat')
                  OPEN(unit=13,file='76/average_thermal.dat')
            END IF 
            
             IF (i==77) THEN
                  OPEN(unit=12,file='77/value_P_at_To.dat')
                  OPEN(unit=13,file='77/average_thermal.dat')
            END IF

            IF (i==78) THEN
                  OPEN(unit=12,file='78/value_P_at_To.dat')
                  OPEN(unit=13,file='78/average_thermal.dat')
            END IF
                             
            IF (i==79) THEN
                  OPEN(unit=12,file='79/value_P_at_To.dat')
                  OPEN(unit=13,file='79/average_thermal.dat')
            END IF 
            
            IF (i==80) THEN
                  OPEN(unit=12,file='80/value_P_at_To.dat')
                  OPEN(unit=13,file='80/average_thermal.dat')
            END IF 
            DO j=1,nE
            READ(12,*)T,E(i,j),ABC3,ABC4,ABC5,ABC6,ABC7,P(i,j)
            !WRITE(*,*)T,E(i,j),P(i,j1)
            END DO

            
            READ(13,*)T2,E_moy(i),M_moy(i),Cv,Ksi,E_2_moy(i),M_2_moy(i)

            
      END DO 

      CLOSE(12) 
                 

!!! Tinh gia tri trung binh P_MP

      OPEN(unit=21,file='value_P_at_To.dat')  
       
      P_MP(:)=0.
          
      DO j=1,nE
            DO i=1,nF
            P_MP(j)=P_MP(j)+P(i,j)
            
            END DO
            
            P_MP(j)=P_MP(j)/real(nF)
            
            WRITE(21,*)T,E(1,j),ABC3,ABC4,ABC5,ABC6,ABC7,P_MP(j)
      END DO      

      CLOSE(21) 

!!! Tinh gia tri trung binh E,M,Cv,Ksi

      OPEN(unit=22,file='average_thermal.dat')      
       
      E_MP=0.
      M_MP=0.
      E_2_MP=0.
      M_2_MP=0.
 
      DO i=1,nF
            E_MP=E_MP+E_moy(i)
            M_MP=M_MP+M_moy(i)
            E_2_MP=E_2_MP+E_2_moy(i)
            M_2_MP=M_2_MP+M_2_moy(i)         
      
      END DO
      E_MP=E_MP/real(nF)
      M_MP=M_MP/real(nF)
      E_2_MP=E_2_MP/real(nF)
      M_2_MP=M_2_MP/real(nF)
      Cv_MP=real(natx*naty*natz)*(E_2_MP-E_MP**2.)/(T2**2.)
      Ksi_MP=real(natx*naty*natz)*(M_2_MP-M_MP**2.)/T2

      WRITE(22,*)T2,E_MP,M_MP,Cv_MP,Ksi_MP

      CLOSE(22) 


!!! ================================================================================================



      
      END PROGRAM


