!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!     HOANG Danh Tai - Laboratoire de Physique Théorique et ModéliStion 
!     UMR 8089 CNRS-Université de Cergy-Pontoise
!     2, Avenue Adolphe Chauvin, 95032 Cergy-Pontoise Cedex, France
!----------------------------------------------------------------------------------------------------!
!     PROGRAMME: MONTE CARLO TRANSPORT OF SPIN ON DIMER 
!!!
!!!    02.01.2012: Sua lai cach tinh khoang cach giua cac Dimer (lay trung diem)
!!!   09.01.2012: Sua lai chi modern Spin (S= +1,-1)
!!!   18.01.2012: Sua lai cach tinh M (tinh nhu Potts, khong phan biet S=+1 hay S=-1)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!  

      PROGRAM main_transport_dimer
      IMPLICIT NONE

      CHARACTER (LEN=150):: CONFIG_INI
      CHARACTER (LEN=150):: TEST_NUMBER_LINE,SEARCH_GS
      CHARACTER (LEN=150):: TEST_UPDATE
      CHARACTER (LEN=50) :: name
      CHARACTER (LEN=15) :: tmp
      CHARACTER (LEN=3)  :: spinAt
      CHARACTER (256)    :: Ligne20,Ligne126
      CHARACTER (256)    :: Ligne14

      REAL    (KIND=8),PARAMETER :: aaa=5.
      REAL    (KIND=8),PARAMETER :: nul=0.


      INTEGER (KIND=4) :: i,j,k,ip,jp,kp,ipp,jpp,kpp,im,jm,km,natx,naty,natz,number,natx_y
      INTEGER (KIND=4) :: natx_p,naty_p,natz_p,i_gs
      INTEGER (KIND=4) :: iT,nT,i0_etat,i_etat,n_etat,number_line_ini,compt_update   
      INTEGER (KIND=4) :: i_loop,i_loop2,n_equi_reseau1,n_equi_reseau2,n_average_thermal,i_times,n_times

      REAL    (KIND=8) :: n_line_x1,n_line_x2,n_line_y1,n_line_y2,n_line_z1,n_line_z2,n_line_total
      REAL    (KIND=8) :: T,delT,Tmax,Tmin,rdn_config,P_config,tab_tmp,rdn_mtp,rdn_etat
      REAL    (KIND=8) :: x,y,z,EN_tmpx,EN_tmpy,EN_tmpz,EN_tmp1x,EN_tmp1y,EN_tmp1z,J1_tmp,J2_tmp
      REAL    (KIND=8) :: energy,order_parameter,energy_2,order_parameter_2
      REAL    (KIND=8) :: rdn_configx,rdn_configy,rdn_configz
      REAL    (KIND=8) :: EN_tmp1,EN_tmp2,EN_tmp
      REAL    (KIND=8) :: EN_moy,EN_2_moy,OP_moy,OP_2_moy,Cv,Ksi

      REAL    (KIND=8) :: order_parameter2,order_parameter2_2,OP2_moy,OP2_2_moy,Ksi2
      REAL    (KIND=8) :: x1,y1,z1,spin1,spin2,rdn_spin,xspin,yspin,zspin
      
         
      INTEGER (KIND=8),DIMENSION(3)                   :: clock
      INTEGER (KIND=4),DIMENSION(:,:,:)  ,ALLOCATABLE :: tab_site

      REAL    (KIND=8),DIMENSION(:,:,:,:),ALLOCATABLE :: S
      REAL    (KIND=8),DIMENSION(3)                   :: n_line
      REAL    (KIND=8),DIMENSION(10)                  :: Hx_tmp2,Hy_tmp2,Hz_tmp2,name_etat
      REAL    (KIND=8),DIMENSION(10)                  :: ENx_D_tmp2,ENy_D_tmp2,ENz_D_tmp2

      INTEGER (KIND=4) :: i1,j1,k1,deli,delj,delk,r01,i0,j0,k0,compt


      REAL    (KIND=8) :: r0,r,D_tmp,A_tmp,EN_D_tmp1,EN_D_tmp,EN_gs,Emin
      REAL    (KIND=8) :: M11,M12,M13,M14


      REAL    (KIND=8) :: n1_line_x1,n1_line_x2,n1_line_y1,n1_line_y2,n1_line_z1,n1_line_z2
      REAL    (KIND=8) :: n2_line_x1,n2_line_x2,n2_line_y1,n2_line_y2,n2_line_z1,n2_line_z2
      REAL    (KIND=8) :: n3_line_x1,n3_line_x2,n3_line_y1,n3_line_y2,n3_line_z1,n3_line_z2
      REAL    (KIND=8) :: n4_line_x1,n4_line_x2,n4_line_y1,n4_line_y2,n4_line_z1,n4_line_z2

      REAL    (KIND=8) :: M1,M11x,M12x,M13x,M14x,M11y,M12y,M13y,M14y,M11z,M12z,M13z,M14z


      !!! Khai bao phan Histogram
      INTEGER (KIND=4) :: i_histo,n_EN_histo
      REAL    (KIND=8) :: EN_min_histo,EN_max_histo,del_EN_histo

      REAL    (KIND=8),DIMENSION(:),ALLOCATABLE :: EN_histo,H_histo,P_histo,Mz_histo
      REAL    (KIND=8),DIMENSION(:),ALLOCATABLE :: EN_2_histo,Mz_2_histo,Mz_EN_histo,Mz_2_EN_histo
      REAL    (KIND=8),DIMENSION(:),ALLOCATABLE :: M0x,M0y,M0z

!!!==================================================================================================
!!!==================================================================================================
      CALL system('rm -r config_ini_3D')
      CALL system('mkdir config_ini_3D')
      CALL system('rm -r config_3D')
      CALL system('mkdir config_3D')
      CALL system('rm *.dat*')


      CALL ini_rdm_number()

      CALL read_input_parameter_file()

      natx_p=natx+1
      naty_p=naty+1
      natz_p=natz+1

      r01=nint(r0)

      natx_y=natx+naty

      ALLOCATE(S(natx_p,naty_p,natz_p,3))
      ALLOCATE(tab_site(natx_p,naty_p,natz_p))
      ALLOCATE(EN_histo(n_EN_histo))
      ALLOCATE(H_histo(n_EN_histo))
      ALLOCATE(P_histo(n_EN_histo))
      ALLOCATE(Mz_histo(n_EN_histo))
      ALLOCATE(EN_2_histo(n_EN_histo))
      ALLOCATE(Mz_2_histo(n_EN_histo))
      ALLOCATE(Mz_EN_histo(n_EN_histo))
      ALLOCATE(Mz_2_EN_histo(n_EN_histo))
      ALLOCATE(M0x(natx_y),M0y(natx_y),M0z(natx_y))

      IF (nT==1) THEN
            delT=0.
      ELSE
            delT=(Tmax-Tmin)/real(nT-1)
      END IF


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! ====== MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM =====
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      OPEN(unit=14,file='1number_line.dat')
      OPEN(unit=20,file='average_thermal.dat')

      CALL load_config_ini()
      CALL write_config_ini_3D()

      IF (SEARCH_GS=='YES') THEN
           
            WRITE(*,*)'Thuc hien chuong trinh tim GS'
            
            DO iT=1,nT
                  WRITE(*,*)'iT = ', iT
                  !CALL ground_state()
                  CALL write_config_3D()
                  CALL energy_gs()

            END DO

            STOP
      END IF


      CALL calcul_time_run()

      DO iT=1,nT
            WRITE(*,*)'iT = ', iT                        

            T=Tmin+delT*real(iT-1)

           ! CALL value_thermal()


            DO i_times=1,n_times
                  CALL equi_reseau1()      
                  CALL average_thermal()
            END DO

            CALL write_config_3D()


            IF (TEST_NUMBER_LINE== 'YES') THEN 
                  CALL cacul_number_line()
            END IF

            IF (TEST_UPDATE=='YES') THEN
                  CALL test_update_line()
            END IF      
            

      END DO

      CLOSE(14)
      CLOSE(20)
     

      DEALLOCATE(S)
      DEALLOCATE(tab_site)

     
      CONTAINS

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE init_rdm_number()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE ini_rdm_number()
      IMPLICIT NONE

      INTEGER (KIND=4) :: i_time,i

      CALL ITIME(clock)
      i_time=(clock(1)+1)*(clock(2)+1)*(clock(3)+1)
        
      DO i=1,i_time
            CALL random_number(tab_tmp)
      ENDDO 

      END SUBROUTINE ini_rdm_number

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE read_input_parameter_file() 
!!! OPEN the parameter from file "parameter.in"
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE read_input_parameter_file()
      IMPLICIT NONE

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
      READ(11, '(A30,(A10))')   tamp, CONFIG_INI
      READ(11, '(A30,(A10))')   tamp, TEST_NUMBER_LINE
      READ(11, '(A30,(A10))')   tamp, TEST_UPDATE
      READ(11, '(A30,(A10))')   tamp, SEARCH_GS
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I5))')    tamp, natx
      READ(11, '(A30,(I5))')    tamp, naty
      READ(11, '(A30,(I5))')    tamp, natz
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F7.4))')  tamp, P_config
      READ(11, '(A30,(F7.4))')  tamp, J1_tmp
      READ(11, '(A30,(F7.4))')  tamp, J2_tmp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F7.4))')  tamp, D_tmp
      READ(11, '(A30,(F7.4))')  tamp, A_tmp
      READ(11, '(A30,(F7.4))')  tamp, r0
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp, nT
      READ(11, '(A30,(F7.4))')  tamp, Tmin
      READ(11, '(A30,(F7.4))')  tamp, Tmax
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp,n_equi_reseau1
      READ(11, '(A30,(I8))')    tamp,n_equi_reseau2
      READ(11, '(A30,(I8))')    tamp,n_average_thermal
      READ(11, '(A30,(I8))')    tamp,n_times
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp,n_EN_histo
      READ(11, '(A30,(F7.4))')  tamp,EN_min_histo
      READ(11, '(A30,(F7.4))')  tamp,EN_max_histo
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp

      CLOSE(11) 

      END SUBROUTINE read_input_parameter_file

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE calcul_time_run
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE calcul_time_run()

      IMPLICIT NONE      
      
      CHARACTER (8)  :: date
      CHARACTER (10) :: time
      CHARACTER (5)  :: zone
      INTEGER,DIMENSION(8) :: value_time1,value_time2

      INTEGER (KIND=4) :: n_date_run,n_hour_run,n_minute_run
      REAL    (KIND=8) :: n_time_run_total

      CALL system('rm time_run.dat')

      OPEN (90,file='time_run.dat')

      CALL date_and_time(date,time,zone,value_time1)
      WRITE(90,'(5I6)')value_time1(5),value_time1(6),value_time1(7),value_time1(8)
      
      CALL equi_reseau()
      
      CALL date_and_time(date,time,zone,value_time2)
      WRITE(90,'(5I6)')value_time2(5),value_time2(6),value_time2(7),value_time2(8)

      n_time_run_total = real(nT*(n_equi_reseau1+(n_equi_reseau2+1)*n_average_thermal))&
           *(real(value_time2(5)-value_time1(5))*60.+real(value_time2(6)-value_time1(6))&
            +real(value_time2(7)-value_time1(7))/60.+real(value_time2(8)-value_time1(8))/60000.)

      n_date_run=int(n_time_run_total/1440.)
      n_hour_run=int(n_time_run_total/60.-n_date_run*24.)
      n_minute_run=int(n_time_run_total-n_date_run*1440.-n_hour_run*60.)

      WRITE(90,*)'n_time_run_total:',n_time_run_total,'mm'
      WRITE(90,*)'n_date_run      :',n_date_run
      WRITE(90,*)'n_hour_run      :',n_hour_run
      WRITE(90,*)'n_minute_run    :',n_minute_run

      CLOSE(90)

      END SUBROUTINE calcul_time_run

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! Initial position configuration
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE load_config_ini()
      IMPLICIT NONE      

      S(:,:,:,:)=0.
      tab_site(:,:,:)=0

!!!--------------------------------------------------------------------
!!! Hinh dang bat ky
      IF (CONFIG_INI == 'NO') THEN
      number_line_ini=0
      
      !!! khi so dimer <= 25%, tung bat ki
      DO WHILE (number_line_ini< int(real(natx*naty*natz)*0.2501))      

            i=0
            j=0
            k=0

            DO WHILE(i==0)
            CALL random_number(rdn_configx)
            i=int(rdn_configx*real(natx_p))
            ENDDO

            DO WHILE(j==0)
                  CALL random_number(rdn_configy)
                  j=int(rdn_configy*real(naty_p))
            ENDDO

            DO WHILE(k==0)
                  CALL random_number(rdn_configz)
                  k=int(rdn_configz*real(natz_p))
            ENDDO

            ip=i+1-(i/natx)*natx
            jp=j+1-(j/naty)*naty
            kp=k+1-(k/natz)*natz

            IF (tab_site(i,j,k)==0) THEN      

                  CALL random_number(rdn_config)

                        !!!------------------ 
                        IF ((rdn_config < 1./3.).and.(tab_site(ip,j,k)==0)) THEN
                              S(i,j,k,1)=1.
                              tab_site(i,j,k)=1
                              tab_site(ip,j,k)=1
                              number_line_ini=number_line_ini+1
                        END IF
                        !!!------------------
                        IF (((1./3.<=rdn_config).and.(rdn_config < 2./3.)).and.(tab_site(i,jp,k)==0)) THEN
                              S(i,j,k,2)=1.
                              tab_site(i,j,k)=1
                              tab_site(i,jp,k)=1
                              number_line_ini=number_line_ini+1
                        END IF
                        !!!------------------
                        IF ((2./3.<=rdn_config).and.(tab_site(i,j,kp)==0)) THEN
                              S(i,j,k,3)=1.
                              tab_site(i,j,k)=1
                              tab_site(i,j,kp)=1
                              number_line_ini=number_line_ini+1
                         END IF
                        !!!------------------

            END IF

      END DO

      WRITE(*,*) number_line_ini

      !!! Lap vao cac khoang trong con thieu cho du so dimer can tung
      
      DO k=1,natz
            kp=k+1-(k/natz)*natz
            DO j=1,naty
                  jp=j+1-(j/naty)*naty
                        DO i=1,natx   
                        ip=i+1-(i/natx)*natx

                        IF (number_line_ini< int(real(natx*naty*natz)*P_config)) THEN

                        CALL random_number(rdn_config)

                        !!!------------------------------------------------------
                        IF (rdn_config<1./3.) THEN

                        IF ((tab_site(i,j,k)==0).and.(tab_site(ip,j,k)==0)) THEN
                              S(i,j,k,1)=1.
                              tab_site(i,j,k)=1
                              tab_site(ip,j,k)=1
                              number_line_ini=number_line_ini+1
                        END IF
                  
                        IF ((tab_site(i,j,k)==0).and.(tab_site(i,jp,k)==0)) THEN
                              S(i,j,k,2)=1.
                              tab_site(i,j,k)=1
                              tab_site(i,jp,k)=1
                              number_line_ini=number_line_ini+1
                        END IF
                       
                        IF ((tab_site(i,j,k)==0).and.(tab_site(i,j,kp)==0)) THEN
                              S(i,j,k,3)=1.
                              tab_site(i,j,k)=1
                              tab_site(i,j,kp)=1
                              number_line_ini=number_line_ini+1
                        END IF

                        END IF
                        !!!------------------------------------------------------

                        IF ((1./3.<=rdn_config).and.(rdn_config<2./3.)) THEN
                                          
                        IF ((tab_site(i,j,k)==0).and.(tab_site(i,jp,k)==0)) THEN
                              S(i,j,k,2)=1.
                              tab_site(i,j,k)=1
                              tab_site(i,jp,k)=1
                              number_line_ini=number_line_ini+1
                        END IF
                       
                        IF ((tab_site(i,j,k)==0).and.(tab_site(i,j,kp)==0)) THEN
                              S(i,j,k,3)=1.
                              tab_site(i,j,k)=1
                              tab_site(i,j,kp)=1
                              number_line_ini=number_line_ini+1
                        END IF

                        IF ((tab_site(i,j,k)==0).and.(tab_site(ip,j,k)==0)) THEN
                              S(i,j,k,1)=1.
                              tab_site(i,j,k)=1
                              tab_site(ip,j,k)=1
                              number_line_ini=number_line_ini+1
                        END IF

                        END IF
                        !!!------------------------------------------------------

                        IF (rdn_config>=2./3.) THEN

                        IF ((tab_site(i,j,k)==0).and.(tab_site(i,j,kp)==0)) THEN
                              S(i,j,k,3)=1.
                              tab_site(i,j,k)=1
                              tab_site(i,j,kp)=1
                              number_line_ini=number_line_ini+1
                        END IF


                        IF ((tab_site(i,j,k)==0).and.(tab_site(ip,j,k)==0)) THEN
                              S(i,j,k,1)=1.
                              tab_site(i,j,k)=1
                              tab_site(ip,j,k)=1
                              number_line_ini=number_line_ini+1
                        END IF
                  
                        IF ((tab_site(i,j,k)==0).and.(tab_site(i,jp,k)==0)) THEN
                              S(i,j,k,2)=1.
                              tab_site(i,j,k)=1
                              tab_site(i,jp,k)=1
                              number_line_ini=number_line_ini+1
                        END IF
                       
                        END IF
                        !!!------------------------------------------------------

                        END IF

                  END DO
            END DO

      END DO

      WRITE(*,*) number_line_ini

      END IF

!!!=================================================================================================
!!!! if used GS1: OX

      IF (CONFIG_INI == 'GS1') THEN
      number_line_ini=0
      S(:,:,:,:)=0.
      tab_site(:,:,:)=0
      
      WRITE(*,*)'Used GS1'

      DO i=1,natx   
      ip=i+1-(i/natx)*natx
            DO j=1,naty
            jp=j+1-(j/naty)*naty
                  DO k=1,natz
                  kp=k+1-(k/natz)*natz

                        IF ((tab_site(i,j,k)==0).and.(tab_site(ip,j,k)==0)&
                           .and.(number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                              S(i,j,k,1)=1.
                              tab_site(i,j,k)=1
                              tab_site(ip,j,k)=1
                              number_line_ini=number_line_ini+1
                        END IF

                  END DO
            END DO
      END DO

      WRITE(*,*) 'number_line_ini=',number_line_ini

      END IF
!!!=================================================================================================
!!! GS dimer spin khi D be, r lon (Ex D=0.45, r=3.2)
!!! 3 lop song song, 3 lop nguoc lai cheo nhau

      IF (CONFIG_INI == 'GS2') THEN
      number_line_ini=0
      S(:,:,:,:)=0.
      tab_site(:,:,:)=0
      
      WRITE(*,*)'Used GS2'

      DO i=1,natx   
      ip=i+1-(i/natx)*natx
      ipp=ip+1-(ip/natx)*natx
            DO j=1,naty
            jp=j+1-(j/naty)*naty
                  DO k=1,natz
                  kp=k+1-(k/natz)*natz

                        IF ((tab_site(i,j,k)==0).and.(tab_site(ip,j,k)==0)&
                       .and.((mod(k,6)==1).or.(mod(k,6)==2).or.(mod(k,6)==3))&
                       .and.(number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                              S(i,j,k,1)=1.
                              tab_site(i,j,k)=1
                              tab_site(ip,j,k)=1
                              number_line_ini=number_line_ini+1
                        END IF

                        IF ((tab_site(ip,j,k)==0).and.(tab_site(ipp,j,k)==0)&
                       .and.((mod(k,6)==4).or.(mod(k,6)==5).or.(mod(k,6)==0))&
                       .and.(number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                              S(ip,j,k,1)=-1.
                              tab_site(ip,j,k)=1
                              tab_site(ipp,j,k)=1
                              number_line_ini=number_line_ini+1
                        END IF


                  END DO
            END DO
      END DO

      WRITE(*,*) 'number_line_ini=',number_line_ini

      END IF

!!!=================================================================================================
!!! GS dimer spin khi D lon, r be (Ex D=0.6, r=1.5)
!!! 1 lop song song OX, 1 song song OX nhung cheo nhau

      IF (CONFIG_INI == 'GS3') THEN
      number_line_ini=0
      S(:,:,:,:)=0.
      tab_site(:,:,:)=0
      
      WRITE(*,*)'Used GS3'

      DO i=1,natx   
      ip=i+1-(i/natx)*natx
      ipp=ip+1-(ip/natx)*natx
            DO j=1,naty
            jp=j+1-(j/naty)*naty
                  DO k=1,natz
                  kp=k+1-(k/natz)*natz

                        IF ((tab_site(i,j,k)==0).and.(tab_site(ip,j,k)==0)&
                       .and.(mod(k,2)==1)&
                       .and.(number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                              S(i,j,k,1)=1.
                              tab_site(i,j,k)=1
                              tab_site(ip,j,k)=1
                              number_line_ini=number_line_ini+1
                        END IF

                        IF ((tab_site(ip,j,k)==0).and.(tab_site(ipp,j,k)==0)&
                       .and.(mod(k,2)==0)&
                       .and.(number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                              S(ip,j,k,1)=1.
                              tab_site(ip,j,k)=1
                              tab_site(ipp,j,k)=1
                              number_line_ini=number_line_ini+1
                        END IF


                  END DO
            END DO
      END DO

      WRITE(*,*) 'number_line_ini=',number_line_ini

      END IF


!!!=================================================================================================
!!! GS dimer spin khi D lon, r trung binh (Ex D=0.5; 0.55; D=0.6 va r=2.1 hoac r=2.3 )
!!! Giong GS2 cua he thu nhat, tung lop xyxyxy cheo cheo

      IF (CONFIG_INI == 'GS4') THEN

      WRITE(*,*)'Used GS4'

      number_line_ini=0
      S(:,:,:,:)=0.
      tab_site(:,:,:)=0
      
      DO i=1,natx   
            ip=i+1-(i/natx)*natx
            DO j=1,naty
                  jp=j+1-(j/naty)*naty
                  DO k=1,natz
                  kp=k+1-(k/natz)*natz

                        IF ((tab_site(i,j,k)==0).and.(tab_site(ip,j,k)==0).and.(mod(i+j,4)==2).and.&
                        (number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                              S(i,j,k,1)=1.
                              tab_site(i,j,k)=1
                              tab_site(ip,j,k)=1
                              number_line_ini=number_line_ini+1
                        END IF

                        IF ((tab_site(i,j,k)==0).and.(tab_site(i,jp,k)==0).and.(mod(i+j,4)==0).and.&
                        (number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                              S(i,j,k,2)=1.
                              tab_site(i,j,k)=1
                              tab_site(i,jp,k)=1
                              number_line_ini=number_line_ini+1
                        END IF

                  END DO
            END DO
      END DO

      WRITE(*,*)'number_line_ini:', number_line_ini

      END IF

!!!!================================================================================================
      IF (CONFIG_INI == 'GS40') THEN

      WRITE(*,*)'Used GS40'

      number_line_ini=0
      S(:,:,:,:)=0.
      tab_site(:,:,:)=0
      
      DO i=1,natx   
            ip=i+1-(i/natx)*natx
            DO j=1,naty
                  jp=j+1-(j/naty)*naty
                  DO k=1,natz
                  kp=k+1-(k/natz)*natz

                        IF ((tab_site(i,j,k)==0).and.(tab_site(ip,j,k)==0).and.(mod(i+j,4)==2).and.&
                        (number_line_ini< int(real(natx*naty*natz)*0.5))) THEN

                        IF (((i/=8).or.(j/=6)).and.((i/=9).or.(j/=5)).and.((i/=10).or.(j/=4))&
                      .and.((i/=10).or.(j/=8)).and.((i/=11).or.(j/=7)).and.((i/=12).or.(j/=6))) THEN
                              S(i,j,k,1)=1.
                              tab_site(i,j,k)=1
                              tab_site(ip,j,k)=1
                              number_line_ini=number_line_ini+1
                        END IF
                        END IF

                        IF ((tab_site(i,j,k)==0).and.(tab_site(i,jp,k)==0).and.(mod(i+j,4)==0).and.&
                        (number_line_ini< int(real(natx*naty*natz)*0.5))) THEN

                        IF (((i/=9).or.(j/=7)).and.((i/=10).or.(j/=6)).and.((i/=11).or.(j/=5))&
                      .and.((i/=12).or.(j/=4)).and.((i/=11).or.(j/=9)).and.((i/=12).or.(j/=8))) THEN

                              S(i,j,k,2)=1.
                              tab_site(i,j,k)=1
                              tab_site(i,jp,k)=1
                              number_line_ini=number_line_ini+1
                        END IF

                        END IF

                  END DO
            END DO
      END DO

      WRITE(*,*)'number_line_ini:', number_line_ini

      END IF

!!!=================================================================================================
!!! GS dimer spin khi D lon, r trung binh (Ex D=0.6, r=2.5)
!!! 1 lop song song OX, 2 lop AF tung doi 

      IF (CONFIG_INI == 'GS5') THEN
      number_line_ini=0
      S(:,:,:,:)=0.
      tab_site(:,:,:)=0
      
      WRITE(*,*)'Used GS5'

      DO i=1,natx   
      ip=i+1-(i/natx)*natx
      ipp=ip+1-(ip/natx)*natx
            DO j=1,naty
            jp=j+1-(j/naty)*naty
                  DO k=1,natz
                  kp=k+1-(k/natz)*natz
                  !!!! Doc theo truc oy
                        IF ((tab_site(i,j,k)==0).and.(tab_site(i,jp,k)==0)&
                       .and.(mod(i,6)==1)&
                       .and.(number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                              IF ((mod(j,4)==1)) THEN                  
                              S(i,j,k,2)=1.
                              tab_site(i,j,k)=1
                              tab_site(i,jp,k)=1
                              number_line_ini=number_line_ini+1
                              END IF

                              IF ((mod(j,4)==3)) THEN                  
                              S(i,j,k,2)=-1.
                              tab_site(i,j,k)=1
                              tab_site(i,jp,k)=1
                              number_line_ini=number_line_ini+1
                              END IF

                        END IF

                        IF ((tab_site(i,j,k)==0).and.(tab_site(i,jp,k)==0)&
                       .and.(mod(i,6)==4)&
                       .and.(number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                              IF ((mod(j,4)==1)) THEN                  
                              S(i,j,k,2)=-1.
                              tab_site(i,j,k)=1
                              tab_site(i,jp,k)=1
                              number_line_ini=number_line_ini+1
                              END IF

                              IF ((mod(j,4)==3)) THEN                  
                              S(i,j,k,2)=1.
                              tab_site(i,j,k)=1
                              tab_site(i,jp,k)=1
                              number_line_ini=number_line_ini+1
                              END IF

                        END IF

                  !!! Doc theo truc Ox
                        IF ((tab_site(i,j,k)==0).and.(tab_site(ip,j,k)==0)&
                       .and.(mod(i,6)==2)&
                       .and.(number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                              IF ((mod(j,4)==0).or.(mod(j,4)==1)) THEN                  
                              S(i,j,k,1)=-1.
                              tab_site(i,j,k)=1
                              tab_site(ip,j,k)=1
                              number_line_ini=number_line_ini+1
                              END IF

                              IF ((mod(j,4)==2).or.(mod(j,4)==3)) THEN                  
                              S(i,j,k,1)=1.
                              tab_site(i,j,k)=1
                              tab_site(ip,j,k)=1
                              number_line_ini=number_line_ini+1
                              END IF

                        END IF

                        IF ((tab_site(i,j,k)==0).and.(tab_site(ip,j,k)==0)&
                       .and.(mod(i,6)==5)&
                       .and.(number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                              IF ((mod(j,4)==0).or.(mod(j,4)==1)) THEN                  
                              S(i,j,k,1)=1.
                              tab_site(i,j,k)=1
                              tab_site(ip,j,k)=1
                              number_line_ini=number_line_ini+1
                              END IF

                              IF ((mod(j,4)==2).or.(mod(j,4)==3)) THEN                  
                              S(i,j,k,1)=-1.
                              tab_site(i,j,k)=1
                              tab_site(ip,j,k)=1
                              number_line_ini=number_line_ini+1
                              END IF

                        END IF



                  END DO
            END DO
      END DO

      WRITE(*,*) 'number_line_ini=',number_line_ini

      END IF

!!!=================================================================================================
!!!! if used GS6: 
!!! GS dimer spin khi D lon, r lon (Ex D=0.6, r=3.2 hoac r=3.6)
!!! 

      IF (CONFIG_INI == 'GS6') THEN
      number_line_ini=0
      S(:,:,:,:)=0.
      tab_site(:,:,:)=0
      
      WRITE(*,*)'Used GS6'

      DO i=1,natx   
      ip=i+1-(i/natx)*natx
            DO j=1,naty
            jp=j+1-(j/naty)*naty
                  DO k=1,natz
                  kp=k+1-(k/natz)*natz

                  !!! Doc theo OX
                  IF (((i==1).and.(j==1)).or.((i==7).and.(j==1))&
                 .or.((i==4).and.(j==4)).or.((i==10).and.(j==4))&
                 .or.((i==3).and.(j==5)).or.((i==9).and.(j==5))&
                 .or.((i==6).and.(j==8)).or.((i==12).and.(j==8))&
                 .or.((i==5).and.(j==9)).or.((i==11).and.(j==9))&
                 .or.((i==2).and.(j==12)).or.((i==8).and.(j==12))) THEN

                  IF ((tab_site(i,j,k)==0).and.(tab_site(ip,j,k)==0)&
                  .and.(number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                        S(i,j,k,1)=1.
                        tab_site(i,j,k)=1
                        tab_site(ip,j,k)=1
                        number_line_ini=number_line_ini+1
                  END IF

                  END IF

                  !!! Doc theo -OX
                  IF (((i==4).and.(j==1)).or.((i==10).and.(j==1))&
                 .or.((i==1).and.(j==4)).or.((i==7).and.(j==4))&
                 .or.((i==6).and.(j==5)).or.((i==12).and.(j==5))&
                 .or.((i==3).and.(j==8)).or.((i==9).and.(j==8))&
                 .or.((i==2).and.(j==9)).or.((i==8).and.(j==9))&
                 .or.((i==5).and.(j==12)).or.((i==11).and.(j==12))) THEN

                  IF ((tab_site(i,j,k)==0).and.(tab_site(ip,j,k)==0)&
                  .and.(number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                        S(i,j,k,1)=-1.
                        tab_site(i,j,k)=1
                        tab_site(ip,j,k)=1
                        number_line_ini=number_line_ini+1
                  END IF

                  END IF


                  !!! Doc theo OY (truong hop 1 hang co 4 cot)
                  IF ((((j==1).or.(j==3)).and.(mod(i,6)==3))&
                 .or.(((j==5).or.(j==7)).and.(mod(i,6)==5))&
                 .or.(((j==9).or.(j==11)).and.(mod(i,6)==1))) THEN

                  IF ((tab_site(i,j,k)==0).and.(tab_site(i,jp,k)==0)&
                  .and.(number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                        S(i,j,k,2)=1.
                        tab_site(i,j,k)=1
                        tab_site(i,jp,k)=1
                        number_line_ini=number_line_ini+1
                  END IF

                  END IF

                  !!! Doc theo -OY (truong hop 1 hang co 4 cot)
                  IF ((((j==1).or.(j==3)).and.(mod(i,6)==0))&
                 .or.(((j==5).or.(j==7)).and.(mod(i,6)==2))&
                 .or.(((j==9).or.(j==11)).and.(mod(i,6)==4))) THEN

                  IF ((tab_site(i,j,k)==0).and.(tab_site(i,jp,k)==0)&
                  .and.(number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                        S(i,j,k,2)=-1.
                        tab_site(i,j,k)=1
                        tab_site(i,jp,k)=1
                        number_line_ini=number_line_ini+1
                  END IF

                  END IF

                  !!! Doc theo OY (truong hop 1 hang co cac cot tung cap)
                  IF (((j==2).and.((i==2).or.(i==4).or.(i==8).or.(i==10)))&
                 .or.((j==6).and.((i==4).or.(i==6).or.(i==10).or.(i==12)))&
                 .or.((j==10).and.((i==2).or.(i==6).or.(i==8).or.(i==12)))) THEN

                  IF ((tab_site(i,j,k)==0).and.(tab_site(i,jp,k)==0)&
                  .and.(number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                        S(i,j,k,2)=1.
                        tab_site(i,j,k)=1
                        tab_site(i,jp,k)=1
                        number_line_ini=number_line_ini+1
                  END IF

                  END IF

                  !!! Doc theo - OY (truong hop 1 hang co cac cot tung cap)
                  IF (((j==2).and.((i==1).or.(i==5).or.(i==7).or.(i==11)))&
                 .or.((j==6).and.((i==1).or.(i==3).or.(i==7).or.(i==9)))&
                 .or.((j==10).and.((i==3).or.(i==5).or.(i==9).or.(i==11)))) THEN

                  IF ((tab_site(i,j,k)==0).and.(tab_site(i,jp,k)==0)&
                  .and.(number_line_ini< int(real(natx*naty*natz)*P_config))) THEN
                        S(i,j,k,2)=-1.
                        tab_site(i,j,k)=1
                        tab_site(i,jp,k)=1
                        number_line_ini=number_line_ini+1
                  END IF

                  END IF


                  END DO
            END DO
      END DO

      WRITE(*,*) 'number_line_ini=',number_line_ini

      END IF

      END SUBROUTINE load_config_ini

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! WRITE initial position configuration in 3D
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE write_config_ini_3D()
      IMPLICIT NONE

          
      OPEN(unit=12,file='config_ini_3D/config_ini_3D_DIMER.pdb')
      OPEN(unit=301,file='config_ini_3D/config_ini_python.dat')

      DO k=1,natz 
            DO j=1,naty
                  DO i=1,natx

                        x=real(i-1)*aaa
                        y=real(j-1)*aaa
                        z=real(k-1)*aaa

                        
                        IF (int(S(i,j,k,1))==1) THEN 
                              spinAt='Au'
                              WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                                   spinAt,x,y,z,nul

                        ELSE
                              
                              IF (int(S(i,j,k,2))==1) THEN 
                              spinAt='Cu'
                              WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                                  spinAt,x,y,z,nul

                              ELSE
                              
                              IF (int(S(i,j,k,3))==1) THEN 
                              spinAt='Mn'
                              WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                                   spinAt,x,y,z,nul
                              ELSE

                              spinAt='H'
                              WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                                   spinAt,x,y,z,nul

                              END IF
                              END IF
                        END IF

                        xspin=real(i) ; yspin=real(j) ; zspin=real(k)

                        IF (S(i,j,k,1)<0) xspin=i+1.
                        IF (S(i,j,k,2)<0) yspin=j+1.
                        IF (S(i,j,k,3)<0) zspin=k+1.

                        WRITE(301,*)xspin,yspin,zspin,S(i,j,k,1),0.8*S(i,j,k,2),0.95*S(i,j,k,3)

                  END DO
            END DO
      END DO

      CLOSE(12)
      CLOSE(301)

      END SUBROUTINE write_config_ini_3D


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! WRITE initial position configuration in 3D
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE write_config_3D()
      IMPLICIT NONE

      number=iT+10000000
      WRITE(tmp,'(I8)') number

      name='_DIMER_'//TRIM(tmp)
      
      OPEN(unit=13,file='config_3D/config_3D'//trim(name)//'.pdb')

      OPEN(unit=302,file='config_3D/config_python'//trim(name)//'.dat')
      DO k=1,natz 
            DO j=1,naty
                  DO i=1,natx

                        x=real(i-1)*aaa
                        y=real(j-1)*aaa
                        z=real(k-1)*aaa

                        
                        IF (int(S(i,j,k,1))==1) THEN 
                              spinAt='Au'
                              WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                                   spinAt,x,y,z,nul

                        ELSE
                              
                              IF (int(S(i,j,k,2))==1) THEN 
                              spinAt='Cu'
                              WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                                  spinAt,x,y,z,nul

                              ELSE
                              
                              IF (int(S(i,j,k,3))==1) THEN 
                              spinAt='Mn'
                              WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                                   spinAt,x,y,z,nul
                              ELSE

                              spinAt='H'
                              WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                                   spinAt,x,y,z,nul

                              END IF
                              END IF
                        END IF

                        xspin=real(i) ; yspin=real(j) ; zspin=real(k)

                        IF (S(i,j,k,1)<0) xspin=i+1.
                        IF (S(i,j,k,2)<0) yspin=j+1.
                        IF (S(i,j,k,3)<0) zspin=k+1.

                        WRITE(302,*)xspin,yspin,zspin,S(i,j,k,1),0.8*S(i,j,k,2),0.95*S(i,j,k,3)

                  END DO
            END DO
      END DO


      CLOSE(13)
      CLOSE(302)

      END SUBROUTINE write_config_3D


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE cacul_number_line()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE cacul_number_line()
      IMPLICIT NONE
      
      n_line_x1=0. ; n_line_x2=0.
      n_line_y1=0. ; n_line_y2=0.
      n_line_z1=0. ; n_line_z2=0.
      n_line_total=0.

      DO k=1,natz
            DO j=1,naty
                  DO i=1,natx

                  IF (S(i,j,k,1)==1.)     n_line_x1=n_line_x1+1.
                  IF (S(i,j,k,1)==-1.)    n_line_x2=n_line_x2+1.
                  IF (S(i,j,k,2)==1.)     n_line_y1=n_line_y1+1.
                  IF (S(i,j,k,2)==-1.)    n_line_y2=n_line_y2+1.
                  IF (S(i,j,k,3)==1.)     n_line_z1=n_line_z1+1.
                  IF (S(i,j,k,3)==-1.)    n_line_z2=n_line_z2+1.

                  ENDDO
            ENDDO
       ENDDO

      n_line_total = n_line_x1+n_line_x2+n_line_y1+n_line_y2+n_line_z1+n_line_z2
                        
      WRITE(Ligne14,*) n_line_total,n_line_x1,n_line_x2,n_line_y1,n_line_y2,n_line_z1,n_line_z2

      WRITE(14,'(a)') trim(Ligne14)

      END SUBROUTINE cacul_number_line

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE test_update_line()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE test_update_line()
      IMPLICIT NONE
      
      DO k=1,natz
            DO j=1,naty
                  DO i=1,natx

                  compt_update=0

                  IF (S(i,j,k,1)==1.)     compt_update=compt_update+1
                  IF (S(i,j,k,1)==-1.)    compt_update=compt_update+1
                  IF (S(i,j,k,2)==1.)     compt_update=compt_update+1
                  IF (S(i,j,k,2)==-1.)    compt_update=compt_update+1
                  IF (S(i,j,k,3)==1.)     compt_update=compt_update+1
                  IF (S(i,j,k,3)==-1.)    compt_update=compt_update+1

                 IF (compt_update>1) THEN
                        WRITE(*,*)'ERROS IN UPDATE'
                 END IF

                  END DO
            END DO
      END DO

      END SUBROUTINE test_update_line

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE equi_reseau1()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE equi_reseau1()
      IMPLICIT NONE

      DO i_loop=1,n_equi_reseau1
            CALL equi_reseau()
            !CALL value_thermal()

      END DO

      END SUBROUTINE equi_reseau1

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE equi_reseau()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE equi_reseau()
      IMPLICIT NONE

      EN_tmp=0. ;  EN_tmp1=0. ;  EN_tmp2=0.     
      EN_tmpx=0.;  EN_tmpy=0. ;  EN_tmpz=0.
      EN_tmp1x=0.; EN_tmp1y=0. ; EN_tmp1z=0.
            
      DO k=1,natz
            kp=k+1-(k/natz)*natz
            km=k-1+(1/k)*natz
            kpp=kp+1-(kp/natz)*natz
            DO j=1,naty
                  jp=j+1-(j/naty)*naty
                  jm=j-1+(1/j)*naty
                  jpp=jp+1-(jp/naty)*naty
                  DO i=1,natx
                        ip=i+1-(i/natx)*natx
                        im=i-1+(1/i)*natx
                        ipp=ip+1-(ip/natx)*natx
                        
                        IF (tab_site(i,j,k)==1) THEN

!!!====================================================================================================
                        !!!---------------------------------------------------------------                        
                        IF ((S(i,j,k,1)==1.).or.(S(i,j,k,1)==-1.)) THEN
                                               
                        CALL value_EN_D_tmp1()

                        EN_tmp1x=-J1_tmp*S(i,j,k,1)*(S(i,jp,k,1)+S(i,jm,k,1)+S(i,j,kp,1)+S(i,j,km,1))&
                                 +D_tmp*EN_D_tmp1

                              CALL random_number(rdn_mtp) 
                              IF (exp(2.*EN_tmp1x/T) > rdn_mtp) THEN
                                    S(i,j,k,1)=-S(i,j,k,1)
                                    EN_tmp1x=-EN_tmp1x
                              END IF                    

                             CALL search_site_void_x()
                              
                             IF (n_etat/=0) THEN
                                   
                                    CALL random_number(rdn_etat)

                                    i0_etat=int(rdn_etat*n_etat)+1
                                    i_etat=name_etat(i0_etat)
                              
                                    spin1=S(i,j,k,1)
                                    
                                    CALL random_number(rdn_spin)
                                    IF (rdn_spin<0.5) THEN
                                          spin2=1.
                                    ELSE
                                          spin2=-1.
                                    END IF


                                    CALL value_Hx_tmp2()
                                    CALL value_ENx_D_tmp2()

                                    EN_tmp2= spin2*Hx_tmp2(i_etat)+D_tmp*ENx_D_tmp2(i_etat)

                                    CALL random_number(rdn_mtp) 
                                    IF (exp(-(EN_tmp2-EN_tmp1x)/T) > rdn_mtp) THEN
                                          CALL update_line_x()
                                          EN_tmp1x=EN_tmp2
                                    END IF

                              END IF
                              
                       END IF   
      
    
                       !!!---------------------------------------------------------------
                              
                       IF ((S(i,j,k,2)==1.).or.(S(i,j,k,2)==-1.)) THEN
                       CALL value_EN_D_tmp1()

                       EN_tmp1y=-J1_tmp*S(i,j,k,2)*(S(ip,j,k,2)+S(im,j,k,2)+S(i,j,kp,2)+S(i,j,km,2))&
                               +D_tmp*EN_D_tmp1
                              
                              CALL random_number(rdn_mtp) 
                              IF (exp(2.*EN_tmp1y/T) > rdn_mtp) THEN
                                    S(i,j,k,2)=-S(i,j,k,2)
                                    EN_tmp1y=-EN_tmp1y
                              END IF 

                              CALL search_site_void_y()
                        
                              IF (n_etat/=0) THEN

                                    CALL random_number(rdn_etat)

                                    i0_etat=int(rdn_etat*n_etat)+1
                                    i_etat=name_etat(i0_etat)

                                    spin1=S(i,j,k,2)
                                    
                                    CALL random_number(rdn_spin)
                                    IF (rdn_spin<0.5) THEN
                                          spin2=1.
                                    ELSE
                                          spin2=-1.
                                    END IF

                                    CALL value_Hy_tmp2()
                                    CALL value_ENy_D_tmp2()
                                    EN_tmp2=spin2*Hy_tmp2(i_etat)+D_tmp*ENy_D_tmp2(i_etat)

                                    CALL random_number(rdn_mtp) 
                                    IF (exp(-(EN_tmp2-EN_tmp1y)/T) > rdn_mtp) THEN
                                          CALL update_line_y()
                                          EN_tmp1y=EN_tmp2
                                    END IF                                  
                                    
                              END IF
                        END IF
                                                                  
                        !!!---------------------------------------------------------------                   
                        
                        IF ((S(i,j,k,3)==1.).or.(S(i,j,k,3)==-1.)) THEN
                        CALL value_EN_D_tmp1()
      
                        EN_tmp1z=-J1_tmp*(S(ip,j,k,3)+S(im,j,k,3)+S(i,jp,k,3)+S(i,jm,k,3))&
                                 +D_tmp*EN_D_tmp1

                              CALL random_number(rdn_mtp) 
                              IF (exp(2.*EN_tmp1z/T) > rdn_mtp) THEN
                                    S(i,j,k,3)=-S(i,j,k,3)
                                    EN_tmp1z=-EN_tmp1z
                              END IF 


                              CALL search_site_void_z()
                              IF (n_etat/=0) THEN

                                    CALL random_number(rdn_etat)

                                    i0_etat=int(rdn_etat*n_etat)+1
                                    i_etat=name_etat(i0_etat)

                                    spin1=S(i,j,k,3)
                                    
                                    CALL random_number(rdn_spin)
                                    IF (rdn_spin<0.5) THEN
                                          spin2=1.
                                    ELSE
                                          spin2=-1.
                                    END IF

                                    CALL value_Hz_tmp2()
                                    CALL value_ENz_D_tmp2()
                                    EN_tmp2=spin2*Hz_tmp2(i_etat)+D_tmp*ENz_D_tmp2(i_etat)

                                    CALL random_number(rdn_mtp) 
                                    IF (exp(-(EN_tmp2-EN_tmp1z)/T) > rdn_mtp) THEN
                                          CALL update_line_z()
                                          EN_tmp1z=EN_tmp2
                                    END IF
                                    

                              END IF


                          END IF
!!!====================================================================================================
                       END IF

                          
                   ENDDO
            ENDDO
      ENDDO
     
      END SUBROUTINE equi_reseau

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_thermal()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE value_thermal()
      IMPLICIT NONE

      EN_tmp=0. ;  EN_tmp1=0. ;  EN_tmp2=0.     
      EN_tmpx=0.;  EN_tmpy=0. ;  EN_tmpz=0.
      EN_tmp1x=0.; EN_tmp1y=0. ; EN_tmp1z=0.
            
      n_line(:)=0.
      order_parameter=0.
      order_parameter2=0.

      DO k=1,natz
            kp=k+1-(k/natz)*natz
            km=k-1+(1/k)*natz
            kpp=kp+1-(kp/natz)*natz
            DO j=1,naty
                  jp=j+1-(j/naty)*naty
                  jm=j-1+(1/j)*naty
                  jpp=jp+1-(jp/naty)*naty
                  DO i=1,natx
                        ip=i+1-(i/natx)*natx
                        im=i-1+(1/i)*natx
                        ipp=ip+1-(ip/natx)*natx
                        
                        IF (tab_site(i,j,k)==1) THEN

!!!====================================================================================================
                        !!!---------------------------------------------------------------                        
                        IF ((S(i,j,k,1)==1.).or.(S(i,j,k,1)==-1.)) THEN
                                               
                        CALL value_EN_D_tmp1()

                        EN_tmp1x=-J1_tmp*S(i,j,k,1)*(S(i,jp,k,1)+S(i,jm,k,1)+S(i,j,kp,1)+S(i,j,km,1))&
                                 +D_tmp*EN_D_tmp1

                              CALL random_number(rdn_mtp) 
                              IF (exp(2.*EN_tmp1x/T) > rdn_mtp) THEN
                                    S(i,j,k,1)=-S(i,j,k,1)
                                    EN_tmp1x=-EN_tmp1x
                              END IF                    

                             CALL search_site_void_x()
                              
                             IF (n_etat/=0) THEN
                                   
                                    CALL random_number(rdn_etat)

                                    i0_etat=int(rdn_etat*n_etat)+1
                                    i_etat=name_etat(i0_etat)
                              
                                    spin1=S(i,j,k,1)
                                    
                                    CALL random_number(rdn_spin)
                                    IF (rdn_spin<0.5) THEN
                                          spin2=1.
                                    ELSE
                                          spin2=-1.
                                    END IF


                                    CALL value_Hx_tmp2()
                                    CALL value_ENx_D_tmp2()

                                    EN_tmp2= spin2*Hx_tmp2(i_etat)+D_tmp*ENx_D_tmp2(i_etat)

                                    CALL random_number(rdn_mtp) 
                                    IF (exp(-(EN_tmp2-EN_tmp1x)/T) > rdn_mtp) THEN
                                          CALL update_line_x()
                                          EN_tmp1x=EN_tmp2
                                    END IF

                              END IF

                              EN_tmpx=EN_tmpx+EN_tmp1x
                              
                       END IF   
      
    
                       !!!---------------------------------------------------------------
                              
                       IF ((S(i,j,k,2)==1.).or.(S(i,j,k,2)==-1.)) THEN
                       CALL value_EN_D_tmp1()

                       EN_tmp1y=-J1_tmp*S(i,j,k,2)*(S(ip,j,k,2)+S(im,j,k,2)+S(i,j,kp,2)+S(i,j,km,2))&
                               +D_tmp*EN_D_tmp1
                              
                              CALL random_number(rdn_mtp) 
                              IF (exp(2.*EN_tmp1y/T) > rdn_mtp) THEN
                                    S(i,j,k,2)=-S(i,j,k,2)
                                    EN_tmp1y=-EN_tmp1y
                              END IF 

                              CALL search_site_void_y()
                        
                              IF (n_etat/=0) THEN

                                    CALL random_number(rdn_etat)

                                    i0_etat=int(rdn_etat*n_etat)+1
                                    i_etat=name_etat(i0_etat)

                                    spin1=S(i,j,k,2)
                                    
                                    CALL random_number(rdn_spin)
                                    IF (rdn_spin<0.5) THEN
                                          spin2=1.
                                    ELSE
                                          spin2=-1.
                                    END IF

                                    CALL value_Hy_tmp2()
                                    CALL value_ENy_D_tmp2()
                                    EN_tmp2=spin2*Hy_tmp2(i_etat)+D_tmp*ENy_D_tmp2(i_etat)

                                    CALL random_number(rdn_mtp) 
                                    IF (exp(-(EN_tmp2-EN_tmp1y)/T) > rdn_mtp) THEN
                                          CALL update_line_y()
                                          EN_tmp1y=EN_tmp2
                                    END IF                                  
                                    
                              END IF
                               EN_tmpy=EN_tmpy+EN_tmp1y
                        END IF
                                                                  
                        !!!---------------------------------------------------------------                   
                        
                        IF ((S(i,j,k,3)==1.).or.(S(i,j,k,3)==-1.)) THEN
                        CALL value_EN_D_tmp1()
      
                        EN_tmp1z=-J1_tmp*(S(ip,j,k,3)+S(im,j,k,3)+S(i,jp,k,3)+S(i,jm,k,3))&
                                 +D_tmp*EN_D_tmp1

                              CALL random_number(rdn_mtp) 
                              IF (exp(2.*EN_tmp1z/T) > rdn_mtp) THEN
                                    S(i,j,k,3)=-S(i,j,k,3)
                                    EN_tmp1z=-EN_tmp1z
                              END IF 


                              CALL search_site_void_z()
                              IF (n_etat/=0) THEN

                                    CALL random_number(rdn_etat)

                                    i0_etat=int(rdn_etat*n_etat)+1
                                    i_etat=name_etat(i0_etat)

                                    spin1=S(i,j,k,3)
                                    
                                    CALL random_number(rdn_spin)
                                    IF (rdn_spin<0.5) THEN
                                          spin2=1.
                                    ELSE
                                          spin2=-1.
                                    END IF

                                    CALL value_Hz_tmp2()
                                    CALL value_ENz_D_tmp2()
                                    EN_tmp2=spin2*Hz_tmp2(i_etat)+D_tmp*ENz_D_tmp2(i_etat)

                                    CALL random_number(rdn_mtp) 
                                    IF (exp(-(EN_tmp2-EN_tmp1z)/T) > rdn_mtp) THEN
                                          CALL update_line_z()
                                          EN_tmp1z=EN_tmp2
                                    END IF
                                    

                              END IF

                              EN_tmpz=EN_tmpz+EN_tmp1z

                          END IF
!!!====================================================================================================
                       END IF

                          
                   ENDDO
            ENDDO
      ENDDO


      !IF ((CONFIG_INI == 'GS1').or.(CONFIG_INI == 'GS3')) THEN   
      
      !n_line_x1=0. ; n_line_x2=0. ; n_line_y1=0. ; n_line_y2=0. ; n_line_z1=0.; n_line_z2=0.

      !DO k=1,natz
          !  DO j=1,naty
                !  DO i=1,natx

                 ! IF (S(i,j,k,1)==1.)     n_line_x1=n_line_x1+1.
                 ! IF (S(i,j,k,1)==-1.)    n_line_x2=n_line_x2+1.
                 ! IF (S(i,j,k,2)==1.)     n_line_y1=n_line_y1+1.
                 ! IF (S(i,j,k,2)==-1.)    n_line_y2=n_line_y2+1.
                 ! IF (S(i,j,k,3)==1.)     n_line_z1=n_line_z1+1.
                 ! IF (S(i,j,k,3)==-1.)    n_line_z2=n_line_z2+1.

                 ! ENDDO
            !ENDDO
       !ENDDO

    ! order_parameter=abs(6.*max(n_line_x1,n_line_x2,n_line_y1,n_line_y2,n_line_z1,n_line_z2)&
                      !  /number_line_ini-1.)/5.     
     
    ! END IF


      !!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !!!! Tinh M khi GS1 nhu truong hop Modern Potts, khong phan biet chieu spin +1 hay -1

      IF ((CONFIG_INI == 'GS1').or.(CONFIG_INI == 'GS3')) THEN     
      DO k=1,natz
            DO j=1,naty
                  DO i=1,natx   
       
                  n_line(1)=n_line(1)+abs(S(i,j,k,1))
                  n_line(2)=n_line(2)+abs(S(i,j,k,2))
                  n_line(3)=n_line(3)+abs(S(i,j,k,3))

                  END DO
           END DO
      END DO

     order_parameter=abs(3.*max(n_line(1),n_line(2),n_line(3))/number_line_ini-1.)/2.     
     
    END IF


    !!!===========================================================================================
     IF ((CONFIG_INI == 'GS4').or.(CONFIG_INI == 'GS40')) THEN
            CALL value_Mz42()                
     END IF
      
     ! IF ((CONFIG_INI == 'GS3').or.(CONFIG_INI == 'GS3333')) THEN
       !     CALL value_Mz3()                 
      !      order_parameter=Mz
     ! END IF    
      
      !IF ((CONFIG_INI == 'GS4').or.(CONFIG_INI == 'GS4444')) THEN
      !      CALL value_Mz4()                 
      !      order_parameter=Mz
      !END IF 

      EN_tmp=EN_tmpx+EN_tmpy+EN_tmpz

      energy=EN_tmp/(2.*number_line_ini)

      energy_2=energy**2.
      order_parameter_2=order_parameter**2.
      !order_parameter2_2=order_parameter2**2.
      
      
     ! WRITE(*,*)'energy =',energy,'order_parameter=',order_parameter

      END SUBROUTINE value_thermal

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_EN_D_tmp1()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE value_EN_D_tmp1()
      IMPLICIT NONE

      !!!----------------------------------------------------------
      IF ((S(i,j,k,1)==1.).or.(S(i,j,k,1)==-1.)) THEN
            x=real(i)+0.5 ; y=real(j) ; z=real(k)
      END IF

      IF ((S(i,j,k,2)==1.).or.(S(i,j,k,2)==-1.)) THEN
            x=real(i) ; y=real(j)+0.5 ; z=real(k)
      END IF

      IF ((S(i,j,k,3)==1.).or.(S(i,j,k,3)==-1.)) THEN
            x=real(i) ; y=real(j) ; z=real(k)+0.5
      END IF
      !!!========================

      EN_D_tmp1=0.   

      DO deli=-r01,r01
      DO delj=-r01,r01
      DO delk=-r01,r01

            i1=i+deli ; j1=j+delj ; k1=k+delk
            
            IF (i1<1)    i1=natx+i1
            IF (i1>natx) i1=i1-natx
            IF (j1<1)    j1=naty+j1                        
            IF (j1>naty) j1=j1-naty                        
            IF (k1<1)    k1=natz+k1                       
            IF (k1>natz) k1=k1-natz      
      !!!------------------------------------------------------------------------------------------
      IF ((S(i1,j1,k1,1)==1.).or.(S(i1,j1,k1,1)==-1.)) THEN
            x1=real(i+deli)+0.5 ; y1=real(j+delj) ; z1=real(k+delk)

            r=sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y)+(z1-z)*(z1-z))
                        
            IF ((0.<r).and.(r<=r0)) THEN                

            EN_D_tmp1=EN_D_tmp1+S(i,j,k,1)*S(i1,j1,k1,1)/(r**3.)&     
            -3.*((x1-x)*S(i,j,k,1)+(y1-y)*S(i,j,k,2)+(z1-z)*S(i,j,k,3))&
             *(x1-x)*S(i1,j1,k1,1)/(r**5.)

            END IF
      END IF
      !!!----------------------------------------------------------
      IF ((S(i1,j1,k1,2)==1.).or.(S(i1,j1,k1,2)==-1.)) THEN
            x1=real(i+deli) ; y1=real(j+delj)+0.5 ; z1=real(k+delk)
            
            r=sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y)+(z1-z)*(z1-z))
                        
            IF ((0.<r).and.(r<=r0)) THEN                

            EN_D_tmp1=EN_D_tmp1+S(i,j,k,2)*S(i1,j1,k1,2)/(r**3.)&     
            -3.*((x1-x)*S(i,j,k,1)+(y1-y)*S(i,j,k,2)+(z1-z)*S(i,j,k,3))&
             *(y1-y)*S(i1,j1,k1,2)/(r**5.)

            END IF
      END IF
      !!!----------------------------------------------------------
      IF ((S(i1,j1,k1,3)==1.).or.(S(i1,j1,k1,3)==-1.)) THEN
            x1=real(i+deli) ; y1=real(j+delj) ; z1=real(k+delk)+0.5

            r=sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y)+(z1-z)*(z1-z))
                        
            IF ((0.<r).and.(r<=r0)) THEN                

            EN_D_tmp1=EN_D_tmp1+S(i,j,k,3)*S(i1,j1,k1,3)/(r**3.)&     
            -3.*((x1-x)*S(i,j,k,1)+(y1-y)*S(i,j,k,2)+(z1-z)*S(i,j,k,3))&
             *(z1-z)*S(i1,j1,k1,3)/(r**5.)

            END IF
      END IF
         
      END DO
      END DO
      END DO


      END SUBROUTINE value_EN_D_tmp1
           
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_EN_D_tmp()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE value_EN_D_tmp()
      IMPLICIT NONE

     !!!----------------------------------------------------------
      IF ((S(i0,j0,k0,1)==1.).or.(S(i0,j0,k0,1)==-1.)) THEN
            x=real(i0)+0.5 ; y=real(j0) ; z=real(k0)
      END IF

      IF ((S(i0,j0,k0,2)==1.).or.(S(i0,j0,k0,2)==-1.)) THEN
            x=real(i0) ; y=real(j0)+0.5 ; z=real(k0)
      END IF

      IF ((S(i0,j0,k0,3)==1.).or.(S(i0,j0,k0,3)==-1.)) THEN
            x=real(i0) ; y=real(j0) ; z=real(k0)+0.5
      END IF
      !!!========================

      EN_D_tmp=0.

      DO deli=-r01,r01
      DO delj=-r01,r01
      DO delk=-r01,r01

            i1=i0+deli ; j1=j0+delj ; k1=k0+delk  
            
            IF (i1<1)    i1=natx+i1
            IF (i1>natx) i1=i1-natx
            IF (j1<1)    j1=naty+j1                        
            IF (j1>naty) j1=j1-naty                        
            IF (k1<1)    k1=natz+k1                       
            IF (k1>natz) k1=k1-natz    
      !!!------------------------------------------------------------------------------------------
      IF ((S(i1,j1,k1,1)==1.).or.(S(i1,j1,k1,1)==-1.)) THEN
            x1=real(i0+deli)+0.5 ; y1=real(j0+delj) ; z1=real(k0+delk)

            r=sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y)+(z1-z)*(z1-z))
                        
            IF ((0.<r).and.(r<=r0)) THEN                

            EN_D_tmp=EN_D_tmp+S(i0,j0,k0,1)*S(i1,j1,k1,1)/(r**3.)&     
            -3.*((x1-x)*S(i0,j0,k0,1)+(y1-y)*S(i0,j0,k0,2)+(z1-z)*S(i0,j0,k0,3))&
             *(x1-x)*S(i1,j1,k1,1)/(r**5.)

            END IF
      END IF
      !!!----------------------------------------------------------
      IF ((S(i1,j1,k1,2)==1.).or.(S(i1,j1,k1,2)==-1.)) THEN
            x1=real(i0+deli) ; y1=real(j0+delj)+0.5 ; z1=real(k0+delk)

            r=sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y)+(z1-z)*(z1-z))
                        
            IF ((0.<r).and.(r<=r0)) THEN                

            EN_D_tmp=EN_D_tmp+S(i0,j0,k0,2)*S(i1,j1,k1,2)/(r**3.)&     
            -3.*((x1-x)*S(i0,j0,k0,1)+(y1-y)*S(i0,j0,k0,2)+(z1-z)*S(i0,j0,k0,3))&
             *(y1-y)*S(i1,j1,k1,2)/(r**5.)

            END IF
      END IF
      !!!----------------------------------------------------------
      IF ((S(i1,j1,k1,3)==1.).or.(S(i1,j1,k1,3)==-1.)) THEN
            x1=real(i0+deli) ; y1=real(j0+delj) ; z1=real(k0+delk)+0.5

            r=sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y)+(z1-z)*(z1-z))
                        
            IF ((0.<r).and.(r<=r0)) THEN                

            EN_D_tmp=EN_D_tmp+S(i0,j0,k0,3)*S(i1,j1,k1,3)/(r**3.)&     
            -3.*((x1-x)*S(i0,j0,k0,1)+(y1-y)*S(i0,j0,k0,2)+(z1-z)*S(i0,j0,k0,3))&
             *(z1-z)*S(i1,j1,k1,3)/(r**5.)

            END IF
      END IF

      END DO
      END DO
      END DO


      END SUBROUTINE value_EN_D_tmp
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE ground_state()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE ground_state()
      IMPLICIT NONE


      DO i_gs=1,n_average_thermal

      EN_tmp=0. ;  EN_tmp1=0. ;  EN_tmp2=0.     
      EN_tmpx=0.;  EN_tmpy=0. ;  EN_tmpz=0.
      EN_tmp1x=0.; EN_tmp1y=0. ; EN_tmp1z=0.
            
      n_line(:)=0.

      DO k=1,natz
            kp=k+1-(k/natz)*natz
            km=k-1+(1/k)*natz
            kpp=kp+1-(kp/natz)*natz
            DO j=1,naty
                  jp=j+1-(j/naty)*naty
                  jm=j-1+(1/j)*naty
                  jpp=jp+1-(jp/naty)*naty
                  DO i=1,natx
                        ip=i+1-(i/natx)*natx
                        im=i-1+(1/i)*natx
                        ipp=ip+1-(ip/natx)*natx

                       
                        IF (tab_site(i,j,k)==1) THEN
                        
                        !!!-------------------------------------------------------------------------                     
                        
                        IF (int(S(i,j,k,1))==1) THEN
                        compt=0
                        Emin=0.
                        
                        CALL value_EN_D_tmp1()
                        EN_tmp1x=-J1_tmp*(S(i,jp,k,1)+S(i,jm,k,1)+S(i,j,kp,1)+S(i,j,km,1))&
                                 +D_tmp*EN_D_tmp1
                        
                        
                             CALL search_site_void_x()
                              
                             IF (n_etat/=0) THEN
                             
                                    DO i0_etat=1,n_etat
                                   
                                          i_etat=name_etat(i0_etat)                            

                                          CALL value_Hx_tmp2()
                                          CALL value_ENx_D_tmp2()

                                          EN_tmp2= Hx_tmp2(i_etat)+D_tmp*ENx_D_tmp2(i_etat)
                                    
                                          IF (EN_tmp2<Emin) THEN
                                                Emin=EN_tmp2
                                                compt=i_etat
                                          END IF
                                    
                                    END DO
                                    
                                    IF ((Emin<EN_tmp1x).and.(compt>0)) THEN
                                          i_etat=compt    
                                          CALL update_line_x()
                                    END IF

                              END IF
                              
                       END IF   
                       
                       !!!---------------------------------------------------------------
                              
                       IF (int(S(i,j,k,2))==1) THEN
                       compt=0
                       Emin=0.
                       CALL value_EN_D_tmp1()

                       EN_tmp1y=-J1_tmp*(S(ip,j,k,2)+S(im,j,k,2)+S(i,j,kp,2)+S(i,j,km,2))&
                               +D_tmp*EN_D_tmp1
                              
                              CALL search_site_void_y()
                        
                              IF (n_etat/=0) THEN

                                    DO i0_etat=1,n_etat
                                    i_etat=name_etat(i0_etat)

                                    CALL value_Hy_tmp2()
                                    CALL value_ENy_D_tmp2()
                                    EN_tmp2=Hy_tmp2(i_etat)+D_tmp*ENy_D_tmp2(i_etat)

                                    IF (EN_tmp2<Emin) THEN
                                          Emin=EN_tmp2
                                          compt=i_etat
                                    END IF
                                    
                                    END DO
                                   
                                    IF ((Emin<EN_tmp1y).and.(compt>0)) THEN
                                          i_etat=compt 
                                          CALL update_line_y()
                                    END IF

                              END IF
                            
                        END IF

                        !!!---------------------------------------------------------------
                        
                        IF (int(S(i,j,k,3))==1) THEN
                        compt=0
                        Emin=0.
                        CALL value_EN_D_tmp1()
      
                        EN_tmp1z=-J1_tmp*(S(ip,j,k,3)+S(im,j,k,3)+S(i,jp,k,3)+S(i,jm,k,3))&
                                 +D_tmp*EN_D_tmp1

                              CALL search_site_void_z()
                              IF (n_etat/=0) THEN
      
                                    DO i0_etat=1,n_etat

                                    i_etat=name_etat(i0_etat)

                                    CALL value_Hz_tmp2()
                                    CALL value_ENz_D_tmp2()
                                    EN_tmp2=Hz_tmp2(i_etat)+D_tmp*ENz_D_tmp2(i_etat)

                                    IF (EN_tmp2<Emin) THEN
                                          Emin=EN_tmp2
                                          compt=i_etat                                                                      
                                    END IF
                                    
                                    END DO
                                    
                                    IF ((Emin<EN_tmp1z).and.(compt>0)) THEN
                                          i_etat=compt 
                                          CALL update_line_z()
                                    END IF

                              END IF

                          END IF
!!!====================================================================================================
                          END IF
                   ENDDO
             ENDDO
      ENDDO   

      END DO     

      END SUBROUTINE ground_state
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!  SUBROUTINE energy_gs: 
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE energy_gs()

      IMPLICIT NONE

      EN_tmp=0. ;  EN_tmp1=0. ;  EN_tmp2=0.     
      EN_tmpx=0.;  EN_tmpy=0. ;  EN_tmpz=0.
      EN_tmp1x=0.; EN_tmp1y=0. ; EN_tmp1z=0.
     
      
      OPEN(unit=210,file='E-gs.dat')
      DO k=1,natz
            kp=k+1-(k/natz)*natz
            km=k-1+(1/k)*natz
            kpp=kp+1-(kp/natz)*natz
            DO j=1,naty
                  jp=j+1-(j/naty)*naty
                  jm=j-1+(1/j)*naty
                  jpp=jp+1-(jp/naty)*naty
                  DO i=1,natx
                        ip=i+1-(i/natx)*natx
                        im=i-1+(1/i)*natx
                        ipp=ip+1-(ip/natx)*natx
                        
                        IF (tab_site(i,j,k)==1) THEN

                        !!!---------------------------------------------------------------                        
                        IF ((S(i,j,k,1)==1.).or.(S(i,j,k,1)==-1.)) THEN
                                               
                        CALL value_EN_D_tmp1()

                        EN_tmp1x=-J1_tmp*S(i,j,k,1)*(S(i,jp,k,1)+S(i,jm,k,1)+S(i,j,kp,1)+S(i,j,km,1))&
                                 +D_tmp*EN_D_tmp1

                        EN_tmpx=EN_tmpx+EN_tmp1x
                              
                       END IF   
      
    
                       !!!---------------------------------------------------------------
                              
                       IF ((S(i,j,k,2)==1.).or.(S(i,j,k,2)==-1.)) THEN
                       CALL value_EN_D_tmp1()

                       EN_tmp1y=-J1_tmp*S(i,j,k,2)*(S(ip,j,k,2)+S(im,j,k,2)+S(i,j,kp,2)+S(i,j,km,2))&
                               +D_tmp*EN_D_tmp1
                              
                       EN_tmpy=EN_tmpy+EN_tmp1y
                       END IF
                                                                  
                        !!!---------------------------------------------------------------                   
                        
                        IF ((S(i,j,k,3)==1.).or.(S(i,j,k,3)==-1.)) THEN
                        CALL value_EN_D_tmp1()
      
                        EN_tmp1z=-J1_tmp*(S(ip,j,k,3)+S(im,j,k,3)+S(i,jp,k,3)+S(i,jm,k,3))&
                                 +D_tmp*EN_D_tmp1

                        EN_tmpz=EN_tmpz+EN_tmp1z

                          END IF
!!!====================================================================================================
                       END IF

                          
                   ENDDO
            ENDDO
      ENDDO

      EN_tmp=EN_tmpx+EN_tmpy+EN_tmpz

      EN_gs=EN_tmp/(2.*number_line_ini)
      
      WRITE(*,*)'D=',D_tmp,'r=',r0,'Egs=',EN_gs

      WRITE(210,*)iT,r0,EN_gs
      
      CLOSE(210)


      END SUBROUTINE energy_gs


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE average_thermal()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE average_thermal()

      IMPLICIT NONE

      EN_moy=0.
      EN_2_moy=0.
      OP_moy=0.
      OP_2_moy=0.

      OP2_moy=0.
      OP2_2_moy=0.

!!!------------------------------------------
      OPEN(unit=126,file='value_P_at_To.dat')

      H_histo(:)=0.
      EN_histo(:)=0.
      Mz_histo(:)=0.
      EN_2_histo(:)=0.
      Mz_2_histo(:)=0.
      Mz_EN_histo(:)=0.
      Mz_2_EN_histo(:)=0.

      del_EN_histo = (EN_max_histo-EN_min_histo)/real(n_EN_histo-1)

!! END khai bao bien cho histogram
     

     
     
      DO i_loop=1,n_average_thermal

            DO i_loop2=1,n_equi_reseau2
                  CALL equi_reseau()
            END DO
            
            CALL value_thermal()        
            
            EN_moy=EN_moy+energy
            EN_2_moy=EN_2_moy+energy_2
            OP_moy=OP_moy+order_parameter
            OP_2_moy=OP_2_moy+order_parameter_2
            OP2_moy=OP2_moy+order_parameter2
            OP2_2_moy=OP2_2_moy+order_parameter2_2

            IF ((EN_min_histo <= energy) .and. (energy <= EN_max_histo)) THEN              
                  i_histo=int((energy-EN_min_histo)/del_EN_histo)+1
                  H_histo(i_histo)=H_histo(i_histo)+1.
                  
                  EN_histo(i_histo)=EN_histo(i_histo)+energy
                  Mz_histo(i_histo)=Mz_histo(i_histo)+order_parameter
                  EN_2_histo(i_histo)=EN_2_histo(i_histo)+energy_2
                  Mz_2_histo(i_histo)=Mz_2_histo(i_histo)+order_parameter_2

                  Mz_EN_histo(i_histo)=Mz_EN_histo(i_histo)+order_parameter*energy
                  Mz_2_EN_histo(i_histo)=Mz_2_EN_histo(i_histo)+order_parameter_2*energy                  

            END IF      
            
      END DO

      DO i_histo=1,n_EN_histo
            P_histo(i_histo)=H_histo(i_histo)/real(n_average_thermal)
           
            IF (H_histo(i_histo)==0.) THEN
                  EN_histo(i_histo)=EN_min_histo+(real(i_histo)-0.5)*del_EN_histo
                  Mz_histo(i_histo)=0.
                  EN_2_histo(i_histo)=0.
                  Mz_2_histo(i_histo)=0.
                  Mz_EN_histo(i_histo)=0.
                  Mz_2_EN_histo(i_histo)=0.
                  
            ELSE
                  EN_histo(i_histo)=EN_histo(i_histo)/H_histo(i_histo)
                  Mz_histo(i_histo)=Mz_histo(i_histo)/H_histo(i_histo)
                  EN_2_histo(i_histo)=EN_2_histo(i_histo)/H_histo(i_histo)
                  Mz_2_histo(i_histo)=Mz_2_histo(i_histo)/H_histo(i_histo)
                  Mz_EN_histo(i_histo)=Mz_EN_histo(i_histo)/H_histo(i_histo)
                  Mz_2_EN_histo(i_histo)=Mz_2_EN_histo(i_histo)/H_histo(i_histo)
            
            END IF

            WRITE(Ligne126,*)T,EN_histo(i_histo),Mz_histo(i_histo),EN_2_histo(i_histo),&
                  Mz_2_histo(i_histo),Mz_EN_histo(i_histo),Mz_2_EN_histo(i_histo),P_histo(i_histo)
            WRITE(126,'(a)') trim(Ligne126)
      ENDDO     
 
      CLOSE(126)
!!---- END program histogram  ---------------------------------------

      EN_moy=EN_moy/real(n_average_thermal)
      EN_2_moy=EN_2_moy/real(n_average_thermal)
      OP_moy=OP_moy/real(n_average_thermal)
      OP_2_moy=OP_2_moy/real(n_average_thermal)

      OP2_moy=OP2_moy/real(n_average_thermal)
      OP2_2_moy=OP2_2_moy/real(n_average_thermal)

      Cv = number_line_ini*(EN_2_moy-EN_moy**2.)/(T**2.)
      Ksi=number_line_ini*(OP_2_moy-OP_moy**2.)/T
      Ksi2=number_line_ini*(OP2_2_moy-OP2_moy**2.)/T

      WRITE(Ligne20,*) T,EN_moy,OP_moy,Cv,Ksi,OP2_moy,Ksi2

      WRITE(20,'(a)') trim(Ligne20)

      WRITE(*,*)'i_times=',i_times,'EN_moy=',EN_moy,'OP_moy=',OP_moy,'Cv=',Cv,'Ksi=',Ksi

      END SUBROUTINE average_thermal

      
 
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE search_site_void_x()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      SUBROUTINE search_site_void_x()

      IMPLICIT NONE

      n_etat=0
      name_etat(:)=0

      !!!! No 1 - x -----------------
      IF (tab_site(im,j,k)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=1
      END IF

      !!!! No 2  + y -----------------                     
      IF (tab_site(i,jp,k)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=2                       
      END IF

      !!!! No 3 - y -----------------
      IF (tab_site(i,jm,k)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=3
      END IF

                        
      !!!! No 4 + z -----------------
      IF (tab_site(i,j,kp)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=4
      END IF
                      
      !!!! No 5 - z -----------------
      IF (tab_site(i,j,km)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=5
      END IF 

      !!!! No 6 + x -----------------
      IF (tab_site(ipp,j,k)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=6
      END IF
                                                
      !!!! No 7 +y  -----------------
      IF (tab_site(ip,jp,k)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=7
      END IF

      !!!! No 8 - y  -----------------
      IF (tab_site(ip,jm,k)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=8
      END IF

      !!!! No 9 + z -----------------
      IF (tab_site(ip,j,kp)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=9
      END IF
     
      !!!! No 10 - z -----------------
      IF (tab_site(ip,j,km)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=10
      END IF
     
      END SUBROUTINE search_site_void_x

                       
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_Hx_tmp2()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      SUBROUTINE value_Hx_tmp2()
      IMPLICIT NONE

      SELECT CASE (i_etat)

      !!!! No 1 - x -----------------
      CASE (1)
      Hx_tmp2(1)=-J1_tmp*(S(im,jp,k,1)+S(im,jm,k,1)+S(im,j,kp,1)+S(im,j,km,1))

      !!!! No 2  + y -----------------                     
      CASE (2)
      Hx_tmp2(2)=-J1_tmp*(S(ip,j,k,2)+S(im,j,k,2)+S(i,j,kp,2)+S(i,j,km,2))

      !!!! No 3 - y -----------------
      CASE (3)
      Hx_tmp2(3)=-J1_tmp*(S(ip,jm,k,2)+S(im,jm,k,2)+S(i,jm,kp,2)+S(i,jm,km,2))
                        
      !!!! No 4 + z----------------- 
      CASE (4)
      Hx_tmp2(4)=-J1_tmp*(S(ip,j,k,3)+S(im,j,k,3)+S(i,jp,k,3)+S(i,jm,k,3))

                      
      !!!! No 5 - z -----------------
      CASE (5)
      Hx_tmp2(5)=-J1_tmp*(S(ip,j,km,3)+S(im,j,km,3)+S(i,jp,km,3)+S(i,jm,km,3))
 
      !!!! No 6 + x -----------------
      CASE (6)
      Hx_tmp2(6)=-J1_tmp*(S(ip,jp,k,1)+S(ip,jm,k,1)+S(ip,j,kp,1)+S(ip,j,km,1))
                                                
      !!!! No 7 +y  -----------------
      CASE (7)
      Hx_tmp2(7)=-J1_tmp*(S(ipp,j,k,2)+S(i,j,k,2)+S(ip,j,kp,2)+S(ip,j,km,2))


      !!!! No 8 - y  -----------------
      CASE (8)
      Hx_tmp2(8)=-J1_tmp*(S(ipp,jm,k,2)+S(i,jm,k,2)+S(ip,jm,kp,2)+S(ip,jm,km,2))

      !!!! No 9 + z  ----------------- 
      CASE (9)
      Hx_tmp2(9)=-J1_tmp*(S(ipp,j,k,3)+S(i,j,k,3)+S(ip,jp,k,3)+S(ip,jm,k,3))

      !!!! No 10 - z -----------------
      CASE (10)
      Hx_tmp2(10)=-J1_tmp*(S(ipp,j,km,3)+S(i,j,km,3)+S(ip,jp,km,3)+S(ip,jm,km,3))

      END SELECT
     

      END SUBROUTINE value_Hx_tmp2

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE search_site_void_y()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      SUBROUTINE search_site_void_y()

      IMPLICIT NONE

      n_etat=0
      name_etat(:)=0

      !!!! No 1 -y -----------------
      IF (tab_site(i,jm,k)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=1
      END IF

      !!!! No 2 +x ----------------- 
      IF (tab_site(ip,j,k)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=2
      END IF

      !!!! No 3 -x -----------------  
      IF (tab_site(im,j,k)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=3
      END IF

      !!!! No 4  + z -----------------                      
      IF (tab_site(i,j,kp)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=4                       
      END IF
                            

      !!!! No 5 -z ----------------- 
      IF (tab_site(i,j,km)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=5
      END IF
              
      !!!! No 6 + y ----------------- 
      IF (tab_site(i,jpp,k)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=6
      END IF

      !!!! No 7 +x ----------------- 
      IF (tab_site(ip,jp,k)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=7
      END IF

      !!!! No 8 -x -----------------  
      IF (tab_site(im,jp,k)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=8
      END IF

      !!!! No 9  + z -----------------                      
      IF (tab_site(i,jp,kp)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=9                       
      END IF
                            
      !!!! No 10 -z ----------------- 
      IF (tab_site(i,jp,km)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=10
      END IF

      END SUBROUTINE search_site_void_y

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_Hy_tmp2()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      SUBROUTINE value_Hy_tmp2()

      IMPLICIT NONE

      SELECT CASE (i_etat)

      !!!! No 1 -y -----------------
      CASE (1)
      Hy_tmp2(1)=-J1_tmp*(S(ip,jm,k,2)+S(im,jm,k,2)+S(i,jm,kp,2)+S(i,jm,km,2))

      !!!! No 2 +x -----------------
      CASE (2)
      Hy_tmp2(2)=-J1_tmp*(S(i,jp,k,1)+S(i,jm,k,1)+S(i,j,kp,1)+S(i,j,km,1))

      !!!! No 3 -x ----------------- 
      CASE (3)
      Hy_tmp2(3)=-J1_tmp*(S(im,jp,k,1)+S(im,jm,k,1)+S(im,j,kp,1)+S(im,j,km,1))

      !!!! No 4  + z -----------------                    
      CASE (4)
      Hy_tmp2(4)=-J1_tmp*(S(ip,j,k,3)+S(im,j,k,3)+S(i,jp,k,3)+S(i,jm,k,3))            

      !!!! No 5 -z ----------------- 
      CASE (5)
      Hy_tmp2(5)=-J1_tmp*(S(ip,j,km,3)+S(im,j,km,3)+S(i,jp,km,3)+S(i,jm,km,3))
              
      !!!! No 6 + y -----------------
      CASE (6)
      Hy_tmp2(6)=-J1_tmp*(S(ip,jp,k,2)+S(im,jp,k,2)+S(i,jp,kp,2)+S(i,jp,km,2))

      !!!! No 7 +x -----------------
      CASE (7)
      Hy_tmp2(7)=-J1_tmp*(S(i,jpp,k,1)+S(i,j,k,1)+S(i,jp,kp,1)+S(i,jp,km,1))

      !!!! No 8 -x -----------------  
      CASE (8)
      Hy_tmp2(8)=-J1_tmp*(S(im,jpp,k,1)+S(im,j,k,1)+S(im,jp,kp,1)+S(im,jp,km,1))

      !!!! No 9  + z -----------------                      
      CASE (9)
      Hy_tmp2(9)=-J1_tmp*(S(ip,jp,k,3)+S(im,jp,k,3)+S(i,jpp,k,3)+S(i,j,k,3))
                           
      !!!! No 10 -z -----------------
      CASE (10)
      Hy_tmp2(10)=-J1_tmp*(S(ip,jp,km,3)+S(im,jp,km,3)+S(i,jpp,km,3)+S(i,j,km,3))
            
      END SELECT

      END SUBROUTINE value_Hy_tmp2


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE search_site_void_z()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      SUBROUTINE search_site_void_z()

      IMPLICIT NONE

      n_etat=0
      name_etat(:)=0

      !!!! No 1 -z -----------------
      IF (tab_site(i,j,km)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=1
      END IF

      !!!! No 2 +x -----------------
      IF (tab_site(ip,j,k)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=2
      END IF

      !!!! No 3 -x -----------------
      IF (tab_site(im,j,k)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=3
      END IF

      !!!! No 4 +y -----------------
      IF (tab_site(i,jp,k)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=4
      END IF
                            
      !!!! No 5 -y -----------------
      IF (tab_site(i,jm,k)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=5
      END IF   

      !!!! No 6  + z -----------------                    
      IF (tab_site(i,j,kpp)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=6                       
      END IF
              
      !!!! No 7 +x ----------------- 
      IF (tab_site(ip,j,kp)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=7
      END IF

      !!!! No 8 -x ----------------- 
      IF (tab_site(im,j,kp)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=8
      END IF

      !!!! No 9 + y -----------------
      IF (tab_site(i,jp,kp)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=9
      END IF

      !!!! No 10 - y ----------------- 
      IF (tab_site(i,jm,kp)==0) THEN
            n_etat=n_etat+1
            name_etat(n_etat)=10
      END IF

      END SUBROUTINE search_site_void_z

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_Hz_tmp2()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      SUBROUTINE value_Hz_tmp2()

      IMPLICIT NONE

      SELECT CASE (i_etat)

      !!!! No 1 -z ----------------- 
      CASE (1)
      Hz_tmp2(1)=-J1_tmp*(S(ip,j,km,3)+S(im,j,km,3)+S(i,jp,km,3)+S(i,jm,km,3))

      !!!! No 2 +x -----------------
      CASE (2)
      Hz_tmp2(2)=-J1_tmp*(S(i,jp,k,1)+S(i,jm,k,1)+S(i,j,kp,1)+S(i,j,km,1))

      !!!! No 3 -x -----------------  
      CASE (3)
      Hz_tmp2(3)=-J1_tmp*(S(im,jp,k,1)+S(im,jm,k,1)+S(im,j,kp,1)+S(im,j,km,1))

      !!!! No 4 +y -----------------
      CASE (4)
      Hz_tmp2(4)=-J1_tmp*(S(ip,j,k,2)+S(im,j,k,2)+S(i,j,kp,2)+S(i,j,km,2))
                            
      !!!! No 5 -y ----------------- 
      CASE (5)
      Hz_tmp2(5)=-J1_tmp*(S(ip,jm,k,2)+S(im,jm,k,2)+S(i,jm,kp,2)+S(i,jm,km,2))

      !!!! No 6  + z -----------------                      
      CASE (6)
      Hz_tmp2(6)=-J1_tmp*(S(ip,j,kp,3)+S(im,j,kp,3)+S(i,jp,kp,3)+S(i,jm,kp,3))
              
      !!!! No 7 +x -----------------
      CASE (7)
      Hz_tmp2(7)=-J1_tmp*(S(i,jp,kp,1)+S(i,jm,kp,1)+S(i,j,kpp,1)+S(i,j,k,1))

      !!!! No 8 -x -----------------
      CASE (8)
      Hz_tmp2(8)=-J1_tmp*(S(im,jp,kp,1)+S(im,jm,kp,1)+S(im,j,kpp,1)+S(im,j,k,1))

      !!!! No 9 + y ----------------- 
      CASE (9)
      Hz_tmp2(9)=-J1_tmp*(S(ip,j,kp,2)+S(im,j,kp,2)+S(i,j,kpp,2)+S(i,j,k,2))

      !!!! No 10 - y -----------------
      CASE (10)
      Hz_tmp2(10)=-J1_tmp*(S(ip,jm,kp,2)+S(im,jm,kp,2)+S(i,jm,kpp,2)+S(i,jm,k,2))

      END SELECT


      END SUBROUTINE value_Hz_tmp2


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE update_line_x()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE update_line_x()
      IMPLICIT NONE                                                   


      SELECT CASE(i_etat)
    
      !!!! No 1  - x -----------------
      CASE (1)
            S(i,j,k,1)=0.
            tab_site(ip,j,k)=0
            S(im,j,k,1)=spin2
            tab_site(im,j,k)=1  

      !!!! No 2  + y ----------------- 
      CASE(2)
            S(i,j,k,1)=0.
            tab_site(ip,j,k)=0 
            S(i,j,k,2)=spin2
            tab_site(i,jp,k)=1


      !!!! No 3  - y -----------------
      CASE (3)
            S(i,j,k,1)=0.
            tab_site(ip,j,k)=0
            S(i,jm,k,2)=spin2
            tab_site(i,jm,k)=1
         
      !!!! No 4  + z -----------------
      CASE (4)
            S(i,j,k,1)=0.
            tab_site(ip,j,k)=0
            S(i,j,k,3)=spin2
            tab_site(i,j,kp)=1

  
       !!!! No 5  - z -----------------
      CASE (5)
            S(i,j,k,1)=0.
            tab_site(ip,j,k)=0
            S(i,j,km,3)=spin2
            tab_site(i,j,km)=1

       !!!! No 6 + x -----------------
       CASE (6)
            S(i,j,k,1)=0.
            tab_site(i,j,k)=0
            S(ip,j,k,1)=spin2
            tab_site(ipp,j,k)=1

      !!!! No 7  + y -----------------
      CASE (7)
            S(i,j,k,1)=0.
            tab_site(i,j,k)=0
            S(ip,j,k,2)=spin2
            tab_site(ip,jp,k)=1

      !!!! No 8  - y -----------------
      CASE (8)
            S(i,j,k,1)=0.
            tab_site(i,j,k)=0
            S(ip,jm,k,2)=spin2
            tab_site(ip,jm,k)=1

                              
      !!!! No 9  + z -----------------
      CASE (9)
            S(i,j,k,1)=0.
            tab_site(i,j,k)=0
            S(ip,j,k,3)=spin2
            tab_site(ip,j,kp)=1

       !!!! No 10  - z ----------------- 
       CASE (10)
            S(i,j,k,1)=0.
            tab_site(i,j,k)=0
            S(ip,j,km,3)=spin2
            tab_site(ip,j,km)=1

       END SELECT

       END SUBROUTINE update_line_x

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE update_line_y()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE update_line_y()
      IMPLICIT NONE                                                   

      SELECT CASE (i_etat)
      !!!! No 1 -y -----------------
      CASE (1)
            S(i,j,k,2)=0.
            tab_site(i,jp,k)=0
            S(i,jm,k,2)=spin2
            tab_site(i,jm,k)=1
                              
      !!!! No 2 +x -----------------
      CASE (2)
            S(i,j,k,2)=0.
            tab_site(i,jp,k)=0 
            S(i,j,k,1)=spin2
            tab_site(ip,j,k)=1

      !!!! No 3 -x -----------------
      CASE (3)
            S(i,j,k,2)=0.
            tab_site(i,jp,k)=0 
            S(im,j,k,1)=spin2
            tab_site(im,j,k)=1

      !!!! No 4 +z -----------------
      CASE (4)
            S(i,j,k,2)=0.
            tab_site(i,jp,k)=0 
            S(i,j,k,3)=spin2
            tab_site(i,j,kp)=1

      !!!! No 5 -z -----------------
      CASE (5)
            S(i,j,k,2)=0.
            tab_site(i,jp,k)=0 
            S(i,j,km,3)=spin2
            tab_site(i,j,km)=1

     !!!! No 6  +y -----------------
      CASE (6)
            S(i,j,k,2)=0.
            tab_site(i,j,k)=0
            S(i,jp,k,2)=spin2
            tab_site(i,jpp,k)=1                     

      !!!! No 7 +x -----------------
      CASE (7)
            S(i,j,k,2)=0.
            tab_site(i,j,k)=0
            S(i,jp,k,1)=spin2
            tab_site(ip,jp,k)=1

      !!!! No 8  -x -----------------
      CASE (8)
            S(i,j,k,2)=0.
            tab_site(i,j,k)=0
            S(im,jp,k,1)=spin2
            tab_site(im,jp,k)=1

       !!!! No 9 + z -----------------
       CASE (9)
            S(i,j,k,2)=0.
            tab_site(i,j,k)=0
            S(i,jp,k,3)=spin2
            tab_site(i,jp,kp)=1

       !!!! No 10 -z -----------------
       CASE (10)
            S(i,j,k,2)=0.
            tab_site(i,j,k)=0
            S(i,jp,km,3)=spin2
            tab_site(i,jp,km)=1

       END SELECT

       END SUBROUTINE update_line_y

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE update_line_z()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE update_line_z()
      IMPLICIT NONE    

      SELECT CASE (i_etat)                                               

      !!!! No 1 -z -----------------
      CASE(1)
            S(i,j,k,3)=0.
            tab_site(i,j,kp)=0 
            S(i,j,km,3)=spin2
            tab_site(i,j,km)=1
                              
      !!!! No 2 +x -----------------
      CASE(2)
            S(i,j,k,3)=0.
            tab_site(i,j,kp)=0 
            S(i,j,k,1)=spin2
            tab_site(ip,j,k)=1


      !!!! No 3 -x -----------------
      CASE (3)
            S(i,j,k,3)=0.
            tab_site(i,j,kp)=0 
            S(im,j,k,1)=spin2
            tab_site(im,j,k)=1

      !!!! No 4 + y -----------------
      CASE (4)
            S(i,j,k,3)=0.
            tab_site(i,j,kp)=0 
            S(i,j,k,2)=spin2
            tab_site(i,jp,k)=1

      !!!! No 5 - y -----------------
      CASE (5)
            S(i,j,k,3)=0.
            tab_site(i,j,kp)=0 
            S(i,jm,k,2)=spin2
            tab_site(i,jm,k)=1

      !!!! No 6 +z -----------------
      CASE (6)
            S(i,j,k,3)=0.
            tab_site(i,j,k)=0 
            S(i,j,kp,3)=spin2
            tab_site(i,j,kpp)=1

      !!!! No 7 +x -----------------
      CASE (7)
            S(i,j,k,3)=0.
            tab_site(i,j,k)=0 
            S(i,j,kp,1)=spin2
            tab_site(ip,j,kp)=1

      !!!! No 8 -x -----------------
      CASE (8)
            S(i,j,k,3)=0.
            tab_site(i,j,k)=0 
            S(im,j,kp,1)=spin2
            tab_site(im,j,kp)=1

      !!!! No 9 + y -----------------
      CASE (9)
            S(i,j,k,3)=0.
            tab_site(i,j,k)=0 
            S(i,j,kp,2)=spin2
            tab_site(i,jp,kp)=1

      !!!! No 10 - y ----------------- 
      CASE (10)
            S(i,j,k,3)=0.
            tab_site(i,j,k)=0 
            S(i,jm,kp,2)=spin2
            tab_site(i,jm,kp)=1


      END SELECT

      END SUBROUTINE update_line_z


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_ENx_D_tmp2()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE value_ENx_D_tmp2()
      IMPLICIT NONE

      SELECT CASE (i_etat)

      !!!! No 1 - x -----------------
      CASE (1)
            i0=im; j0=j; k0=k
            S(i,j,k,1)=0.
            S(i0,j0,k0,1)=spin2
            CALL value_EN_D_tmp()
            ENx_D_tmp2(1)=EN_D_tmp
            S(i,j,k,1)=spin1
            S(i0,j0,k0,1)=0.

      !!!! No 2  + y -----------------                     
      CASE (2)
            i0=i; j0=j; k0=k
            S(i,j,k,1)=0.
            S(i0,j0,k0,2)=spin2
            CALL value_EN_D_tmp()
            ENx_D_tmp2(2)=EN_D_tmp
            S(i,j,k,1)=spin1
            S(i0,j0,k0,2)=0.

      !!!! No 3 - y -----------------
      CASE (3)
            i0=i; j0=jm; k0=k
            S(i,j,k,1)=0.
            S(i0,j0,k0,2)=spin2
            CALL value_EN_D_tmp()
            ENx_D_tmp2(3)=EN_D_tmp
            S(i,j,k,1)=spin1
            S(i0,j0,k0,2)=0.
                        
      !!!! No 4 + z----------------- 
      CASE (4)
            i0=i; j0=j; k0=k
            S(i,j,k,1)=0.
            S(i0,j0,k0,3)=spin2
            CALL value_EN_D_tmp()
            ENx_D_tmp2(4)=EN_D_tmp
            S(i,j,k,1)=spin1
            S(i0,j0,k0,3)=0.
                      
      !!!! No 5 - z -----------------
      CASE (5)
            i0=i; j0=j; k0=km
            S(i,j,k,1)=0.
            S(i0,j0,k0,3)=spin2
            CALL value_EN_D_tmp()
            ENx_D_tmp2(5)=EN_D_tmp
            S(i,j,k,1)=spin1
            S(i0,j0,k0,3)=0.
 
      !!!! No 6 + x -----------------
      CASE (6)
            i0=ip; j0=j; k0=k
            S(i,j,k,1)=0.
            S(i0,j0,k0,1)=spin2
            CALL value_EN_D_tmp()
            ENx_D_tmp2(6)=EN_D_tmp
            S(i,j,k,1)=spin1
            S(i0,j0,k0,1)=0.
                                                
      !!!! No 7 +y  -----------------
      CASE (7)
            i0=ip; j0=j; k0=k
            S(i,j,k,1)=0.
            S(i0,j0,k0,2)=spin2
            CALL value_EN_D_tmp()
            ENx_D_tmp2(7)=EN_D_tmp
            S(i,j,k,1)=spin1
            S(i0,j0,k0,2)=0.

      !!!! No 8 - y  -----------------
      CASE (8)
            i0=ip; j0=jm; k0=k
            S(i,j,k,1)=0.
            S(i0,j0,k0,2)=spin2
            CALL value_EN_D_tmp()
            ENx_D_tmp2(8)=EN_D_tmp
            S(i,j,k,1)=spin1
            S(i0,j0,k0,2)=0.

      !!!! No 9 + z  ----------------- 
      CASE (9)
            i0=ip; j0=j; k0=k
            S(i,j,k,1)=0.
            S(i0,j0,k0,3)=spin2
            CALL value_EN_D_tmp()
            ENx_D_tmp2(9)=EN_D_tmp
            S(i,j,k,1)=spin1
            S(i0,j0,k0,3)=0.

      !!!! No 10 - z -----------------
      CASE (10)
            i0=ip; j0=j; k0=km
            S(i,j,k,1)=0.
            S(i0,j0,k0,3)=spin2
            CALL value_EN_D_tmp()
            ENx_D_tmp2(10)=EN_D_tmp
            S(i,j,k,1)=spin1
            S(i0,j0,k0,3)=0.

      END SELECT
     

      END SUBROUTINE value_ENx_D_tmp2

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_ENy_D_tmp2()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      SUBROUTINE value_ENy_D_tmp2()

      IMPLICIT NONE

      SELECT CASE (i_etat)

      !!!! No 1 -y -----------------
      CASE (1)
            i0=i; j0=jm; k0=k
            S(i,j,k,2)=0.
            S(i0,j0,k0,2)=spin2
            CALL value_EN_D_tmp()
            ENy_D_tmp2(1)=EN_D_tmp
            S(i,j,k,2)=spin1
            S(i0,j0,k0,2)=0.

      !!!! No 2 +x -----------------
      CASE (2)
            i0=i; j0=j; k0=k
            S(i,j,k,2)=0.
            S(i0,j0,k0,1)=spin2
            CALL value_EN_D_tmp()
            ENy_D_tmp2(2)=EN_D_tmp
            S(i,j,k,2)=spin1
            S(i0,j0,k0,1)=0.

      !!!! No 3 -x ----------------- 
      CASE (3)
            i0=im; j0=j; k0=k
            S(i,j,k,2)=0.
            S(i0,j0,k0,1)=spin2
            CALL value_EN_D_tmp()
            ENy_D_tmp2(3)=EN_D_tmp
            S(i,j,k,2)=spin1
            S(i0,j0,k0,1)=0.

      !!!! No 4  + z -----------------                    
      CASE (4)
            i0=i; j0=j; k0=k
            S(i,j,k,2)=0.
            S(i0,j0,k0,3)=spin2
            CALL value_EN_D_tmp()
            ENy_D_tmp2(4)=EN_D_tmp
            S(i,j,k,2)=spin1
            S(i0,j0,k0,3)=0.

      !!!! No 5 -z ----------------- 
      CASE (5)
            i0=i; j0=j; k0=km
            S(i,j,k,2)=0.
            S(i0,j0,k0,3)=spin2
            CALL value_EN_D_tmp()
            ENy_D_tmp2(5)=EN_D_tmp
            S(i,j,k,2)=spin1
            S(i0,j0,k0,3)=0.
              
      !!!! No 6 + y -----------------
      CASE (6)
            i0=i; j0=jp; k0=k
            S(i,j,k,2)=0.
            S(i0,j0,k0,2)=spin2
            CALL value_EN_D_tmp()
            ENy_D_tmp2(6)=EN_D_tmp
            S(i,j,k,2)=spin1
            S(i0,j0,k0,2)=0.

      !!!! No 7 +x -----------------
      CASE (7)
            i0=i; j0=jp; k0=k
            S(i,j,k,2)=0.
            S(i0,j0,k0,1)=spin2
            CALL value_EN_D_tmp()
            ENy_D_tmp2(7)=EN_D_tmp
            S(i,j,k,2)=spin1
            S(i0,j0,k0,1)=0.

      !!!! No 8 -x -----------------  
      CASE (8)
            i0=im; j0=jp; k0=k
            S(i,j,k,2)=0.
            S(i0,j0,k0,1)=spin2
            CALL value_EN_D_tmp()
            ENy_D_tmp2(8)=EN_D_tmp
            S(i,j,k,2)=spin1
            S(i0,j0,k0,1)=0.

      !!!! No 9  + z -----------------                      
      CASE (9)
            i0=i; j0=jp; k0=k
            S(i,j,k,2)=0.
            S(i0,j0,k0,3)=spin2
            CALL value_EN_D_tmp()
            ENy_D_tmp2(9)=EN_D_tmp
            S(i,j,k,2)=spin1
            S(i0,j0,k0,3)=0.
          
      !!!! No 10 -z -----------------
      CASE (10)
            i0=i; j0=jp; k0=km
            S(i,j,k,2)=0.
            S(i0,j0,k0,3)=spin2
            CALL value_EN_D_tmp()
            ENy_D_tmp2(10)=EN_D_tmp
            S(i,j,k,2)=spin1
            S(i0,j0,k0,3)=0.
            
      END SELECT

      END SUBROUTINE value_ENy_D_tmp2

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_ENz_D_tmp2()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE value_ENz_D_tmp2()

      IMPLICIT NONE

      SELECT CASE (i_etat)

      !!!! No 1 -z ----------------- 
      CASE (1)
            i0=i; j0=j; k0=km
            S(i,j,k,3)=0.
            S(i0,j0,k0,3)=spin2
            CALL value_EN_D_tmp()
            ENz_D_tmp2(1)=EN_D_tmp
            S(i,j,k,3)=spin1
            S(i0,j0,k0,3)=0.

      !!!! No 2 +x -----------------
      CASE (2)
            i0=i; j0=j; k0=k
            S(i,j,k,3)=0.
            S(i0,j0,k0,1)=spin2
            CALL value_EN_D_tmp()
            ENz_D_tmp2(2)=EN_D_tmp
            S(i,j,k,3)=spin1
            S(i0,j0,k0,1)=0.

      !!!! No 3 -x -----------------  
      CASE (3)
            i0=im; j0=j; k0=k
            S(i,j,k,3)=0.
            S(i0,j0,k0,1)=spin2
            CALL value_EN_D_tmp()
            ENz_D_tmp2(3)=EN_D_tmp
            S(i,j,k,3)=spin1
            S(i0,j0,k0,1)=0.

      !!!! No 4 +y -----------------
      CASE (4)
            i0=i; j0=j; k0=k
            S(i,j,k,3)=0.
            S(i0,j0,k0,2)=spin2
            CALL value_EN_D_tmp()
            ENz_D_tmp2(4)=EN_D_tmp
            S(i,j,k,3)=spin1
            S(i0,j0,k0,2)=0.
                            
      !!!! No 5 -y ----------------- 
      CASE (5)
            i0=i; j0=jm; k0=k
            S(i,j,k,3)=0.
            S(i0,j0,k0,2)=spin2
            CALL value_EN_D_tmp()
            ENz_D_tmp2(5)=EN_D_tmp
            S(i,j,k,3)=spin1
            S(i0,j0,k0,2)=0.

      !!!! No 6  + z -----------------                      
      CASE (6)
            i0=i; j0=j; k0=kp
            S(i,j,k,3)=0.
            S(i0,j0,k0,3)=spin2
            CALL value_EN_D_tmp()
            ENz_D_tmp2(6)=EN_D_tmp
            S(i,j,k,3)=spin1
            S(i0,j0,k0,3)=0.
              
      !!!! No 7 +x -----------------
      CASE (7)
            i0=i; j0=j; k0=kp
            S(i,j,k,3)=0.
            S(i0,j0,k0,1)=spin2
            CALL value_EN_D_tmp()
            ENz_D_tmp2(7)=EN_D_tmp
            S(i,j,k,3)=spin1
            S(i0,j0,k0,1)=0.

      !!!! No 8 -x -----------------
      CASE (8)
            i0=im; j0=j; k0=kp
            S(i,j,k,3)=0.
            S(i0,j0,k0,1)=spin2
            CALL value_EN_D_tmp()
            ENz_D_tmp2(8)=EN_D_tmp
            S(i,j,k,3)=spin1
            S(i0,j0,k0,1)=0.

      !!!! No 9 + y ----------------- 
      CASE (9)
            i0=i; j0=j; k0=kp
            S(i,j,k,3)=0.
            S(i0,j0,k0,2)=spin2
            CALL value_EN_D_tmp()
            ENz_D_tmp2(9)=EN_D_tmp
            S(i,j,k,3)=spin1
            S(i0,j0,k0,2)=0.

      !!!! No 10 - y -----------------
      CASE (10)
            i0=i; j0=jm; k0=kp
            S(i,j,k,3)=0.
            S(i0,j0,k0,2)=spin2
            CALL value_EN_D_tmp()
            ENz_D_tmp2(10)=EN_D_tmp
            S(i,j,k,3)=spin1
            S(i0,j0,k0,2)=0.

      END SELECT


      END SUBROUTINE value_ENz_D_tmp2

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_Mz4() 13.01.2012
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE value_Mz4()
      IMPLICIT NONE


!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      M11=0. ; M12=0. ; M13=0. ; M14=0. 
      
      n1_line_x1=0. ; n1_line_x2=0. ; n1_line_y1=0. ; n1_line_y2=0. ; n1_line_z1=0.; n1_line_z2=0.
      n2_line_x1=0. ; n2_line_x2=0. ; n2_line_y1=0. ; n2_line_y2=0. ; n2_line_z1=0.; n2_line_z2=0.
      n3_line_x1=0. ; n3_line_x2=0. ; n3_line_y1=0. ; n3_line_y2=0. ; n3_line_z1=0.; n3_line_z2=0.
      n4_line_x1=0. ; n4_line_x2=0. ; n4_line_y1=0. ; n4_line_y2=0. ; n4_line_z1=0.; n4_line_z2=0.

      DO j=1,naty
            DO i=1,natx
                  SELECT CASE(mod(i+j,4))
                        CASE(1)
                        DO k=1,natz
                        IF (S(i,j,k,1)==1.)     n1_line_x1=n1_line_x1+1.
                        IF (S(i,j,k,1)==-1.)    n1_line_x2=n1_line_x2+1.
                        IF (S(i,j,k,2)==1.)     n1_line_y1=n1_line_y1+1.
                        IF (S(i,j,k,2)==-1.)    n1_line_y2=n1_line_y2+1.
                        IF (S(i,j,k,3)==1.)     n1_line_z1=n1_line_z1+1.
                        IF (S(i,j,k,3)==-1.)    n1_line_z2=n1_line_z2+1.       
                        END DO

                        CASE(2)
                        DO k=1,natz
                        IF (S(i,j,k,1)==1.)     n2_line_x1=n2_line_x1+1.
                        IF (S(i,j,k,1)==-1.)    n2_line_x2=n2_line_x2+1.
                        IF (S(i,j,k,2)==1.)     n2_line_y1=n2_line_y1+1.
                        IF (S(i,j,k,2)==-1.)    n2_line_y2=n2_line_y2+1.
                        IF (S(i,j,k,3)==1.)     n2_line_z1=n2_line_z1+1.
                        IF (S(i,j,k,3)==-1.)    n2_line_z2=n2_line_z2+1.  
                        END DO           

                        CASE(3)
                        DO k=1,natz
                        IF (S(i,j,k,1)==1.)     n3_line_x1=n3_line_x1+1.
                        IF (S(i,j,k,1)==-1.)    n3_line_x2=n3_line_x2+1.
                        IF (S(i,j,k,2)==1.)     n3_line_y1=n3_line_y1+1.
                        IF (S(i,j,k,2)==-1.)    n3_line_y2=n3_line_y2+1.
                        IF (S(i,j,k,3)==1.)     n3_line_z1=n3_line_z1+1.
                        IF (S(i,j,k,3)==-1.)    n3_line_z2=n3_line_z2+1.  
                        END DO           

                        CASE(0)
                        DO k=1,natz
                        IF (S(i,j,k,1)==1.)     n4_line_x1=n4_line_x1+1.
                        IF (S(i,j,k,1)==-1.)    n4_line_x2=n4_line_x2+1.
                        IF (S(i,j,k,2)==1.)     n4_line_y1=n4_line_y1+1.
                        IF (S(i,j,k,2)==-1.)    n4_line_y2=n4_line_y2+1.
                        IF (S(i,j,k,3)==1.)     n4_line_z1=n4_line_z1+1.
                        IF (S(i,j,k,3)==-1.)    n4_line_z2=n4_line_z2+1.  
                        END DO          

                  END SELECT   
            END DO
      END DO


      IF ((n1_line_x1+n1_line_x2+n1_line_y1+n1_line_y2+n1_line_z1+n1_line_z2)/=0.) THEN
      M11=abs(6.*max(n1_line_x1,n1_line_x2,n1_line_y1,n1_line_y2,n1_line_z1,n1_line_z2)&
            /(n1_line_x1+n1_line_x2+n1_line_y1+n1_line_y2+n1_line_z1+n1_line_z2)-1.)/5. 
      END IF

      IF ((n2_line_x1+n2_line_x2+n2_line_y1+n2_line_y2+n2_line_z1+n2_line_z2)/=0.) THEN
      M12=abs(6.*max(n2_line_x1,n2_line_x2,n2_line_y1,n2_line_y2,n2_line_z1,n2_line_z2)&
            /(n2_line_x1+n2_line_x2+n2_line_y1+n2_line_y2+n2_line_z1+n2_line_z2)-1.)/5. 
      END IF

      IF ((n3_line_x1+n3_line_x2+n3_line_y1+n3_line_y2+n3_line_z1+n3_line_z2)/=0.) THEN
      M13=abs(6.*max(n3_line_x1,n3_line_x2,n3_line_y1,n3_line_y2,n3_line_z1,n3_line_z2)&
            /(n3_line_x1+n3_line_x2+n3_line_y1+n3_line_y2+n3_line_z1+n3_line_z2)-1.)/5. 
      END IF

      IF ((n4_line_x1+n4_line_x2+n4_line_y1+n4_line_y2+n4_line_z1+n4_line_z2)/=0.) THEN
      M14=abs(6.*max(n4_line_x1,n4_line_x2,n4_line_y1,n4_line_y2,n4_line_z1,n4_line_z2)&
            /(n4_line_x1+n4_line_x2+n4_line_y1+n4_line_y2+n4_line_z1+n4_line_z2)-1.)/5. 
      END IF
      
      order_parameter=max((M11+M13),(M12+M14))/2.
      
      END SUBROUTINE value_Mz4

!!!!================================================================================================
!!!! Tinh M cua GS4, khong phan biet Spin+1 hay -1
!!!! Nhu tinh GS2 cua system non polarite

      SUBROUTINE value_Mz42()
      IMPLICIT NONE

      M1=0. ; M11x=0.; M12x=0.; M13x=0.; M14x=0; M11y=0.; M12y=0.; M13y=0.; M14y=0.
      M11z=0.; M12z=0.; M13z=0.; M14z=0.
      
      M11=0. ; M12=0. ; M13=0. ; M14=0. 

      DO j=1,naty
            DO i=1,natx
                  SELECT CASE(mod(i+j,4))
                        CASE(1)
                        DO k=1,natz
                        M11x=M11x+abs(S(i,j,k,1))
                        M11y=M11y+abs(S(i,j,k,2))
                        M11z=M11z+abs(S(i,j,k,3))
                        END DO

                        CASE(2)
                        DO k=1,natz
                        M12x=M12x+abs(S(i,j,k,1))
                        M12y=M12y+abs(S(i,j,k,2))
                        M12z=M12z+abs(S(i,j,k,3))
                        END DO           

                        CASE(3)
                        DO k=1,natz
                        M13x=M13x+abs(S(i,j,k,1))
                        M13y=M13y+abs(S(i,j,k,2))
                        M13z=M13z+abs(S(i,j,k,3))
                        END DO           

                        CASE(0)
                        DO k=1,natz
                        M14x=M14x+abs(S(i,j,k,1))
                        M14y=M14y+abs(S(i,j,k,2))
                        M14z=M14z+abs(S(i,j,k,3))
                        END DO          

                  END SELECT   
            END DO
      END DO
      
      IF ((M11x+M11y+M11z)/=0)    M11=abs(3.*max(M11x,M11y,M11z)/(M11x+M11y+M11z)-1.)/2.      
      IF ((M12x+M12y+M12z)/=0)    M12=abs(3.*max(M12x,M12y,M12z)/(M12x+M12y+M12z)-1.)/2.
      IF ((M13x+M13y+M13z)/=0)    M13=abs(3.*max(M13x,M13y,M13z)/(M13x+M13y+M13z)-1.)/2.
      IF ((M14x+M14y+M14z)/=0)    M14=abs(3.*max(M14x,M14y,M14z)/(M14x+M14y+M14z)-1.)/2.
      
      order_parameter=max((M11+M13),(M12+M14))/2.
      
      !WRITE(*,*)M11,M12,M13,M14
      
      END SUBROUTINE value_Mz42

!!!!================================================================================================


      END PROGRAM main_transport_dimer
      

      
