      SUBROUTINE LIMA_PY   (K1, K2, K3, K4, K5,                                  &
!
                            LWARM, LACTIT_1BYTE, LSEDC_1BYTE, LRAIN_1BYTE,                  &
                            LCOLD, LHHONI_1BYTE, LSEDI_1BYTE,                                &
!                            
                            PTSTEP, KRR,                                         &
!                            
                            KSPLITR, KSPLITG,                                    &
!                            
                            PDZZ, PRHODJ, PRHODREF, PEXNREF,                     &
!                            
                            PPABST, PW_NU,                                       &
!                            
                            PRT, PSVT, PTHT,                                     &
!                            
                            PRS, PSVS, PTHS,                                     &
!                            
                            PINPRC,                                     &
                            PINPRR, PINPRR3D, PEVAP3D,                           &
                            PINPRS, PINPRG, PINPRH,PINDEP )
!
USE MODI_LIMA_WARM
USE MODI_LIMA_COLD
USE MODI_LIMA_MIXED
USE MODI_LIMA_ADJUST
!
USE MODD_IO_ll,   ONLY: TFILEDATA
!
IMPLICIT NONE
INTEGER,                     INTENT(IN)    :: K1      ! Nx
INTEGER,                     INTENT(IN)    :: K2      ! Ny
INTEGER,                     INTENT(IN)    :: K3      ! Nz
INTEGER,                     INTENT(IN)    :: K4      ! size(PRT)
INTEGER,                     INTENT(IN)    :: K5      ! size(PSVT)
!
LOGICAL(KIND=1),             INTENT(IN)    :: LWARM
LOGICAL(KIND=1),             INTENT(IN)    :: LACTIT_1BYTE
LOGICAL(KIND=1),             INTENT(IN)    :: LSEDC_1BYTE
LOGICAL(KIND=1),             INTENT(IN)    :: LRAIN_1BYTE
LOGICAL(KIND=1),             INTENT(IN)    :: LCOLD
LOGICAL(KIND=1),             INTENT(IN)    :: LHHONI_1BYTE
LOGICAL(KIND=1),             INTENT(IN)    :: LSEDI_1BYTE
!
REAL,                        INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                     INTENT(IN)    :: KRR     ! Number of moist variable
!
INTEGER,                     INTENT(IN)    :: KSPLITR ! Number of small time step integration for  rain sedimendation
INTEGER,                     INTENT(IN)    :: KSPLITG ! Number of small time step integration for  hail sedimendation
!
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PEXNREF ! Reference Exner function
!
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PPABST  ! absolute pressure at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PW_NU
!
REAL, DIMENSION(K1,K2,K3,K4),INTENT(IN)    :: PRT
REAL, DIMENSION(K1,K2,K3,K5),INTENT(IN)    :: PSVT
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PTHT
!
REAL, DIMENSION(K1,K2,K3,K4),INTENT(INOUT) :: PRS
REAL, DIMENSION(K1,K2,K3,K5),INTENT(INOUT) :: PSVS
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PTHS
!
REAL, DIMENSION(K1,K2),      INTENT(OUT)   :: PINPRC! Cloud instant precip
REAL, DIMENSION(K1,K2),      INTENT(OUT)   :: PINPRR! Rain instant precip
REAL, DIMENSION(K1,K2,K3),   INTENT(OUT)   :: PINPRR3D! Rain inst precip 3D
REAL, DIMENSION(K1,K2,K3),   INTENT(OUT)   :: PEVAP3D! Rain evap profile
REAL, DIMENSION(K1,K2),      INTENT(OUT)   :: PINPRS! Snow instant precip
REAL, DIMENSION(K1,K2),      INTENT(OUT)   :: PINPRG! Graupel instant precip
REAL, DIMENSION(K1,K2),      INTENT(OUT)   :: PINPRH! Graupel instant precip
REAL, DIMENSION(K1,K2),      INTENT(OUT)   :: PINDEP!Cloud droplets deposition
!
REAL, DIMENSION(K1,K2,K3) :: ZWORK1, ZWORK2, ZWORK3
LOGICAL :: LSEDI, LACTIT, LSEDC, LRAIN, LHHONI
TYPE(TFILEDATA) :: TPFILE
!
LSEDI=LSEDI_1BYTE
LACTIT=LACTIT_1BYTE
LSEDC=LSEDC_1BYTE
LRAIN=LRAIN_1BYTE
LHHONI=LHHONI_1BYTE
!
IF (LWARM) CALL LIMA_WARM(OACTIT=LACTIT, OSEDC=LSEDC, ORAIN=LRAIN, KSPLITR=KSPLITR, PTSTEP=PTSTEP, KMI=1, &
                          TPFILE=TPFILE, OCLOSE_OUT=.FALSE., KRR=KRR, PZZ=PDZZ, PRHODJ=PRHODJ,     &
                          PRHODREF=PRHODREF, PEXNREF=PEXNREF, PW_NU=PW_NU, PPABSM=PPABST, PPABST=PPABST,         &
                          PTHM=PTHT, PRCM=PRT(:,:,:,2),                  &
                          PTHT=PTHT, PRT=PRT, PSVT=PSVT,                                   &
                          PTHS=PTHS, PRS=PRS, PSVS=PSVS,                                   &
                          PINPRC=PINPRC, PINPRR=PINPRR, PINDEP=PINDEP, PINPRR3D=PINPRR3D, PEVAP3D=PEVAP3D         )
!
IF (LCOLD) CALL LIMA_COLD(OSEDI=LSEDI, OHHONI=LHHONI, KSPLITG=KSPLITG, PTSTEP=PTSTEP, KMI=1,               &
                          KRR=KRR, PZZ=PDZZ, PRHODJ=PRHODJ,     &
                          PRHODREF=PRHODREF, PEXNREF=PEXNREF, PPABST=PPABST, PW_NU=PW_NU,                 &
                          PTHM=PTHT, PPABSM=PPABST,                                      &
                          PTHT=PTHT, PRT=PRT, PSVT=PSVT,                                   &
                          PTHS=PTHS, PRS=PRS, PSVS=PSVS,                                   &
                          PINPRS=PINPRS, PINPRG=PINPRG, PINPRH=PINPRH)
!
IF (LWARM .AND. LCOLD) CALL LIMA_MIXED(OSEDI=LSEDI, OHHONI=LHHONI, KSPLITG=KSPLITG, PTSTEP=PTSTEP, KMI=1,     &
                                 KRR=KRR, PZZ=PDZZ, PRHODJ=PRHODJ, &
                                 PRHODREF=PRHODREF, PEXNREF=PEXNREF, PPABST=PPABST, PW_NU=PW_NU,             &
                                 PTHM=PTHT, PPABSM=PPABST,                                  &
                                 PTHT=PTHT, PRT=PRT, PSVT=PSVT,                               &
                                 PTHS=PTHS, PRS=PRS, PSVS=PSVS                                )
!
CALL LIMA_ADJUST(KRR=KRR, KMI=1, TPFILE=TPFILE, HRAD='DUMMY',                  &
                 HTURBDIM='DUMMY', OCLOSE_OUT=.FALSE., OSUBG_COND=.FALSE., PTSTEP=PTSTEP,         &
                 PRHODREF=PRHODREF, PRHODJ=PRHODJ, PEXNREF=PEXNREF, PPABSM=PPABST, PSIGS=ZWORK1, PPABST=PPABST, &
                 PRT=PRT, PRS=PRS, PSVT=PSVT, PSVS=PSVS,                             &
                 PTHS=PTHS, PSRCS=ZWORK2, PCLDFR=ZWORK3                               )
!add ZINPRC in PINPRR
PINPRR=PINPRR+PINPRC
!
END SUBROUTINE LIMA_PY
