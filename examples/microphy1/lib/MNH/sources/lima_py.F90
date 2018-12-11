      SUBROUTINE LIMA_PY   (K1, K2, K3, K4, K5,                                  &
!
                            LWARM, LACTIT, LSEDC, LRAIN,                  &
                            LCOLD, LHHONI, LSEDI,                                &
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
                            PINPRS, PINPRG, PINPRH )
!
USE MODI_LIMA_WARM
USE MODI_LIMA_COLD
USE MODI_LIMA_MIXED
USE MODI_LIMA_ADJUST
!
IMPLICIT NONE
INTEGER,                     INTENT(IN)    :: K1      ! Nx
INTEGER,                     INTENT(IN)    :: K2      ! Ny
INTEGER,                     INTENT(IN)    :: K3      ! Nz
INTEGER,                     INTENT(IN)    :: K4      ! size(PRT)
INTEGER,                     INTENT(IN)    :: K5      ! size(PSVT)
!
LOGICAL,                     INTENT(IN)    :: LWARM
LOGICAL,                     INTENT(IN)    :: LACTIT
LOGICAL,                     INTENT(IN)    :: LSEDC
LOGICAL,                     INTENT(IN)    :: LRAIN
LOGICAL,                     INTENT(IN)    :: LCOLD
LOGICAL,                     INTENT(IN)    :: LHHONI
LOGICAL,                     INTENT(IN)    :: LSEDI
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
!
REAL, DIMENSION(K1,K2,K3) :: ZWORK1, ZWORK2, ZWORK3
!
IF (LWARM) CALL LIMA_WARM(OACTIT=LACTIT, OSEDC=LSEDC, ORAIN=LRAIN, KSPLITR=KSPLITR, PTSTEP=PTSTEP, KMI=1, &
                          HFMFILE='DUMMY', HLUOUT='DUMMY', OCLOSE_OUT=.FALSE., KRR=KRR, PZZ=PDZZ, PRHODJ=PRHODJ,     &
                          PRHODREF=PRHODREF, PEXNREF=PEXNREF, PW_NU=PW_NU, PPABSM=PPABST, PPABST=PPABST,         &
                          PTHM=PTHT, PRCM=PRT(:,:,:,2),                  &
                          PTHT=PTHT, PRT=PRT, PSVT=PSVT,                                   &
                          PTHS=PTHS, PRS=PRS, PSVS=PSVS,                                   &
                          PINPRC=PINPRC, PINPRR=PINPRR, PINPRR3D=PINPRR3D, PEVAP3D=PEVAP3D         )
!
IF (LCOLD) CALL LIMA_COLD(OSEDI=LSEDI, OHHONI=LHHONI, KSPLITG=KSPLITG, PTSTEP=PTSTEP, KMI=1,               &
                          HFMFILE='DUMMY', HLUOUT='DUMMY', OCLOSE_OUT=.FALSE., KRR=KRR, PZZ=PDZZ, PRHODJ=PRHODJ,     &
                          PRHODREF=PRHODREF, PEXNREF=PEXNREF, PPABST=PPABST, PW_NU=PW_NU,                 &
                          PTHM=PTHT, PPABSM=PPABST,                                      &
                          PTHT=PTHT, PRT=PRT, PSVT=PSVT,                                   &
                          PTHS=PTHS, PRS=PRS, PSVS=PSVS,                                   &
                          PINPRS=PINPRS, PINPRG=PINPRG, PINPRH=PINPRH)
!
IF (LWARM .AND. LCOLD) CALL LIMA_MIXED(OSEDI=LSEDI, OHHONI=LHHONI, KSPLITG=KSPLITG, PTSTEP=PTSTEP, KMI=1,     &
                                 HFMFILE='DUMMY', HLUOUT='DUMMY', OCLOSE_OUT=.FALSE., KRR=KRR, PZZ=PDZZ, PRHODJ=PRHODJ, &
                                 PRHODREF=PRHODREF, PEXNREF=PEXNREF, PPABST=PPABST, PW_NU=PW_NU,             &
                                 PTHM=PTHT, PPABSM=PPABST,                                  &
                                 PTHT=PTHT, PRT=PRT, PSVT=PSVT,                               &
                                 PTHS=PTHS, PRS=PRS, PSVS=PSVS                                )
!
CALL LIMA_ADJUST(KRR=KRR, KMI=1, HFMFILE='DUMMY', HLUOUT='DUMMY', HRAD='DUMMY',                  &
                 HTURBDIM='DUMMY', OCLOSE_OUT=.FALSE., OSUBG_COND=.FALSE., PTSTEP=PTSTEP,         &
                 PRHODREF=PRHODREF, PRHODJ=PRHODJ, PEXNREF=PEXNREF, PPABSM=PPABST, PSIGS=ZWORK1, PPABST=PPABST, &
                 PRT=PRT, PRS=PRS, PSVT=PSVT, PSVS=PSVS,                             &
                 PTHS=PTHS, PSRCS=ZWORK2, PCLDFR=ZWORK3                               )
!add ZINPRC in PINPRR
PINPRR=PINPRR+PINPRC
!
END SUBROUTINE LIMA_PY
