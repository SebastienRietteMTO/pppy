SUBROUTINE RAIN_ICE_RED_PY_4(K1, K2, K3, OSEDIC, HSEDIM, HSUBG_AUCV_RC, OWARM,KKA,KKU,KKL,        &
                             PTSTEP, KRR,              &
                             PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
                             PTHT, PRVT, PRCT, PRRT, PRIT, PRST, &
                             PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS, &
                             PINPRC, PINPRR, PINPRR3D, PEVAP3D,                    &
                             PINPRS, PINPRG, PSIGS, PINDEP, PSEA, PTOWN,                   &
                             PRHT, PRHS, PINPRH                        )
USE MODI_RAIN_ICE_RED
IMPLICIT NONE
INTEGER, INTENT(IN) :: K1,K2,K3
LOGICAL,                  INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
CHARACTER(LEN=4),         INTENT(IN)    :: HSEDIM ! Sedimentation scheme
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV_RC ! Kind of Subgrid autoconversion method
LOGICAL,                  INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to
                                                   !   form by warm processes
                                                   !      (Kessler scheme)
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
!
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PCLDFR  ! Convective Mass Flux Cloud fraction
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRC! Cloud instant precip
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRR! Rain instant precip
REAL, DIMENSION(K1,K2,K3), INTENT(OUT)   :: PINPRR3D! Rain inst precip 3D
REAL, DIMENSION(K1,K2,K3), INTENT(OUT)     :: PEVAP3D! Rain evap profile
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRS! Snow instant precip
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRG! Graupel instant precip
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
REAL, DIMENSION(K1,K2), INTENT(OUT)     :: PINDEP  ! Cloud instant deposition
!
REAL, DIMENSION(K1,K2), INTENT(IN)        :: PSEA ! Sea Mask
REAL, DIMENSION(K1,K2), INTENT(IN)        :: PTOWN! Fraction that is town
REAL, DIMENSION(K1,K2,K3), INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(K1,K2,K3), INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRH! Hail instant precip

LOGICAL, DIMENSION(K1, K2, K3) :: LLMICRO
LLMICRO(:,:,:)=.TRUE.

CALL RAIN_ICE_RED(OSEDIC, HSEDIM, HSUBG_AUCV_RC, OWARM,KKA,KKU,KKL,        &
                  PTSTEP, KRR, LLMICRO, PEXNREF,             &
                  PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
                  PTHT, PRVT, PRCT, PRRT, PRIT, PRST, &
                  PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS, &
                  PINPRC, PINPRR, PINPRR3D, PEVAP3D,                    &
                  PINPRS, PINPRG, PSIGS, PINDEP, PSEA, PTOWN,                   &
                  PRHT, PRHS, PINPRH)
END SUBROUTINE RAIN_ICE_RED_PY_4

SUBROUTINE RAIN_ICE_RED_PY_3(K1, K2, K3, OSEDIC, HSEDIM, HSUBG_AUCV_RC, OWARM,KKA,KKU,KKL,        &
                             PTSTEP, KRR,              &
                             PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
                             PTHT, PRVT, PRCT, PRRT, PRIT, PRST, &
                             PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS, &
                             PINPRC, PINPRR, PINPRR3D, PEVAP3D,                    &
                             PINPRS, PINPRG, PSIGS, PINDEP, PSEA, PTOWN)
USE MODI_RAIN_ICE_RED
IMPLICIT NONE
INTEGER, INTENT(IN) :: K1,K2,K3
LOGICAL,                  INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
CHARACTER(LEN=4),         INTENT(IN)    :: HSEDIM ! Sedimentation scheme
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV_RC ! Kind of Subgrid autoconversion method
LOGICAL,                  INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to
                                                   !   form by warm processes
                                                   !      (Kessler scheme)
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
!
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PCLDFR  ! Convective Mass Flux Cloud fraction
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRC! Cloud instant precip
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRR! Rain instant precip
REAL, DIMENSION(K1,K2,K3), INTENT(OUT)   :: PINPRR3D! Rain inst precip 3D
REAL, DIMENSION(K1,K2,K3), INTENT(OUT)     :: PEVAP3D! Rain evap profile
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRS! Snow instant precip
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRG! Graupel instant precip
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
REAL, DIMENSION(K1,K2), INTENT(OUT)     :: PINDEP  ! Cloud instant deposition
!
REAL, DIMENSION(K1,K2), INTENT(IN)        :: PSEA ! Sea Mask
REAL, DIMENSION(K1,K2), INTENT(IN)        :: PTOWN! Fraction that is town

LOGICAL, DIMENSION(K1, K2, K3) :: LLMICRO
LLMICRO(:,:,:)=.TRUE.

CALL RAIN_ICE_RED(OSEDIC, HSEDIM, HSUBG_AUCV_RC, OWARM,KKA,KKU,KKL,        &
                  PTSTEP, KRR, LLMICRO, PEXNREF,             &
                  PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
                  PTHT, PRVT, PRCT, PRRT, PRIT, PRST, &
                  PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS, &
                  PINPRC, PINPRR, PINPRR3D, PEVAP3D,                    &
                  PINPRS, PINPRG, PSIGS, PINDEP, PSEA, PTOWN)
END SUBROUTINE RAIN_ICE_RED_PY_3

SUBROUTINE RAIN_ICE_PY_4(K1, K2, K3, OSEDIC, HSEDIM, HSUBG_AUCV_RC, OWARM,KKA,KKU,KKL,        &
                         KSPLITR, PTSTEP, KRR,              &
                         PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
                         PTHT, PRVT, PRCT, PRRT, PRIT, PRST, &
                         PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS, &
                         PINPRC, PINPRR, PINPRR3D, PEVAP3D,                    &
                         PINPRS, PINPRG, PSIGS, PINDEP, PSEA, PTOWN,                   &
                         PRHT, PRHS, PINPRH                        )
USE MODI_RAIN_ICE
IMPLICIT NONE
INTEGER, INTENT(IN) :: K1,K2,K3
LOGICAL,                  INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
CHARACTER(LEN=4),         INTENT(IN)    :: HSEDIM ! Sedimentation scheme
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV_RC ! Kind of Subgrid autoconversion method
LOGICAL,                  INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to
                                                   !   form by warm processes
                                                   !      (Kessler scheme)
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)    :: KSPLITR ! Number of small time step integration for  rain sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
!
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PCLDFR  ! Convective Mass Flux Cloud fraction
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRC! Cloud instant precip
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRR! Rain instant precip
REAL, DIMENSION(K1,K2,K3), INTENT(OUT)   :: PINPRR3D! Rain inst precip 3D
REAL, DIMENSION(K1,K2,K3), INTENT(OUT)     :: PEVAP3D! Rain evap profile
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRS! Snow instant precip
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRG! Graupel instant precip
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
REAL, DIMENSION(K1,K2), INTENT(OUT)     :: PINDEP  ! Cloud instant deposition
!
REAL, DIMENSION(K1,K2), INTENT(IN)        :: PSEA ! Sea Mask
REAL, DIMENSION(K1,K2), INTENT(IN)        :: PTOWN! Fraction that is town
REAL, DIMENSION(K1,K2,K3), INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(K1,K2,K3), INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRH! Hail instant precip

CALL RAIN_ICE(OSEDIC, HSEDIM, HSUBG_AUCV_RC, OWARM,KKA,KKU,KKL,        &
              KSPLITR, PTSTEP, KRR,              &
              PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
              PTHT, PRVT, PRCT, PRRT, PRIT, PRST, &
              PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS, &
              PINPRC, PINPRR, PINPRR3D, PEVAP3D,                    &
              PINPRS, PINPRG, PSIGS, PINDEP, PSEA, PTOWN,                   &
              PRHT, PRHS, PINPRH)
END SUBROUTINE RAIN_ICE_PY_4

SUBROUTINE RAIN_ICE_PY_3(K1, K2, K3, OSEDIC, HSEDIM, HSUBG_AUCV_RC, OWARM,KKA,KKU,KKL,        &
                         KSPLITR, PTSTEP, KRR,              &
                         PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
                         PTHT, PRVT, PRCT, PRRT, PRIT, PRST, &
                         PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS, &
                         PINPRC, PINPRR, PINPRR3D, PEVAP3D,                    &
                         PINPRS, PINPRG, PSIGS, PINDEP, PSEA, PTOWN)
USE MODI_RAIN_ICE
IMPLICIT NONE
INTEGER, INTENT(IN) :: K1,K2,K3
LOGICAL,                  INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
CHARACTER(LEN=4),         INTENT(IN)    :: HSEDIM ! Sedimentation scheme
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV_RC ! Kind of Subgrid autoconversion method
LOGICAL,                  INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to
                                                   !   form by warm processes
                                                   !      (Kessler scheme)
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)    :: KSPLITR ! Number of small time step integration for  rain sedimendation
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
!
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PCLDFR  ! Convective Mass Flux Cloud fraction
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRC! Cloud instant precip
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRR! Rain instant precip
REAL, DIMENSION(K1,K2,K3), INTENT(OUT)   :: PINPRR3D! Rain inst precip 3D
REAL, DIMENSION(K1,K2,K3), INTENT(OUT)     :: PEVAP3D! Rain evap profile
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRS! Snow instant precip
REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRG! Graupel instant precip
REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
REAL, DIMENSION(K1,K2), INTENT(OUT)     :: PINDEP  ! Cloud instant deposition
!
REAL, DIMENSION(K1,K2), INTENT(IN)        :: PSEA ! Sea Mask
REAL, DIMENSION(K1,K2), INTENT(IN)        :: PTOWN! Fraction that is town

CALL RAIN_ICE(OSEDIC, HSEDIM, HSUBG_AUCV_RC, OWARM,KKA,KKU,KKL,        &
              KSPLITR, PTSTEP, KRR,              &
              PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
              PTHT, PRVT, PRCT, PRRT, PRIT, PRST, &
              PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS, &
              PINPRC, PINPRR, PINPRR3D, PEVAP3D,                    &
              PINPRS, PINPRG, PSIGS, PINDEP, PSEA, PTOWN)
END SUBROUTINE RAIN_ICE_PY_3

SUBROUTINE ICE_ADJUST_PY(K1, K2, K3, PTSTEP, PSIGQSAT, &
                         PRHODJ, PEXNREF, PSIGS, PPABST, PTH, &
                         PRV, PRC, PRI, PRR, PRS, PRG, PRH)
USE MODI_ICE_ADJUST
USE MODD_PARAM_ICE, ONLY : CFRAC_ICE_ADJUST
IMPLICIT NONE
INTEGER, INTENT(IN) :: K1,K2,K3
REAL, INTENT(IN) :: PTSTEP, PSIGQSAT
REAL, DIMENSION(K1, K2, K3), INTENT(INOUT) :: PRV, PRC, PRI, PTH
REAL, DIMENSION(K1, K2, K3), INTENT(IN) :: PRR, PRS, PRG, PRHODJ, &
                                           PEXNREF, PSIGS, PPABST
REAL, DIMENSION(K1, K2, K3), INTENT(IN) :: PRH
!
REAL, DIMENSION(K1, K2, K3) :: ZMFCONV, ZZZ, ZCF_MF,ZRC_MF,ZRI_MF, &
                             & ZCLDFR, ZRVS, ZRCS, ZRIS, &
                             & ZRH, ZSRCS, ZTHS, ZRVOUT, ZRCOUT, &
                             & ZRIOUT, ZTHOUT
!
ZMFCONV = 0.
ZZZ = 0.
ZCF_MF = 0.
ZRC_MF = 0.
ZRI_MF = 0.
ZRVS=PRV/PTSTEP
ZRCS=PRC/PTSTEP
ZRIS=PRI/PTSTEP
ZTHS=PTH/PTSTEP
!
CALL ICE_ADJUST (1, K3, 1, 7, CFRAC_ICE_ADJUST,                    &
                 'XXXX', .TRUE., .TRUE., PTSTEP, PSIGQSAT,  &
                 PRHODJ, PEXNREF, PSIGS, ZMFCONV, PPABST, ZZZ, PEXNREF, &
                 ZCF_MF, ZRC_MF, ZRI_MF,                             &
                 PRV, PRC, ZRVS, ZRCS, PTH, ZTHS, ZSRCS, ZCLDFR,     &
                 PRR, PRI, ZRIS, PRS, PRG,    &
                             PRH, ZRVOUT, ZRCOUT, ZRIOUT, ZTHOUT)
PRV=ZRVOUT
PRC=ZRCOUT
PRI=ZRIOUT
PTH=ZTHOUT
END SUBROUTINE ICE_ADJUST_PY
