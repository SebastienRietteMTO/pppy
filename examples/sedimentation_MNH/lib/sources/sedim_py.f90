SUBROUTINE SEDIM_PY(HSEDIM, LHAIL, IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                   &PTSTEP, KRR, KSPLITR, OMOMENTUM, &
                   &PSEA, PTOWN, PDZZ, &
                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                   &PRCT, PRRT, PRIT, PRST, PRGT, PRHT, &
                   &PRCS, PRRS, PRIS, PRSS, PRGS, PRHS, &
                   &PINPRC, PINPRR, PINPRI, PINPRS, PINPRG, PINPRH)
!
USE MODI_ICE4_SEDIMENTATION_STAT
USE MODI_ICE4_SEDIMENTATION_SPLIT
USE MODI_ICE4_SEDIMENTATION_SPLIT_MOMENTUM
USE MODI_ICE4_SEDIMENTATION_SPLIT_OLD
!
IMPLICIT NONE
!
CHARACTER(LEN=4), INTENT(IN) :: HSEDIM
LOGICAL(KIND=8), INTENT(IN) :: LHAIL
INTEGER, INTENT(IN) :: IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, KRR, KSPLITR
LOGICAL(KIND=8), INTENT(IN) :: OMOMENTUM
REAL, INTENT(IN) :: PTSTEP
REAL, DIMENSION(IIT, IJT), INTENT(IN) :: PSEA, PTOWN
REAL, DIMENSION(IIT, IJT, IKT), INTENT(IN) :: PDZZ, PRHODREF, PPABST, PTHT, PRHODJ, &
                                            & PRCT, PRRT, PRIT, PRST, PRGT, PRHT
REAL, DIMENSION(IIT, IJT, IKT), INTENT(INOUT) :: PRCS, PRRS, PRIS, PRSS, PRGS, PRHS
REAL, DIMENSION(IIT, IJT), INTENT(OUT) :: PINPRC, PINPRR, PINPRI, PINPRS, PINPRG, PINPRH
!
LOGICAL :: LLMOMENTUM
REAL, DIMENSION(IIT, IJT) :: ZINDEP
!
IF(HSEDIM=='STAT') THEN
   IF(LHAIL) THEN
      CALL ICE4_SEDIMENTATION_STAT(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                  &PTSTEP, KRR, .TRUE., .FALSE., 0., &
                                  &PDZZ, &
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT,&
                                  &PRSS, PRST, PRGS, PRGT,&
                                  &PINPRC, ZINDEP, PINPRR, PINPRI, PINPRS, PINPRG, &
                                  &PSEA=PSEA, PTOWN=PTOWN, PINPRH=PINPRH, PRHT=PRHT, PRHS=PRHS)
   ELSE
      CALL ICE4_SEDIMENTATION_STAT(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                  &PTSTEP, KRR, .TRUE., .FALSE., 0., &
                                  &PDZZ, &
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT,&
                                  &PRSS, PRST, PRGS, PRGT,&
                                  &PINPRC, ZINDEP, PINPRR, PINPRI, PINPRS, PINPRG, &
                                  &PSEA=PSEA, PTOWN=PTOWN)
   ENDIF
ELSEIF(HSEDIM=='SPLN') THEN
   IF(LHAIL) THEN
      CALL ICE4_SEDIMENTATION_SPLIT(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                   &PTSTEP, KRR, .TRUE., .FALSE., 0., &
                                   &PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, ZINDEP, PINPRR, PINPRI, PINPRS, PINPRG, &
                                   &PSEA=PSEA, PTOWN=PTOWN, PINPRH=PINPRH, PRHT=PRHT, PRHS=PRHS)
   ELSE
      CALL ICE4_SEDIMENTATION_SPLIT(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                   &PTSTEP, KRR, .TRUE., .FALSE., 0., &
                                   &PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, ZINDEP, PINPRR, PINPRI, PINPRS, PINPRG, &
                                   &PSEA=PSEA, PTOWN=PTOWN)
   ENDIF
ELSEIF(HSEDIM=='SPLI') THEN
   IF(LHAIL) THEN
      CALL ICE4_SEDIMENTATION_SPLIT_OLD(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                   &PTSTEP, KRR, .TRUE., KSPLITR, &
                                   &PSEA, PTOWN, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINPRR, PINPRI, PINPRS, PINPRG, &
                                   &PINPRH=PINPRH, PRHT=PRHT, PRHS=PRHS)
   ELSE
      CALL ICE4_SEDIMENTATION_SPLIT_OLD(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                   &PTSTEP, KRR, .TRUE., KSPLITR, &
                                   &PSEA, PTOWN, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINPRR, PINPRI, PINPRS, PINPRG)
   ENDIF
ELSEIF(HSEDIM=='SPL2') THEN
   LLMOMENTUM=OMOMENTUM
   IF(LHAIL) THEN
      CALL ICE4_SEDIMENTATION_SPLIT_MOMENTUM(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                   &PTSTEP, KRR, .TRUE., LLMOMENTUM, &
                                   &PSEA, PTOWN, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINPRR, PINPRI, PINPRS, PINPRG, &
                                   &PINPRH=PINPRH, PRHT=PRHT, PRHS=PRHS)
   ELSE
      CALL ICE4_SEDIMENTATION_SPLIT_MOMENTUM(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                   &PTSTEP, KRR, .TRUE., LLMOMENTUM, &
                                   &PSEA, PTOWN, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINPRR, PINPRI, PINPRS, PINPRG)
   ENDIF
ELSE
      PRINT*, "Wrong HSEDIM value"
ENDIF
!
END SUBROUTINE SEDIM_PY
