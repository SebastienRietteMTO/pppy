!     ######spl
SUBROUTINE INIT_LIMA_PY(PTSTEP, LDWARM, CMICRO, KSPLITR, KSPLITG, CCSEDIM,    &
                             LDCRIAUTI, PCRIAUTI, PT0CRIAUTI, PCRIAUTC,                    &
!
                             LNUCL_IN, LCOLD_IN, LSNOW_IN, LHAIL_IN, LHHONI_IN, LMEYERS_IN,&
                             NMOD_IFN_IN, CIFN1_IN, CIFN2_IN, IFN_HOM_IN,                  &
                             IFN_SPECIES_IN, INT_MIXING_IN, NMOD_IMM_IN, IND_SPECIE_IN,    &
                             CPRISTINE_ICE_LIMA_IN, CHEVRIMED_ICE_LIMA_IN,                 &
                             XALPHAI_IN, XNUI_IN, XALPHAS_IN, XNUS_IN, XALPHAG_IN, XNUG_IN,&
                             XFACTNUC_DEP_IN, XFACTNUC_CON_IN, NPHILLIPS_IN,               &
!
                             LWARM_IN, LACTI_IN, LRAIN_IN, LSEDC_IN, LACTIT_IN, LBOUND_IN, &
                             NMOD_CCN_IN, CCCN1_IN, CCCN2_IN, CCCN3_IN, CCCN4_IN,          &
                             CCN_HOM_IN, CCN_MODES_IN, HINI_CCN_IN, HTYPE_CCN_IN,          &
                             XALPHAC_IN, XNUC_IN, XALPHAR_IN, XNUR_IN,                     &
                             XFSOLUB_CCN_IN, XACTEMP_CCN_IN, XAERDIFF_IN, XAERHEIGHT_IN,   &
                             LSCAV_IN, LAERO_MASS_IN                                       )

!**** *INI_MICRO*   - Initialize common meso_NH MODD_ used in microphysics for AROME

!     Purpose.
!     --------
!           Initialize 
!           MODD_RAIN_ICE_DESCR, MODD_RAIN_ICE_PARAM and MODD_PARAM_ICE  
!           parameters used in AROME microphysics 

!**   Interface.
!     ----------
!        *CALL* *INI_MICRO (KULOUT,KSTEP,KSPLITR)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output
!        PTSTEP  : Time step
!        KSPLITR : Number of small time step interation for rain sedimentation 
!        LDWARM : value assigned to LWARM       

!        Implicit arguments :
!        --------------------
!        

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Documentation AROME 

!     Author.
!     -------
!        Y. Seity 

!     Modifications.
!     --------------
!        Original : 03-12-12
!        05-08-25 Kovacic  added LDWARM
!     ------------------------------------------------------------------

USE YOMLUN   , ONLY : NULNAM

USE MODD_NSV
USE MODD_LIMA_PRECIP_SCAVENGING_n
USE MODD_PARAM_LIMA_COLD
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_MIXED
USE MODD_PARAM_LIMA_WARM
 
USE MODI_INI_LIMA
USE MODI_INIT_AEROSOL_PROPERTIES

IMPLICIT NONE
! -----------------------------------------------------------------------
!     DUMMY INTEGER SCALARS
REAL, INTENT (IN) :: PTSTEP
LOGICAL, INTENT (IN) :: LDWARM
CHARACTER(4), INTENT (IN) :: CMICRO 
CHARACTER(4), INTENT (IN) :: CCSEDIM
INTEGER, INTENT (OUT) :: KSPLITR
INTEGER, INTENT (OUT) :: KSPLITG
LOGICAL, INTENT (IN) :: LDCRIAUTI
REAL, INTENT (IN) :: PCRIAUTI
REAL, INTENT (IN) :: PT0CRIAUTI
REAL, INTENT (IN) :: PCRIAUTC
!
LOGICAL, INTENT (IN)       :: LNUCL_IN
LOGICAL, INTENT (IN)       :: LCOLD_IN
LOGICAL, INTENT (IN)       :: LSNOW_IN
LOGICAL, INTENT (IN)       :: LHAIL_IN
LOGICAL, INTENT (IN)       :: LHHONI_IN
LOGICAL, INTENT (IN)       :: LMEYERS_IN
INTEGER, INTENT (IN)       :: NMOD_IFN_IN
REAL, INTENT (IN)          :: CIFN1_IN
REAL, INTENT (IN)          :: CIFN2_IN
LOGICAL, INTENT (IN)       :: IFN_HOM_IN
CHARACTER(6), INTENT (IN) :: IFN_SPECIES_IN
CHARACTER(4), INTENT (IN) :: INT_MIXING_IN
INTEGER, INTENT (IN)       :: NMOD_IMM_IN
INTEGER, INTENT (IN)       :: IND_SPECIE_IN
CHARACTER(4), INTENT (IN) :: CPRISTINE_ICE_LIMA_IN
CHARACTER(4), INTENT (IN) :: CHEVRIMED_ICE_LIMA_IN
REAL, INTENT (IN)          :: XALPHAI_IN
REAL, INTENT (IN)          :: XNUI_IN
REAL, INTENT (IN)          :: XALPHAS_IN
REAL, INTENT (IN)          :: XNUS_IN
REAL, INTENT (IN)          :: XALPHAG_IN
REAL, INTENT (IN)          :: XNUG_IN
REAL, INTENT (IN)          :: XFACTNUC_DEP_IN
REAL, INTENT (IN)          :: XFACTNUC_CON_IN
INTEGER, INTENT (IN)       :: NPHILLIPS_IN
!
LOGICAL, INTENT (IN)       :: LWARM_IN
LOGICAL, INTENT (IN)       :: LACTI_IN
LOGICAL, INTENT (IN)       :: LRAIN_IN
LOGICAL, INTENT (IN)       :: LSEDC_IN
LOGICAL, INTENT (IN)       :: LACTIT_IN
LOGICAL, INTENT (IN)       :: LBOUND_IN
INTEGER, INTENT (IN)       :: NMOD_CCN_IN
REAL, INTENT (IN)          :: CCCN1_IN
REAL, INTENT (IN)          :: CCCN2_IN
REAL, INTENT (IN)          :: CCCN3_IN
REAL, INTENT (IN)          :: CCCN4_IN
LOGICAL, INTENT (IN)       :: CCN_HOM_IN
CHARACTER(8), INTENT (IN) :: CCN_MODES_IN
CHARACTER(3), INTENT (IN) :: HINI_CCN_IN
CHARACTER(1), INTENT (IN) :: HTYPE_CCN_IN
REAL, INTENT (IN)          :: XALPHAC_IN
REAL, INTENT (IN)          :: XNUC_IN
REAL, INTENT (IN)          :: XALPHAR_IN
REAL, INTENT (IN)          :: XNUR_IN
REAL, INTENT (IN)          :: XFSOLUB_CCN_IN
REAL, INTENT (IN)          :: XACTEMP_CCN_IN
REAL, INTENT (IN)          :: XAERDIFF_IN
REAL, INTENT (IN)          :: XAERHEIGHT_IN
LOGICAL, INTENT (IN)       :: LSCAV_IN
LOGICAL, INTENT (IN)       :: LAERO_MASS_IN
                             
!-----------------------------------------------------------------------
!    LOCAL VARIABLES
REAL :: ZCRI0, ZTCRI0   
INTEGER :: ISV 
! -----------------------------------------------------------------------
!
!
! -----------------------------------------------------------------------
CALL INI_CST
! -----------------------------------------------------------------------
! lecture Valeurs par défaut pour les paramètres de la namelist LIMA
LNUCL              = LNUCL_IN
LCOLD              = LCOLD_IN
LSNOW              = LSNOW_IN
LHAIL              = LHAIL_IN
LHHONI             = LHHONI_IN
LMEYERS            = LMEYERS_IN
NMOD_IFN           = NMOD_IFN_IN
XIFN_CONC(1)       = CIFN1_IN
XIFN_CONC(2)       = CIFN2_IN
LIFN_HOM           = IFN_HOM_IN
CIFN_SPECIES       = TRIM(ADJUSTL(IFN_SPECIES_IN))
CINT_MIXING        = TRIM(ADJUSTL(INT_MIXING_IN))
NMOD_IMM           = NMOD_IMM_IN
NIND_SPECIE        = IND_SPECIE_IN
CPRISTINE_ICE_LIMA = TRIM(ADJUSTL(CPRISTINE_ICE_LIMA_IN))
CHEVRIMED_ICE_LIMA = TRIM(ADJUSTL(CHEVRIMED_ICE_LIMA_IN))
XALPHAI            = XALPHAI_IN
XNUI               = XNUI_IN
XALPHAS            = XALPHAS_IN
XNUS               = XNUS_IN
XALPHAG            = XALPHAG_IN
XNUG               = XNUG_IN
XFACTNUC_DEP       = XFACTNUC_DEP_IN
XFACTNUC_CON       = XFACTNUC_CON_IN
NPHILLIPS          = NPHILLIPS_IN
!
LWARM              = LWARM_IN
LACTI              = LACTI_IN
LRAIN              = LRAIN_IN
LSEDC              = LSEDC_IN
LACTIT             = LACTIT_IN
LBOUND             = LBOUND_IN
NMOD_CCN           = NMOD_CCN_IN
XCCN_CONC(1)       = CCCN1_IN
XCCN_CONC(2)       = CCCN2_IN
XCCN_CONC(3)       = CCCN3_IN
XCCN_CONC(4)       = CCCN4_IN
LCCN_HOM           = CCN_HOM_IN
CCCN_MODES         = TRIM(ADJUSTL(CCN_MODES_IN))
HINI_CCN           = TRIM(ADJUSTL(HINI_CCN_IN))
HTYPE_CCN          = TRIM(ADJUSTL(HTYPE_CCN_IN))
XALPHAC            = XALPHAC_IN
XNUC               = XNUC_IN
XALPHAR            = XALPHAR_IN
XNUR               = XNUR_IN
XFSOLUB_CCN        = XFSOLUB_CCN_IN
XACTEMP_CCN        = XACTEMP_CCN_IN
XAERDIFF           = XAERDIFF_IN
XAERHEIGHT         = XAERHEIGHT_IN
LSCAV              = LSCAV_IN
LAERO_MASS         = LAERO_MASS_IN
!
!!$LNUCL              = .TRUE.
LSEDI              = .FALSE.
!!$LSNOW              = .TRUE.
!!$LHAIL              = .FALSE.
!!$LHHONI             = .FALSE.
!!$LMEYERS            = .FALSE.
!!$NMOD_IFN           = 1
!!$CIFN1              = 1000
!!$CIFN2              = 0.001
!!$IFN_HOM            = .TRUE.
!!$IFN_SPECIES        = 'PHILLIPS'
!!$INT_MIXING         = ''
!!$NMOD_IMM           = 0
!!$IND_SPECIE         = 1
!!$CPRISTINE_ICE_LIMA = 'PLAT'
!!$CHEVRIMED_ICE_LIMA = 'GRAU'
!!$XALPHAI            = 0.
!!$XNUI               = 0.
!!$XALPHAS            = 0.
!!$XNUS               = 0.
!!$XALPHAG            = 0.
!!$XNUG               = 0.
!!$XFACTNUC_DEP       = 1.
!!$XFACTNUC_CON       = 1.
!!$NPHILLIPS          = 8
!!$!
!!$LWARM              = .TRUE.
!!$LACTI              = .TRUE.
!!$LRAIN              = .TRUE.
!!$LSEDC              = .FALSE.
!!$LACTIT             = .FALSE.
!!$LBOUND             = .FALSE.
!!$NMOD_CCN           = 1
!!$CCCN1              = 350.
!!$CCCN2              = 100.
!!$CCCN3              = 35.
!!$CCCN4              = 1.
!!$CCN_HOM            = .TRUE.
!!$CCN_MODES          = ''
!!$HINI_CCN           = 'XXX'
!!$HTYPE_CCN          = 'X'
!!$XALPHAC            = 3.
!!$XNUC               = 1.
!!$XALPHAR            = ALPHAR
!!$XNUR               = NUR
!!$XFSOLUB_CCN        = 1.
!!$XACTEMP_CCN        = 280.
!!$XAERDIFF           = 0.
!!$XAERHEIGHT         = 2000.
!!$LSCAV              = .FALSE.
!!$LAERO_MASS         = .FALSE.
! -----------------------------------------------------------------------
! lecture de la namelist LIMA
!  CALL POSNAM(NULNAM,'NAMLIMA')
!  READ(NULNAM,NAMLIMA)
! -----------------------------------------------------------------------
! initialisation des NSV
  ISV = 1

  NSV_LIMA_BEG = ISV
  IF (LWARM) THEN
! Nc
     NSV_LIMA_NC = ISV
     ISV = ISV+1
! Nr
     IF (LRAIN) THEN
        NSV_LIMA_NR = ISV
        ISV = ISV+1
     END IF
  END IF ! LWARM
! CCN
  IF (NMOD_CCN .GT. 0) THEN
     NSV_LIMA_CCN_FREE = ISV
     ISV = ISV + NMOD_CCN
     NSV_LIMA_CCN_ACTI = ISV
     ISV = ISV + NMOD_CCN
  END IF
! Scavenging
  IF (LSCAV .AND. LAERO_MASS) THEN
     NSV_LIMA_SCAVMASS = ISV
     ISV = ISV+1
  END IF ! LSCAV
! 
  IF (LCOLD) THEN
! Ni
     NSV_LIMA_NI = ISV
     ISV = ISV+1
  END IF ! LCOLD
! IFN
  IF (NMOD_IFN .GT. 0) THEN
     NSV_LIMA_IFN_FREE = ISV
     ISV = ISV + NMOD_IFN
     NSV_LIMA_IFN_NUCL = ISV
     ISV = ISV + NMOD_IFN
  END IF
! IMM
  IF (NMOD_IMM .GT. 0) THEN
     NSV_LIMA_IMM_NUCL = ISV
     ISV = ISV + MAX(1,NMOD_IMM)
  END IF
! Homogeneous freezing of CCN
  IF (LCOLD .AND. LHHONI) THEN
     NSV_LIMA_HOM_HAZE = ISV
     ISV = ISV + 1
  END IF
! End and total variables
  ISV = ISV - 1
  NSV_LIMA_END = ISV
  NSV_LIMA = NSV_LIMA_END - NSV_LIMA_BEG + 1

NSV=NSV_LIMA

! -----------------------------------------------------------------------
! initialisation de LIMA
CALL INIT_AEROSOL_PROPERTIES
! PDZMIN = 20 comme dans l'appel à INI_RAIN_ICE !
CALL INI_LIMA(PTSTEP, 20., KSPLITR, KSPLITG)




RETURN
END SUBROUTINE INIT_LIMA_PY
