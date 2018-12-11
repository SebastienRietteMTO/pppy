MODULE module_mp_python
   USE module_physics_init
   USE module_microphysics_driver

CONTAINS

SUBROUTINE other_init_py()
   USE module_state_description, ONLY: P_QNWFA, P_QNIFA
   P_QNWFA=1
   P_QNIFA=2
END SUBROUTINE other_init_py

SUBROUTINE mp_init_py(mp_physics,cycling_1byte, RAINNC,SNOWNC,GRAUPELNC,restart_1byte,warm_rain_1byte,      &
                      adv_moist_cond_1byte,                             &
                      MPDT, DT, DX, DY, LOWLYR,                   & ! for eta mp
                      F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY,           & ! for eta mp
                      mp_restart_state,tbpvs_state,tbpvs0_state,   & ! eta mp
                      allowed_to_read_1byte, start_of_simulation_1byte,       &
!CAMMGMP specific variables
                      ixcldliq, ixcldice, ixnumliq, ixnumice,     &       
                      nssl_cccn, nssl_alphah, nssl_alphahl,       &
                      nssl_ipelec, nssl_isaund,                  &
                         nssl_cnoh, nssl_cnohl,                  &
                         nssl_cnor, nssl_cnos,                   &
                         nssl_rho_qh, nssl_rho_qhl,              &
                         nssl_rho_qs,                            &
                      ccn_conc,                                  & ! RAS
                      z_at_q, qnwfa2d, scalar, num_sc,        &  ! G. Thompson
                      ids, ide, jds, jde, kds, kde,               &
                      ims, ime, jms, jme, kms, kme,               &
                      its, ite, jts, jte, kts, kte                )

IMPLICIT NONE
   LOGICAL(KIND=1) , INTENT(IN)        :: restart_1byte
   LOGICAL , INTENT(OUT)       :: warm_rain_1byte,adv_moist_cond_1byte
   REAL    , INTENT(IN)        :: MPDT, DT, DX, DY
   REAL, INTENT(IN)  :: nssl_cccn, nssl_alphah, nssl_alphahl, &
                                  nssl_cnoh, nssl_cnohl,                  &
                                  nssl_cnor, nssl_cnos,                   &
                                  nssl_rho_qh, nssl_rho_qhl,              &
                                  nssl_rho_qs
   INTEGER, INTENT(IN)  :: nssl_ipelec, nssl_isaund

   LOGICAL(KIND=1) , INTENT(IN)        :: start_of_simulation_1byte
   INTEGER , INTENT(IN)        :: ixcldliq, ixcldice, ixnumliq, ixnumice ! CAMMGMP specific variables

   INTEGER , INTENT(IN)        :: ids, ide, jds, jde, kds, kde,   &
                                  ims, ime, jms, jme, kms, kme,   &
                                  its, ite, jts, jte, kts, kte

   INTEGER , DIMENSION( ims:ime , jms:jme ) ,INTENT(INOUT)  :: LOWLYR
   REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(INOUT) :: RAINNC,SNOWNC,GRAUPELNC
   REAL,     DIMENSION( ims:ime , kms:kme, jms:jme ) , INTENT(INOUT) :: &
                                  F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY
   REAL , DIMENSION(43) ,INTENT(INOUT)  :: mp_restart_state !dimension found in eta_new
   REAL , DIMENSION(7501) ,INTENT(INOUT)  :: tbpvs_state,tbpvs0_state !dimension found in etanew
   LOGICAL(KIND=1) , INTENT(IN)  :: allowed_to_read_1byte
   REAL, INTENT(INOUT) :: ccn_conc ! RAS
   REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT):: qnwfa2d                       !  G. Thompson
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN):: z_at_q                   !  G. Thompson
   INTEGER, INTENT(IN) :: num_sc                                               !  G. Thompson
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme, num_sc), INTENT(INOUT):: scalar    !  G. Thompson

   INTEGER, INTENT(IN) :: mp_physics
   LOGICAL(KIND=1), INTENT(IN) :: cycling_1byte

  TYPE (grid_config_rec_type) ::     config_flags
  LOGICAL :: restart, warm_rain,adv_moist_cond, start_of_simulation, allowed_to_read, cycling

  !Logicals exchanged by interface must be 1-byte long but we need default length logical
  !to give to mp_init subroutine
  restart=restart_1byte
  start_of_simulation=start_of_simulation_1byte
  allowed_to_read=allowed_to_read_1byte
  cycling=cycling_1byte

  !config_flags is a fortran structured that cannot be passed through fortran/python interface
  !so, we pass usefull components of the structure
  config_flags%mp_physics = mp_physics
  config_flags%cycling = cycling

  CALL INIT_GFORTRAN_BIG_ENDIAN()

  CALL mp_init(RAINNC,SNOWNC,GRAUPELNC,config_flags,restart,warm_rain,      &
                      adv_moist_cond,                             &
                      MPDT, DT, DX, DY, LOWLYR,                   & 
                      F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY,           & 
                      mp_restart_state,tbpvs_state,tbpvs0_state,   & 
                      allowed_to_read, start_of_simulation,       &
                      ixcldliq, ixcldice, ixnumliq, ixnumice,     &       
                      nssl_cccn, nssl_alphah, nssl_alphahl,       &
                      nssl_ipelec, nssl_isaund,                  &
                         nssl_cnoh, nssl_cnohl,                  &
                         nssl_cnor, nssl_cnos,                   &
                         nssl_rho_qh, nssl_rho_qhl,              &
                         nssl_rho_qs,                            &
                      ccn_conc,                                  & 
                      z_at_q, qnwfa2d, scalar, num_sc,        &  
                      ids, ide, jds, jde, kds, kde,               &
                      ims, ime, jms, jme, kms, kme,               &
                      its, ite, jts, jte, kts, kte                )
  !Logicals exchanged with mp_init have the default compiler length but we need 1-byte long logical
  !to exchange through the interface
  adv_moist_cond_1byte=adv_moist_cond
  warm_rain_1byte=warm_rain

  CALL INIT_GFORTRAN_NATIVE_ENDIAN()

END SUBROUTINE mp_init_py

SUBROUTINE microphysics_driver_py(                                &
                       th, rho, pi_phy, p                                &
                      ,ht, dz8w, p8w, dt,dx,dy                           &
                      ,mp_physics, spec_zone                             &
                      ,specified_1byte, channel_switch_1byte                         &
                      ,warm_rain_1byte                                         &
                      ,t8w                                               &
                      ,chem_opt, progn                                   &
                      ,cldfra, cldfra_old, exch_h, nsource               &
                      ,qlsink, precr, preci, precs, precg                &
                      ,xland,snowh,itimestep                             &
                      ,f_ice_phy,f_rain_phy,f_rimef_phy                  &
                      ,lowlyr,sr, id                                     &
                      ,ids,ide, jds,jde, kds,kde                         &
                      ,ims,ime, jms,jme, kms,kme                         &
                      ,ips,ipe, jps,jpe, kps,kpe                         &
                      ,i_start,i_end,j_start,j_end,kts,kte               &
                      ,num_tiles, naer                                   &
!======================
                      !Variables required for CAMMGMP Scheme
                      ,dlf,dlf2,t_phy,p_hyd,p8w_hyd,tke_pbl,z_at_w,qfx   &
                      ,rliq,turbtype3d,smaw3d,wsedl3d,cldfra_old_mp      &
                      ,cldfra_mp,cldfra_mp_all,lradius,iradius           &
                      ,cldfrai,cldfral,cldfra_conv                       &
                      ,alt                                               &
                      ,accum_mode,aitken_mode,coarse_mode                &
                      ,icwmrsh3d,icwmrdp3d,shfrc3d,cmfmc3d,cmfmc2_3d     &
                      ,cycling_1byte ,fnm,fnp,rh_old_mp,lcd_old_mp         &
#if ( WRF_CHEM == 1 )
                      ,chem                                              &! For CAMMGMP scheme Prognostic aerosols
                      ,qme3d,prain3d,nevapr3d,rate1ord_cw2pr_st3d        &
!                      ,dgnum4D,dgnumwet4D                                &
#endif
!======================                                   
                      ,qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr   &
                      ,qic_curr,qip_curr,qid_curr &
                      ,qnic_curr,qnip_curr,qnid_curr &
                      ,qndrop_curr,qni_curr,qh_curr,qnh_curr             &
                      ,qzr_curr,qzi_curr,qzs_curr,qzg_curr,qzh_curr      &
                      ,qns_curr,qnr_curr,qng_curr,qnn_curr,qnc_curr      &
                      ,qnwfa_curr,qnifa_curr                             & ! for water/ice-friendly aerosols
                      ,f_qnwfa_1byte,f_qnifa_1byte                                   & ! for water/ice-friendly aerosols
                      ,qvolg_curr,qvolh_curr                             &
                      ,qir_curr,qib_curr                                 & ! for P3
                      ,effr_curr,ice_effr_curr,tot_effr_curr             &
                       ,qic_effr_curr,qip_effr_curr,qid_effr_curr        &             
                      ,f_qv_1byte,f_qc_1byte,f_qr_1byte,f_qi_1byte,f_qs_1byte,f_qg_1byte,f_qndrop_1byte,f_qni_1byte      &
                      ,f_qns_1byte,f_qnr_1byte,f_qng_1byte,f_qnc_1byte,f_qnn_1byte,f_qh_1byte,f_qnh_1byte          &
                      ,            f_qzr_1byte,f_qzi_1byte,f_qzs_1byte,f_qzg_1byte,f_qzh_1byte         &
                      ,f_qvolg_1byte,f_qvolh_1byte                                   &
                      ,f_qic_1byte,f_qip_1byte,f_qid_1byte &
                      ,f_qnic_1byte,f_qnip_1byte,f_qnid_1byte &
                      ,f_qir_1byte,f_qib_1byte                                       & ! for P3
                      ,f_effr_1byte,f_ice_effr_1byte,f_tot_effr_1byte                      &
                      ,f_qic_effr_1byte,f_qip_effr_1byte,f_qid_effr_1byte                  &                 
                      ,cu_used                                           &
                      ,qrcuten, qscuten, qicuten, qccuten                &
                      ,qt_curr,f_qt_1byte                                      &
                      ,mp_restart_state,tbpvs_state,tbpvs0_state         & ! for etampnew or fer_mp_hires
                      ,hail,ice2                                         & ! for mp_gsfcgce
!                     ,ccntype                                           & ! for mp_milbrandt2mom
                      ,u,v,w,z                                          &   
                      ,rainnc,    rainncv                                &
                      ,snownc,    snowncv                                &
                      ,hailnc,    hailncv                                &
                      ,graupelnc, graupelncv                             &
#if ( WRF_CHEM == 1 )
                      ,rainprod, evapprod                                &
                      ,qv_b4mp, qc_b4mp, qi_b4mp, qs_b4mp                &
#endif
                      ,qnwfa2d                                           & ! for water/ice-friendly aerosols
                      ,refl_10cm                                         & ! HM, 9/22/09, add for refl
                      ,vmi3d                                             & ! for P3 
                      ,di3d                                              & ! for P3 
                      ,rhopo3d                                           & ! for P3 
! YLIN
! Added the RI_CURR array to the call
                      ,ri_curr                                           &
                      ,diagflag_1byte,   do_radar_ref                          &
                      ,re_cloud, re_ice, re_snow                         & ! G. Thompson
                      ,has_reqc, has_reqi, has_reqs                      & ! G. Thompson
                      ,ccn_conc                                          & ! RAS
                      ,scalar,num_scalar                                   &
                      ,kext_ql,kext_qs,kext_qg            &
                      ,kext_qh,kext_qa                         &
                      ,kext_qic,kext_qid,kext_qip         &
                      ,kext_ft_qic,kext_ft_qid,kext_ft_qip         &
                      ,kext_ft_qs,kext_ft_qg            &
                      ,height,tempc &
                      ,TH_OLD                                            &
                      ,QV_OLD                                            &
                      ,xlat,xlong,ivgtyp                                 &
                      ,qrimef_curr,f_qrimef_1byte                              &
                                                   )


! Model Layer
  USE module_configure, only: grid_config_rec_type


   IMPLICIT NONE

   
   INTEGER, INTENT(IN) :: mp_physics
!-----------------------------------------------------------------------------------------------------
   LOGICAL, INTENT(IN) :: specified_1byte, cycling_1byte
!-----------------------------------------------------------------------------------------------------
   INTEGER, INTENT(IN) :: chem_opt, progn
   INTEGER, INTENT(IN) :: hail, ice2 !, ccntype
!
   INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde
   INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme,num_scalar
   INTEGER, INTENT(IN) :: ips,ipe, jps,jpe, kps,kpe
   INTEGER, INTENT(IN) :: kts,kte
   INTEGER, INTENT(IN) :: itimestep,num_tiles,spec_zone
   INTEGER, DIMENSION(num_tiles), INTENT(IN) ::                       &
     &           i_start,i_end,j_start,j_end

   LOGICAL, INTENT(IN) ::   warm_rain_1byte
!
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                    &
         INTENT(INOUT) ::                                         th
!

!
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme),                    &
         INTENT(IN) ::                                             &
                                                                 rho, &
                                                                dz8w, &
                                                                 p8w, &
                                                              pi_phy, &
                                                                   p
    REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),INTENT(INOUT) :: th_old,qv_old
    REAL,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar
    INTEGER, DIMENSION( ims:ime , jms:jme ), INTENT(IN)::   IVGTYP
    REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN)    :: XLAT, XLONG
    
!=================
!Data for CAMMGMP scheme
   REAL,INTENT(IN) ::accum_mode,aitken_mode,coarse_mode  
!1D variables required for CAMMGMP scheme
   REAL , DIMENSION( kms:kme ) ,                                      &
        INTENT(IN   ) ::                                        fnm,  & !Factors for interpolation at "w" grid (interfaces)
                                                                fnp     
!2D variables required for CAMMGMP scheme
   REAL, DIMENSION( ims:ime, jms:jme ),                               &
        INTENT(IN) ::                                                 &
                                                                 qfx, &    !Moisture flux at surface (kg m-2 s-1)
                                                                 rliq      !Vertically-integrated reserved cloud condensate(m/s)
 
 !3D variables required for CAMMGMP scheme
 REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                        &
      INTENT(IN) ::                                                   &
                                                                 dlf, &    !Detraining cloud water tendendcy
                                                                dlf2, &    !dq/dt due to export of cloud water into environment by shallow convection(kg/kg/s)
                                                               t_phy, &    !Temprature at the mid points (K)
                                                               p_hyd, &    !Hydrostatic pressure(Pa)
                                                             p8w_hyd, &    !Hydrostatic Pressure at level interface (Pa)
                                                              z_at_w, &    !Height above sea level at layer interfaces (m) 
                                                             tke_pbl, &    !Turbulence kinetic energy
                                                          turbtype3d, &    !Turbulent interface types [ no unit ]
                                                              smaw3d, &    !Normalized Galperin instability function for momentum  ( 0<= <=4.964 and 1 at neutral ) [no units]
                                                                 alt, &    !inverse density(m3/kg)
                                                           icwmrsh3d, &    !Shallow cumulus in-cloud water mixing ratio (kg/m2)
                                                           icwmrdp3d, &    !Deep Convection in-cloud water mixing ratio (kg/m2)
                                                             shfrc3d, &    !Shallow cloud fraction
                                                             cmfmc3d, &    !Deep + Shallow Convective mass flux [ kg /s/m^2 ]
                                                           cmfmc2_3d       !Shallow convective mass flux [ kg/s/m^2 ]
!#if ( WRF_CHEM == 1 )
!4D variables required for CAMMGMP scheme
! REAL, DIMENSION( ims:ime, kms:kme, jms:jme,ntot_amode_cam_mam ),     &
!        INTENT(IN) ::                                                 &
!                                                             dgnum4D, &
!                                                          dgnumwet4D 
!#endif
!In-outs
 REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                        &
      INTENT(INOUT) ::                                                &
                                                       cldfra_old_mp, &    !Old Cloud fraction for CAMMGMP microphysics only
                                                           rh_old_mp, &    !Old RH
                                                          lcd_old_mp       !Old liquid cloud fraction
!In-outs -optional
#if ( WRF_CHEM == 1 )
 REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem),     &
      INTENT(INOUT) ::                                                &
                                                                 chem      !Chem array for CAMMGMP scheme Prognostic aerosols      
#endif
!outs
REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                        &
      INTENT(INOUT)::                                                 &
                                                            wsedl3d, &    !Sedimentation velocity of stratiform liquid cloud droplet (m/s) 
                                                          cldfra_mp, &    !Old Cloud fraction for CAMMGMP microphysics only
                                                      cldfra_mp_all, &    !Old Cloud fraction for CAMMGMP microphysics only
                                                            cldfrai, &    !Old Cloud fraction for CAMMGMP microphysics only
                                                            cldfral, &    !Old Cloud fraction for CAMMGMP microphysics only
                                                            lradius, &    !Old Cloud fraction for CAMMGMP microphysics only
                                                            iradius, &    !Old Cloud fraction for CAMMGMP microphysics only                                                            
                                                        cldfra_conv 
#if ( WRF_CHEM == 1 )
REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                        &
      INTENT(INOUT) ::                                                 &
                                                              qme3d, &     !Net condensation rate (kg/kg/s)
                                                            prain3d, &     !Rate of conversion of condensate to precipitation (kg/kg/s)
                                                           nevapr3d, &     !Evaporation rate of rain + snow (kg/kg/s)
                                                rate1ord_cw2pr_st3d        !1st order rate for direct conversion of strat. cloud water to precip (1/s)
#endif

   REAL, INTENT(INOUT),  DIMENSION(ims:ime, kms:kme, jms:jme ) ::     &
                                     F_ICE_PHY,F_RAIN_PHY,F_RIMEF_PHY
!!$#if ( WRF_CHEM == 1 )
!  REAL, INTENT(OUT), DIMENSION(ims:ime, kms:kme, jms:jme ) ::     &
   REAL, INTENT(OUT), DIMENSION(ims:ime, kms:kme, jms:jme ) ::     &
!!$#else
!!$  REAL, DIMENSION(ims:ime, kms:kme, jms:jme ) ::     &
!!$#endif
         qlsink, & ! cloud water sink (/s)
         precr, & ! rain precipitation rate at all levels (kg/m2/s)
         preci, & ! ice precipitation rate at all levels (kg/m2/s)
         precs, & ! snow precipitation rate at all levels (kg/m2/s)
         precg    ! graupel precipitation rate at all levels (kg/m2/s)

!

   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN)   :: XLAND
   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN)   :: SNOWH

   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(OUT)   :: SR

   REAL, INTENT(IN   ) :: dt,dx,dy

   INTEGER, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT) :: LOWLYR

!
! Optional
!
   REAL, DIMENSION( ims:ime , kms:kme, jms:jme ) , INTENT(OUT) :: refl_10cm
   REAL, DIMENSION( ims:ime , kms:kme, jms:jme ) , INTENT(OUT) :: vmi3d,di3d,rhopo3d ! for P3

   LOGICAL,   INTENT(IN   )    :: channel_switch_1byte
   REAL,  INTENT(INOUT   ) :: naer  ! aerosol number concentration (/kg)
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN) :: qnwfa2d      ! Added by G. Thompson
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(INOUT ) ::                                        &
                  u,v,w, z, t8w                                       &
                 ,cldfra, cldfra_old, exch_h                      &
                 ,qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr &
                 ,qt_curr,qndrop_curr,qni_curr,qh_curr,qnh_curr   &
                 ,qns_curr,qnr_curr,qng_curr,qnn_curr,qnc_curr    &
                 ,qic_curr,qip_curr,qid_curr &
                 ,qnic_curr,qnip_curr,qnid_curr &
                 ,qzr_curr,qzi_curr,qzs_curr,qzg_curr,qzh_curr    &
                 ,qir_curr,qib_curr                               & ! for P3
                 ,effr_curr,ice_effr_curr,tot_effr_curr           &
                 ,qic_effr_curr,qip_effr_curr,qid_effr_curr           &
                  ,kext_ql,kext_qs,kext_qg          &
                 ,kext_qh,kext_qa                       &
                 ,kext_qic,kext_qip,kext_qid,tempc,height      &
                 ,kext_ft_qic,kext_ft_qip,kext_ft_qid &
                 ,kext_ft_qs,kext_ft_qg                           &
                 ,qnwfa_curr,qnifa_curr                           & ! Added by G. Thompson
                 ,qvolg_curr,qvolh_curr, qrimef_curr



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(IN) :: qrcuten, qscuten, qicuten, qccuten
   INTEGER, INTENT(IN), optional ::     cu_used
#if ( WRF_CHEM == 1 )
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(INOUT) :: rainprod, evapprod
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(INOUT) :: qv_b4mp, qc_b4mp, qi_b4mp, qs_b4mp
#endif
! YLIN
! Added RI_CURR similar to microphysics fields above
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                  &
         INTENT(INOUT) :: ri_curr


   REAL, DIMENSION(ims:ime, kms:kme, jms:jme ),                   &
         INTENT(OUT ) ::                                          &
                  nsource

!
   REAL, DIMENSION( ims:ime , jms:jme ),                          &
         INTENT(INOUT) ::                                            &
                                                           RAINNC &
                                                         ,RAINNCV &
                                                          ,SNOWNC &
                                                         ,SNOWNCV &
                                                       ,GRAUPELNC &
                                                      ,GRAUPELNCV &
                                                          ,HAILNC &
                                                         ,HAILNCV
   INTEGER,INTENT(IN   )    ::                        id

   REAL , DIMENSION( ims:ime , jms:jme )  ,             &
         INTENT(IN)   ::                                       ht

   REAL, DIMENSION (43), INTENT(INOUT) :: mp_restart_state !dimension found in eta_new
   REAL, DIMENSION (7501), INTENT(INOUT) :: tbpvs_state,tbpvs0_state !dimension found in eta_new
!

   LOGICAL, INTENT(IN) :: f_qv_1byte,f_qc_1byte,f_qr_1byte,f_qi_1byte,f_qs_1byte,f_qg_1byte, &
                       f_qndrop_1byte,f_qni_1byte,f_qt_1byte    &
                       ,f_qns_1byte,f_qnr_1byte,f_qng_1byte,f_qnn_1byte,f_qnc_1byte,f_qh_1byte,&
                       f_qnh_1byte,f_qzr_1byte       &
                      ,f_effr_1byte,f_ice_effr_1byte,f_tot_effr_1byte &
                       ,f_qic_effr_1byte,f_qip_effr_1byte,f_qid_effr_1byte &
                      ,f_qic_1byte,f_qip_1byte,f_qid_1byte &
                      ,f_qnic_1byte,f_qnip_1byte,f_qnid_1byte                                  &
                       ,f_qzi_1byte,f_qzs_1byte,f_qzg_1byte,f_qzh_1byte,f_qvolg_1byte,f_qvolh_1byte              &
                       ,f_qrimef_1byte                                             &
                       ,f_qir_1byte,f_qib_1byte                                          & ! for P3
                       ,f_qnwfa_1byte, f_qnifa_1byte                         ! Added by G. Thompson


   LOGICAL, INTENT(IN) :: diagflag_1byte
   REAL, INTENT(IN) :: ccn_conc ! RAS
   INTEGER, INTENT(IN) :: do_radar_ref
   REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) ::  & ! G. Thompson
                 re_cloud, re_ice, re_snow
   INTEGER, INTENT(IN):: has_reqc, has_reqi, has_reqs
!  REAL , DIMENSION( ims:ime , jms:jme ), INTENT(INOUT) :: lwp

 
   TYPE(grid_config_rec_type):: config_flags
   LOGICAL :: specified, cycling, warm_rain, channel_switch, f_qv,f_qc,f_qr,&
              f_qi,f_qs,f_qg,f_qndrop,f_qni,f_qt    &
                       ,f_qns,f_qnr,f_qng,f_qnn,f_qnc,f_qh,f_qnh,f_qzr       &
                      ,f_effr,f_ice_effr,f_tot_effr &
                       ,f_qic_effr,f_qip_effr,f_qid_effr &
                      ,f_qic,f_qip,f_qid &
                      ,f_qnic,f_qnip,f_qnid                                  &
                       ,f_qzi,f_qzs,f_qzg,f_qzh,f_qvolg,f_qvolh              &
                       ,f_qrimef                                             &
                       ,f_qir,f_qib                                          & ! for P3
                       ,f_qnwfa, f_qnifa, diagflag

  !Logicals exchanged by interface must be 1-byte long but we need default length logical
  !to give to mp_init subroutine
  specified=specified_1byte
  cycling=cycling_1byte
  warm_rain=warm_rain_1byte
  channel_switch=channel_switch_1byte
  f_qv=f_qv_1byte
  f_qc=f_qc_1byte
  f_qr=f_qr_1byte
  f_qi=f_qi_1byte
  f_qs=f_qs_1byte
  f_qg=f_qg_1byte
  f_qndrop=f_qndrop_1byte
  f_qni=f_qni_1byte
  f_qt=f_qt_1byte
  f_qns=f_qns_1byte
  f_qnr=f_qnr_1byte
  f_qng=f_qng_1byte
  f_qnn=f_qnn_1byte
  f_qnc=f_qnc_1byte
  f_qh=f_qh_1byte
  f_qnh=f_qnh_1byte
  f_qzr=f_qzr_1byte
  f_effr=f_effr_1byte
  f_ice_effr=f_ice_effr_1byte
  f_tot_effr=f_tot_effr_1byte
  f_qic_effr=f_qic_effr_1byte
  f_qip_effr=f_qip_effr_1byte
  f_qid_effr=f_qid_effr_1byte
  f_qic=f_qic_1byte
  f_qip=f_qip_1byte
  f_qid=f_qid_1byte
  f_qnic=f_qnic_1byte
  f_qnip=f_qnip_1byte
  f_qnid=f_qnid_1byte
  f_qzi=f_qzi_1byte
  f_qzs=f_qzs_1byte
  f_qzg=f_qzg_1byte
  f_qzh=f_qzh_1byte
  f_qvolg=f_qvolg_1byte
  f_qvolh=f_qvolh_1byte
  f_qrimef=f_qrimef_1byte
  f_qir=f_qir_1byte
  f_qib=f_qib_1byte
  f_qnwfa=f_qnwfa_1byte
  f_qnifa=f_qnifa_1byte
  diagflag=diagflag_1byte

   !config_flags is a fortran structured that cannot be passed through fortran/python interface
   !so, we pass usefull components of the structure
   config_flags%mp_physics = mp_physics
   config_flags%cycling = cycling

   CALL INIT_GFORTRAN_BIG_ENDIAN()


   CALL microphysics_driver(                                          &
                       th, rho, pi_phy, p                                &
                      ,ht, dz8w, p8w, dt,dx,dy                           &
                      ,mp_physics, spec_zone                             &
                      ,specified, channel_switch                         &
                      ,warm_rain                                         &
                      ,t8w                                               &
                      ,chem_opt, progn                                   &
                      ,cldfra, cldfra_old, exch_h, nsource               &
                      ,qlsink, precr, preci, precs, precg                &
                      ,xland,snowh,itimestep                             &
                      ,f_ice_phy,f_rain_phy,f_rimef_phy                  &
                      ,lowlyr,sr, id                                     &
                      ,ids,ide, jds,jde, kds,kde                         &
                      ,ims,ime, jms,jme, kms,kme                         &
                      ,ips,ipe, jps,jpe, kps,kpe                         &
                      ,i_start,i_end,j_start,j_end,kts,kte               &
                      ,num_tiles, naer                                   &
!======================
                      !Variables required for CAMMGMP Scheme
                      ,dlf,dlf2,t_phy,p_hyd,p8w_hyd,tke_pbl,z_at_w,qfx   &
                      ,rliq,turbtype3d,smaw3d,wsedl3d,cldfra_old_mp      &
                      ,cldfra_mp,cldfra_mp_all,lradius,iradius           &
                      ,cldfrai,cldfral,cldfra_conv                       &
                      ,alt                                               &
                      ,accum_mode,aitken_mode,coarse_mode                &
                      ,icwmrsh3d,icwmrdp3d,shfrc3d,cmfmc3d,cmfmc2_3d     &
                      ,config_flags,fnm,fnp,rh_old_mp,lcd_old_mp         &
#if ( WRF_CHEM == 1 )
                      ,chem                                              &! For CAMMGMP scheme Prognostic aerosols
                      ,qme3d,prain3d,nevapr3d,rate1ord_cw2pr_st3d        &
!                      ,dgnum4D,dgnumwet4D                                &
#endif
!======================                                   
                      ,qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr   &
                      ,qic_curr,qip_curr,qid_curr &
                      ,qnic_curr,qnip_curr,qnid_curr &
                      ,qndrop_curr,qni_curr,qh_curr,qnh_curr             &
                      ,qzr_curr,qzi_curr,qzs_curr,qzg_curr,qzh_curr      &
                      ,qns_curr,qnr_curr,qng_curr,qnn_curr,qnc_curr      &
                      ,qnwfa_curr,qnifa_curr                             & ! for water/ice-friendly aerosols
                      ,f_qnwfa,f_qnifa                                   & ! for water/ice-friendly aerosols
                      ,qvolg_curr,qvolh_curr                             &
                      ,qir_curr,qib_curr                                 & ! for P3
                      ,effr_curr,ice_effr_curr,tot_effr_curr             &
                       ,qic_effr_curr,qip_effr_curr,qid_effr_curr        &             
                      ,f_qv,f_qc,f_qr,f_qi,f_qs,f_qg,f_qndrop,f_qni      &
                      ,f_qns,f_qnr,f_qng,f_qnc,f_qnn,f_qh,f_qnh          &
                      ,            f_qzr,f_qzi,f_qzs,f_qzg,f_qzh         &
                      ,f_qvolg,f_qvolh                                   &
                      ,f_qic,f_qip,f_qid &
                      ,f_qnic,f_qnip,f_qnid &
                      ,f_qir,f_qib                                       & ! for P3
                      ,f_effr,f_ice_effr,f_tot_effr                      &
                      ,f_qic_effr,f_qip_effr,f_qid_effr                  &                 
                      ,cu_used                                           &
                      ,qrcuten, qscuten, qicuten, qccuten                &
                      ,qt_curr,f_qt                                      &
                      ,mp_restart_state,tbpvs_state,tbpvs0_state         & ! for etampnew or fer_mp_hires
                      ,hail,ice2                                         & ! for mp_gsfcgce
!                     ,ccntype                                           & ! for mp_milbrandt2mom
                      ,u,v,w,z                                          &   
                      ,rainnc,    rainncv                                &
                      ,snownc,    snowncv                                &
                      ,hailnc,    hailncv                                &
                      ,graupelnc, graupelncv                             &
#if ( WRF_CHEM == 1 )
                      ,rainprod, evapprod                                &
                      ,qv_b4mp, qc_b4mp, qi_b4mp, qs_b4mp                &
#endif
                      ,qnwfa2d                                           & ! for water/ice-friendly aerosols
                      ,refl_10cm                                         & ! HM, 9/22/09, add for refl
                      ,vmi3d                                             & ! for P3 
                      ,di3d                                              & ! for P3 
                      ,rhopo3d                                           & ! for P3 
! YLIN
! Added the RI_CURR array to the call
                      ,ri_curr                                           &
                      ,diagflag,   do_radar_ref                          &
                      ,re_cloud, re_ice, re_snow                         & ! G. Thompson
                      ,has_reqc, has_reqi, has_reqs                      & ! G. Thompson
                      ,ccn_conc                                          & ! RAS
                      ,scalar,num_scalar                                   &
                      ,kext_ql,kext_qs,kext_qg            &
                      ,kext_qh,kext_qa                         &
                      ,kext_qic,kext_qid,kext_qip         &
                      ,kext_ft_qic,kext_ft_qid,kext_ft_qip         &
                      ,kext_ft_qs,kext_ft_qg            &
                      ,height,tempc &
                      ,TH_OLD                                            &
                      ,QV_OLD                                            &
                      ,xlat,xlong,ivgtyp                                 &
                      ,qrimef_curr,f_qrimef                              &
                                                   )


CALL INIT_GFORTRAN_NATIVE_ENDIAN()

END SUBROUTINE microphysics_driver_py

END MODULE module_mp_python
