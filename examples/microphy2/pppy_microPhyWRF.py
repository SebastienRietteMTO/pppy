#!/usr/bin/env python 3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

#Written by A. Riandet, L. Richecoeur and M.L. Roussel
#Modified by S. Riette

"""
This module contains an implementation of PPPY suitable to call microphysical
schemes from the WRF model. This implementation call the microphysical driver
of the model, hence, it is virtually possible to call any of the WRF scheme.
In reality, some adjustments are necessary (like suppression of sedimentation,
initialisation of specific variables...) and only a subset of schemes have
been tested.

Here are some notes about the WRF schemes that are implemented:

option keys:
- enableCCNsource for THOMPSONAERO
- dx and dy for THOMPSONAERO and FULL_KHAIN_LYNN
- enableCCNinit for FULL_KHAIN_LYNN

THOMPSONAERO, FAST_KHAIN_LYNN, FULL_KHAIN_LYNN
In WRF, all number concentrations are stored in a variable called scalar.
This variable is 'unpacked' before calling mp_init and microphysics_driver.
Then the different microphysics schemes does not have to worry about the scalar variable
except THOMPSONAERO, FAST_KHAIN_LYNN and FULL_KHAIN_LYNN.
For the initialization of THOMPSONAERO, CCN and IN number must be in scalar variable.
For FAST_KHAIN_LYNN and FULL_KHAIN_LYNN, content of each bin is stored in the scalar variable,
then this variable is the pronostic variable of these schemes.

THOMPSONAERO
* Variable qnwfa2d is set during initialization of THOMPSONAERO to represent a ground tendency
of CCN. This tendency is, then, applied during integration to refill nwfa variable near ground.
CCN concentration takes into account a source if option enableCCNsource is True and discards
this source otherwise.
after having called self._mp_init_py
* If nwfa and/or nifa are null, they are filled with a profile computed from height

FULL_KHAIN_LYNN
* Scheme computes derivatives. Remember that, presently, halo is filled with the mean.
* In normal usage, scheme initializes CCN concentration when:
  - this is the first timestep
  - this is not the first timestep and dx is greater than 7500.
  If option enableCCNinit is False, a timestep different from 1 is always given to scheme
  to switch off first case. Second case can be switched off by providing a small dx (scheme option).
  In case we go through initialization, altitude is considered to be of some meters (depending of
  the number of points chosen). To change this behavior, code must be updated to enable the usage
  of a user's defined dz8w variable.
* Scheme needs u, v and w to extrapolate in time tendency of T and qv. We need to provide
  not null dz8w for derivative computation.
  Normally this extrapolation gives zero in a 0D case because all points are identical.
* scheme crashes for initial condition P=100000., T=290., rv=rc=1.E-2, rr=1.E-4,
  ri=ric=rid=rip=rs=rg=rh=0., nc=3.E8, nr=2000., ni=ns=ng=nh=nic=nid=nip=0.,
  ccn1ft=ccn1at=1.E8, ifn1ft=ifn1at=0., u=v=w=0., xland=1.
  for dt=20s. Iterations of breakup subroutine leads to infinite values.
"""

import numpy
import logging
import glob
import os
import shutil
from pppy import ctypesForFortran
import pppy

class pppy_microPhyWRF(pppy.PPPY):
    """
    PPPY implementation for calling microphysical schemes of the WRF model.
    """

    #Constants used to identify schemes
    KESSLERSCHEME = 1
    LINSCHEME = 2
    WSM3SCHEME = 3
    WSM5SCHEME = 4
    WSM6SCHEME = 6
    FER_MP_HIRES = 5
    FER_MP_HIRES_ADVECT = 15
    GSFCGCESCHEME = 7
    THOMPSON = 8
    THOMPSONAERO = 28
    MILBRANDT2MOM = 9
    MORR_TWO_MOMENT = 10
    CAMMGMPSCHEME = 11
    SBU_YLINSCHEME = 13
    WDM5SCHEME = 14
    WDM6SCHEME = 16
    NSSL_2MOM = 17
    NSSL_2MOMCCN = 18
    NSSL_1MOM = 19
    NSSL_1MOMLFO = 21
    NSSL_2MOMG = 22
    FAST_KHAIN_LYNN = 30
    FULL_KHAIN_LYNN = 32
    P3_1CATEGORY = 50
    P3_1CATEGORY_NC = 51
    ETAMPNEW = 95

    def __init__(self, dt, method, name, tag, solib,
                 mp_physics,
                 enableCCNsource=None, dx=None, dy=None,
                 enableCCNinit=None):
        """
        In addition to dt, method, name and tag parameters
        defined in the PPPY class, this parameterization
        needs the following parameters:
        solib           : path to the shared library to use
        mp_physics      : number of the scheme to run
        enableCCNsource : True or False to enable CCN source term for THOMPSONAERO scheme
        dx, dy          : size of the grid box for THOMPSONAERO and FULL_KHAIN_LYNN schemes
        enableCCNinit   : True or False to compute initial values for CCN for FULL_KHAIN_LYNN

        This methods also deals with the resource files.
        """
        WRF_options = dict(mp_physics=mp_physics)
        self._mp_physics = mp_physics
        if mp_physics == self.THOMPSONAERO:
            WRF_options['enableCCNsource'] = enableCCNsource
            self._enableCCNsource = enableCCNsource
        if mp_physics in [self.THOMPSONAERO, self.FULL_KHAIN_LYNN]:
            WRF_options['dx'] = dx
            WRF_options['dy'] = dy
            self._dx, self._dy = dx, dy
        if mp_physics == self.FULL_KHAIN_LYNN:
            WRF_options['enableCCNinit'] = enableCCNinit
            self._enableCCNinit = enableCCNinit

        super().__init__(dt, method, name, tag,
                         solib=solib, **WRF_options)
        self._solib = solib

        self._handle = None
        self._mp_init_py = None
        self._nl_set_force_read_thompson = None
        self._nl_set_write_thompson_tables = None
        self._microphysics_driver_py = None
        self._other_init_py = None
        self._resources_files = []
        self._new_resources_files = []

    def setup(self, init_state, duration):
        """
        This method opens the shared library and do the init stuff.
        """
        super().setup(init_state, duration)

        #Opening of shared lib is done here to avoid loading it at initialization
        #Because if shared lib is laod at init, all declared schemes will have
        #their shared lib loaded simultaneously which can be a source of error.
        #This way, shared lib is loaded in a sub-process (if class instance is used
        #through a PPPYComp instance).

        IN = ctypesForFortran.IN
        OUT = ctypesForFortran.OUT
        INOUT = ctypesForFortran.INOUT

        ctypesFF, self._handle = ctypesForFortran.ctypesForFortranFactory(self._solib)

        @ctypesFF(prefix='__module_mp_python_MOD_', suffix='')
        def mp_init_py(mp_physics, cycling,
                       RAINNC, SNOWNC, GRAUPELNC,
                       restart,
                       MPDT, DT, DX, DY, LOWLYR, F_ICE_PHY, F_RAIN_PHY, F_RIMEF_PHY,
                       mp_restart_state, tbpvs_state,tbpvs0_state, allowed_to_read,
                       start_of_simulation, ixcldliq, ixcldice, ixnumliq, ixnumice,
                       nssl_cccn, nssl_alphah, nssl_alphahl, nssl_ipelec, nssl_isaund,
                       nssl_cnoh, nssl_cnohl, nssl_cnor, nssl_cnos, nssl_rho_qh,
                       nssl_rho_qhl, nssl_rho_qs, ccn_conc,z_at_q, qnwfa2d, scalar, num_sc,
                       ids, ide, jds, jde, kds, kde, ims, ime,
                       jms, jme, kms, kme, its, ite, jts, jte, kts, kte):
            "This function calls the mp_init_py fortran subroutine"
            shape2Dhalo = (ime-ims+1, jme-jms+1) #Shape of a 2D field with halo
            shape3Dhalo = (ime-ims+1, kme-kms+1, jme-jms+1) #Shape of a 3D field with halo
            return ([mp_physics, cycling,
                     RAINNC, SNOWNC, GRAUPELNC,
                     restart,
                     MPDT, DT, DX, DY, LOWLYR, F_ICE_PHY, F_RAIN_PHY, F_RIMEF_PHY,
                     mp_restart_state, tbpvs_state,tbpvs0_state, allowed_to_read,
                     start_of_simulation, ixcldliq, ixcldice, ixnumliq, ixnumice,
                     nssl_cccn, nssl_alphah, nssl_alphahl, nssl_ipelec, nssl_isaund,
                     nssl_cnoh, nssl_cnohl, nssl_cnor, nssl_cnos, nssl_rho_qh,
                     nssl_rho_qhl, nssl_rho_qs, ccn_conc,z_at_q, qnwfa2d, scalar, num_sc,
                     ids, ide, jds, jde, kds, kde, ims, ime,
                     jms, jme, kms, kme, its, ite, jts, jte, kts, kte],
                    [(numpy.int32, None, IN), #INTEGER, INTENT(IN) :: mp_physics
                     (numpy.bool, None, IN), #cycling :
                     (numpy.float32, shape2Dhalo, INOUT), #REAL, DIMENSION(ims:ime,jms:jme) , INTENT(INOUT) :: RAINNC
                     (numpy.float32, shape2Dhalo, INOUT), #REAL, DIMENSION(ims:ime,jms:jme) , INTENT(INOUT) :: SNOWNC
                     (numpy.float32, shape2Dhalo, INOUT), #REAL, DIMENSION(ims:ime,jms:jme) , INTENT(INOUT) :: GRAUPELNC
                     (numpy.bool, None, IN), #LOGICAL, INTENT(IN) :: restart
                     (numpy.bool, None, OUT), #LOGICAL , INTENT(OUT) :: warm_rain
                     (numpy.bool, None, OUT), #LOGICAL , INTENT(OUT) :: adv_moist_cond
                     (numpy.float32, None, IN), #REAL, INTENT(IN) :: MPDT
                     (numpy.float32, None, IN), #REAL, INTENT(IN) :: DT
                     (numpy.float32, None, IN), #REAL, INTENT(IN) :: DX
                     (numpy.float32, None, IN), #REAL, INTENT(IN) :: DY
                     (numpy.int32, shape2Dhalo, INOUT), #INTEGER , DIMENSION( ims:ime , jms:jme ) ,INTENT(INOUT)  :: LOWLYR
                     (numpy.float32, shape3Dhalo, INOUT), #REAL,     DIMENSION( ims:ime , kms:kme, jms:jme ) , INTENT(INOUT) :: F_ICE_PHY
                     (numpy.float32, shape3Dhalo, INOUT), #REAL,     DIMENSION( ims:ime , kms:kme, jms:jme ) , INTENT(INOUT) :: F_RAIN_PHY
                     (numpy.float32, shape3Dhalo, INOUT), #REAL,     DIMENSION( ims:ime , kms:kme, jms:jme ) , INTENT(INOUT) :: F_RIMEF_PHY
                     (numpy.float32, (43, ), INOUT), #REAL , DIMENSION(:) ,INTENT(INOUT)  :: mp_restart_state
                     (numpy.float32, (7501, ), INOUT), #REAL , DIMENSION(:) ,INTENT(INOUT)  :: tbpvs_state
                     (numpy.float32, (7501, ), INOUT), #REAL , DIMENSION(:) ,INTENT(INOUT)  :: tbpvs0_state
                     (numpy.bool, None, IN), #LOGICAL , INTENT(IN)  :: allowed_to_read
                     (numpy.bool, None, IN), #LOGICAL , INTENT(IN) :: start_of_simulation
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: ixcldliq1
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: ixcldice1
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: ixnumliq1
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: ixnumice1
                     (numpy.float32, None, IN), #REAL, INTENT(IN), OPTIONAL  :: nssl_cccn
                     (numpy.float32, None, IN), #REAL, INTENT(IN), OPTIONAL  :: nssl_alphah
                     (numpy.float32, None, IN), #REAL, INTENT(IN), OPTIONAL  :: nssl_alphahl
                     (numpy.int32, None, IN), #INTEGER, INTENT(IN), OPTIONAL  :: nssl_ipelec
                     (numpy.int32, None, IN), #INTEGER, INTENT(IN), OPTIONAL  :: nssl_isaund
                     (numpy.float32, None, IN), #REAL, INTENT(IN), OPTIONAL  :: nssl_cnoh
                     (numpy.float32, None, IN), #REAL, INTENT(IN), OPTIONAL  :: nssl_cnohl
                     (numpy.float32, None, IN), #REAL, INTENT(IN), OPTIONAL  :: nssl_cnor
                     (numpy.float32, None, IN), #REAL, INTENT(IN), OPTIONAL  :: nssl_cnos
                     (numpy.float32, None, IN), #REAL, INTENT(IN), OPTIONAL  :: nssl_rho_qh
                     (numpy.float32, None, IN), #REAL, INTENT(IN), OPTIONAL  :: nssl_rho_qhl
                     (numpy.float32, None, IN), #REAL, INTENT(IN), OPTIONAL  :: nssl_rho_qs
                     (numpy.float32, None, INOUT), #REAL, INTENT(INOUT) :: ccn_conc
                     (numpy.float32, shape3Dhalo, IN), #REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN):: z_at_q   
                     (numpy.float32, shape2Dhalo, INOUT), #REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT):: qnwfa2d 
                     (numpy.float32, tuple(list(shape3Dhalo) + [num_sc]), INOUT), #REAL, DIMENSION(ims:ime,kms:kme,jms:jme, num_sc), INTENT(INOUT):: scalar  
                     (numpy.int32, None, IN), #INTEGER, INTENT(IN) :: num_sc 
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: ids
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: ide
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: jds
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: jde
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: kds
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: kde
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: ims
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: ime
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: jms
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: jme
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: kms
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: kme
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: its
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: ite
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: jts
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: jte
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: kts
                     (numpy.int32, None, IN), #INTEGER , INTENT(IN)        :: kte
                    ],
                    None)
        self._mp_init_py = mp_init_py

        @ctypesFF(prefix='', suffix='_')
        def nl_set_force_read_thompson(id_id, force_read_thompson):
            "This function calls the nl_set_force_read_thompson fortran subroutine"
            return ([id_id, force_read_thompson],
                    [(numpy.int32, None, IN),
                     (numpy.bool, None, IN)],
                    None)
        self._nl_set_force_read_thompson = nl_set_force_read_thompson

        @ctypesFF(prefix='', suffix='_')
        def nl_set_write_thompson_tables(id_id, write_thompson_tables):
            "This function calls the nl_set_write_thompson_tables fortran subroutine"
            return ([id_id, write_thompson_tables],
                    [(numpy.int32, None, IN),
                     (numpy.bool, None, IN)],
                    None)
        self._nl_set_write_thompson_tables = nl_set_write_thompson_tables

        @ctypesFF(prefix='__module_mp_python_MOD_', suffix='')
        def microphysics_driver_py(th, rho, pi_phy, p,
                                   ht, dz8w, p8w, dt,dx,dy,
                                   mp_physics, spec_zone,
                                   specified, channel_switch,
                                   warm_rain,
                                   t8w,
                                   chem_opt, progn,
                                   cldfra, cldfra_old, exch_h,
                                   xland,snowh,itimestep,
                                   f_ice_phy,f_rain_phy,f_rimef_phy,
                                   lowlyr, gid,
                                   ids,ide, jds,jde, kds,kde,
                                   ims,ime, jms,jme, kms,kme,
                                   ips,ipe, jps,jpe, kps,kpe,
                                   i_start,i_end,j_start,j_end,kts,kte,
                                   num_tiles, naer,
                                   # Variables required for CAMMGMP Scheme
                                   dlf,dlf2,t_phy,p_hyd,p8w_hyd,tke_pbl,z_at_w,qfx,
                                   rliq,turbtype3d,smaw3d,wsedl3d,cldfra_old_mp,
                                   cldfra_mp,cldfra_mp_all,lradius,iradius,
                                   cldfrai,cldfral,cldfra_conv,
                                   alt,
                                   accum_mode,aitken_mode,coarse_mode,
                                   icwmrsh3d,icwmrdp3d,shfrc3d,cmfmc3d,cmfmc2_3d,
                                   cycling,
                                   fnm,fnp,rh_old_mp,lcd_old_mp,
                                   #if ( WRF_CHEM == 1 )
                                   #chem,
                                   #For CAMMGMP scheme Prognostic aerosols
                                   #qme3d,prain3d,nevapr3d,rate1ord_cw2pr_st3d,
                                   #dgnum4D,dgnumwet4D,
                                   #endif
                                   qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr,
                                   qic_curr,qip_curr,qid_curr,
                                   qnic_curr,qnip_curr,qnid_curr,
                                   qndrop_curr,qni_curr,qh_curr,qnh_curr,
                                   qzr_curr,qzi_curr,qzs_curr,qzg_curr,qzh_curr,
                                   qns_curr,qnr_curr,qng_curr,qnn_curr,qnc_curr,
                                   qnwfa_curr,qnifa_curr,
                                   # for water/ice-friendly aerosols
                                   f_qnwfa,f_qnifa,
                                   # for water/ice-friendly aerosols
                                   qvolg_curr,qvolh_curr,
                                   qir_curr,qib_curr,
                                   # for P3
                                   effr_curr,ice_effr_curr,tot_effr_curr,
                                   qic_effr_curr,qip_effr_curr,qid_effr_curr,       
                                   f_qv,f_qc,f_qr,f_qi,f_qs,f_qg,f_qndrop,f_qni,
                                   f_qns,f_qnr,f_qng,f_qnc,f_qnn,f_qh,f_qnh,
                                   f_qzr,f_qzi,f_qzs,f_qzg,f_qzh,
                                   f_qvolg,f_qvolh,
                                   f_qic,f_qip,f_qid,
                                   f_qnic,f_qnip,f_qnid,
                                   f_qir,f_qib,
                                   # for P3
                                   f_effr,f_ice_effr,f_tot_effr,
                                   f_qic_effr,f_qip_effr,f_qid_effr,           
                                   cu_used,
                                   qrcuten, qscuten, qicuten, qccuten,
                                   qt_curr,f_qt,
                                   mp_restart_state,tbpvs_state,tbpvs0_state ,
                                   # for etampnew or fer_mp_hires
                                   hail,ice2,
                                   # for mp_gsfcgce
                                   #ccntype ,
                                   # for mp_milbrandt2mom
                                   u,v,w,z,
                                   rainnc,    rainncv,
                                   snownc,    snowncv,
                                   hailnc,    hailncv,
                                   graupelnc, graupelncv,
                                   #if ( WRF_CHEM == 1 )
                                   #rainprod, evapprod,
                                   #qv_b4mp, qc_b4mp, qi_b4mp, qs_b4mp,
                                   #endif
                                   qnwfa2d,
                                   # for water/ice-friendly aerosols
                                   # HM, 9/22/09, add for refl
                                   # for P3 
                                   # for P3 
                                   # for P3 
                                   # YLIN
                                   # Added the RI_CURR ariy to the call
                                   ri_curr,
                                   diagflag,   do_radar_ref ,
                                   re_cloud, re_ice, re_snow,
                                   # G. Thompson
                                   has_reqc, has_reqi, has_reqs,
                                   # G. Thompson
                                   ccn_conc,
                                   # RAS
                                   scalar,num_scalar,
                                   kext_ql,kext_qs,kext_qg,
                                   kext_qh,kext_qa,
                                   kext_qic,kext_qid,kext_qip,
                                   kext_ft_qic,kext_ft_qid,kext_ft_qip,
                                   kext_ft_qs,kext_ft_qg,
                                   height,tempc,
                                   TH_OLD ,
                                   QV_OLD,
                                   xlat,xlong,ivgtyp,
                                   qrimef_curr,f_qrimef):
            "This function calls the microphysics_driver_py fortran subroutine"
            shape3Dhalo = (ime-ims+1, kme-kms+1, jme-jms+1) #Shape of a 3D field with halo
            shape2Dhalo = (ime-ims+1, jme-jms+1) #Shape of a 2D field with halo
            shape1Dhalo = (kme-kms+1, ) #Shape of a 1D field with halo
            varList = [th, rho, pi_phy, p,
                      ht, dz8w, p8w, dt,dx,dy,
                      mp_physics, spec_zone,
                      specified, channel_switch,
                      warm_rain,
                      t8w,
                      chem_opt, progn,
                      cldfra, cldfra_old, exch_h,
                      xland,snowh,itimestep,
                      f_ice_phy,f_rain_phy,f_rimef_phy,
                      lowlyr, gid,
                      ids,ide, jds,jde, kds,kde,
                      ims,ime, jms,jme, kms,kme,
                      ips,ipe, jps,jpe, kps,kpe,
                      i_start,i_end,j_start,j_end,kts,kte,
                      num_tiles, naer,
                      dlf,dlf2,t_phy,p_hyd,p8w_hyd,tke_pbl,z_at_w,qfx,
                      rliq,turbtype3d,smaw3d,wsedl3d,cldfra_old_mp,
                      cldfra_mp,cldfra_mp_all,lradius,iradius,
                      cldfrai,cldfral,cldfra_conv,
                      alt,
                      accum_mode,aitken_mode,coarse_mode,
                      icwmrsh3d,icwmrdp3d,shfrc3d,cmfmc3d,cmfmc2_3d,
                      cycling,
                      fnm,fnp,rh_old_mp,lcd_old_mp,
                      ##if ( WRF_CHEM == 1 )
                      #chem,
                      #qme3d,prain3d,nevapr3d,rate1ord_cw2pr_st3d,
                      #!dgnum4D,dgnumwet4D           ,
                      ##endif
                      qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr,
                      qic_curr,qip_curr,qid_curr,
                      qnic_curr,qnip_curr,qnid_curr,
                      qndrop_curr,qni_curr,qh_curr,qnh_curr,
                      qzr_curr,qzi_curr,qzs_curr,qzg_curr,qzh_curr,
                      qns_curr,qnr_curr,qng_curr,qnn_curr,qnc_curr,
                      qnwfa_curr,qnifa_curr,
                      f_qnwfa,f_qnifa,
                      qvolg_curr,qvolh_curr,
                      qir_curr,qib_curr,
                      effr_curr,ice_effr_curr,tot_effr_curr,
                      qic_effr_curr,qip_effr_curr,qid_effr_curr            ,
                      f_qv,f_qc,f_qr,f_qi,f_qs,f_qg,f_qndrop,f_qni,
                      f_qns,f_qnr,f_qng,f_qnc,f_qnn,f_qh,f_qnh,
                      f_qzr,f_qzi,f_qzs,f_qzg,f_qzh,
                      f_qvolg,f_qvolh,
                      f_qic,f_qip,f_qid,
                      f_qnic,f_qnip,f_qnid,
                      f_qir,f_qib,
                      f_effr,f_ice_effr,f_tot_effr,
                      f_qic_effr,f_qip_effr,f_qid_effr               ,
                      cu_used,
                      qrcuten, qscuten, qicuten, qccuten,
                      qt_curr,f_qt,
                      mp_restart_state,tbpvs_state,tbpvs0_state ,
                      hail,ice2,
                      #ccntype,
                      u,v,w,z   ,
                      rainnc,    rainncv,
                      snownc,    snowncv,
                      hailnc,    hailncv,
                      graupelnc, graupelncv,
                      ##if ( WRF_CHEM == 1 )
                      #rainprod, evapprod,
                      #qv_b4mp, qc_b4mp, qi_b4mp, qs_b4mp,
                      ##endif
                      qnwfa2d,
                      ri_curr,
                      diagflag,   do_radar_ref ,
                      re_cloud, re_ice, re_snow,
                      has_reqc, has_reqi, has_reqs,
                      ccn_conc,
                      scalar,num_scalar,
                      kext_ql,kext_qs,kext_qg,
                      kext_qh,kext_qa,
                      kext_qic,kext_qid,kext_qip,
                      kext_ft_qic,kext_ft_qid,kext_ft_qip,
                      kext_ft_qs,kext_ft_qg,
                      height,tempc,
                      TH_OLD ,
                      QV_OLD,
                      xlat,xlong,ivgtyp,
                      qrimef_curr,f_qrimef]

            signature = [(numpy.float32, shape3Dhalo, INOUT), #th : potential temperature   
                         (numpy.float32, shape3Dhalo, IN), #rho : density of air 
                         (numpy.float32, shape3Dhalo, IN), #pi_phy : exner function 
                         (numpy.float32, shape3Dhalo, IN), #p : pressure
                         (numpy.float32, shape2Dhalo, IN), #ht
                         (numpy.float32, shape3Dhalo, IN), #dz8w : dz between full levels (m)
                         (numpy.float32, shape3Dhalo, IN), #p8w : pressure at full levels (Pa)
                         (numpy.float32, None, IN), #dt : time step
                         (numpy.float32, None, IN), #dx : space variation
                         (numpy.float32, None, IN), #dy : space variation
                         (numpy.int32, None, IN), #mp_physics : scheme index
                         (numpy.int32, None, IN), #spec_zone : number of points in specified zone
                         (numpy.bool, None, IN), #specified : 
                         (numpy.bool, None, IN), #channel_switch : 
                         (numpy.bool, None, IN), #warm_rain : 
                         (numpy.float32, shape3Dhalo, INOUT), #t8w : temperature at layer interfaces
                         (numpy.int32, None, IN), #chem_opt : 
                         (numpy.int32, None, IN), #progn : 
                         (numpy.float32, shape3Dhalo, INOUT), # cldfra : current cloud fraction
                         (numpy.float32, shape3Dhalo, INOUT), # cldfra_old : previous cloud fraction
                         (numpy.float32, shape3Dhalo, INOUT), # exch_h : vertical diffusivity
                         (numpy.float32, shape3Dhalo, OUT), #nsource : OUT
                         (numpy.float32, shape3Dhalo, OUT), #qlsink : fractional cloud water sink OUT
                         (numpy.float32, shape3Dhalo, OUT), #precr : rain precipitation rate at all levels OUT
                         (numpy.float32, shape3Dhalo, OUT), #preci : ice precipitation rate at all levels OUT
                         (numpy.float32, shape3Dhalo, OUT), #precs : snow precipitation rate at all levels OUT
                         (numpy.float32, shape3Dhalo, OUT), #precg: graupel precipitation rate at all levels OUT
                         (numpy.float32, shape2Dhalo, IN), #xland :
                         (numpy.float32, shape2Dhalo, IN), #snowh : 
                         (numpy.int32, None, IN), #itimestep : if 1 use ice processes in full sbm scheme
                         (numpy.float32, shape3Dhalo, INOUT), #f_ice_phy : fraction of ice
                         (numpy.float32, shape3Dhalo, INOUT), #f_rain_phy : fraction of rain
                         (numpy.float32, shape3Dhalo, INOUT), #f_rimef_phy : mass ratio of rimed ice
                         (numpy.int32, shape2Dhalo, INOUT), #lowlyr : 
                         (numpy.float32, shape2Dhalo, OUT), #sr : one time step mass ratio of snow to total precipitation OUT
                         (numpy.int32, None, IN), #id : grid id number
                         (numpy.int32, None, IN), #ids : start index for i in domain
                         (numpy.int32, None, IN), #ide : end index for i in domain
                         (numpy.int32, None, IN), #jds : start index for j in domain
                         (numpy.int32, None, IN), #jde : end index for j in domain
                         (numpy.int32, None, IN), #kds : start index for k in domain
                         (numpy.int32, None, IN), #kde : end index for k in domain
                         (numpy.int32, None, IN), #ims : start index for i in memory
                         (numpy.int32, None, IN), #ime : end index for i in memory
                         (numpy.int32, None, IN), #jms : start index for j in memory
                         (numpy.int32, None, IN), #jme : end index for j in memory
                         (numpy.int32, None, IN), #kms : start index for k in memory
                         (numpy.int32, None, IN), #kme : end index for k in memory
                         (numpy.int32, None, IN), #ips : start index for i
                         (numpy.int32, None, IN), #ipe : end index for i
                         (numpy.int32, None, IN), #jps : start index for j
                         (numpy.int32, None, IN), #jpe : end index for j
                         (numpy.int32, None, IN), #kps : start index for k
                         (numpy.int32, None, IN), #kpe : end index for k
                         (numpy.int32, (num_tiles,), IN), #i_start : start indice for i in tile
                         (numpy.int32, (num_tiles,), IN), #i_end : end indice for i in tile 
                         (numpy.int32, (num_tiles,), IN), #j_start : start indice for j in tile
                         (numpy.int32, (num_tiles,), IN), #j_end : end indice for j in tile
                         (numpy.int32, None, IN), #kts : start index for k in tile
                         (numpy.int32, None, IN), #kte : end index for k in tile
                         (numpy.int32, None, IN), #num_tiles : number of tiles
                         (numpy.float32, None, INOUT), #naer : aerosol number concentration
                         (numpy.float32, shape3Dhalo, IN), #dlf : detraining cloud water tendency
                         (numpy.float32, shape3Dhalo, IN), #dlf2 : dq/dt due to export of cloud water into environment by shallow convection
                         (numpy.float32, shape3Dhalo, IN), #t_phy : temperature at the mid points
                         (numpy.float32, shape3Dhalo, IN), #p_hyd : hydrostatic pressure
                         (numpy.float32, shape3Dhalo, IN), #p8w_hyd : hydrostatic pressure at level interface
                         (numpy.float32, shape3Dhalo, IN), #tke_pbl : turbulence kinetic energy
                         (numpy.float32, shape3Dhalo, IN), #z_at_w : height above sea level at layer interfaces 
                         (numpy.float32, shape2Dhalo,IN), #qfx : moisture flux at surface
                         (numpy.float32, shape2Dhalo,IN), #rliq : vertically integrated reserved cloud condensate
                         (numpy.float32, shape3Dhalo, IN), #turbtype3d : turbulence interface type 
                         (numpy.float32, shape3Dhalo, IN), #smaw3d : Normalized Galperin instability function for momentum 
                         (numpy.float32, shape3Dhalo, INOUT), #wsedl3d : sedimentation velocity of stratiform liquid cloud dropplets
                         (numpy.float32, shape3Dhalo, INOUT), #cldfra_old_mp : old cloud fraction for cammgmp microphysics only
                         (numpy.float32, shape3Dhalo, INOUT), #cldfra_mp : old cloud fraction for cammgmp microphysics only
                         (numpy.float32, shape3Dhalo, INOUT), #cldfra_mp_all : old cloud fraction for cammgmp microphysics only
                         (numpy.float32, shape3Dhalo, INOUT), #lradius : old cloud fraction for cammgmp microphysics only
                         (numpy.float32, shape3Dhalo, INOUT), #iradius : old cloud fraction for cammgmp microphysics only
                         (numpy.float32, shape3Dhalo, INOUT), #cldfrai : old cloud fraction for cammgmp microphysics only  
                         (numpy.float32, shape3Dhalo, INOUT), #cldfral : old cloud fraction for cammgmp microphysics only
                         (numpy.float32, shape3Dhalo, INOUT), #cldfra_conv : 
                         (numpy.float32, shape3Dhalo, IN), #alt : inverse density
                         (numpy.float32, None, IN), #accum_mode : 
                         (numpy.float32, None, IN), #aitken_mode :
                         (numpy.float32, None, IN), #coarse_mode : 
                         (numpy.float32, shape3Dhalo, IN), #icwmrsh3d : shallow cumulus in cloud water mixing ratio
                         (numpy.float32, shape3Dhalo, IN), #icwmrdp3d : deep convection in cloud water mixing ratio
                         (numpy.float32, shape3Dhalo, IN), #shfrc3d : shallow cloud fraction
                         (numpy.float32, shape3Dhalo, IN), #cmfmc3d : deep + shallow convective mass flux
                         (numpy.float32, shape3Dhalo, IN), #cmfmc2_3d : shallow convective mass flux
                         (numpy.int32, None, IN), #cycling
                         (numpy.float32, shape1Dhalo, IN), #fnm : factors for interpolation at WGRID
                         (numpy.float32, shape1Dhalo, IN), #fnp : factors for interpolation at WGRID
                         (numpy.float32, shape3Dhalo, INOUT), #rh_old_mp : old rh
                         (numpy.float32, shape3Dhalo, INOUT), #lcd_old_mp :  old liquid cloud fraction
                         #(numpy.float32, tuple(list(shape3Dhalo) + [1]), INOUT), #chem : scheme array for cammgmp scheme pronostic aerosol (dim : num_chem)
                         #(numpy.float32, shape3Dhalo, INOUT), #qme3d : net condensation rate
                         #(numpy.float32, shape3Dhalo, INOUT), #prain3d : rate of conversion of condensate to precipitation
                         #(numpy.float32, shape3Dhalo, INOUT), #nevapr3d : evaporation rate of rain + snow
                         #(numpy.float32, shape3Dhalo, INOUT), #rate1ord_cw2pr_st3d : first order rate for direct conversion strat cloud water to precip
                         #(numpy.float32, tuple(list(shape3Dhalo) + [3]), IN), #dgnum4D (dim : ntot_amode_cam_mam)
                         #(numpy.float32, tuple(list(shape3Dhalo) + [3]), IN), #dgnumwet4D (dim : ntot_amode_cam_mam)
                         (numpy.float32, shape3Dhalo, INOUT), #qv_curr : vapor mixing ratio
                         (numpy.float32, shape3Dhalo, INOUT), #qc_curr : cloud mixing ratio
                         (numpy.float32, shape3Dhalo, INOUT), #qr_curr : rain mixing ratio
                         (numpy.float32, shape3Dhalo, INOUT), #qi_curr : ice mixing ratio
                         (numpy.float32, shape3Dhalo, INOUT), #qs_curr : snow mixing ratio
                         (numpy.float32, shape3Dhalo, INOUT), #qg_curr : graupel mixing ratio
                         (numpy.float32, shape3Dhalo, INOUT), #qic_curr : ice column mixing ratio
                         (numpy.float32, shape3Dhalo, INOUT), #qip_curr : ice plates mixing ratio
                         (numpy.float32, shape3Dhalo, INOUT), #qid_curr : ice dendrites mixing ratio
                         (numpy.float32, shape3Dhalo, INOUT), #qnic_curr : number of concentration of ice column
                         (numpy.float32, shape3Dhalo, INOUT), #qnip_curr : number of concentration of ice plates
                         (numpy.float32, shape3Dhalo, INOUT), #qnid_curr : number of concentration of ice dendrites
                         (numpy.float32, shape3Dhalo, INOUT), #qndrop_curr 
                         (numpy.float32, shape3Dhalo, INOUT), #qni_curr : number of concentration of ice 
                         (numpy.float32, shape3Dhalo, INOUT), #qh_curr : hail mixing ratio
                         (numpy.float32, shape3Dhalo, INOUT), #qnh_curr : number of concentration of hail
                         (numpy.float32, shape3Dhalo, INOUT), #qzr_curr
                         (numpy.float32, shape3Dhalo, INOUT), #qzi_curr
                         (numpy.float32, shape3Dhalo, INOUT), #qzs_curr
                         (numpy.float32, shape3Dhalo, INOUT), #qzg_curr
                         (numpy.float32, shape3Dhalo, INOUT), #qzh_curr
                         (numpy.float32, shape3Dhalo, INOUT), #qns_curr : number of concentration of snow
                         (numpy.float32, shape3Dhalo, INOUT), #qnr_curr : number of concentration of rain
                         (numpy.float32, shape3Dhalo, INOUT), #qng_curr : number of concentration of graupel
                         (numpy.float32, shape3Dhalo, INOUT), #qnn_curr
                         (numpy.float32, shape3Dhalo, INOUT), #qnc_curr
                         (numpy.float32, shape3Dhalo, INOUT), #qnwfa_curr
                         (numpy.float32, shape3Dhalo, INOUT), #qnifa_curr
                         (numpy.bool, None, IN), #f_qnwfa : 
                         (numpy.bool, None, IN), #f_qnifa :
                         (numpy.float32, shape3Dhalo, INOUT), #qvolg_curr :
                         (numpy.float32, shape3Dhalo, INOUT), #qvolh_curr
                         (numpy.float32, shape3Dhalo, INOUT), #qir_curr :
                         (numpy.float32, shape3Dhalo, INOUT), #qib_curr :
                         (numpy.float32, shape3Dhalo, INOUT), #effr_curr : 
                         (numpy.float32, shape3Dhalo, INOUT), #ice_effr_curr :
                         (numpy.float32, shape3Dhalo, INOUT), #tot_effr_curr :
                         (numpy.float32, shape3Dhalo, INOUT), #qic_effr_curr :
                         (numpy.float32, shape3Dhalo, INOUT), #qip_effr_curr :
                         (numpy.float32, shape3Dhalo, INOUT), #qid_effr_curr :
                         (numpy.bool, None, IN), #f_qv :
                         (numpy.bool, None, IN), #f_qc :
                         (numpy.bool, None, IN), #f_qr :
                         (numpy.bool, None, IN), #f_qi :
                         (numpy.bool, None, IN), #f_qs :
                         (numpy.bool, None, IN), #f_qg :
                         (numpy.bool, None, IN), #f_qndrop :
                         (numpy.bool, None, IN), #f_qni :
                         (numpy.bool, None, IN), #f_qns :
                         (numpy.bool, None, IN), #f_qnr :
                         (numpy.bool, None, IN), #f_qng :
                         (numpy.bool, None, IN), #f_qnc :
                         (numpy.bool, None, IN), #f_qnn :
                         (numpy.bool, None, IN), #f_qh :
                         (numpy.bool, None, IN), #f_qnh :
                         (numpy.bool, None, IN), #f_qzr :
                         (numpy.bool, None, IN), #f_qzi :
                         (numpy.bool, None, IN), #f_qzs :
                         (numpy.bool, None, IN), #f_qzg :
                         (numpy.bool, None, IN), #f_qzh :
                         (numpy.bool, None, IN), #f_qvolg :
                         (numpy.bool, None, IN), #f_qvolh : 
                         (numpy.bool, None, IN), #f_qic :
                         (numpy.bool, None, IN), #f_qip :
                         (numpy.bool, None, IN), #f_qid :
                         (numpy.bool, None, IN), #f_qnic :
                         (numpy.bool, None, IN), #f_qnip :
                         (numpy.bool, None, IN), #f_qnid : 
                         (numpy.bool, None, IN), #f_qir :
                         (numpy.bool, None, IN), #f_qib :
                         (numpy.bool, None, IN), #f_effr :
                         (numpy.bool, None, IN), #f_ice_effr :
                         (numpy.bool, None, IN), #f_tot_effr :
                         (numpy.bool, None, IN), #f_qic_effr :
                         (numpy.bool, None, IN), #f_qip_effr : 
                         (numpy.bool, None, IN), #f_qid_effr :
                         (numpy.int32, None, IN), #cu_used : 
                         (numpy.float32,shape3Dhalo, IN), #qrcuten : 
                         (numpy.float32,shape3Dhalo, IN), #qscuten : 
                         (numpy.float32,shape3Dhalo, IN), #qicuten :
                         (numpy.float32,shape3Dhalo, IN), #qccuten : 
                         (numpy.float32, shape3Dhalo, INOUT), #qt_curr :
                         (numpy.bool, None, IN), #f_qt :
                         (numpy.float32, (43, ), INOUT), #mp_restart_state :
                         (numpy.float32, (7501, ), INOUT), #tbpvs_state :
                         (numpy.float32, (7501, ), INOUT), #tbpvs0_state : 
                         (numpy.int32, None, IN), #hail : ccn type
                         (numpy.int32, None, IN), #ice2 : ccn type
                         #(), #ccntype :
                         (numpy.float32, shape3Dhalo, INOUT), #u :
                         (numpy.float32, shape3Dhalo, INOUT), #v :
                         (numpy.float32, shape3Dhalo, INOUT), #w :
                         (numpy.float32, shape3Dhalo, INOUT), #z :
                         (numpy.float32, shape2Dhalo, INOUT ), #rainnc :
                         (numpy.float32, shape2Dhalo, INOUT ), #rainncv :
                         (numpy.float32, shape2Dhalo, INOUT ), #snownc :
                         (numpy.float32, shape2Dhalo, INOUT ), #snowncv :
                         (numpy.float32, shape2Dhalo, INOUT ), #hailnc :
                         (numpy.float32, shape2Dhalo, INOUT ), #hailncv :
                         (numpy.float32, shape2Dhalo, INOUT ), #graupelnc :
                         (numpy.float32, shape2Dhalo, INOUT ), #graupelncv : 
                         #(numpy.float32, shape3Dhalo, INOUT), #rainprod :
                         #(numpy.float32, shape3Dhalo, INOUT), #evapprod :
                         #(numpy.float32, shape3Dhalo, INOUT), #qv_b4mp :
                         #(numpy.float32, shape3Dhalo, INOUT), #qc_b4mp
                         #(numpy.float32, shape3Dhalo, INOUT), #qi_b4mp :
                         #(numpy.float32, shape3Dhalo, INOUT), #qs_b4mp :
                         (numpy.float32, shape2Dhalo, IN), #qnwfa2d : 
                         (numpy.float32, shape3Dhalo, OUT), #refl_10cm : OUT
                         (numpy.float32, shape3Dhalo, OUT), #vmi3d :OUT
                         (numpy.float32, shape3Dhalo, OUT), #di3d :OUT
                         (numpy.float32, shape3Dhalo, OUT), #rhopo3d : OUT
                         (numpy.float32, shape3Dhalo, INOUT), #ri_curr : OUT
                         (numpy.bool, None, IN), #diagflag : logical to tell us when to produce diagnostic for history or restart
                         (numpy.int32, None, IN), #do_radar_ref 
                         (numpy.float32, shape3Dhalo, INOUT), #re_cloud : 
                         (numpy.float32, shape3Dhalo, INOUT), #re_ice :
                         (numpy.float32, shape3Dhalo, INOUT), #re_snow :
                         (numpy.int32, None, IN), #has_reqc : 
                         (numpy.int32, None, IN), #has_reqi :
                         (numpy.int32, None, IN), #has_reqs :
                         (numpy.float32, None, IN), #ccn_conc :
                         (numpy.float32, tuple(list(shape3Dhalo) + [num_scalar]), INOUT), #scalar : 
                         (numpy.int32, None, IN), #num_scalar :
                         (numpy.float32, shape3Dhalo, INOUT), #kext_ql :
                         (numpy.float32, shape3Dhalo, INOUT), #kext_qs :
                         (numpy.float32, shape3Dhalo, INOUT), #kext_qg :
                         (numpy.float32, shape3Dhalo, INOUT), #kext_qh :
                         (numpy.float32, shape3Dhalo, INOUT), #kext_qa :
                         (numpy.float32, shape3Dhalo, INOUT), #kext_qic :
                         (numpy.float32, shape3Dhalo, INOUT), #kext_qid :
                         (numpy.float32, shape3Dhalo, INOUT), #kext_qip :
                         (numpy.float32, shape3Dhalo, INOUT), #kext_ft_qic :
                         (numpy.float32, shape3Dhalo, INOUT), #kext_ft_qid :
                         (numpy.float32, shape3Dhalo, INOUT), #kext_ft_qip :
                         (numpy.float32, shape3Dhalo, INOUT), #kext_ft_qs :
                         (numpy.float32, shape3Dhalo, INOUT), #kext_ft_qg : 
                         (numpy.float32, shape3Dhalo, INOUT), #height :
                         (numpy.float32, shape3Dhalo, INOUT), #tempc :
                         (numpy.float32, shape3Dhalo, INOUT), #TH_OLD :
                         (numpy.float32, shape3Dhalo, INOUT), #QV_OLD :
                         (numpy.float32, shape2Dhalo, IN), #xlat :
                         (numpy.float32, shape2Dhalo, IN), #xlong : 
                         (numpy.int32, shape2Dhalo, IN), #ivgtyp :
                         (numpy.float32,shape3Dhalo, INOUT), #qrimef_curr
                         (numpy.bool, None, IN), #f_qrimef
                        ]
            return (varList, signature, None)
        self._microphysics_driver_py = microphysics_driver_py

        @ctypesFF(prefix='__module_mp_python_MOD_', suffix='')
        def other_init_py():
            "This function calls the other_init_py fortran subroutine"
            return ([], [], None)
        self._other_init_py = other_init_py

        #Performs the init
        self._other_init_py()

        #Copy resource files to current working directory
        self._resources_files = []
        path = os.path.dirname(os.path.abspath(__file__))
        for filename in glob.glob(path + "/resources_WRF/*"):
            dest = os.path.basename(filename)
            self._resources_files.append(dest)
            shutil.copyfile(filename, dest)

        #Some resource files can be missing for Thompson scheme
        #In this case, these files will be computed during the first use
        #of the scheme. This section of code asks the scheme to write
        #these missing files to disk and fill the _new_resources_files
        #variable to allow the finalize method to store them with the other
        #WRF resource files.
        self._nl_set_write_thompson_tables(1, False)
        self._nl_set_force_read_thompson(1, False)
        self._new_resources_files = []
        if self._mp_physics in [self.THOMPSON, self.THOMPSONAERO]:
            missing = False
            for filename in ['freezeH2O.dat', 'qr_acr_qg.dat', 'qr_acr_qs.dat']:
                if not os.path.exists(filename):
                    self._new_resources_files.append((filename,
                                                      path + "/resources_WRF/" + filename))
                    missing = True
            if missing:
                #With this option, table will be saved on disk in current directory
                #Saved file will be moved in resources_WRF directory in method finalize
                self._nl_set_write_thompson_tables(1, True)
            else:
                #With this option, an error is issued if read fails
                self._nl_set_force_read_thompson(1, True)

    def finalize(self):
        """
        The method closes the shared library so we are able to load
        a new one with the same symbol names.
        Moreover missing resource files for Thompson (see end of setup method) are moved
        to the resource folder and other resource files are removed.
        """
        super().finalize()
        try:
            del self._mp_init_py, self._microphysics_driver_py
            del self._nl_set_force_read_thompson, self._nl_set_write_thompson_tables
            del self._other_init_py
            ctypesForFortran.dlclose(self._handle)
            for src, dest in self._new_resources_files:
                shutil.copyfile(src, dest)
                os.remove(src)
            for filename in self._resources_files:
                os.remove(filename)
        except KeyboardInterrupt:
            raise
        except:
            pass

    def build_init_state(self, state):
        """
        Needed variables are strongly scheme-dependent.
        This method computes the missing variables to insure a comparable
        initial state for all schemes
        """
        state = super().build_init_state(state)

        needed = ['T', 'P', 'rv', 'rc', 'rr']
        if self._mp_physics in [self.THOMPSON, self.THOMPSONAERO]:
            needed.extend(['w', 'ri', 'rh', 'rs', 'rg',
                           'nc', 'nr', 'ni', 'ns', 'ng', 'nh', 'ccn1ft'])
        if self._mp_physics == self.THOMPSONAERO:
            needed.extend(['ccn1ft', 'ifn1ft'])
        if self._mp_physics in [self.FULL_KHAIN_LYNN]:
            needed.extend(['u', 'v', 'w', 'xland'])
        for var in needed:
            if var not in state:
                raise ValueError(var + " must be in state")
        if all([item in state for item in ['ri', 'ric', 'rip', 'rid']]):
            if not numpy.allclose(state['ri'], state['ric'] + state['rip'] + state['rid']):
                raise ValueError("ri and ric+rip+rid are not consistent")

        if self._mp_physics == self.KESSLERSCHEME:
            pass
        elif self._mp_physics in [self.THOMPSON, self.THOMPSONAERO]:
            if 'nwfa2d' not in state:
                #Values will be set during mp_init
                state['nwfa2d'] = numpy.zeros(tuple(list(state['T'].shape)))
        elif self._mp_physics == self.FULL_KHAIN_LYNN:
            #Old values
            P0 = 100000.
            Boltz = 1.380658E-23
            Avogadro = 6.0221367E+23
            Md = 28.9644E-3
            Rd = Avogadro * Boltz / Md
            Cpd = 7. * Rd / 2.
            state['qv_old'] = state['rv'].copy()
            pi_phy = (state['P'] / P0) ** (Rd / Cpd) #exner function
            state['th_old'] = state['T'] / pi_phy #potential temperature

            #Bin values
            if not 'scalar' in state:
                logging.warning("Initialization of scalar must be done according to the " +
                                "distribution of LIMA (if number available) or ICE " +
                                "(how to use nic, nid and nip if available?)")
                scalar = numpy.zeros(tuple(list(state['T'].shape) + [265]))
                scalar[..., 1:18] = state['rc'] / 17.
                scalar[..., 18:34] = state['rr'] / 16.
                scalar[..., 133:166] = state['ric'] / 33.
                scalar[..., 166:199] = state['rip'] / 33.
                scalar[..., 199:232] = state['rid'] / 33.
                scalar[..., 67:100] = state['rg'] / 33.
                scalar[..., 232:265] = state['rh'] / 33.
                scalar[..., 100:133] = state['ccn1ft'] / 33.
                scalar[..., 34:67] = state['rs'] / 33.
                state['scalar'] = scalar
            else:
                if state['scalar'].shape[-1] != 265:
                    raise ValueError("Size of scalar variable is not OK")
                if 'rc' in state and 'rr' in state:
                    #1:18 for rc, 18:34 for rr
                    if not numpy.allclose(scalar[..., 1:34].sum(axis=-1),
                                          state['rc'] + state['rr']):
                        raise ValueError("scalar and rc+rr are not consistent.")
                if 'ric' in state:
                    if not numpy.allclose(scalar[..., 133:166].sum(axis=-1), state['ric']):
                        raise ValueError("scalar and ric are not consistent.")
                if 'rip' in state:
                    if not numpy.allclose(scalar[..., 166:199].sum(axis=-1), state['rip']):
                        raise ValueError("scalar and rip are not consistent.")
                if 'rid' in state:
                    if not numpy.allclose(scalar[..., 199:232].sum(axis=-1), state['rid']):
                        raise ValueError("scalar and rid are not consistent.")
                if 'ri' in state:
                    if not numpy.allclose(scalar[..., 133:231].sum(axis=-1), state['ri']):
                        raise ValueError("scalar and ri are not consistent.")
                if 'rs' in state:
                    if not numpy.allclose(scalar[..., 34:67].sum(axis=-1), state['rs']):
                        raise ValueError("scalar and rs are not consistent.")
                if 'rg' in state:
                    if not numpy.allclose(scalar[..., 67:100].sum(axis=-1), state['rg']):
                        raise ValueError("scalar and rg are not consistent.")
                if 'rh' in state:
                    if not numpy.allclose(scalar[..., 232:265].sum(axis=-1), state['rh']):
                        raise ValueError("scalar and rh are not consistent.")
                if 'ccn1ft' in state:
                    if not numpy.allclose(scalar[..., 100:133].sum(axis=-1), state['ccn1ft']):
                        raise ValueError("scalar and ccn1ft are not consistent.")
        else:
            raise NotImplementedError("Initialization have not been checked for this scheme")

        return state

    @staticmethod
    def _to32(array):
        "return array expressed on 32 bits"
        return array.astype(dtype={numpy.dtype(numpy.float64):numpy.float32,
                                   numpy.dtype(numpy.int64):numpy.int32}[array.dtype])

    @staticmethod
    def _to64(array):
        "return array expressed on 64 bits"
        return array.astype(dtype={numpy.dtype(numpy.float32):numpy.float64,
                                   numpy.dtype(numpy.int32):numpy.int64}[array.dtype])

    @staticmethod
    def _add_halo(array):
        "Adds halo to array"
        #Three first dimensions are geographical dimensions
        #return array with array 1-point halo
        if len(array.shape) >= 3:
            new_shape = tuple([array.shape[0] + 2,
                               array.shape[1] + 2,
                               array.shape[2] + 2] + list(array.shape[3:]))
            result = numpy.zeros(new_shape, dtype=array.dtype)
            result[...] = array.mean()
            result[1:-1, 1:-1, 1:-1, ...] = array
        elif len(array.shape) == 2:
            new_shape = (array.shape[0] + 2, array.shape[1] + 2)
            result = numpy.zeros(new_shape, dtype=array.dtype)
            result[...] = array.mean()
            result[1:-1, 1:-1] = array
        elif len(array.shape) == 1:
            new_shape = (array.shape[0] + 2, )
            result = numpy.zeros(new_shape, dtype=array.dtype)
            result[...] = array.mean()
            result[1:-1] = array
        return result

    @staticmethod
    def _del_halo(array):
        "Removes halo to array"
        #Three first dimensions are geographical dimensions
        #return array without the 1-point halo
        if len(array.shape) >= 3:
            return array[1:-1, 1:-1, 1:-1, ...]
        elif len(array.shape) == 2:
            return array[1:-1, 1:-1]
        elif len(array.shape) == 1:
            return array[1:-1]

    def execute(self, previous_state, dt, timestep_number):
        super().execute(previous_state, dt, timestep_number)

        #Dimensions
        shapeOri3D = previous_state['T'].shape
        if len(shapeOri3D) > 3:
            raise ValueError("Maximum shape length is 3")
        if len(shapeOri3D) == 1:
            shape3D = tuple(list(shapeOri3D) + [1, 1])
            shapeOri2D = shape3D
        elif len(shapeOri3D) == 2:
            shape3D = tuple(list(shapeOri3D) + [1])
            shapeOri2D = shapeOri3D
        elif len(shapeOri3D) == 3:
            shape3D = shapeOri3D
            shapeOri2D = shapeOri3D[:2]

        if shape3D[1:] != (1, 1):
            raise NotImplementedError("It seems that we need to change dimension order. " +
                                      "We receive vertical dimension in last index but it seems " +
                                      "that WRF expect vertical dim to be in second position.")
        shape2D = shape3D[0:2]
        shape3Dhalo = (shape3D[0] + 2, shape3D[2] + 2, shape3D[1] + 2)
        shape2Dhalo = (shape3D[0] + 2, shape3D[1] + 2)
        shape1Dhalo = (shape3D[2] + 2, )
        #shapeOri3D: shape of an input state variables which can be 3D
        #shapeOri2D: shape of an input state variables which can be, at most, 2D
        #shape3D: 3D shape of input state variables, vertical dim is in third position
        #shape2D: 2D shape of input state variables
        #shape3DHalo: 3D shape with halo, vertical dim is in second position
        #shape2DHalo: 2D shape with halo
        #shape1Dhalo: 1D (vertical) shape with halo
        zeros2D = numpy.zeros(shape2Dhalo, numpy.float32)
        zeros3D = numpy.zeros(shape3Dhalo, numpy.float32)

        #Index
        ids, jds, kds = 1, 1, 1 # start domain
        ide, jde, kde = shape3D[0] + 4, shape3D[1] + 4, shape3D[2] + 4 # end domain
        ims, jms, kms = 2, 2, 2 # start memory
        ime, jme, kme = shape3D[0] + 3, shape3D[1] + 3, shape3D[2] + 3 # end memory
        its, jts, kts = 3, 3, 3 # start tile
        ite, jte, kte = shape3D[0] + 2, shape3D[1] + 2, shape3D[2] + 2 # end tile
        ips, jps, kps = 1, 1, 1 # start ?
        ipe, jpe, kpe = 1, 1, 1 # end ?
        num_tiles = 1
        i_start = numpy.ones((num_tiles, ), numpy.int32) * its
        i_end = numpy.ones((num_tiles, ), numpy.int32) * ite
        j_start = numpy.ones((num_tiles, ), numpy.int32) * jts
        j_end = numpy.ones((num_tiles, ), numpy.int32) * jte

        #Constants
        P0 = 100000.
        Boltz = 1.380658E-23
        Avogadro = 6.0221367E+23
        Md = 28.9644E-3
        Rd = Avogadro * Boltz / Md
        Rv = 461.6
        Cpd = 7. * Rd / 2.

        #Previous state values
        add_halo = self._add_halo
        to32 = self._to32
        T = add_halo(to32(previous_state['T'].reshape(shape3D))) #Temperature (K)
        p = add_halo(to32(previous_state['P'].reshape(shape3D))) #pressupre (Pa)
        qv_curr = add_halo(to32(previous_state['rv'].reshape(shape3D))) #mixing ratio (kg/kg)
        qc_curr = add_halo(to32(previous_state['rc'].reshape(shape3D))) #mixing ratio (kg/kg)
        qr_curr = add_halo(to32(previous_state['rr'].reshape(shape3D))) #mixing ratio (kg/kg)
        qi_curr = add_halo(to32(previous_state['ri'].reshape(shape3D))) #mixing ratio (kg/kg)
        qic_curr = add_halo(to32(previous_state['ric'].reshape(shape3D))) if 'ric' in previous_state else zeros3D.copy() #mixing ratio (kg/kg)
        qip_curr = add_halo(to32(previous_state['rip'].reshape(shape3D))) if 'rip' in previous_state else  zeros3D.copy() #mixing ratio (kg/kg)
        qid_curr = add_halo(to32(previous_state['rid'].reshape(shape3D))) if 'rid' in previous_state else zeros3D.copy() #mixing ratio (kg/kg)
        qs_curr = add_halo(to32(previous_state['rs'].reshape(shape3D))) if 'rs' in previous_state else zeros3D.copy() #mixing ratio (kg/kg)
        qg_curr = add_halo(to32(previous_state['rg'].reshape(shape3D))) if 'rg' in previous_state else zeros3D.copy() #mixing ratio (kg/kg)
        qh_curr = add_halo(to32(previous_state['rh'].reshape(shape3D))) if 'rh' in previous_state else zeros3D.copy() #mixing ratio (kg/kg)
        qnc_curr = add_halo(to32(previous_state['nc'].reshape(shape3D))) if 'nc' in previous_state else zeros3D.copy() #number concentration (#/kg)
        qnr_curr = add_halo(to32(previous_state['nr'].reshape(shape3D))) if 'nr' in previous_state else zeros3D.copy() #number concentration (#/kg)
        qni_curr = add_halo(to32(previous_state['ni'].reshape(shape3D))) if 'ni' in previous_state else zeros3D.copy()
        qnip_curr = add_halo(to32(previous_state['nip'].reshape(shape3D))) if 'nip' in previous_state else zeros3D.copy() #number concentration (#/kg)
        qnic_curr = add_halo(to32(previous_state['nic'].reshape(shape3D))) if 'nic' in previous_state else zeros3D.copy() #number concentration (#/kg)
        qnid_curr = add_halo(to32(previous_state['nid'].reshape(shape3D))) if 'nid' in previous_state else zeros3D.copy() #number concentration (#/kg)
        qns_curr = add_halo(to32(previous_state['ns'].reshape(shape3D))) if 'ns' in previous_state else zeros3D.copy() #number concentration (#/kg)
        qng_curr = add_halo(to32(previous_state['ng'].reshape(shape3D))) if 'ng' in previous_state else zeros3D.copy() #number concentration (#/kg)
        qnh_curr = add_halo(to32(previous_state['nh'].reshape(shape3D))) if 'ng' in previous_state else zeros3D.copy() #number concentration (#/kg)
        ccn1ft = add_halo(to32(previous_state['ccn1ft'].reshape(shape3D))) if 'ccn1ft' in previous_state else zeros3D.copy() #Free CCN concentration number (#/kg)
        ifn1ft = add_halo(to32(previous_state['ifn1ft'].reshape(shape3D))) if 'ifn1ft' in previous_state else zeros3D.copy() #Free IFN concentration number (#/kg)
        u = add_halo(to32(previous_state['u'].reshape(shape3D))) if 'u' in previous_state else zeros3D.copy() #u wind component
        v = add_halo(to32(previous_state['v'].reshape(shape3D))) if 'v' in previous_state else zeros3D.copy() #v wind component
        w = add_halo(to32(previous_state['w'].reshape(shape3D))) if 'w' in previous_state else zeros3D.copy() #w wind component
        TH_OLD = add_halo(to32(previous_state['th_old'].reshape(shape3D))) if 'th_old' in previous_state else zeros3D.copy() #old potential temperature
        QV_OLD = add_halo(to32(previous_state['qv_old'].reshape(shape3D))) if 'qv_old' in previous_state else zeros3D.copy() #old water vapor mixing ratio
        xland = add_halo(to32(previous_state['xland'].reshape(shape2D))) if 'xland' in previous_state else zeros2D.copy() #1 for land
        qnn_curr = ccn1ft
        qnwfa_curr = ccn1ft
        qnifa_curr = ifn1ft
        qnwfa2d = add_halo(to32(previous_state['nwfa2d'].reshape(shape2D))) if 'nwfa2d' in previous_state else zeros2D.copy()
        if self._mp_physics == self.FULL_KHAIN_LYNN:
            num_scalar = 265
            scalar = add_halo(to32(previous_state['scalar'].reshape(tuple(list(shape3D) + [num_scalar]))))
        elif self._mp_physics == self.FAST_KHAIN_LYNN:
            num_scalar = 133
            scalar = add_halo(to32(previous_state['scalar'].reshape(tuple(list(shape3D) + [num_scalar]))))
        elif self._mp_physics == self.THOMPSONAERO:
            num_scalar = 2
            scalar = numpy.zeros(tuple(list(shape3Dhalo) + [num_scalar]), numpy.float32)
            scalar[..., 0] = qnwfa_curr
            scalar[..., 1] = qnifa_curr
        else:
            num_scalar = 0
            scalar = numpy.zeros(tuple(list(shape3Dhalo) + [num_scalar]), numpy.float32)

        #Thermo
        rho = p / (T * (Rd + qv_curr * Rv))
        pi_phy = (p / P0) ** (Rd / Cpd) #exner function
        th = T / pi_phy #potential temperature

        #Old values
        th_save = th.copy()
        qv_save = qv_curr.copy()

        #Space and time
        if self._mp_physics in [self.ETAMPNEW, self.FER_MP_HIRES, self.FER_MP_HIRES_ADVECT,
                          self.THOMPSONAERO,
                          self.FAST_KHAIN_LYNN, self.FULL_KHAIN_LYNN]:
            dx, dy = self._dx, self._dy
        else:
            dx, dy = 1., 1.
        itimestep = timestep_number
        MPDT = 1. #ETAMPNEW, FER_MP_HIRES and FER_MP_HIRES_ADVECT
        z_at_q = numpy.zeros(shape3Dhalo, numpy.float32) #height of levels, used for init; is-it the same as z used for micro_driver?
        if self._mp_physics in [self.THOMPSONAERO] and \
           (qnwfa_curr.max() == 0. or qnifa_curr.max() == 0.):
            logging.warning("CCN/IN initialization in Thompson will consider that " +
                            "all points are at the surface (because of z_at_q values)")
            if shape3D[2] != 1:
                raise NotImplementedError("If we have several points on the vertical, " +
                                          "we must set a profile for z_at_q and " +
                                          "insure consistency with dz8w")
        z = zeros3D.copy() #height above sea level, used for micro_driver; is-it the same as z_at_q used for init?
        z_at_w = zeros3D.copy() #height above sea level at layer interfaces
        dz8w = zeros3D.copy() #layer thickness (between full?, half?), in hm?
        for k in range(dz8w.shape[1]):
            dz8w[:, k, :] = k - 1 #To have 0 on the first physical level
        ht = zeros2D.copy() #Terrain height

        #Start strategy
        cycling = False
        restart = False
        start_of_simulation = timestep_number == 1

        #External tendencies from parameterized cumulus used
        #by MORR_TWO_MOMENT, NSSL_2MOM, NSSL_2MOMCCN and NSSL_2MOMG
        qrcuten = zeros3D.copy()
        qscuten = zeros3D.copy()
        qicuten = zeros3D.copy()
        qccuten = zeros3D.copy()

        #Ground precipitations (c for cumulative, cv for cumulative over current timestep)
        rainnc = zeros2D.copy()
        rainncv = zeros2D.copy()
        snownc = zeros2D.copy()
        snowncv = zeros2D.copy()
        hailnc = zeros2D.copy()
        hailncv = zeros2D.copy()
        graupelnc = zeros2D.copy()
        graupelncv = zeros2D.copy()

        #Diagnostics
        diagflag = False #to enable diagnostics computation
        do_radar_ref = 0 #1 to compute reflectivity (valid value for THOMPSON, at least)
        has_reqc = 0 #!=0 to enable effective radius computation (valid value for THOMPSON, at least)
        has_reqi = 0 #!=0 to enable effective radius computation (valid value for THOMPSON, at least)
        has_reqs = 0 #!=0 to enable effective radius computation (valid value for THOMPSON, at least)
        re_cloud = zeros3D.copy() #cloud effective radius, OUT array but declared INOUT
        re_ice = zeros3D.copy() #ice effective radius, OUT array but declared INOUT
        re_snow = zeros3D.copy() #snow effective radius, OUT array but declared INOUT
        effr_curr = zeros3D.copy() #cloud effective radius?, OUT array but declared INOUT
        ice_effr_curr = zeros3D.copy() #ice effective radius?, OUT array but declared INOUT
        tot_effr_curr = zeros3D.copy() #unused variable (enters FULL_KHAIN_LYNN but is not used inside)
        qic_effr_curr = zeros3D.copy() #ice effective radius?, OUT array but declared INOUT
        qip_effr_curr = zeros3D.copy() #ice effective radius?, OUT array but declared INOUT
        qid_effr_curr = zeros3D.copy() #ice effective radius?, OUT array but declared INOUT
        tempc = zeros3D.copy() #input temperature in Celsius degrees, OUT array but declared INOUT
        kext_ql = zeros3D.copy() #??, OUT array but declared INOUT
        kext_qs = zeros3D.copy() #??, OUT array but declared INOUT
        kext_qg = zeros3D.copy() #??, OUT array but declared INOUT
        kext_qh = zeros3D.copy() #??, OUT array but declared INOUT
        kext_qa = zeros3D.copy() #??, OUT array but declared INOUT
        kext_qic = zeros3D.copy() #??, OUT array but declared INOUT
        kext_qid = zeros3D.copy() #??, OUT array but declared INOUT
        kext_qip = zeros3D.copy() #??, OUT array but declared INOUT
        kext_ft_qic = zeros3D.copy() #??, OUT array but declared INOUT
        kext_ft_qid = zeros3D.copy() #??, OUT array but declared INOUT
        kext_ft_qip = zeros3D.copy() #??, OUT array but declared INOUT
        kext_ft_qs = zeros3D.copy() #??, OUT array but declared INOUT
        kext_ft_qg = zeros3D.copy() #??, OUT array but declared INOUT

        #Not very usefull variables
        f_qvolg = True #To enable NNSL schemes call
        f_qvolh = True #To enable NNSL schemes call

        #Unused variables: these variable is not used at all in microphysics_driver
        ivgtyp = numpy.zeros(shape2Dhalo, numpy.int32) #variable enters FAST_KHAIN_LYNN and FULL_KHAIN_LYNN but is not used inside
        xlong = zeros2D.copy() #variable enters FULL_KHAIN_LYNN but is not used inside
        xlat = zeros2D.copy() #variable enters FULL_KHAIN_LYNN but is not used inside
        height = zeros3D.copy() #variable enters FULL_KHAIN_LYNN but is not used inside
        fnm = numpy.ones(shape1Dhalo, numpy.float32) #Factors for interpolation at "w" grid (interfaces)
        fnp = numpy.ones(shape1Dhalo, numpy.float32) #Factors for interpolation at "w" grid (interfaces)
        qzr_curr = zeros3D.copy()
        qzi_curr = zeros3D.copy()
        qzs_curr = zeros3D.copy()
        qzg_curr = zeros3D.copy()
        qzh_curr = zeros3D.copy()
        f_effr = True
        f_ice_effr = True
        f_tot_effr = True
        f_qic_effr = True
        f_qip_effr = True
        f_qid_effr = True
        f_qic = True
        f_qip = True
        f_qid = True
        f_qnic = True
        f_qnip = True
        f_qnid = True
        f_qv = True
        f_qr = True
        f_qs = True
        f_qh = True
        f_qnc = True
        f_qnr = True
        f_qni = True
        f_qns = True
        f_qng = True
        f_qnh = True
        f_qnn = True
        f_qzr = True
        f_qzi = True
        f_qzs = True
        f_qzg = True
        f_qzh = True
        f_qnwfa = True
        f_qnifa = True
        f_qir = True
        f_qib = True
        f_qt = True
        f_qrimef = True

        #Currently only some schemes are totally plugged:
        #KESSLERSCHEME, THOMPSON, THOMPSONAERO and FULL_KHAIN_LYNN
        #For these schemes, all necessary variables entering mp_init and/or microphysics_driver are
        #correctly initialized (at least, are supposed to be...)
        #The variable which are not correctly initialized are shared in two categories:
        # * first are the immediately following variables:
        #       the schemes that need these variables are identified,
        #       scheme name is written as comment on the following lines
        #       and a tentative to use one of these schemes will
        #       raise a python error. Please double-check before using this
        #       information.
        # * then, other variables (labeled "Other Variables initialization") follow:
        #   schemes that need them are not identified

        #Variables that must be correctly initialized if using some currently untested schemes
        f_ice_phy = zeros3D.copy() #fraction of ice for FER_MP_HIRES, FER_MP_HIRES_ADVECT, ETAMPNEW, CAMMGMPSCHEME
        f_rain_phy = zeros3D.copy() #fraction of rain for FER_MP_HIRES, FER_MP_HIRES_ADVECT, ETAMPNEW, CAMMGMPSCHEME
        f_rimef_phy = zeros3D.copy() #rimed fraction for FER_MP_HIRES, FER_MP_HIRES_ADVECT, ETAMPNEW, CAMMGMPSCHEME
        lowlyr = numpy.zeros(shape2Dhalo, numpy.int32) #ETAMPNEW, FER_MP_HIRES and FER_MP_HIRES_ADVECT
        mp_restart_state = numpy.zeros((43, ), numpy.float32) #ETAMPNEW, FER_MP_HIRES and FER_MP_HIRES_ADVECT
        tbpvs_state = numpy.zeros((7501, ), numpy.float32) #ETAMPNEW, FER_MP_HIRES and FER_MP_HIRES_ADVECT
        tbpvs0_state = numpy.zeros((7501, ), numpy.float32) #ETAMPNEW, FER_MP_HIRES and FER_MP_HIRES_ADVECT
        ixcldliq = 0 #CAMMGMP
        ixcldice = 0 #CAMMGMP
        ixnumliq = 0 #CAMMGMP
        ixnumice = 0 #CAMMGMP
        nssl_cccn = 0. #NSSL_1MOMLFO, NSSL_1MOM, NSSL_2MOM, NSSL_2MOMG, NSSL_2MOMCCN
        nssl_alphah = 0. #NSSL_1MOMLFO, NSSL_1MOM, NSSL_2MOM, NSSL_2MOMG, NSSL_2MOMCCN
        nssl_alphahl = 0. #NSSL_1MOMLFO, NSSL_1MOM, NSSL_2MOM, NSSL_2MOMG, NSSL_2MOMCCN
        nssl_ipelec = 0 #NSSL_1MOMLFO, NSSL_1MOM, NSSL_2MOM, NSSL_2MOMG, NSSL_2MOMCCN
        nssl_isaund = 0 #NSSL_1MOMLFO, NSSL_1MOM, NSSL_2MOM, NSSL_2MOMG, NSSL_2MOMCCN
        nssl_cnoh = 0. #NSSL_1MOMLFO, NSSL_1MOM, NSSL_2MOM, NSSL_2MOMG, NSSL_2MOMCCN
        nssl_cnohl = 0. #NSSL_1MOMLFO, NSSL_1MOM, NSSL_2MOM, NSSL_2MOMG, NSSL_2MOMCCN
        nssl_cnor = 0. #NSSL_1MOMLFO, NSSL_1MOM, NSSL_2MOM, NSSL_2MOMG, NSSL_2MOMCCN
        nssl_cnos = 0. #NSSL_1MOMLFO, NSSL_1MOM, NSSL_2MOM, NSSL_2MOMG, NSSL_2MOMCCN
        nssl_rho_qh = 0. #NSSL_1MOMLFO, NSSL_1MOM, NSSL_2MOM, NSSL_2MOMG, NSSL_2MOMCCN
        nssl_rho_qhl = 0. #NSSL_1MOMLFO, NSSL_1MOM, NSSL_2MOM, NSSL_2MOMG, NSSL_2MOMCCN
        nssl_rho_qs = 0. #NSSL_1MOMLFO, NSSL_1MOM, NSSL_2MOM, NSSL_2MOMG, NSSL_2MOMCCN
        ccn_conc = 0. #WDM5SCHEME, WDM6SCHEME, LINSCHEME, MORR_TWO_MOMENT
        allowed_to_read = True #WSM3SCHEME, WSM5SCHEME, WSM6SCHEME, ETAMPNEW, FER_MP_HIRES, FER_MP_HIRES_ADVECT
        spec_zone = 0 #WSM5SCHEME
        specified = False #WSM5SCHEME
        channel_switch = False #WSM5SCHEME
        accum_mode, aitken_mode, coarse_mode = 0., 0., 0. #CAMMGMP
        ice2, hail = 0, 0 #GSFCGCESCHEME
        snowh = zeros2D.copy() #CAMMGMP
        qfx = zeros2D.copy() #CAMMGMP (Moisture flux at surface (kg m-2 s-1))
        rliq = numpy.zeros(shape2Dhalo, numpy.float32) #CAMMGMP (Vertically-integrated reserved cloud condensate(m/s))
        f_qndrop = True #LINSCHEME and MORR_TWO_MOMENT
        qrimef_curr = zeros3D.copy() #FER_MP_HIRES
        qt_curr = zeros3D.copy() #ETAMPNEW, FER_MP_HIRES
        lradius = zeros3D.copy() #CAMMGMP
        iradius = zeros3D.copy() #CAMMGMP
        cldfra_old_mp = zeros3D.copy() #CAMMGMP
        cldfra_mp = zeros3D.copy() #CAMMGMP
        cldfra_mp_all = zeros3D.copy() #CAMMGMP
        cldfrai = zeros3D.copy() #CAMMGMP
        cldfral = zeros3D.copy() #CAMMGMP
        cldfra_conv = zeros3D.copy() #CAMMGMP
        wsedl3d = zeros3D.copy() #CAMMGMP
        dlf = zeros3D.copy() #CAMMGMP (Detraining cloud water tendency)
        dlf2 = zeros3D.copy() #CAMMGMP (dq/dt due to export of cloud water into environment by shallow convection (kg/kg/s))
        t_phy = zeros3D.copy() #CAMMGMP (Temperature at the mid points (K))
        p_hyd = zeros3D.copy() #CAMMGMP (Hydrostatic pressure(Pa))
        p8w_hyd = zeros3D.copy() #CAMMGMP (Hydrostatic Pressure at level interface (Pa))
        tke_pbl = zeros3D.copy() #CAMMGMP (Turbulence kinetic energy)
        turbtype3d = zeros3D.copy() #CAMMGMP (Turbulent interface types [ no unit ])
        smaw3d = zeros3D.copy() #CAMMGMP (Normalized Galperin instability function for momentum  ( 0<= <=4.964 and 1 at neutral ) [no units])
        alt = zeros3D.copy() #CAMMGMP (inverse density(m3/kg))
        icwmrsh3d = zeros3D.copy() #CAMMGMP (Shallow cumulus in-cloud water mixing ratio (kg/m2))
        icwmrdp3d = zeros3D.copy() #CAMMGMP (Deep Convection in-cloud water mixing ratio (kg/m2))
        shfrc3d = zeros3D.copy() #CAMMGMP (Shallow cloud fraction)
        cmfmc3d = zeros3D.copy() #CAMMGMP (Deep + Shallow Convective mass flux [ kg /s/m^2 ])
        cmfmc2_3d = zeros3D.copy() #CAMMGMP (Shallow convective mass flux [ kg/s/m^2 ])
        rh_old_mp = zeros3D.copy() #CAMMGMP (Old RH)
        lcd_old_mp = zeros3D.copy() #CAMMGMP (Old liquid cloud fraction)
        qvolg_curr  = zeros3D.copy() #NSSL_1MOM, NSSL_2MOM, NSSL_2MOMG, NSSL_2MOMCCN
        qvolh_curr = zeros3D.copy() #NSSL_2MOM, NSSL_2MOMCCN
        qir_curr = zeros3D.copy() #P3_1CATEGORY, P3_1CATEGORY_NC
        qib_curr = zeros3D.copy() #P3_1CATEGORY, P3_1CATEGORY_NC
        ri_curr = zeros3D.copy() #SBU_YLINSCHEME
        cu_used = 0 #NSSL_2MOM, NSSL_2MOMG, NSSL_2MOMCCN

        ###Other Variables initialization: boolean
        f_qc = True
        f_qi = True
        f_qg = True
        ###Other Variables initialization: Integer 0D
        chem_opt = 0
        progn = 0
        gid = 0 #grid id
        ###Other Variables initialization: Float 0D
        naer = 0.
        ###Other Variables initialization: Float 3D
        t8w = numpy.zeros(shape3Dhalo, numpy.float32)
        p8w = numpy.zeros(shape3Dhalo, numpy.float32)
        cldfra = numpy.zeros(shape3Dhalo, numpy.float32)
        cldfra_old = numpy.zeros(shape3Dhalo, numpy.float32)
        exch_h = numpy.zeros(shape3Dhalo, numpy.float32)
        qndrop_curr = numpy.zeros(shape3Dhalo, numpy.float32)

        #Known limitations, there are other
        messages = ""
        if self._mp_physics in [self.NSSL_2MOM, self.NSSL_2MOMCCN]:
            messages += "We must set config_flags%elec_physics\n"
            messages += "qvolh_curr must be initialized correctly\n"
        if self._mp_physics in [self.NSSL_1MOM, self.NSSL_2MOM, self.NSSL_2MOMG, self.NSSL_2MOMCCN]:
            messages += "qvolg_curr must be initialized correctly\n"
        if self._mp_physics in [self.NSSL_2MOM, self.NSSL_2MOMG, self.NSSL_2MOMCCN]:
            messages += "cu_used must be initialized correctly\n"
        if self._mp_physics in [self.WSM6SCHEME, self.MORR_TWO_MOMENT, self.WDM6SCHEME]:
            messages += "We must set config_flags%hail_opt\n"
        if self._mp_physics in [self.CAMMGMPSCHEME]:
            messages += "We must set config_flags%chem_opt, config_flags%shcu_physics, " + \
                        "config_flags%cu_physics and config_flags%CAM_MP_MAM_cpled\n"
            messages += "accum_mode, aitken_mode and coarse_mode must be initialized correctly\n"
            messages += "snowh, qfx and rliq must be initialized correctly\n"
            messages += "lradius, iradius must be initialized correctly\n"
            messages += "cldfra_old_mp, cldfra_mp, cldfra_mp_all, " + \
                        "cldfrai, cldfral and cldfra_conv must be initialized correctly\n"
            messages += "wsedl3d must be initialized correctly\n"
            messages += "ixcldliq, ixcldice, ixnumliq, ixnumice must be initialized correctly\n"
            messages += "dlf, dlf2, icwmrsh3d, icwmrdp3d, shfrc3d, cmfmc3d, " + \
                        "cmfmc2_3d must be initialized correctly\n"
            messages += "t_phy, p_hyd, p8w_hyd, tke_pbl, turbtype3d, smaw3d, alt, " + \
                        "rh_old_mp, lcd_old_mp must be initialized correctly\n"
        if self._mp_physics in [self.ETAMPNEW]:
            messages += "mp_restart_state,tbpvs_state,tbpvs0_state must be initialized\n"
        if self._mp_physics in [self.ETAMPNEW, self.FER_MP_HIRES, self.FER_MP_HIRES_ADVECT]:
            messages += "MPDT and lowlyr must be initialized correctly\n"
        if self._mp_physics in [self.FER_MP_HIRES, self.FER_MP_HIRES_ADVECT,
                          self.ETAMPNEW, self.CAMMGMPSCHEME]:
            messages += "f_ice_phy, f_rain_phy must be initialized correctly\n"
        if self._mp_physics in [self.NSSL_1MOMLFO, self.NSSL_1MOM, self.NSSL_2MOM,
                          self.NSSL_2MOMG, self.NSSL_2MOMCCN]:
            messages += "nssl_* variables must be initialized correctly\n"
        if self._mp_physics in [self.WDM5SCHEME, self.WDM6SCHEME,
                                self.LINSCHEME, self.MORR_TWO_MOMENT]:
            messages += "ccn_conc must be initialized correctly\n"
        if self._mp_physics == self.WSM5SCHEME:
            messages += "ips, jps, ipe and jpe must be initialized correctly\n"
            messages += "spec_zone, specified and channel_switch must be initialized correctly\n"
        if self._mp_physics in [self.WSM3SCHEME, self.WSM5SCHEME, self.WSM6SCHEME, self.ETAMPNEW,
                          self.FER_MP_HIRES, self.FER_MP_HIRES_ADVECT]:
            messages += "allowed_to_read must be initialized correctly\n"
        if self._mp_physics in [self.LINSCHEME, self.MORR_TWO_MOMENT, self.GSFCGCESCHEME,
                          self.SBU_YLINSCHEME, self.CAMMGMPSCHEME]:
            messages += "z must be initialized correctly " + \
                        "(except if it is only used for the sedimentation)\n"
        if self._mp_physics == self.GSFCGCESCHEME:
            messages += "ice2 and hail must be used to select mode\n"
        if self._mp_physics in [self.LINSCHEME, self.MORR_TWO_MOMENT]:
            messages += "f_qndrop must be initialized correctly\n"
        if self._mp_physics == self.FER_MP_HIRES:
            messages += "qrimef_curr must be initialized correctly\n"
        if self._mp_physics in [self.ETAMPNEW, self.FER_MP_HIRES]:
            messages += "qt_curr must be initialized correctly\n"
        if self._mp_physics in [self.P3_1CATEGORY, self.P3_1CATEGORY_NC]:
            messages += "qir_curr and qib_curr must be initialized correctly\n"
        if self._mp_physics == self.SBU_YLINSCHEME:
            messages += "ri_curr must be initialized correctly\n"
        if messages != "":
            raise NotImplementedError(messages)
        if self._mp_physics not in [self.KESSLERSCHEME, self.THOMPSON,
                              self.THOMPSONAERO, self.FULL_KHAIN_LYNN]:
            raise NotImplementedError("Even if no identified issue has been found, " +
                                      "this scheme has not been tested")

        #Calling Initialization
        result = self._mp_init_py(self._mp_physics, cycling,
                                  rainnc, snownc, graupelnc,
                                  restart,
                                  MPDT, dt, dx, dy, lowlyr, f_ice_phy, f_rain_phy, f_rimef_phy,
                                  mp_restart_state, tbpvs_state, tbpvs0_state, allowed_to_read,
                                  start_of_simulation, ixcldliq, ixcldice, ixnumliq, ixnumice,
                                  nssl_cccn, nssl_alphah, nssl_alphahl, nssl_ipelec, nssl_isaund,
                                  nssl_cnoh, nssl_cnohl, nssl_cnor, nssl_cnos, nssl_rho_qh,
                                  nssl_rho_qhl, nssl_rho_qs, ccn_conc, z_at_q,
                                  qnwfa2d, scalar, num_scalar,
                                  ids, ide, jds, jde, kds, kde, ims, ime,
                                  jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
        (rainnc, snownc, graupelnc, warm_rain,
         adv_moist_cond, lowly, f_ice_phy, f_rain_phy, f_rimef_phy,
         mp_restart_state, tbpvs_state, tbpvs0_state, ccn_conc, qnwfa2d, scalar) = result

        if self._mp_physics == self.THOMPSONAERO:
            qnwfa_curr = scalar[..., 0]
            qnifa_curr = scalar[..., 1]
            if not self._enableCCNsource:
                qnwfa2d = qnwfa2d * 0.

        if self._mp_physics == self.FULL_KHAIN_LYNN:
            if (self._enableCCNinit or self._dx > 7500.) and itimestep == 1:
                logging.warning("CCN initialization of FULL_KHAIN_LYNN will consider that " +
                                "all points are near the surface (because of dz8w values)")
            if self._dx > 7500. and itimestep == 1:
                logging.info("CCN of FULL_KHAIN_LYNN will be reset at each timestep")
            if not self._enableCCNinit:
                itimestep = max(itimestep, 2)

        #Calling Microphysics driver
        result = self._microphysics_driver_py(th, rho, pi_phy, p,
                                              ht, dz8w, p8w, dt,dx,dy,
                                              self._mp_physics, spec_zone,
                                              specified, channel_switch,
                                              warm_rain,
                                              t8w,
                                              chem_opt, progn,
                                              cldfra, cldfra_old, exch_h,
                                              xland,snowh,itimestep,
                                              f_ice_phy,f_rain_phy,f_rimef_phy,
                                              lowlyr, gid,
                                              ids,ide, jds,jde, kds,kde,
                                              ims,ime, jms,jme, kms,kme,
                                              ips,ipe, jps,jpe, kps,kpe,
                                              i_start,i_end,j_start,j_end,kts,kte,
                                              num_tiles, naer,
                                              dlf,dlf2,t_phy,p_hyd,p8w_hyd,tke_pbl,z_at_w,qfx,
                                              rliq,turbtype3d,smaw3d,wsedl3d,cldfra_old_mp,
                                              cldfra_mp,cldfra_mp_all,lradius,iradius,
                                              cldfrai,cldfral,cldfra_conv,
                                              alt,
                                              accum_mode,aitken_mode,coarse_mode,
                                              icwmrsh3d,icwmrdp3d,shfrc3d,cmfmc3d,cmfmc2_3d,
                                              cycling,
                                              fnm,fnp,rh_old_mp,lcd_old_mp,
                                              #chem,
                                              #qme3d,prain3d,nevapr3d,rate1ord_cw2pr_st3d,
                                              #dgnum4D,dgnumwet4D           ,
                                              qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr,
                                              qic_curr,qip_curr,qid_curr,
                                              qnic_curr,qnip_curr,qnid_curr,
                                              qndrop_curr,qni_curr,qh_curr,qnh_curr,
                                              qzr_curr,qzi_curr,qzs_curr,qzg_curr,qzh_curr,
                                              qns_curr,qnr_curr,qng_curr,qnn_curr,qnc_curr,
                                              qnwfa_curr,qnifa_curr,
                                              f_qnwfa,f_qnifa,
                                              qvolg_curr,qvolh_curr,
                                              qir_curr,qib_curr,
                                              effr_curr,ice_effr_curr,tot_effr_curr,
                                              qic_effr_curr,qip_effr_curr,qid_effr_curr            ,
                                              f_qv,f_qc,f_qr,f_qi,f_qs,f_qg,f_qndrop,f_qni,
                                              f_qns,f_qnr,f_qng,f_qnc,f_qnn,f_qh,f_qnh,
                                              f_qzr,f_qzi,f_qzs,f_qzg,f_qzh,
                                              f_qvolg,f_qvolh,
                                              f_qic,f_qip,f_qid,
                                              f_qnic,f_qnip,f_qnid,
                                              f_qir,f_qib,
                                              f_effr,f_ice_effr,f_tot_effr,
                                              f_qic_effr,f_qip_effr,f_qid_effr               ,
                                              cu_used,
                                              qrcuten, qscuten, qicuten, qccuten,
                                              qt_curr,f_qt,
                                              mp_restart_state,tbpvs_state,tbpvs0_state ,
                                              hail,ice2, #ccntype,
                                              u,v,w,z   ,
                                              rainnc,    rainncv,
                                              snownc,    snowncv,
                                              hailnc,    hailncv,
                                              graupelnc, graupelncv,
                                              #rainprod, evapprod,
                                              #qv_b4mp, qc_b4mp, qi_b4mp, qs_b4mp,
                                              qnwfa2d,
                                              ri_curr,
                                              diagflag,   do_radar_ref ,
                                              re_cloud, re_ice, re_snow,
                                              has_reqc, has_reqi, has_reqs,
                                              ccn_conc,
                                              scalar,num_scalar,
                                              kext_ql,kext_qs,kext_qg,
                                              kext_qh,kext_qa,
                                              kext_qic,kext_qid,kext_qip,
                                              kext_ft_qic,kext_ft_qid,kext_ft_qip,
                                              kext_ft_qs,kext_ft_qg,
                                              height,tempc,
                                              TH_OLD ,
                                              QV_OLD,
                                              xlat,xlong,ivgtyp,
                                              qrimef_curr,f_qrimef)
        #Result
        (th, t8w, cldfra, cldfra_old ,exch_h, nsource, qlsink, precr,
        preci, precs, precg, f_ice_phy, f_rain_phy, f_rimef_phy, 
        lowlyr, sr, naer, wsedl3d, cldfra_old_mp, cldfra_mp, cldfra_mp_all,
        lradius, iradius, cldfrai, cldfral, cldfra_conv, rh_old_mp, lcd_old_mp,
        qv_curr, qc_curr,qr_curr, qi_curr, qs_curr, qg_curr, qic_curr, qip_curr, qid_curr,
        qnic_curr, qnip_curr, qnid_curr, qndrop_curr, qni_curr, qh_curr, qnh_curr,
        qzr_curr, qzi_curr, qzs_curr, qzg_curr, qzh_curr,
        qns_curr, qnr_curr, qng_curr, qnn_curr, qnc_curr,
        qnwfa_curr, qnifa_curr, qvolg_curr, qvolh_curr, qir_curr, qib_curr,
        effr_curr, ice_effr_curr, tot_effr_curr, qic_effr_curr, qip_effr_curr, qid_effr_curr,
        qt_curr, mp_restart_state, tbpvs_state, tbpvs0_state, u, v, w, z,
        rainnc, rainncv, snownc, snowncv, hailnc, hailncv, graupelnc, graupelncv,
        refl_10cm, vmi3d, di3d, rhopo3d, ri_curr, re_cloud, re_ice, re_snow,
        scalar, kext_ql, kext_qs, kext_qg, kext_qh, kext_qa,
        kext_qic, kext_qid, kext_qip, kext_ft_qic, kext_ft_qid, kext_ft_qip, kext_ft_qs, kext_ft_qg,
        height, tempc, TH_OLD, QV_OLD, qrimef_curr) = result

        next_state = {}

        #Old values
        to64 = self._to64
        del_halo = self._del_halo
        if 'th_old' in previous_state:
            next_state['th_old'] = to64(del_halo(th_save)).reshape(shapeOri3D)
        if 'qv_old' in previous_state:
            next_state['qv_old'] = to64(del_halo(qv_save)).reshape(shapeOri3D)

        #Output
        next_state['T'] = to64(del_halo(th * pi_phy)).reshape(shapeOri3D)
        next_state['rv'] = to64(del_halo(qv_curr)).reshape(shapeOri3D)
        next_state['rc'] = to64(del_halo(qc_curr)).reshape(shapeOri3D)
        next_state['rr'] = to64(del_halo(qr_curr)).reshape(shapeOri3D)

        if self._mp_physics == self.FULL_KHAIN_LYNN:
            next_state['ric'] = to64(del_halo(qic_curr)).reshape(shapeOri3D)
            next_state['rip'] = to64(del_halo(qip_curr)).reshape(shapeOri3D)
            next_state['rid'] = to64(del_halo(qid_curr)).reshape(shapeOri3D)
            next_state['ri'] = next_state['ric'] + next_state['rip'] + next_state['rid']
            next_state['nip'] = to64(del_halo(qnip_curr)).reshape(shapeOri3D)
            next_state['nic'] = to64(del_halo(qnic_curr)).reshape(shapeOri3D)
            next_state['nid'] = to64(del_halo(qnid_curr)).reshape(shapeOri3D)
            next_state['ni'] = next_state['nip'] + next_state['nic'] + next_state['nid']
            next_state['scalar'] = to64(del_halo(scalar)).reshape(tuple(list(shapeOri3D) + [265]))
            next_state['ccn1ft'] = to64(del_halo(qnn_curr)).reshape(shapeOri3D)
        if self._mp_physics in [self.THOMPSON, self.THOMPSONAERO]:
            next_state['ri'] = to64(del_halo(qi_curr)).reshape(shapeOri3D)
            next_state['ni'] = to64(del_halo(qni_curr)).reshape(shapeOri3D)
        if self._mp_physics == self.THOMPSONAERO:
            next_state['ccn1ft'] = to64(del_halo(qnwfa_curr)).reshape(shapeOri3D)
            next_state['ifn1ft'] = to64(del_halo(qnifa_curr)).reshape(shapeOri3D)
            next_state['nwfa2d'] = to64(del_halo(qnwfa2d)).reshape(shapeOri2D)
        if self._mp_physics in [self.THOMPSON, self.THOMPSONAERO, self.FULL_KHAIN_LYNN]:
            next_state['rs'] = to64(del_halo(qs_curr)).reshape(shapeOri3D)
            next_state['rg'] = to64(del_halo(qg_curr)).reshape(shapeOri3D)
            next_state['rh'] = to64(del_halo(qh_curr)).reshape(shapeOri3D)
            next_state['nc'] = to64(del_halo(qnc_curr)).reshape(shapeOri3D)
            next_state['nr'] = to64(del_halo(qnr_curr)).reshape(shapeOri3D)
            next_state['ns'] = to64(del_halo(qns_curr)).reshape(shapeOri3D)
            next_state['ng'] = to64(del_halo(qng_curr)).reshape(shapeOri3D)
            next_state['nh'] = to64(del_halo(qnh_curr)).reshape(shapeOri3D)

        return next_state
