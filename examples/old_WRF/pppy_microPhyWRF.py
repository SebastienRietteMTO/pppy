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
- hail_opt for WSM6 and WDM6
- ccn0 for WDM5 and WDM6
- gid for FER_MP_HIRES_ADVECT

There are several variables related to ccn/ifn and it is possible that intialization and
usage of these variables is wrong:
- qnn_curr
- ccn0 (ccn_conc in WRF source code)
- naer
- nssl_cccn
- qnwfa2d
- qnwfa_curr
- qnifa_curr
- scalar
- enableProgNc

THOMPSONAERO
* For the initialization, CCN and IN number must be in scalar variable.
* Variable qnwfa2d is set during initialization of THOMPSONAERO to represent a ground tendency
  of CCN. This tendency is, then, applied during integration to refill nwfa variable near ground.
  CCN concentration takes into account a source if option enableCCNsource is True and discards
  this source otherwise.
* After having called self._mp_init_py, if nwfa and/or nifa are null, they are filled with a
  profile computed from height
* A water conservation appears with this scheme (is it an error here? in the scheme?)

THOMPSON
* A water conservation appears with this scheme (is it an error here? in the scheme?)

FULL_KHAIN_LYNN
* The content of each bin is stored in the scalar variable, then this variable is the
  pronostic variable of this scheme for the following iterations.
* Scheme computes derivatives. Remember that, presently, halo is filled with the mean.
* In normal usage, scheme initializes CCN concentration when:
  - this is the first timestep
  - this is not the first timestep and dx is greater than 7500.
  If option enableCCNinit is False, a timestep different from 1 is always given to scheme
  to switch off first case. Second case can be switched off by providing a small dx (scheme option).
  In case we go through initialization, altitude is considered to be of some meters (depending of
  the number of points chosen). To change this behavior, code must be updated to enable the usage
  of a user's defined dz8w variable.
* The scheme needs u, v and w to extrapolate in time tendency of T and qv. We need to provide
  not null dz8w for derivative computation.
  Normally this extrapolation gives zero in a 0D case because all points are identical.
* The scheme crashes for initial condition P=100000., T=290., rv=rc=1.E-2, rr=1.E-4,
  ri=ric=rid=rip=rs=rg=rh=0., nc=3.E8, nr=2000., ni=ns=ng=nh=nic=nid=nip=0.,
  ccn1ft=ccn1at=1.E8, ifn1ft=ifn1at=0., u=v=w=0., xland=1.
  for dt=20s. Iterations of breakup subroutine leads to infinite values.
* A water conservation appears with this scheme (is it an error here? in the scheme?)

WSM6 and WDM6
* The hail_opt option is used internally to define the characteristics of the rimed specie.
  The number of species does not change with this option.

WDM5 and WDM6
* At first timestep, the value of ccn0 (scalar) is used to initialiaze the 3D array qnn_curr

FER_MP_HIRES, FER_MP_HIRES_ADVECT
* These two options are for the same scheme. The 'ADVECT' version share the same set of
  pronostic variables with the other schemes whereas the non 'ADVECT' scheme uses the qt
  variable and several fractions to split this qt content into internal variables.
  It is decided to only implement the 'ADVECT' version in this script
* The gid argument impacts the choice of the critic relative humidity threshold

MORR_TWO_MOMENT
* WRF normally allows to deal with activation with this scheme. The behaviour is
  not enabled here because it would need more intialisation (to be able to call
  the prescribe_aerosol_mixiactivate subroutine). Nonetheless, it is technically
  possible to call the scheme with prognostic Nc, using the enableProgNc option.
* hail_opt selects the set of parameters to use for the graupel category.
"""

import numpy
import logging
import glob
import os
import shutil
import ctypesForFortran
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
                 enableCCNinit=None, hail_opt=None,
                 ccn0=None, gid=None, enableProgNc=None):
        """
        In addition to dt, method, name and tag parameters
        defined in the PPPY class, this parametrisation
        needs the following parameters:
        solib           : path to the shared library to use
        mp_physics      : number of the scheme to run
        enableCCNsource : True or False to enable CCN source term for THOMPSONAERO scheme
        dx, dy          : size of the grid box for THOMPSONAERO and FULL_KHAIN_LYNN schemes
        enableCCNinit   : True or False to compute initial values for CCN for FULL_KHAIN_LYNN
        hail_opt        : for WSM6 and WDM6 schemes: 1 for hail
        ccn0            : initial ccn value for the WDM5 and WDM6 schemes
        gid             : grid id used by FER_MP_HIRES_ADVECT

        This methods also deals with the resource files.
        """
        WRF_options = dict(mp_physics=mp_physics)
        if mp_physics == self.THOMPSONAERO:
            assert enableCCNsource is not None, "enableCCNsource must be defined for THOMPSONAERO"
            WRF_options['enableCCNsource'] = enableCCNsource
        if mp_physics in [self.THOMPSONAERO, self.FULL_KHAIN_LYNN, self.FER_MP_HIRES_ADVECT]:
            #seems to be useless for FER_MP_HIRES_ADVECT but enter the scheme
            assert dx is not None and dy is not None, "dx and dy must be defined for THOMPSONAERO, FULL_KHAIN_LYNN and FER_MP_HIRES_ADVECT"
            WRF_options['dx'] = dx
            WRF_options['dy'] = dy
        if mp_physics == self.FULL_KHAIN_LYNN:
            assert enableCCNinit is not None, "enableCCNinit must be defined for FULL_KHAIN_LYNN"
            WRF_options['enableCCNinit'] = enableCCNinit
        if mp_physics in [self.WSM6SCHEME, self.MORR_TWO_MOMENT, self.WDM6SCHEME]:
            assert hail_opt is not None, "hail_opt must be defined for WSM6SCHEME, MORR_TWO_MOMENT and WDM6SCHEME"
            WRF_options['hail_opt'] = hail_opt
        if mp_physics in [self.WDM5SCHEME, self.WDM6SCHEME]:
            assert ccn0 is not None, "ccn0 must be defined for WDM5SCHEME and WDM6SCHEME"
            WRF_options['ccn0'] = ccn0
        if mp_physics == self.FER_MP_HIRES:
            raise ValueError("It is decided to not implement this version of " + \
                             "scheme and use only the 'ADVECT' version of it.")
        if mp_physics == self.FER_MP_HIRES_ADVECT:
            assert gid is not None, "gid option must be defined for FER_MP_HIRES_ADVECT"
            WRF_options['gid'] = gid
        if mp_physics == self.MORR_TWO_MOMENT:
            assert enableProgNc is not None, "enableProgNc must be defined for MORR_TWO_MOMENT"
            WRF_options['enableProgNc'] = enableProgNc

        super().__init__(dt, method, name, tag,
                         solib=solib, **WRF_options)

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
        #Because if shared lib is load at init, all declared schemes will have
        #their shared lib loaded simultaneously which can be a source of error.
        #This way, shared lib is loaded in a sub-process (if class instance is used
        #through a PPPYComp instance).

        IN = ctypesForFortran.IN
        OUT = ctypesForFortran.OUT
        INOUT = ctypesForFortran.INOUT

        ctypesFF, self._handle = ctypesForFortran.ctypesForFortranFactory(self._options['solib'])

        @ctypesFF(prefix='__module_mp_python_MOD_', suffix='')
        def mp_init_py(mp_physics, cycling, hail_opt,
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
            return ([mp_physics, cycling, hail_opt,
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
                     (numpy.int32, None, IN), #hail_opt
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
            shape3D = (ime-ims+1, kme-kms+1, jme-jms+1) #Shape of a 3D field
            shape2D = (ime-ims+1, jme-jms+1) #Shape of a 2D field
            shape1D = (kme-kms+1, ) #Shape of a 1D field
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

            signature = [(numpy.float32, shape3D, INOUT), #th : potential temperature   
                         (numpy.float32, shape3D, IN), #rho : density of air 
                         (numpy.float32, shape3D, IN), #pi_phy : exner function 
                         (numpy.float32, shape3D, IN), #p : pressure
                         (numpy.float32, shape2D, IN), #ht
                         (numpy.float32, shape3D, IN), #dz8w : dz between full levels (m)
                         (numpy.float32, shape3D, IN), #p8w : pressure at full levels (Pa)
                         (numpy.float32, None, IN), #dt : time step
                         (numpy.float32, None, IN), #dx : space variation
                         (numpy.float32, None, IN), #dy : space variation
                         (numpy.int32, None, IN), #mp_physics : scheme index
                         (numpy.int32, None, IN), #spec_zone : number of points in specified zone
                         (numpy.bool, None, IN), #specified : 
                         (numpy.bool, None, IN), #channel_switch : 
                         (numpy.bool, None, IN), #warm_rain : 
                         (numpy.float32, shape3D, INOUT), #t8w : temperature at layer interfaces
                         (numpy.int32, None, IN), #chem_opt : 
                         (numpy.int32, None, IN), #progn : 
                         (numpy.float32, shape3D, INOUT), # cldfra : current cloud fraction
                         (numpy.float32, shape3D, INOUT), # cldfra_old : previous cloud fraction
                         (numpy.float32, shape3D, INOUT), # exch_h : vertical diffusivity
                         (numpy.float32, shape3D, OUT), #nsource : OUT
                         (numpy.float32, shape3D, OUT), #qlsink : fractional cloud water sink OUT
                         (numpy.float32, shape3D, OUT), #precr : rain precipitation rate at all levels OUT
                         (numpy.float32, shape3D, OUT), #preci : ice precipitation rate at all levels OUT
                         (numpy.float32, shape3D, OUT), #precs : snow precipitation rate at all levels OUT
                         (numpy.float32, shape3D, OUT), #precg: graupel precipitation rate at all levels OUT
                         (numpy.float32, shape2D, IN), #xland :
                         (numpy.float32, shape2D, IN), #snowh : 
                         (numpy.int32, None, IN), #itimestep : if 1 use ice processes in full sbm scheme
                         (numpy.float32, shape3D, INOUT), #f_ice_phy : fraction of ice
                         (numpy.float32, shape3D, INOUT), #f_rain_phy : fraction of rain
                         (numpy.float32, shape3D, INOUT), #f_rimef_phy : mass ratio of rimed ice
                         (numpy.int32, shape2D, INOUT), #lowlyr : 
                         (numpy.float32, shape2D, OUT), #sr : one time step mass ratio of snow to total precipitation OUT
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
                         (numpy.float32, shape3D, IN), #dlf : detraining cloud water tendency
                         (numpy.float32, shape3D, IN), #dlf2 : dq/dt due to export of cloud water into environment by shallow convection
                         (numpy.float32, shape3D, IN), #t_phy : temperature at the mid points
                         (numpy.float32, shape3D, IN), #p_hyd : hydrostatic pressure
                         (numpy.float32, shape3D, IN), #p8w_hyd : hydrostatic pressure at level interface
                         (numpy.float32, shape3D, IN), #tke_pbl : turbulence kinetic energy
                         (numpy.float32, shape3D, IN), #z_at_w : height above sea level at layer interfaces 
                         (numpy.float32, shape2D,IN), #qfx : moisture flux at surface
                         (numpy.float32, shape2D,IN), #rliq : vertically integrated reserved cloud condensate
                         (numpy.float32, shape3D, IN), #turbtype3d : turbulence interface type 
                         (numpy.float32, shape3D, IN), #smaw3d : Normalized Galperin instability function for momentum 
                         (numpy.float32, shape3D, INOUT), #wsedl3d : sedimentation velocity of stratiform liquid cloud dropplets
                         (numpy.float32, shape3D, INOUT), #cldfra_old_mp : old cloud fraction for cammgmp microphysics only
                         (numpy.float32, shape3D, INOUT), #cldfra_mp : old cloud fraction for cammgmp microphysics only
                         (numpy.float32, shape3D, INOUT), #cldfra_mp_all : old cloud fraction for cammgmp microphysics only
                         (numpy.float32, shape3D, INOUT), #lradius : old cloud fraction for cammgmp microphysics only
                         (numpy.float32, shape3D, INOUT), #iradius : old cloud fraction for cammgmp microphysics only
                         (numpy.float32, shape3D, INOUT), #cldfrai : old cloud fraction for cammgmp microphysics only  
                         (numpy.float32, shape3D, INOUT), #cldfral : old cloud fraction for cammgmp microphysics only
                         (numpy.float32, shape3D, INOUT), #cldfra_conv : 
                         (numpy.float32, shape3D, IN), #alt : inverse density
                         (numpy.float32, None, IN), #accum_mode : 
                         (numpy.float32, None, IN), #aitken_mode :
                         (numpy.float32, None, IN), #coarse_mode : 
                         (numpy.float32, shape3D, IN), #icwmrsh3d : shallow cumulus in cloud water mixing ratio
                         (numpy.float32, shape3D, IN), #icwmrdp3d : deep convection in cloud water mixing ratio
                         (numpy.float32, shape3D, IN), #shfrc3d : shallow cloud fraction
                         (numpy.float32, shape3D, IN), #cmfmc3d : deep + shallow convective mass flux
                         (numpy.float32, shape3D, IN), #cmfmc2_3d : shallow convective mass flux
                         (numpy.int32, None, IN), #cycling
                         (numpy.float32, shape1D, IN), #fnm : factors for interpolation at WGRID
                         (numpy.float32, shape1D, IN), #fnp : factors for interpolation at WGRID
                         (numpy.float32, shape3D, INOUT), #rh_old_mp : old rh
                         (numpy.float32, shape3D, INOUT), #lcd_old_mp :  old liquid cloud fraction
                         #(numpy.float32, tuple(list(shape3D) + [1]), INOUT), #chem : scheme array for cammgmp scheme pronostic aerosol (dim : num_chem)
                         #(numpy.float32, shape3D, INOUT), #qme3d : net condensation rate
                         #(numpy.float32, shape3D, INOUT), #prain3d : rate of conversion of condensate to precipitation
                         #(numpy.float32, shape3D, INOUT), #nevapr3d : evaporation rate of rain + snow
                         #(numpy.float32, shape3D, INOUT), #rate1ord_cw2pr_st3d : first order rate for direct conversion strat cloud water to precip
                         #(numpy.float32, tuple(list(shape3D) + [3]), IN), #dgnum4D (dim : ntot_amode_cam_mam)
                         #(numpy.float32, tuple(list(shape3D) + [3]), IN), #dgnumwet4D (dim : ntot_amode_cam_mam)
                         (numpy.float32, shape3D, INOUT), #qv_curr : vapor mixing ratio
                         (numpy.float32, shape3D, INOUT), #qc_curr : cloud mixing ratio
                         (numpy.float32, shape3D, INOUT), #qr_curr : rain mixing ratio
                         (numpy.float32, shape3D, INOUT), #qi_curr : ice mixing ratio
                         (numpy.float32, shape3D, INOUT), #qs_curr : snow mixing ratio
                         (numpy.float32, shape3D, INOUT), #qg_curr : graupel mixing ratio
                         (numpy.float32, shape3D, INOUT), #qic_curr : ice column mixing ratio
                         (numpy.float32, shape3D, INOUT), #qip_curr : ice plates mixing ratio
                         (numpy.float32, shape3D, INOUT), #qid_curr : ice dendrites mixing ratio
                         (numpy.float32, shape3D, INOUT), #qnic_curr : number of concentration of ice column
                         (numpy.float32, shape3D, INOUT), #qnip_curr : number of concentration of ice plates
                         (numpy.float32, shape3D, INOUT), #qnid_curr : number of concentration of ice dendrites
                         (numpy.float32, shape3D, INOUT), #qndrop_curr 
                         (numpy.float32, shape3D, INOUT), #qni_curr : number of concentration of ice 
                         (numpy.float32, shape3D, INOUT), #qh_curr : hail mixing ratio
                         (numpy.float32, shape3D, INOUT), #qnh_curr : number of concentration of hail
                         (numpy.float32, shape3D, INOUT), #qzr_curr
                         (numpy.float32, shape3D, INOUT), #qzi_curr
                         (numpy.float32, shape3D, INOUT), #qzs_curr
                         (numpy.float32, shape3D, INOUT), #qzg_curr
                         (numpy.float32, shape3D, INOUT), #qzh_curr
                         (numpy.float32, shape3D, INOUT), #qns_curr : number of concentration of snow
                         (numpy.float32, shape3D, INOUT), #qnr_curr : number of concentration of rain
                         (numpy.float32, shape3D, INOUT), #qng_curr : number of concentration of graupel
                         (numpy.float32, shape3D, INOUT), #qnn_curr
                         (numpy.float32, shape3D, INOUT), #qnc_curr
                         (numpy.float32, shape3D, INOUT), #qnwfa_curr
                         (numpy.float32, shape3D, INOUT), #qnifa_curr
                         (numpy.bool, None, IN), #f_qnwfa : 
                         (numpy.bool, None, IN), #f_qnifa :
                         (numpy.float32, shape3D, INOUT), #qvolg_curr :
                         (numpy.float32, shape3D, INOUT), #qvolh_curr
                         (numpy.float32, shape3D, INOUT), #qir_curr :
                         (numpy.float32, shape3D, INOUT), #qib_curr :
                         (numpy.float32, shape3D, INOUT), #effr_curr : 
                         (numpy.float32, shape3D, INOUT), #ice_effr_curr :
                         (numpy.float32, shape3D, INOUT), #tot_effr_curr :
                         (numpy.float32, shape3D, INOUT), #qic_effr_curr :
                         (numpy.float32, shape3D, INOUT), #qip_effr_curr :
                         (numpy.float32, shape3D, INOUT), #qid_effr_curr :
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
                         (numpy.float32,shape3D, IN), #qrcuten : 
                         (numpy.float32,shape3D, IN), #qscuten : 
                         (numpy.float32,shape3D, IN), #qicuten :
                         (numpy.float32,shape3D, IN), #qccuten : 
                         (numpy.float32, shape3D, INOUT), #qt_curr :
                         (numpy.bool, None, IN), #f_qt :
                         (numpy.float32, (43, ), INOUT), #mp_restart_state :
                         (numpy.float32, (7501, ), INOUT), #tbpvs_state :
                         (numpy.float32, (7501, ), INOUT), #tbpvs0_state : 
                         (numpy.int32, None, IN), #hail : ccn type
                         (numpy.int32, None, IN), #ice2 : ccn type
                         #(), #ccntype :
                         (numpy.float32, shape3D, INOUT), #u :
                         (numpy.float32, shape3D, INOUT), #v :
                         (numpy.float32, shape3D, INOUT), #w :
                         (numpy.float32, shape3D, INOUT), #z :
                         (numpy.float32, shape2D, INOUT ), #rainnc :
                         (numpy.float32, shape2D, INOUT ), #rainncv :
                         (numpy.float32, shape2D, INOUT ), #snownc :
                         (numpy.float32, shape2D, INOUT ), #snowncv :
                         (numpy.float32, shape2D, INOUT ), #hailnc :
                         (numpy.float32, shape2D, INOUT ), #hailncv :
                         (numpy.float32, shape2D, INOUT ), #graupelnc :
                         (numpy.float32, shape2D, INOUT ), #graupelncv : 
                         #(numpy.float32, shape3D, INOUT), #rainprod :
                         #(numpy.float32, shape3D, INOUT), #evapprod :
                         #(numpy.float32, shape3D, INOUT), #qv_b4mp :
                         #(numpy.float32, shape3D, INOUT), #qc_b4mp
                         #(numpy.float32, shape3D, INOUT), #qi_b4mp :
                         #(numpy.float32, shape3D, INOUT), #qs_b4mp :
                         (numpy.float32, shape2D, IN), #qnwfa2d : 
                         (numpy.float32, shape3D, OUT), #refl_10cm : OUT
                         (numpy.float32, shape3D, OUT), #vmi3d :OUT
                         (numpy.float32, shape3D, OUT), #di3d :OUT
                         (numpy.float32, shape3D, OUT), #rhopo3d : OUT
                         (numpy.float32, shape3D, INOUT), #ri_curr : OUT
                         (numpy.bool, None, IN), #diagflag : logical to tell us when to produce diagnostic for history or restart
                         (numpy.int32, None, IN), #do_radar_ref 
                         (numpy.float32, shape3D, INOUT), #re_cloud : 
                         (numpy.float32, shape3D, INOUT), #re_ice :
                         (numpy.float32, shape3D, INOUT), #re_snow :
                         (numpy.int32, None, IN), #has_reqc : 
                         (numpy.int32, None, IN), #has_reqi :
                         (numpy.int32, None, IN), #has_reqs :
                         (numpy.float32, None, IN), #ccn_conc :
                         (numpy.float32, tuple(list(shape3D) + [num_scalar]), INOUT), #scalar : 
                         (numpy.int32, None, IN), #num_scalar :
                         (numpy.float32, shape3D, INOUT), #kext_ql :
                         (numpy.float32, shape3D, INOUT), #kext_qs :
                         (numpy.float32, shape3D, INOUT), #kext_qg :
                         (numpy.float32, shape3D, INOUT), #kext_qh :
                         (numpy.float32, shape3D, INOUT), #kext_qa :
                         (numpy.float32, shape3D, INOUT), #kext_qic :
                         (numpy.float32, shape3D, INOUT), #kext_qid :
                         (numpy.float32, shape3D, INOUT), #kext_qip :
                         (numpy.float32, shape3D, INOUT), #kext_ft_qic :
                         (numpy.float32, shape3D, INOUT), #kext_ft_qid :
                         (numpy.float32, shape3D, INOUT), #kext_ft_qip :
                         (numpy.float32, shape3D, INOUT), #kext_ft_qs :
                         (numpy.float32, shape3D, INOUT), #kext_ft_qg : 
                         (numpy.float32, shape3D, INOUT), #height :
                         (numpy.float32, shape3D, INOUT), #tempc :
                         (numpy.float32, shape3D, INOUT), #TH_OLD :
                         (numpy.float32, shape3D, INOUT), #QV_OLD :
                         (numpy.float32, shape2D, IN), #xlat :
                         (numpy.float32, shape2D, IN), #xlong : 
                         (numpy.int32, shape2D, IN), #ivgtyp :
                         (numpy.float32,shape3D, INOUT), #qrimef_curr
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
        if self._options['mp_physics'] in [self.THOMPSON, self.THOMPSONAERO]:
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
        mp_physics = self._options['mp_physics']

        #Absolutely required arrays
        needed = ['T', 'P', 'rv', 'rc', 'rr']
        if mp_physics in [self.THOMPSON, self.THOMPSONAERO]:
            needed.extend(['w', 'ri', 'rs', 'rg',
                           'nr', 'ni'])
            if mp_physics == self.THOMPSONAERO:
                needed.extend(['nc', 'ccn1ft', 'ifn1ft'])
        elif mp_physics in [self.FULL_KHAIN_LYNN]:
            needed = ['T', 'P', 'rv', 'u', 'v', 'w', 'xland']
        elif mp_physics == self.WSM6SCHEME:
            needed.extend(['ri', 'rs', 'rg'])
        elif mp_physics == self.WSM5SCHEME:
            needed.extend(['ri', 'rs'])
        elif mp_physics == self.WDM6SCHEME:
            needed.extend(['ccn1ft', 'nc', 'nr', 'ri', 'rs', 'rg'])
        elif mp_physics == self.WDM5SCHEME:
            needed.extend(['ccn1ft', 'nc', 'nr', 'ri', 'rs'])
        elif mp_physics == self.MILBRANDT2MOM:
            needed.extend(['w', 'ri', 'rs', 'rg', 'rh',
                           'nc', 'nr', 'ni', 'ns', 'ng', 'nh'])
        elif mp_physics == self.FER_MP_HIRES_ADVECT:
            needed.append('ri')
        elif mp_physics in [self.KESSLERSCHEME, self.WSM3SCHEME]:
            pass
        elif mp_physics == self.MORR_TWO_MOMENT:
            needed.extend(['ri', 'rs', 'rg'])
            if self._options['enableProgNc']:
                needed.append('nc')
        else:
            raise NotImplementedError("Required arrays have not been defined for this scheme: " + str(mp_physics))
        for var in needed:
            if var not in state:
                raise ValueError(var + " must be in state")
        if all([item in state for item in ['ri', 'ric', 'rip', 'rid']]):
            if not numpy.allclose(state['ri'], state['ric'] + state['rip'] + state['rid']):
                raise ValueError("ri and ric+rip+rid are not consistent")

        #Building missing arrays and fusionning categories
        if mp_physics == self.KESSLERSCHEME:
            if 'ri' in state:
                state['rc'] += state['ri']
                state['ri'] -= state['ri']
            for item in ['rs', 'rh', 'rg']:
                if item  in state:
                    state['rr'] += state[item]
                    state[item] -= state[item]
        elif mp_physics == self.MILBRANDT2MOM:
            pass
        elif mp_physics == self.FER_MP_HIRES_ADVECT:
            for item in ['rs', 'rh', 'rg']:
                if item  in state:
                    state['rr'] += state[item]
                    state[item] -= state[item]
            if 'rimef' not in state:
                state['rimef'] = state['ri'].copy()
        elif mp_physics in [self.WSM6SCHEME, self.WDM6SCHEME, self.MORR_TWO_MOMENT]:
            #only one category for graupel and hail, option hail_opt
            #selects the set of constant to use inside the scheme
            if 'rh' in state:
                state['rg'] += state['rh']
                state['rh'] -= state['rh']
        elif mp_physics in [self.WSM5SCHEME, self.WDM5SCHEME]:
            #Only one icy precipitating category (rs) to represent rs+rg+rh
            for item in ['rh', 'rg']:
                if item  in state:
                    state['rs'] += state[item]
                    state[item] -= state[item]
        elif mp_physics == self.WSM3SCHEME:
            #only one cloud category to represent rc+ri
            if 'ri' in state:
                state['rc'] += state['ri']
                state['ri'] -= state['ri']
            #only one precipitating category to represent rr+rs+rg+rh
            for item in ['rh', 'rg', 'rs']:
                if item  in state:
                    state['rr'] += state[item]
                    state[item] -= state[item]
        elif mp_physics in [self.THOMPSON, self.THOMPSONAERO]:
            if 'rh' in state:
                state['rg'] += state['rh']
                state['rh'] -= state['rh']
            if mp_physics == self.THOMPSONAERO and 'nwfa2d' not in state:
                #Values will be set during mp_init
                state['nwfa2d'] = numpy.zeros(tuple(list(state['T'].shape)))
        elif mp_physics == self.FULL_KHAIN_LYNN:
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
                for var in ['rc', 'rr', 'ri', 'rs', 'rg', 'rh', 'ccn1ft']:
                   if var not in state:
                       raise ValueError(var + " must be in state (when scalar is missing)")
                logging.warning("Initialization of scalar must be done according to a " + \
                                "given distribution")
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
        #return array with a 1-point halo on the horizontal
        #and a minimum of 2 values on the vertical
        #method move the position of the k axis
        if len(array.shape) >= 3:
            #input shape is (i, j, k [, n]), output shape is (i + 2, max(k, 2), j + 2 [, n])
            new_shape = tuple([array.shape[0] + 2,
                               max(array.shape[2], 2),
                               array.shape[1] + 2] + list(array.shape[3:]))
            result = numpy.zeros(new_shape, dtype=array.dtype)
            result[...] = array.mean() #to have physical values in the halo
            result[1:-1, :array.shape[2], 1:-1, ...] = numpy.swapaxes(array, 1, 2)
            if array.shape[2] == 1:
                result[:, 1, :] = result[:, 0, :]
        elif len(array.shape) == 2:
            #input shape is (i, j), output shape is (i + 2, j+ 2)
            new_shape = (array.shape[0] + 2, array.shape[1] + 2)
            result = numpy.zeros(new_shape, dtype=array.dtype)
            result[...] = array.mean() #to have physical values in the halo
            result[1:-1, 1:-1] = array
        elif len(array.shape) == 1:
            #input shape is (k, ), output shape is (max(k, 2), )
            new_shape = (max(array.shape[0], 2), )
            result = numpy.zeros(new_shape, dtype=array.dtype)
            result[:array.shape[0]] = array
            if array.shape[0] == 1:
                result[1] = result[0]
        return result

    @staticmethod
    def _del_halo(array, halo_on_vert):
        "Removes halo to array"
        #Three first dimensions are geographical dimensions
        #return array without the 1-point halo
        #and with vertical axes moved to third position
        if len(array.shape) >= 3:
            result = numpy.swapaxes(array[1:-1, :, 1:-1, ...], 1, 2)
            return result[:, :, :-1, ...] if halo_on_vert else result
        elif len(array.shape) == 2:
            return array[1:-1, 1:-1]
        elif len(array.shape) == 1:
            return array[:-1] if halo_on_vert else array

    def execute(self, previous_state, dt, timestep_number):
        super().execute(previous_state, dt, timestep_number)
        mp_physics = self._options['mp_physics']

        #To be compatible with other microphysical pppy already coded
        #it is decided to hypothesised that incoming variables are
        #64bits and axes are in the (i,j,k) order whereas variables
        #intering the WRF driver are 32bits and axes is are in the (i,k,j) order

        #Dimensions
        #shapeOri3D: shape of an input state variables which can be 3D in (i, j, k) order
        #shapeOri2D: shape of an input state variables which can be, at most, 2D
        #shape3D: 3D shape of input state variables in (i, j, k) order
        #shape2D: 2D shape of input state variables
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
        shape2D = shape3D[0:2]

        #Mecanism around domain/memory/tile is not entirely understood
        #To be sure, a 1-point halo is added, on the horizontal, to build
        #the memory array around the input array and another 1-point
        #halo is added around the memory array to get the domain array.
        #On the vertical, 1-point is added on the top if the original
        #array contains only one point.
        #The size of the halo is hard coded here and in _add_halo and _del_halo
        #Array entering the WRF subroutines are defined on the memory shape
        #shape3DHalo: 3D shape of memory array in (i, k, j) order
        #shape2DHalo: 2D shape of memory array
        #shape1Dhalo: 1D (vertical) shape of memory array
        ids, jds, kds = 1, 1, 1 # start domain
        ide, jde, kde = shape3D[0] + 4, shape3D[1] + 4, max(shape3D[2], 2) # end domain
        ims, jms, kms = 2, 2, 1 # start memory
        ime, jme, kme = shape3D[0] + 3, shape3D[1] + 3, max(shape3D[2], 2) # end memory
        its, jts, kts = 3, 3, 1 # start tile
        ite, jte, kte = shape3D[0] + 2, shape3D[1] + 2, max(shape3D[2], 2) # end tile
        shape3Dhalo = (ime - ims + 1, kme - kms + 1, jme - jms + 1)
        shape2Dhalo = (ime - ims + 1, jme - jms + 1)
        shape1Dhalo = (kme - kms + 1, )

        spec_zone = 0 #is somewhat related to dimensions
        specified = False #is somewhat related to dimensions
        channel_switch = False #is somewhat related to dimensions
        gid = self._options.get('gid', 0) #grid id, relevant for FER_MP_HIRES_ADVECT schemes
        ips, jps, kps = its, jts, kts # start ?
        ipe, jpe, kpe = ite, jte, kte # end ?
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

        #Previous state values, pronostic arrays
        add_halo = self._add_halo #adds halo and invert axes
        to32 = self._to32
        zeros2D = numpy.zeros(shape2Dhalo, numpy.float32)
        zeros3D = numpy.zeros(shape3Dhalo, numpy.float32)
        T = add_halo(to32(previous_state['T'].reshape(shape3D))) #Temperature (K)
        p = add_halo(to32(previous_state['P'].reshape(shape3D))) #pressupre (Pa)
        p8w = p #p8w is normally the pressure at interfaces
        qv_curr = add_halo(to32(previous_state['rv'].reshape(shape3D))) #mixing ratio (kg/kg)
        qc_curr = add_halo(to32(previous_state['rc'].reshape(shape3D))) #mixing ratio (kg/kg)
        qr_curr = add_halo(to32(previous_state['rr'].reshape(shape3D))) #mixing ratio (kg/kg)
        qi_curr = add_halo(to32(previous_state['ri'].reshape(shape3D)))  if 'ri' in previous_state else zeros3D.copy()#mixing ratio (kg/kg)
        qrimef_curr = add_halo(to32(previous_state['rimef'].reshape(shape3D)))  if 'rimef' in previous_state else qi_curr #for FER_MP_HIRES_ADVECT
        qic_curr = add_halo(to32(previous_state['ric'].reshape(shape3D))) if 'ric' in previous_state else zeros3D.copy() #mixing ratio (kg/kg)
        qip_curr = add_halo(to32(previous_state['rip'].reshape(shape3D))) if 'rip' in previous_state else  zeros3D.copy() #mixing ratio (kg/kg)
        qid_curr = add_halo(to32(previous_state['rid'].reshape(shape3D))) if 'rid' in previous_state else zeros3D.copy() #mixing ratio (kg/kg)
        qs_curr = add_halo(to32(previous_state['rs'].reshape(shape3D))) if 'rs' in previous_state else zeros3D.copy() #mixing ratio (kg/kg)
        qg_curr = add_halo(to32(previous_state['rg'].reshape(shape3D))) if 'rg' in previous_state else zeros3D.copy() #mixing ratio (kg/kg)
        qh_curr = add_halo(to32(previous_state['rh'].reshape(shape3D))) if 'rh' in previous_state else zeros3D.copy() #mixing ratio (kg/kg)
        qnc_curr = add_halo(to32(previous_state['nc'].reshape(shape3D))) if 'nc' in previous_state else zeros3D.copy() #number concentration (#/kg)
        qndrop_curr = qnc_curr
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
        QV_OLD = add_halo(to32(previous_state['qv_old'].reshape(shape3D))) if 'qv_old' in previous_state else qv_curr #old water vapor mixing ratio
        xland = add_halo(to32(previous_state['xland'].reshape(shape2D))) if 'xland' in previous_state else zeros2D.copy() #1 for land
        qnn_curr = ccn1ft
        qnwfa_curr = ccn1ft
        qnifa_curr = ifn1ft
        qnwfa2d = add_halo(to32(previous_state['nwfa2d'].reshape(shape2D))) if 'nwfa2d' in previous_state else zeros2D.copy()
        if mp_physics == self.FULL_KHAIN_LYNN:
            num_scalar = 265
            scalar = add_halo(to32(previous_state['scalar'].reshape(tuple(list(shape3D) + [num_scalar]))))
        elif mp_physics == self.FAST_KHAIN_LYNN:
            num_scalar = 133
            scalar = add_halo(to32(previous_state['scalar'].reshape(tuple(list(shape3D) + [num_scalar]))))
        elif mp_physics == self.THOMPSONAERO:
            num_scalar = 2
            scalar = numpy.zeros(tuple(list(shape3Dhalo) + [num_scalar]), numpy.float32)
            scalar[..., 0] = qnwfa_curr
            scalar[..., 1] = qnifa_curr
        else:
            num_scalar = 0
            scalar = numpy.zeros(tuple(list(shape3Dhalo) + [num_scalar]), numpy.float32)
        ccn_conc = self._options.get('ccn0', 0)
        naer = ccn_conc #LINSCHEME, MORR_TWO_MOMENT
        nssl_cccn = 1.225 * ccn_conc #NSSL_1MOMLFO, NSSL_1MOM, NSSL_2MOM, NSSL_2MOMG, NSSL_2MOMCCN

        #Thermo
        rho = p / (T * (Rd + qv_curr * Rv))
        pi_phy = (p / P0) ** (Rd / Cpd) #exner function
        th = T / pi_phy #potential temperature
        TH_OLD = add_halo(to32(previous_state['th_old'].reshape(shape3D))) if 'th_old' in previous_state else th #old potential temperature

        #Old values
        th_save = th.copy()
        qv_save = qv_curr.copy()

        #Space and time
        dx = self._options.get('dx', 1.)
        dy = self._options.get('dy', 1.)
        itimestep = timestep_number
        allowed_to_read = timestep_number == 1 #ETAMPNEW, FER_MP_HIRES_ADVECT (enters WSM3SCHEME, WSM5SCHEME but unused)
        MPDT = 0. #time step in minutes of the microphysics (ETAMPNEW and FER_MP_HIRES_ADVECT)
                  #a zero value forces to use the general time step
        z_at_q = numpy.zeros(shape3Dhalo, numpy.float32) #height of levels, used for init; is-it the same as z used for micro_driver?
        if mp_physics in [self.THOMPSONAERO] and \
           (qnwfa_curr.max() == 0. or qnifa_curr.max() == 0.):
            logging.warning("CCN/IN initialization in Thompson will consider that " +
                            "all points are at the surface (because of z_at_q values)")
            if shape3D[2] != 1:
                raise NotImplementedError("If we have several points on the vertical, " +
                                          "we must set a profile for z_at_q and " +
                                          "insure consistency with dz8w")
        z = zeros3D.copy() #height above sea level, used for micro_driver; is-it the same as z_at_q used for init?
        z_at_w = zeros3D.copy() #height above sea level at layer interfaces
        dz8w = zeros3D.copy() + 100. #layer thickness (between full?, half?)
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
        lowlyr = numpy.zeros(shape2Dhalo, numpy.int32) #ETAMPNEW and FER_MP_HIRES_ADVECT (inout but overwritten in init)

        #Unused variables: these variable are not used at all in microphysics_driver
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

        #Currently only some schemes are totally plugged
        #For these schemes, all necessary variables entering mp_init and/or microphysics_driver are
        #correctly initialized (at least, are supposed to be...)
        #The variable which are not correctly initialized receive default values here
        #Each initisation is followed by a comment indicatint the schemes which need the variable
        #and a tentative to use one of these schemes will raise a python error.
        #Please double-check before using this information.

        #Variables that must be correctly initialized if using some currently untested schemes
        f_ice_phy = zeros3D.copy() #fraction of ice for ETAMPNEW, CAMMGMPSCHEME (enters but is not used by FER_MP_HIRES_ADVECT)
        f_rain_phy = zeros3D.copy() #fraction of rain for ETAMPNEW, CAMMGMPSCHEME (enters but is not used by FER_MP_HIRES_ADVECT)
        f_rimef_phy = zeros3D.copy() + 1. #rimed fraction for ETAMPNEW, CAMMGMPSCHEME (enters but is not used by FER_MP_HIRES_ADVECT)
        mp_restart_state = numpy.zeros((43, ), numpy.float32) #ETAMPNEW and FER_MP_HIRES_ADVECT
        tbpvs_state = numpy.zeros((7501, ), numpy.float32) #ETAMPNEW and FER_MP_HIRES_ADVECT
        tbpvs0_state = numpy.zeros((7501, ), numpy.float32) #ETAMPNEW and FER_MP_HIRES_ADVECT
        ixcldliq = 0 #CAMMGMP
        ixcldice = 0 #CAMMGMP
        ixnumliq = 0 #CAMMGMP
        ixnumice = 0 #CAMMGMP
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
        accum_mode, aitken_mode, coarse_mode = 0., 0., 0. #CAMMGMP
        ice2, hail = 0, 0 #GSFCGCESCHEME
        snowh = zeros2D.copy() #CAMMGMP
        qfx = zeros2D.copy() #CAMMGMP (Moisture flux at surface (kg m-2 s-1))
        rliq = numpy.zeros(shape2Dhalo, numpy.float32) #CAMMGMP (Vertically-integrated reserved cloud condensate(m/s))
        qt_curr = zeros3D.copy() #ETAMPNEW
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
        f_qg = True #GSFCGCESCHEME, LINSCHEME
        #The following variables are related to orognostic ccn
        f_qc = True #LINSSHEME and MORR_TWO_MOMENT under some circonstances (call to prescribe_aerosol_mixactivate controlled by chem_opt and progn)
        f_qi = True #LINSCHEME and MORR_TWO_MOMENT under some circonstances (call to prescribe_aerosol_mixactivate controlled by chem_opt and progn)
        cldfra = numpy.zeros(shape3Dhalo, numpy.float32) #LINSSHEME and MORR_TWO_MOMENT under some circonstances (call to prescribe_aerosol_mixactivate controlled by chem_opt and progn)
        cldfra_old = numpy.zeros(shape3Dhalo, numpy.float32) #LINSSHEME and MORR_TWO_MOMENT under some circonstances (call to prescribe_aerosol_mixactivate controlled by chem_opt and progn)
        chem_opt = 0 #CAMMGMP, LINSCHEME, MORR_TWO_MOMENT, NSSL_2MOM, NSSL_2MOMCCN, NSSL_2MOMG
        progn = 0 #CAMMGMP, LINSCHEME, MORR_TWO_MOMENT, NSSL_2MOM, NSSL_2MOMCCN, NSSL_2MOMG
        exch_h = numpy.zeros(shape3Dhalo, numpy.float32) #LINSSHEME and MORR_TWO_MOMENT under some circonstances (call to prescribe_aerosol_mixactivate controlled by chem_opt and progn) and CAMMGMP
        t8w = numpy.zeros(shape3Dhalo, numpy.float32) #LINSSHEME and MORR_TWO_MOMENT under some circonstances (call to prescribe_aerosol_mixactivate controlled by chem_opt and progn)

        if mp_physics not in [self.KESSLERSCHEME, self.THOMPSON,
                              self.THOMPSONAERO, self.FULL_KHAIN_LYNN,
                              self.WSM3SCHEME, self.WSM5SCHEME, self.WSM6SCHEME,
                              self.WDM5SCHEME, self.WDM6SCHEME,
                              self.MILBRANDT2MOM, self.MORR_TWO_MOMENT,
                              self.FER_MP_HIRES_ADVECT]:
            raise NotImplementedError("Even if no identified issue has been found, " +
                                      "this scheme has not been tested")

        #Calling Initialization
        result = self._mp_init_py(mp_physics, cycling, self._options.get('hail_opt', 0),
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
         adv_moist_cond, lowlyr, f_ice_phy, f_rain_phy, f_rimef_phy,
         mp_restart_state, tbpvs_state, tbpvs0_state, ccn_conc, qnwfa2d, scalar) = result

        if mp_physics == self.THOMPSONAERO:
            qnwfa_curr = scalar[..., 0]
            qnifa_curr = scalar[..., 1]
            if not self._options['enableCCNsource']:
                qnwfa2d = qnwfa2d * 0.

        if mp_physics == self.FULL_KHAIN_LYNN:
            if (self._options['enableCCNinit'] or self._options['dx'] > 7500.) and itimestep == 1:
                logging.warning("CCN initialization of FULL_KHAIN_LYNN will consider that " +
                                "all points are near the surface (because of dz8w values)")
            if self._options['dx'] > 7500. and itimestep == 1:
                logging.info("CCN of FULL_KHAIN_LYNN will be reset at each timestep")
            if not self._options['enableCCNinit']:
                itimestep = max(itimestep, 2)

        #Calling Microphysics driver
        if chem_opt in [0, 401] and progn == 1 and mp_physics in [self.LINSCHEME, self.MORR_TWO_MOMENT]:
            #At least, f_qc, f_qi, exch_h, cldfra, clfra_old and t8w must be initialised correctly
            raise NotImplementedError("In this case (chem_opt, progn and mp_physics combination), " + \
                                      "prescribe_aerosol_mixactivate will be called. Code is not ready for that.")
        f_qndrop = self._options.get('enableProgNc', False)
        result = self._microphysics_driver_py(th, rho, pi_phy, p,
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
        del_halo = lambda a: self._del_halo(a,  shape3Dhalo[1] != shape3D[2]) #revert axes and del halo
        if 'th_old' in previous_state:
            next_state['th_old'] = to64(del_halo(th_save)).reshape(shapeOri3D)
        if 'qv_old' in previous_state:
            next_state['qv_old'] = to64(del_halo(qv_save)).reshape(shapeOri3D)

        #Output
        next_state['T'] = to64(del_halo(th * pi_phy)).reshape(shapeOri3D)
        next_state['rv'] = to64(del_halo(qv_curr)).reshape(shapeOri3D)
        next_state['rc'] = to64(del_halo(qc_curr)).reshape(shapeOri3D)
        next_state['rr'] = to64(del_halo(qr_curr)).reshape(shapeOri3D)

        if mp_physics == self.FULL_KHAIN_LYNN:
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
            next_state['rs'] = to64(del_halo(qs_curr)).reshape(shapeOri3D)
            next_state['rg'] = to64(del_halo(qg_curr)).reshape(shapeOri3D)
            next_state['rh'] = to64(del_halo(qh_curr)).reshape(shapeOri3D)
            next_state['nc'] = to64(del_halo(qnc_curr)).reshape(shapeOri3D)
            next_state['nr'] = to64(del_halo(qnr_curr)).reshape(shapeOri3D)
            next_state['ns'] = to64(del_halo(qns_curr)).reshape(shapeOri3D)
            next_state['ng'] = to64(del_halo(qng_curr)).reshape(shapeOri3D)
            next_state['nh'] = to64(del_halo(qnh_curr)).reshape(shapeOri3D)
        elif mp_physics in [self.THOMPSON, self.THOMPSONAERO]:
            next_state['rs'] = to64(del_halo(qs_curr)).reshape(shapeOri3D)
            next_state['rg'] = to64(del_halo(qg_curr)).reshape(shapeOri3D)
            next_state['ri'] = to64(del_halo(qi_curr)).reshape(shapeOri3D)
            next_state['ni'] = to64(del_halo(qni_curr)).reshape(shapeOri3D)
            next_state['nr'] = to64(del_halo(qnr_curr)).reshape(shapeOri3D)
            if mp_physics == self.THOMPSONAERO:
                next_state['ccn1ft'] = to64(del_halo(qnwfa_curr)).reshape(shapeOri3D)
                next_state['ifn1ft'] = to64(del_halo(qnifa_curr)).reshape(shapeOri3D)
                next_state['nwfa2d'] = to64(del_halo(qnwfa2d)).reshape(shapeOri2D)
                next_state['nc'] = to64(del_halo(qnc_curr)).reshape(shapeOri3D)
        elif mp_physics in [self.WSM6SCHEME, self.WDM6SCHEME]:
            next_state['ri'] = to64(del_halo(qi_curr)).reshape(shapeOri3D)
            next_state['rs'] = to64(del_halo(qs_curr)).reshape(shapeOri3D)
            next_state['rg'] = to64(del_halo(qg_curr)).reshape(shapeOri3D)
            if mp_physics == self.WDM6SCHEME:
                next_state['ccn1ft'] = to64(del_halo(qnn_curr)).reshape(shapeOri3D)
                next_state['nc'] = to64(del_halo(qnc_curr)).reshape(shapeOri3D)
                next_state['nr'] = to64(del_halo(qnr_curr)).reshape(shapeOri3D)
        elif mp_physics in [self.WSM5SCHEME, self.WDM5SCHEME]:
            next_state['ri'] = to64(del_halo(qi_curr)).reshape(shapeOri3D)
            next_state['rs'] = to64(del_halo(qs_curr)).reshape(shapeOri3D)
            if mp_physics == self.WDM5SCHEME:
                next_state['ccn1ft'] = to64(del_halo(qnn_curr)).reshape(shapeOri3D)
                next_state['nc'] = to64(del_halo(qnc_curr)).reshape(shapeOri3D)
                next_state['nr'] = to64(del_halo(qnr_curr)).reshape(shapeOri3D)
        elif mp_physics == self.MILBRANDT2MOM:
            next_state['ri'] = to64(del_halo(qi_curr)).reshape(shapeOri3D)
            next_state['rs'] = to64(del_halo(qs_curr)).reshape(shapeOri3D)
            next_state['rg'] = to64(del_halo(qg_curr)).reshape(shapeOri3D)
            next_state['rh'] = to64(del_halo(qh_curr)).reshape(shapeOri3D)
            next_state['nc'] = to64(del_halo(qnc_curr)).reshape(shapeOri3D)
            next_state['nr'] = to64(del_halo(qnr_curr)).reshape(shapeOri3D)
            next_state['ni'] = to64(del_halo(qni_curr)).reshape(shapeOri3D)
            next_state['ns'] = to64(del_halo(qns_curr)).reshape(shapeOri3D)
            next_state['ng'] = to64(del_halo(qng_curr)).reshape(shapeOri3D)
            next_state['nh'] = to64(del_halo(qnh_curr)).reshape(shapeOri3D)
        elif mp_physics == self.FER_MP_HIRES_ADVECT:
            next_state['ri'] = to64(del_halo(qi_curr)).reshape(shapeOri3D)
            next_state['rimef'] = to64(del_halo(qrimef_curr)).reshape(shapeOri3D)
        elif mp_physics in [self.KESSLERSCHEME, self.WSM3SCHEME]:
            pass
        elif mp_physics == self.MORR_TWO_MOMENT:
            next_state['ri'] = to64(del_halo(qi_curr)).reshape(shapeOri3D)
            next_state['rs'] = to64(del_halo(qs_curr)).reshape(shapeOri3D)
            next_state['rg'] = to64(del_halo(qg_curr)).reshape(shapeOri3D)
            next_state['nr'] = to64(del_halo(qnr_curr)).reshape(shapeOri3D)
            next_state['ni'] = to64(del_halo(qni_curr)).reshape(shapeOri3D)
            next_state['ns'] = to64(del_halo(qns_curr)).reshape(shapeOri3D)
            next_state['ng'] = to64(del_halo(qng_curr)).reshape(shapeOri3D)
            if self._options['enableProgNc']:
                next_state['nc'] = to64(del_halo(qndrop_curr)).reshape(shapeOri3D)
        else:
            raise NotImplementedError("Output arrays are not saved for this scheme: " + str(mp_physics))

        return next_state
