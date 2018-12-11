#!/usr/bin/env python 3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

#Written by B. Vié
#Modified by S. Riette

"""
This module contains an implementation of PPPY suitable to call microphysical
LIMA scheme from the Meso-NH model.
"""

import numpy
import os
from pppy import ctypesForFortran
import pppy

class pppy_microPhyLima_MNH54(pppy.PPPY):
    """
    PPPY implementation for calling the microphysical LIMA scheme of the Meso-NH model.
    """

    def __init__(self, dt, method, name, tag,
                 solib, LDWARM=True, LNUCL=True, LCOLD=True, LSNOW=True,
                 LHAIL=False, LHHONI=False, LMEYERS=False, NMOD_IFN=1,
                 CIFN1=1000., CIFN2=0., IFN_HOM=True, IFN_SPECIES='MOCAGE',
                 INT_MIXING='DM1', NMOD_IMM=0, IND_SPECIE=1, CPRISTINE_ICE_LIMA='PLAT',
                 CHEVRIMED_ICE_LIMA='GRAU', XALPHAI=0., XNUI=3., XALPHAS=1., XNUS=1.,
                 XALPHAG=1., XNUG=1., XFACTNUC_DEP=1., XFACTNUC_CON=1., NPHILLIPS=8,
                 LWARM=True, LACTI=True, LRAIN=True, LSEDC=False, LACTIT=False,
                 NMOD_CCN=1, CCCN1=350., CCCN2=0., CCCN3=0., CCCN4=0., CCN_HOM=True,
                 CCN_MODES='MOCAGE', HINI_CCN='AER', HTYPE_CCN='C', XALPHAC=3.,
                 XNUC=1., XALPHAR=1., XNUR=2., XFSOLUB_CCN=1., XACTEMP_CCN=280.,
                 XAERDIFF=0., XAERHEIGHT=2000., LSCAV=False, LAERO_MASS=False):
        """
        In addition to dt, method, name and tag parameters
        defined in the PPPY class, this parameterization
        needs or accept the following parameters:
        solib              : path to the shared lib for the scheme
        LDWARM             : (optional)
        LNUCL              : (optional)
        LCOLD              : (optional)
        LSNOW              : (optional)
        LHAIL              : (optional)
        LHHONI             : (optional)
        LMEYERS            : (optional)
        NMOD_IFN           : (optional)
        CIFN1              : (optional)
        CIFN2              : (optional)
        IFN_HOM            : (optional)
        IFN_SPECIES        : (optional)
        INT_MIXING         : (optional)
        NMOD_IMM           : (optional)
        IND_SPECIE         : (optional)
        CPRISTINE_ICE_LIMA : (optional)
        CHEVRIMED_ICE_LIMA : (optional)
        XALPHAI            : (optional)
        XNUI               : (optional)
        XALPHAS            : (optional)
        XNUS               : (optional)
        XALPHAG            : (optional)
        XNUG               : (optional)
        XFACTNUC_DEP       : (optional)
        XFACTNUC_CON       : (optional)
        NPHILLIPS          : (optional)
        LWARM              : (optional)
        LACTI              : (optional)
        LRAIN              : (optional)
        LSEDC              : (optional)
        LACTIT             : (optional)
        NMOD_CCN           : (optional)
        CCCN1              : (optional)
        CCCN2              : (optional)
        CCCN3              : (optional)
        CCCN4              : (optional)
        CCN_HOM            : (optional)
        CCN_MODES          : (optional)
        HINI_CCN           : (optional)
        HTYPE_CCN          : (optional)
        XALPHAC            : (optional)
        XNUC               : (optional)
        XALPHAR            : (optional)
        XNUR               : (optional)
        XFSOLUB_CCN        : (optional)
        XACTEMP_CCN        : (optional)
        XAERDIFF           : (optional)
        XAERHEIGHT         : (optional)
        LSCAV              : (optional)
        LAERO_MASS         : (optional)
        """

        lima_options = dict(LDWARM=LDWARM, LNUCL=LNUCL, LCOLD=LCOLD, LSNOW=LSNOW,
                            LHAIL=LHAIL, LHHONI=LHHONI, LMEYERS=LMEYERS, NMOD_IFN=NMOD_IFN,
                            CIFN1=CIFN1, CIFN2=CIFN2, IFN_HOM=IFN_HOM, IFN_SPECIES=IFN_SPECIES,
                            INT_MIXING=INT_MIXING, NMOD_IMM=NMOD_IMM, IND_SPECIE=IND_SPECIE,
                            CPRISTINE_ICE_LIMA=CPRISTINE_ICE_LIMA,
                            CHEVRIMED_ICE_LIMA=CHEVRIMED_ICE_LIMA,
                            XALPHAI=XALPHAI, XNUI=XNUI, XALPHAS=XALPHAS, XNUS=XNUS,
                            XALPHAG=XALPHAG, XNUG=XNUG, XFACTNUC_DEP=XFACTNUC_DEP,
                            XFACTNUC_CON=XFACTNUC_CON, NPHILLIPS=NPHILLIPS,
                            LWARM=LWARM, LACTI=LACTI, LRAIN=LRAIN, LSEDC=LSEDC, LACTIT=LACTIT,
                            NMOD_CCN=NMOD_CCN, CCCN1=CCCN1, CCCN2=CCCN2,
                            CCCN3=CCCN3, CCCN4=CCCN4, CCN_HOM=CCN_HOM,
                            CCN_MODES=CCN_MODES, HINI_CCN=HINI_CCN,
                            HTYPE_CCN=HTYPE_CCN, XALPHAC=XALPHAC,
                            XNUC=XNUC, XALPHAR=XALPHAR, XNUR=XNUR,
                            XFSOLUB_CCN=XFSOLUB_CCN, XACTEMP_CCN=XACTEMP_CCN,
                            XAERDIFF=XAERDIFF, XAERHEIGHT=XAERHEIGHT,
                            LSCAV=LSCAV, LAERO_MASS=LAERO_MASS)
        self._lima_options = lima_options
        self._solib = solib
        super().__init__(dt, method, name, tag, solib=solib, **lima_options)

        self._handle = None
        self._init_lima_py = None
        self._lima_py = None
        self._KSPLITR = None
        self._KSPLITG = None

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

        @ctypesFF()
        def init_lima_py(PTSTEP, LDWARM, CMICRO, CCSEDIM,
                         LDCRIAUTI, PCRIAUTI, PT0CRIAUTI, PCRIAUTC,
                         LNUCL_IN, LCOLD_IN, LSNOW_IN, LHAIL_IN, LHHONI_IN, LMEYERS_IN,
                         NMOD_IFN_IN, CIFN1_IN, CIFN2_IN, IFN_HOM_IN,
                         IFN_SPECIES_IN, INT_MIXING_IN, NMOD_IMM_IN, IND_SPECIE_IN,
                         CPRISTINE_ICE_LIMA_IN, CHEVRIMED_ICE_LIMA_IN,
                         XALPHAI_IN, XNUI_IN, XALPHAS_IN, XNUS_IN, XALPHAG_IN, XNUG_IN,
                         XFACTNUC_DEP_IN, XFACTNUC_CON_IN, NPHILLIPS_IN,
                         LWARM_IN, LACTI_IN, LRAIN_IN, LSEDC_IN, LACTIT_IN, LBOUND_IN,
                         NMOD_CCN_IN, CCCN1_IN, CCCN2_IN, CCCN3_IN, CCCN4_IN,
                         CCN_HOM_IN, CCN_MODES_IN, HINI_CCN_IN, HTYPE_CCN_IN,
                         XALPHAC_IN, XNUC_IN, XALPHAR_IN, XNUR_IN,
                         XFSOLUB_CCN_IN, XACTEMP_CCN_IN, XAERDIFF_IN, XAERHEIGHT_IN,
                         LSCAV_IN, LAERO_MASS_IN):
            "This function calls the init_lima_py fortran subroutine"
            return ([PTSTEP, LDWARM, CMICRO, CCSEDIM,
                     LDCRIAUTI, PCRIAUTI, PT0CRIAUTI, PCRIAUTC,
                     LNUCL_IN, LCOLD_IN, LSNOW_IN, LHAIL_IN, LHHONI_IN, LMEYERS_IN,
                     NMOD_IFN_IN, CIFN1_IN, CIFN2_IN, IFN_HOM_IN,
                     IFN_SPECIES_IN, INT_MIXING_IN, NMOD_IMM_IN, IND_SPECIE_IN,
                     CPRISTINE_ICE_LIMA_IN, CHEVRIMED_ICE_LIMA_IN,
                     XALPHAI_IN, XNUI_IN, XALPHAS_IN, XNUS_IN, XALPHAG_IN, XNUG_IN,
                     XFACTNUC_DEP_IN, XFACTNUC_CON_IN, NPHILLIPS_IN,
                     LWARM_IN, LACTI_IN, LRAIN_IN, LSEDC_IN, LACTIT_IN, LBOUND_IN,
                     NMOD_CCN_IN, CCCN1_IN, CCCN2_IN, CCCN3_IN, CCCN4_IN,
                     CCN_HOM_IN, CCN_MODES_IN, HINI_CCN_IN, HTYPE_CCN_IN,
                     XALPHAC_IN, XNUC_IN, XALPHAR_IN, XNUR_IN,
                     XFSOLUB_CCN_IN, XACTEMP_CCN_IN, XAERDIFF_IN, XAERHEIGHT_IN,
                     LSCAV_IN, LAERO_MASS_IN],
                    [(numpy.float64, None, IN), #REAL, INTENT (IN) :: PTSTEP
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN) :: LDWARM
                     (numpy.str, (4, ), IN), #CHARACTER(4), INTENT (IN) :: CMICRO 
                     (numpy.str, (4, ), IN), #CHARACTER(4), INTENT (IN) :: CCSEDIM
                     (numpy.int64, None, OUT), #INTEGER, INTENT (OUT) :: KSPLITR
                     (numpy.int64, None, OUT), #INTEGER, INTENT (OUT) :: KSPLITG
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN) :: LDCRIAUTI
                     (numpy.float64, None, IN), #REAL, INTENT (IN) :: PCRIAUTI
                     (numpy.float64, None, IN), #REAL, INTENT (IN) :: PT0CRIAUTI
                     (numpy.float64, None, IN), #REAL, INTENT (IN) :: PCRIAUTC
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN)       :: LNUCL_IN
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN)       :: LCOLD_IN
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN)       :: LSNOW_IN
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN)       :: LHAIL_IN
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN)       :: LHHONI_IN
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN)       :: LMEYERS_IN
                     (numpy.int64, None, IN), #INTEGER, INTENT (IN)       :: NMOD_IFN_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: CIFN1_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: CIFN2_IN
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN)       :: IFN_HOM_IN
                     (numpy.str, (6, ), IN), #CHARACTER(6), INTENT (IN) :: IFN_SPECIES_IN
                     (numpy.str, (4, ), IN), #CHARACTER(4), INTENT (IN) :: INT_MIXING_IN
                     (numpy.int64, None, IN), #INTEGER, INTENT (IN)       :: NMOD_IMM_IN
                     (numpy.int64, None, IN), #INTEGER, INTENT (IN)       :: IND_SPECIE_IN
                     (numpy.str, (4, ), IN), #CHARACTER(4), INTENT (IN) :: CPRISTINE_ICE_LIMA_IN
                     (numpy.str, (4, ), IN), #CHARACTER(4), INTENT (IN) :: CHEVRIMED_ICE_LIMA_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: XALPHAI_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: XNUI_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: XALPHAS_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: XNUS_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: XALPHAG_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: XNUG_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: XFACTNUC_DEP_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: XFACTNUC_CON_IN
                     (numpy.int64, None, IN), #INTEGER, INTENT (IN)       :: NPHILLIPS_IN
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN)       :: LWARM_IN
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN)       :: LACTI_IN
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN)       :: LRAIN_IN
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN)       :: LSEDC_IN
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN)       :: LACTIT_IN
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN)       :: LBOUND_IN
                     (numpy.int64, None, IN), #INTEGER, INTENT (IN)       :: NMOD_CCN_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: CCCN1_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: CCCN2_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: CCCN3_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: CCCN4_IN
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN)       :: CCN_HOM_IN
                     (numpy.str, (8, ), IN), #CHARACTER(8), INTENT (IN) :: CCN_MODES_IN
                     (numpy.str, (3, ), IN), #CHARACTER(3), INTENT (IN) :: HINI_CCN_IN
                     (numpy.str, (1, ), IN), #CHARACTER(1), INTENT (IN) :: HTYPE_CCN_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: XALPHAC_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: XNUC_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: XALPHAR_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: XNUR_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: XFSOLUB_CCN_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: XACTEMP_CCN_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: XAERDIFF_IN
                     (numpy.float64, None, IN), #REAL, INTENT (IN)          :: XAERHEIGHT_IN
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN)       :: LSCAV_IN
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN)       :: LAERO_MASS_IN
                    ],
                    None)
        self._init_lima_py = init_lima_py

        @ctypesFF()
        def lima_py(K1, K2, K3, K4, K5,
                    LWARM, LACTIT, LSEDC, LRAIN,
                    LCOLD, LHHONI, LSEDI,
                    PTSTEP, KRR,
                    KSPLITR, KSPLITG,
                    PDZZ, PRHODJ, PRHODREF, PEXNREF,
                    PPABST, PW_NU,
                    PRT, PSVT, PTHT,
                    PRS, PSVS, PTHS):
            "This function calls the lima_py fortran subroutine"
            return ([K1, K2, K3, K4, K5,
                     LWARM, LACTIT, LSEDC, LRAIN,
                     LCOLD, LHHONI, LSEDI,
                     PTSTEP, KRR,
                     KSPLITR, KSPLITG,
                     PDZZ, PRHODJ, PRHODREF, PEXNREF,
                     PPABST, PW_NU,
                     PRT, PSVT, PTHT,
                     PRS, PSVS, PTHS],
                    [
                     (numpy.int64, None, IN), #INTEGER,                     INTENT(IN)    :: K1      ! Nx
                     (numpy.int64, None, IN), #INTEGER,                     INTENT(IN)    :: K2      ! Ny
                     (numpy.int64, None, IN), #INTEGER,                     INTENT(IN)    :: K3      ! Nz
                     (numpy.int64, None, IN), #INTEGER,                     INTENT(IN)    :: K4      ! size(PRT)
                     (numpy.int64, None, IN), #INTEGER,                     INTENT(IN)    :: K5      ! size(PSVT)
                     (numpy.bool, None, IN), #LOGICAL,                     INTENT(IN)    :: LWARM
                     (numpy.bool, None, IN), #LOGICAL,                     INTENT(IN)    :: LACTIT
                     (numpy.bool, None, IN), #LOGICAL,                     INTENT(IN)    :: LSEDC
                     (numpy.bool, None, IN), #LOGICAL,                     INTENT(IN)    :: LRAIN
                     (numpy.bool, None, IN), #LOGICAL,                     INTENT(IN)    :: LCOLD
                     (numpy.bool, None, IN), #LOGICAL,                     INTENT(IN)    :: LHHONI
                     (numpy.bool, None, IN), #LOGICAL,                     INTENT(IN)    :: LSEDI
                     (numpy.float64, None, IN), #REAL,                        INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
                     (numpy.int64, None, IN), #INTEGER,                     INTENT(IN)    :: KRR     ! Number of moist variable
                     (numpy.int64, None, IN), #INTEGER,                     INTENT(IN)    :: KSPLITR ! Number of small time step integration for  rain sedimendation
                     (numpy.int64, None, IN), #INTEGER,                     INTENT(IN)    :: KSPLITG ! Number of small time step integration for  hail sedimendation
                     (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
                     (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
                     (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PRHODREF! Reference density
                     (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PEXNREF ! Reference Exner function
                     (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PPABST  ! absolute pressure at t
                     (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PW_NU
                     (numpy.float64, (K1, K2, K3, K4), IN), #REAL, DIMENSION(K1,K2,K3,K4),INTENT(IN)    :: PRT
                     (numpy.float64, (K1, K2, K3, K5), IN), #REAL, DIMENSION(K1,K2,K3,K5),INTENT(IN)    :: PSVT
                     (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(K1,K2,K3),   INTENT(IN)    :: PTHT
                     (numpy.float64, (K1, K2, K3, K4), INOUT), #REAL, DIMENSION(K1,K2,K3,K4),INTENT(INOUT) :: PRS
                     (numpy.float64, (K1, K2, K3, K5), INOUT), #REAL, DIMENSION(K1,K2,K3,K5),INTENT(INOUT) :: PSVS
                     (numpy.float64, (K1, K2, K3), INOUT), #REAL, DIMENSION(K1,K2,K3),   INTENT(INOUT) :: PTHS
                     (numpy.float64, (K1, K2), OUT), #REAL, DIMENSION(K1,K2),      INTENT(OUT)   :: PINPRC! Cloud instant precip
                     (numpy.float64, (K1, K2), OUT), #REAL, DIMENSION(K1,K2),      INTENT(OUT)   :: PINPRR! Rain instant precip
                     (numpy.float64, (K1, K2, K3), OUT), #REAL, DIMENSION(K1,K2,K3),   INTENT(OUT)   :: PINPRR3D! Rain inst precip 3D
                     (numpy.float64, (K1, K2, K3), OUT), #REAL, DIMENSION(K1,K2,K3),   INTENT(OUT)   :: PEVAP3D! Rain evap profile
                     (numpy.float64, (K1, K2), OUT), #REAL, DIMENSION(K1,K2),      INTENT(OUT)   :: PINPRS! Snow instant precip
                     (numpy.float64, (K1, K2), OUT), #REAL, DIMENSION(K1,K2),      INTENT(OUT)   :: PINPRG! Graupel instant precip
                     (numpy.float64, (K1, K2), OUT), #REAL, DIMENSION(K1,K2),      INTENT(OUT)   :: PINPRH! Graupel instant precip
                     (numpy.float64, (K1, K2), OUT), #REAL, DIMENSION(K1,K2),      INTENT(OUT)   :: PINDEP!Cloud droplets deposition
                    ],
                    None)
        self._lima_py = lima_py

        CCSEDIM = 'SPLI'
        LDCRIAUTI = False
        PCRIAUTI, PT0CRIAUTI, PCRIAUTC = 0., 0., 0.

        LDWARM = self._lima_options['LDWARM']
        LNUCL_IN = self._lima_options['LNUCL']
        LCOLD_IN = self._lima_options['LCOLD']
        LSNOW_IN = self._lima_options['LSNOW']
        LHAIL_IN = self._lima_options['LHAIL']
        LHHONI_IN = self._lima_options['LHHONI']
        LMEYERS_IN = self._lima_options['LMEYERS']
        NMOD_IFN_IN = self._lima_options['NMOD_IFN']
        CIFN1_IN = self._lima_options['CIFN1']
        CIFN2_IN = self._lima_options['CIFN2']
        IFN_HOM_IN = self._lima_options['IFN_HOM']
        IFN_SPECIES_IN = self._lima_options['IFN_SPECIES']
        INT_MIXING_IN = self._lima_options['INT_MIXING']
        NMOD_IMM_IN = self._lima_options['NMOD_IMM']
        IND_SPECIE_IN = self._lima_options['IND_SPECIE']
        CPRISTINE_ICE_LIMA_IN = self._lima_options['CPRISTINE_ICE_LIMA']
        CHEVRIMED_ICE_LIMA_IN = self._lima_options['CHEVRIMED_ICE_LIMA']
        XALPHAI_IN = self._lima_options['XALPHAI']
        XNUI_IN = self._lima_options['XNUI']
        XALPHAS_IN = self._lima_options['XALPHAS']
        XNUS_IN = self._lima_options['XNUS']
        XALPHAG_IN = self._lima_options['XALPHAG']
        XNUG_IN = self._lima_options['XNUG']
        XFACTNUC_DEP_IN = self._lima_options['XFACTNUC_DEP']
        XFACTNUC_CON_IN = self._lima_options['XFACTNUC_CON']
        NPHILLIPS_IN = self._lima_options['NPHILLIPS']
        LWARM_IN = self._lima_options['LWARM']
        LACTI_IN = self._lima_options['LACTI']
        LRAIN_IN = self._lima_options['LRAIN']
        LSEDC_IN = self._lima_options['LSEDC']
        LACTIT_IN = self._lima_options['LACTIT']
        NMOD_CCN_IN = self._lima_options['NMOD_CCN']
        CCCN1_IN = self._lima_options['CCCN1']
        CCCN2_IN = self._lima_options['CCCN2']
        CCCN3_IN = self._lima_options['CCCN3']
        CCCN4_IN = self._lima_options['CCCN4']
        CCN_HOM_IN = self._lima_options['CCN_HOM']
        CCN_MODES_IN = self._lima_options['CCN_MODES']
        HINI_CCN_IN = self._lima_options['HINI_CCN']
        HTYPE_CCN_IN = self._lima_options['HTYPE_CCN']
        XALPHAC_IN = self._lima_options['XALPHAC']
        XNUC_IN = self._lima_options['XNUC']
        XALPHAR_IN = self._lima_options['XALPHAR']
        XNUR_IN = self._lima_options['XNUR']
        XFSOLUB_CCN_IN = self._lima_options['XFSOLUB_CCN']
        XACTEMP_CCN_IN = self._lima_options['XACTEMP_CCN']
        XAERDIFF_IN = self._lima_options['XAERDIFF']
        XAERHEIGHT_IN = self._lima_options['XAERHEIGHT']
        LSCAV_IN = self._lima_options['LSCAV']
        LAERO_MASS_IN = self._lima_options['LAERO_MASS']

        self._KSPLITR, self._KSPLITG = self._init_lima_py(
                         self._dt, LDWARM, 'LIMA', CCSEDIM,
                         LDCRIAUTI, PCRIAUTI, PT0CRIAUTI, PCRIAUTC,
                         LNUCL_IN, LCOLD_IN, LSNOW_IN, LHAIL_IN, LHHONI_IN, LMEYERS_IN,
                         NMOD_IFN_IN, CIFN1_IN, CIFN2_IN, IFN_HOM_IN,
                         IFN_SPECIES_IN, INT_MIXING_IN, NMOD_IMM_IN, IND_SPECIE_IN,
                         CPRISTINE_ICE_LIMA_IN, CHEVRIMED_ICE_LIMA_IN,
                         XALPHAI_IN, XNUI_IN, XALPHAS_IN, XNUS_IN, XALPHAG_IN, XNUG_IN,
                         XFACTNUC_DEP_IN, XFACTNUC_CON_IN, NPHILLIPS_IN,
                         LWARM_IN, LACTI_IN, LRAIN_IN, LSEDC_IN, LACTIT_IN, False,
                         NMOD_CCN_IN, CCCN1_IN, CCCN2_IN, CCCN3_IN, CCCN4_IN,
                         CCN_HOM_IN, CCN_MODES_IN, HINI_CCN_IN, HTYPE_CCN_IN,
                         XALPHAC_IN, XNUC_IN, XALPHAR_IN, XNUR_IN,
                         XFSOLUB_CCN_IN, XACTEMP_CCN_IN, XAERDIFF_IN, XAERHEIGHT_IN,
                         LSCAV_IN, LAERO_MASS_IN)

    def finalize(self):
        """
        The method closes the shared library so we are able
        to load a new one with the same symbol names.
        The method removes temporary files too.
        """
        super().finalize()
        del self._init_lima_py, self._lima_py
        ctypesForFortran.dlclose(self._handle)
        if os.path.exists('fort.1'):
            os.remove('fort.1')

    def build_init_state(self, state):
        state = super().build_init_state(state)
        needed = ['T', 'P', 'rv', 'rc', 'rr', 'ri', 'rs', 'rg',
                  'ccn1ft', 'ccn1at', 'ifn1ft', 'ifn1at']
        if self._lima_options['LHAIL']:
            needed += ['rh']
        for var in needed:
            if var not in state:
                raise ValueError(var + " must be in state")
        if (not self._lima_options['LHAIL']) and 'rh' in state:
            state['rg'] += state['rh']
            state['rh'] = state['rh'] * 0.
        for var in ['ni', 'nc', 'nr']:
            if not var in state:
                state[var] = numpy.zeros(state['ri'].shape)
        if not 'SIGS' in state:
            state['SIGS'] = numpy.zeros(state['ri'].shape)
        return state

    def execute(self, previous_state, timestep, timestep_number):
        """
        We add an extension zone of one point all around.
        We compute input variables, then we call the microphysics.
        We come back to the original variables without the extension zone.
        """
        super().execute(previous_state, timestep, timestep_number)

        #Options
        LSEDI=False

        #Dimensions
        PT = previous_state['T']
        shape = PT.shape
        shape_ori = shape
        if len(shape) > 3:
            raise ValueError("Maximum shape length is 3")
        if len(shape) == 1:
            shape = tuple(list(shape) + [1, 1])
        elif len(shape) == 2:
            shape = tuple(list(shape) + [1])
        assert len(shape) == 3, "Internal error"
        KRR = 7 if self._lima_options['LHAIL'] else 6
        shapeExt = tuple([shape[i] + 2 for i in range(3)])

        #Auxiliary fields (not really needed)
        PDZZ = numpy.ones(shapeExt) #values used only for sedimentation
        PRHODJ = numpy.ones(shapeExt)
        PRHODREF = numpy.ones(shapeExt)

        #Temperature
        P0 = 100000.
        Boltz = 1.380658E-23
        Avogadro = 6.0221367E+23
        Md = 28.9644E-3
        Rd = Avogadro * Boltz / Md
        Cpd = 7. * Rd / 2.
        P = numpy.ndarray(shape=shapeExt, dtype=numpy.float64, order='F')
        P[...] = 100000.
        P[1:-1, 1:-1, 1:-1] = previous_state['P'].reshape(shape)
        PEXN = (P / P0) ** (Rd / Cpd)
        PTHT = PT.reshape(shape) / PEXN

        #State variables
        PRT = numpy.ndarray(shape=tuple(list(shapeExt) + [KRR]), dtype=numpy.float64, order='F')
        PRT[...] = 0.
        PRT[1:-1, 1:-1, 1:-1, 0] = previous_state['rv'].reshape(shape)
        PRT[1:-1, 1:-1, 1:-1, 1] = previous_state['rc'].reshape(shape)
        PRT[1:-1, 1:-1, 1:-1, 2] = previous_state['rr'].reshape(shape)
        PRT[1:-1, 1:-1, 1:-1, 3] = previous_state['ri'].reshape(shape)
        PRT[1:-1, 1:-1, 1:-1, 4] = previous_state['rs'].reshape(shape)
        PRT[1:-1, 1:-1, 1:-1, 5] = previous_state['rg'].reshape(shape)
        if self._lima_options['LHAIL']:
            PRT[1:-1, 1:-1, 1:-1, 6] = previous_state['rh'].reshape(shape)
        PSVT = numpy.ndarray(shape=tuple(list(shapeExt) + [8]), dtype=numpy.float64, order='F')
        PSVT[...] = 0.
        PSVT[1:-1, 1:-1, 1:-1, 0] = previous_state['nc']
        PSVT[1:-1, 1:-1, 1:-1, 1] = previous_state['nr']
        PSVT[1:-1, 1:-1, 1:-1, 2] = previous_state['ccn1ft']
        PSVT[1:-1, 1:-1, 1:-1, 3] = previous_state['ccn1at']
        PSVT[1:-1, 1:-1, 1:-1, 4] = previous_state['ni']
        PSVT[1:-1, 1:-1, 1:-1, 5] = previous_state['ifn1ft']
        PSVT[1:-1, 1:-1, 1:-1, 6] = previous_state['ifn1at']
        PSVT[1:-1, 1:-1, 1:-1, 7] = 0.

        #Tendencies
        PTHS = PTHT / timestep
        PRS = PRT / timestep
        PSVS = PSVT / timestep

        PW_NU = numpy.ndarray(shapeExt, dtype=numpy.float64, order='F')
        PW_NU[...] = 5.

        result = self._lima_py(shapeExt[0], shapeExt[1], shapeExt[2], KRR, 8,
                               self._lima_options['LWARM'], self._lima_options['LACTIT'],
                               self._lima_options['LSEDC'], self._lima_options['LRAIN'],
                               self._lima_options['LCOLD'], self._lima_options['LHHONI'], LSEDI,
                               timestep, KRR,
                               self._KSPLITR, self._KSPLITG,
                               PDZZ, PRHODJ, PRHODREF, PEXN,
                               P, PW_NU,
                               PRT, PSVT, PTHT,
                               PRS, PSVS, PTHS)
        (PRS, PSVS, PTHS, PINPRC, PINPRR, PINPRR3D, PEVAP3D, PINPRS, PINPRG, PINPRH, PINDEP) = result
        next_state = {}
        next_state['T'] = (PTHS[1:-1, 1:-1, 1:-1] * timestep * PEXN[1:-1, 1:-1, 1:-1]).reshape(shape_ori)
        next_state['rv'] = PRS[1:-1, 1:-1, 1:-1, 0].reshape(shape_ori) * timestep
        next_state['rc'] = PRS[1:-1, 1:-1, 1:-1, 1].reshape(shape_ori) * timestep
        next_state['rr'] = PRS[1:-1, 1:-1, 1:-1, 2].reshape(shape_ori) * timestep
        next_state['ri'] = PRS[1:-1, 1:-1, 1:-1, 3].reshape(shape_ori) * timestep
        next_state['rs'] = PRS[1:-1, 1:-1, 1:-1, 4].reshape(shape_ori) * timestep
        next_state['rg'] = PRS[1:-1, 1:-1, 1:-1, 5].reshape(shape_ori) * timestep
        if self._lima_options['LHAIL']:
            next_state['rh'] = PRS[1:-1, 1:-1, 1:-1, 6].reshape(shape_ori) * timestep
        next_state['nc'] = PSVS[1:-1, 1:-1, 1:-1, 0].reshape(shape_ori) * timestep
        next_state['nr'] = PSVS[1:-1, 1:-1, 1:-1, 1].reshape(shape_ori) * timestep
        next_state['ccn1ft'] = PSVS[1:-1, 1:-1, 1:-1, 2].reshape(shape_ori) * timestep
        next_state['ccn1at'] = PSVS[1:-1, 1:-1, 1:-1, 3].reshape(shape_ori) * timestep
        next_state['ni'] = PSVS[1:-1, 1:-1, 1:-1, 4].reshape(shape_ori) * timestep
        next_state['ifn1ft'] = PSVS[1:-1, 1:-1, 1:-1, 5].reshape(shape_ori) * timestep
        next_state['ifn1at'] = PSVS[1:-1, 1:-1, 1:-1, 6].reshape(shape_ori) * timestep
        return next_state


