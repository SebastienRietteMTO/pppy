#!/usr/bin/env python 3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

#Written by S. Riette

"""
This module contains the PPPY implementation for calling ICE microphysical
scheme of Meso-NH model.
"""

import numpy
import os
from pppy import ctypesForFortran
import pppy

class pppy_microPhyIce_MNH53(pppy.PPPY):
    """
    PPPY implementation for calling ICE microphysical scheme of Meso-NH model.
    """
    def __init__(self, dt, method, name, tag, solib,
                 ice_version, HSUBG_AUCV_RC, adjustment=True):
        """
        In addition to dt, method, name and tag parameters
        defined in the PPPY class, this parameterization
        needs the following parameters:
        solib            : path to the shared library to use
        ice_version      :
        HSUBG_AUCV_RC    :
        adjustment       : (optional)
        """
        super().__init__(dt, method, name, tag, solib=solib,
                         ice_version=ice_version, HSUBG_AUCV_RC=HSUBG_AUCV_RC,
                         adjustment=adjustment)
        self._solib = solib
        self._ice_version = ice_version
        self._HSUBG_AUCV_RC = HSUBG_AUCV_RC
        self._adjustment = adjustment

        self._handle = None
        self._init_py = None
        self._rain_ice_py_3 = None
        self._rain_ice_py_4 = None
        self._ice_adjust_py = None
        self._KSPLITR = None

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

        #Define fortran routine to use through ctypesForFortran
        IN = ctypesForFortran.IN
        OUT = ctypesForFortran.OUT
        INOUT = ctypesForFortran.INOUT

        ctypesFF, self._handle = ctypesForFortran.ctypesForFortranFactory(self._solib)

        @ctypesFF()
        def init_py(KULOUT, PTSTEP, LDWARM, CMICRO, CCSEDIM, LDCRIAUTI, PCRIAUTI,
                    PT0CRIAUTI, PCRIAUTC):
            "This function calls the init_py fortran subroutine"
            return ([KULOUT, PTSTEP, LDWARM, CMICRO, CCSEDIM, LDCRIAUTI, PCRIAUTI,
                     PT0CRIAUTI, PCRIAUTC],
                    [(numpy.int64, None, IN), #INTEGER, INTENT (IN) :: KULOUT
                     (numpy.float64, None, IN), #REAL, INTENT (IN) :: PTSTEP
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN) :: LDWARM
                     (numpy.str, (4, ), IN), #CHARACTER(4), INTENT (IN) :: CMICRO
                     (numpy.int64, None, OUT), #INTEGER, INTENT (OUT) :: KSPLITR
                     (numpy.str, (4, ), IN), #CHARACTER(4), INTENT (IN) :: CCSEDIM
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN) :: LDCRIAUTI
                     (numpy.float64, None, IN), #REAL, INTENT (IN) :: PCRIAUTI
                     (numpy.float64, None, IN), #REAL, INTENT (IN) :: PT0CRIAUTI
                     (numpy.float64, None, IN), #REAL, INTENT (IN) :: PCRIAUTC
                    ],
                    None)
        self._init_py = init_py

        #SUBROUTINE RAIN_ICE_PY (K1, K2, K3, OSEDIC, HSEDIM, HSUBG_AUCV_RC, OWARM,KKA,KKU,KKL,  &
        #                        KSPLITR, PTSTEP, KMI, KRR,              &
        #                        PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR, &
        #                        PTHT, PRVT, PRCT, PRRT, PRIT, PRST, &
        #                        PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS, &
        #                        PINPRC, PINPRR, PINPRR3D, PEVAP3D,                    &
        #                        PINPRS, PINPRG, PSIGS, PINDEP, PSEA, PTOWN,                   &
        #                        PRHT, PRHS, PINPRH, OCONVHG
        def rain_ice_py_34(ice_version,
                           K1, K2, K3, OSEDIC, HSEDIM, HSUBG_AUCV_RC, OWARM,KKA,KKU,KKL,
                           KSPLITR, PTSTEP, KMI, KRR,
                           PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,
                           PTHT, PRVT, PRCT, PRRT, PRIT, PRST,
                           PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,
                           PSIGS, PSEA, PTOWN,
                           PRHT=None, PRHS=None, OCONVHG=None):
            """
            This function is the body of rain_ice_py_3, rain_ice_py_4,
            rain_ice_red_py_3, rain_ice_red_py_4
            """
            varList = [K1, K2, K3, OSEDIC, HSEDIM, HSUBG_AUCV_RC, OWARM,KKA,KKU,KKL,
                       KSPLITR, PTSTEP, KMI, KRR,
                       PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,
                       PTHT, PRVT, PRCT, PRRT, PRIT, PRST,
                       PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,
                       PSIGS, PSEA, PTOWN]
            if ice_version == 4:
                varList.extend([PRHT, PRHS, OCONVHG])

            signature = [(numpy.int64, None, IN), (numpy.int64, None, IN), (numpy.int64, None, IN), #K1,K2,K3
                         (numpy.bool, None, IN), #LOGICAL,                  INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
                         (numpy.str, (4, ), IN), #CHARACTER(LEN=4),         INTENT(IN)    :: HSEDIM ! Sedimentation scheme
                         (numpy.str, (4, ), IN), #CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV_RC ! Kind of Subgrid autoconversion method
                         (numpy.bool, None, IN), #LOGICAL,                  INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to form by warm processes (Kessler scheme)
                         (numpy.int64, None, IN), #INTEGER,                  INTENT(IN)    :: KKA   !near ground array index
                         (numpy.int64, None, IN), #INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
                         (numpy.int64, None, IN), #INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
                         (numpy.int64, None, IN), #INTEGER,                  INTENT(IN)    :: KSPLITR ! Number of small time step integration for  rain sedimendation
                         (numpy.float64, None, IN), #REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
                         (numpy.int64, None, IN), #INTEGER,                  INTENT(IN)    :: KMI     ! Model index
                         (numpy.int64, None, IN), #INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
                         (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
                         (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
                         (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
                         (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
                         (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! absolute pressure at t
                         (numpy.float64, (K1, K2, K3), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
                         (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCLDFR  ! Convective Mass Flux Cloud fraction
                         (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
                         (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
                         (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
                         (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
                         (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
                         (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
                         (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
                         (numpy.float64, (K1, K2, K3), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
                         (numpy.float64, (K1, K2, K3), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
                         (numpy.float64, (K1, K2, K3), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
                         (numpy.float64, (K1, K2, K3), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
                         (numpy.float64, (K1, K2, K3), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
                         (numpy.float64, (K1, K2, K3), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
                         (numpy.float64, (K1, K2, K3), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source
                         (numpy.float64, (K1, K2), OUT), #REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRC! Cloud instant precip        ####Over dimensionned
                         (numpy.float64, (K1, K2), OUT), #REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRR! Rain instant precip        ####Over dimensionned
                         (numpy.float64, (K1, K2, K3), OUT), #REAL, DIMENSION(K1,K2,K3), INTENT(OUT)   :: PINPRR3D! Rain inst precip 3D
                         (numpy.float64, (K1, K2, K3), OUT), #REAL, DIMENSION(:,:,:), INTENT(OUT)     :: PEVAP3D! Rain evap profile
                         (numpy.float64, (K1, K2), OUT), #REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRS! Snow instant precip        ####Over dimensionned
                         (numpy.float64, (K1, K2), OUT), #REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRG! Graupel instant precip        ####Over dimensionned
                         (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
                         (numpy.float64, (K1, K2), OUT), #REAL, DIMENSION(K1,K2), INTENT(OUT)     :: PINDEP  ! Cloud instant deposition
                         (numpy.float64, (K1, K2), IN), #REAL, DIMENSION(:,:), INTENT(IN)        :: PSEA ! Sea Mask
                         (numpy.float64, (K1, K2), IN), #REAL, DIMENSION(:,:), INTENT(IN)        :: PTOWN! Fraction that is town
                        ]
            if ice_version == 4:
                signature.extend([(numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:), OPTIONAL,  INTENT(IN)    :: PRHT    ! Hail m.r. at t
                                  (numpy.float64, (K1, K2, K3), INOUT), #REAL, DIMENSION(:,:,:), OPTIONAL,  INTENT(INOUT) :: PRHS    ! Hail m.r. source
                                  (numpy.float64, (K1, K2), OUT), #REAL, DIMENSION(K1,K2), INTENT(OUT)       :: PINPRH! Hail instant precip
                                  (numpy.bool, None, IN), #LOGICAL, OPTIONAL,                 INTENT(IN)    :: OCONVHG! Switch for conversion from hail to graupel
                                 ])
            return (varList, signature, None)

        @ctypesFF()
        def rain_ice_py_3(K1, K2, K3, OSEDIC, HSEDIM, HSUBG_AUCV_RC, OWARM,KKA,KKU,KKL,
                          KSPLITR, PTSTEP, KMI, KRR,
                          PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,
                          PTHT, PRVT, PRCT, PRRT, PRIT, PRST,
                          PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,
                          PSIGS, PSEA, PTOWN):
            "This function calls the rain_ice_py_3 fortran subroutine"
            return rain_ice_py_34(3,
                                  K1, K2, K3, OSEDIC, HSEDIM, HSUBG_AUCV_RC, OWARM,KKA,KKU,KKL,
                                  KSPLITR, PTSTEP, KMI, KRR,
                                  PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,
                                  PTHT, PRVT, PRCT, PRRT, PRIT, PRST,
                                  PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,
                                  PSIGS, PSEA, PTOWN)
        self._rain_ice_py_3 = rain_ice_py_3

        @ctypesFF()
        def rain_ice_py_4(K1, K2, K3, OSEDIC, HSEDIM, HSUBG_AUCV_RC, OWARM,KKA,KKU,KKL,
                          KSPLITR, PTSTEP, KMI, KRR,
                          PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,
                          PTHT, PRVT, PRCT, PRRT, PRIT, PRST,
                          PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,
                          PSIGS, PSEA, PTOWN,
                          PRHT, PRHS, OCONVHG):
            "This function calls the rain_ice_py_4 fortran subroutine"
            return rain_ice_py_34(4,
                                  K1, K2, K3, OSEDIC, HSEDIM, HSUBG_AUCV_RC, OWARM,KKA,KKU,KKL,
                                  KSPLITR, PTSTEP, KMI, KRR,
                                  PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,
                                  PTHT, PRVT, PRCT, PRRT, PRIT, PRST,
                                  PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,
                                  PSIGS, PSEA, PTOWN,
                                  PRHT, PRHS, OCONVHG)
        self._rain_ice_py_4 = rain_ice_py_4

        @ctypesFF()
        def ice_adjust_py(K1, K2, K3, PTSTEP, PSIGQSAT,
                         PRHODJ, PEXNREF, PSIGS, PPABST, PTH,
                         PRV, PRC, PRI, PRR, PRS, PRG, PRH):
            "This function calls the ice_adjust_py fortran subroutine"
            return ([K1, K2, K3, PTSTEP, PSIGQSAT,
                     PRHODJ, PEXNREF, PSIGS,PPABST, PTH,
                     PRV, PRC, PRI, PRR, PRS, PRG, PRH],
                    [(numpy.int64, None, IN), (numpy.int64, None, IN), (numpy.int64, None, IN), #K1,K2,K3
                     (numpy.float64, None, IN), #REAL, INTENT (IN) :: PTSTEP
                     (numpy.float64, None, IN), #REAL, INTENT (IN) :: PSIGQSAT
                     (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ
                     (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF
                     (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PSIGS
                     (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! absolute pressure at t
                     (numpy.float64, (K1, K2, K3), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTH    ! Theta at time t
                     (numpy.float64, (K1, K2, K3), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRV    ! Water vapor m.r. at t
                     (numpy.float64, (K1, K2, K3), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRC    ! Cloud water m.r. at t
                     (numpy.float64, (K1, K2, K3), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRI    ! Pristine ice m.r. at t
                     (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRR    ! Rain water m.r. at t
                     (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRS    ! Snow/aggregate m.r. at t
                     (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRG    ! Graupel/hail m.r. at t
                     (numpy.float64, (K1, K2, K3), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRH    ! Hail m.r. at t
                    ],
                    None)
        self._ice_adjust_py = ice_adjust_py

        if self._ice_version not in ['ICE3', 'ICE4']:
            raise ValueError("ice_version must be 'ICE3', 'ICE4'")
        self._KSPLITR = self._init_py(1, self._dt, True,
                                      self._ice_version, 'SPLI', False, 0., 0., 0.)

    def finalize(self):
        """
        The method closes the shared library so we are able to load a new one
        with the same symbol names.
        Moreover the fort.1 temporary file is removed.
        """
        super().finalize()
        del self._init_py, self._rain_ice_py_3, self._rain_ice_py_4
        del self._ice_adjust_py
        ctypesForFortran.dlclose(self._handle)
        if os.path.exists('fort.1'):
            os.remove('fort.1')

    def build_init_state(self, state):
        "The method deals with the rh, ni and sigs variables"
        state = super().build_init_state(state)
        needed = ['T', 'P', 'rv', 'rc', 'rr', 'ri', 'rs', 'rg']
        if self._ice_version in ['ICE4', 'OLD4']:
            needed += ['rh']
        for var in needed:
            if var not in state:
                raise ValueError(var + " must be in state")
        if self._ice_version in ['ICE3', 'OLD3'] and 'rh' in state:
            state['rg'] += state['rh']
            state['rh'] = state['rh'] * 0.
        if not 'ni' in state:
            state['ni'] = numpy.zeros(state['ri'].shape)
        if not 'SIGS' in state:
            state['SIGS'] = numpy.zeros(state['ri'].shape)
        return state

    def execute(self, previous_state, timestep, timestep_number):
        "The method performs the computational part of calling the ice scheme"
        super().execute(previous_state, timestep, timestep_number)

        #Options
        OSEDIC = False
        HSEDIM = 'NONE'
        OWARM = True
        hail = self._ice_version in ['ICE4', 'OLD4']

        #Dimensions
        shape = previous_state['T'].shape
        shape_ori = shape
        if len(shape) > 3:
            raise ValueError("Maximum shape length is 3")
        if len(shape) == 1:
            shape = tuple(list(shape) + [1, 1])
        elif len(shape) == 2:
            shape = tuple(list(shape) + [1])
        assert len(shape) == 3, "Internal error"
        KKA = 1
        KKU = shape[2]
        KKL = 1
        KMI = 0

        #Auxiliary fields (not really needed)
        PDZZ = numpy.ones(shape) #values used only for sedimentation
        PRHODJ = numpy.ones(shape)
        PRHODREF = numpy.ones(shape)
        PCLDFR = numpy.ones(shape)
        PSEA = numpy.zeros(tuple(shape[:2]))
        PTOWN = numpy.zeros(tuple(shape[:2]))

        #Other variables than T
        P = previous_state['P'].reshape(shape)
        PRV = previous_state['rv'].reshape(shape)
        PRC = previous_state['rc'].reshape(shape)
        PRR = previous_state['rr'].reshape(shape)
        PRI = previous_state['ri'].reshape(shape)
        PRS = previous_state['rs'].reshape(shape)
        PRG = previous_state['rg'].reshape(shape)
        if hail:
            PRH = previous_state['rh'].reshape(shape)
        else:
            PRH = previous_state['rg'].reshape(shape) * 0.

        #Temperature
        P0 = 100000.
        Boltz = 1.380658E-23
        Avogadro = 6.0221367E+23
        Md = 28.9644E-3
        Rd = Avogadro * Boltz / Md
        Cpd = 7. * Rd / 2.
        PEXN = (P / P0) ** (Rd / Cpd)
        PTHT = previous_state['T'].reshape(shape) / PEXN

        #Adjustment
        if self._adjustment:
            PTHT, PRV, PRC, PRI = self._ice_adjust_py(shape[0], shape[1], shape[2], timestep, 0.,
                                                      PRHODJ, PEXN,
                                                      previous_state['SIGS'].reshape(shape),
                                                      P, PTHT,
                                                      PRV, PRC, PRI, PRR, PRS, PRG, PRH)

        #Initial tendencies
        PTHS = PTHT / timestep
        PRVS = PRV / timestep
        PRCS = PRC / timestep
        PRRS = PRR / timestep
        PRIS = PRI / timestep
        PRSS = PRS / timestep
        PRGS = PRG / timestep
        if hail:
            PRHS = PRH / timestep

        #Microphysics
        rain_ice_py = self._rain_ice_py_4 if hail else self._rain_ice_py_3
        result = rain_ice_py(shape[0], shape[1], shape[2],
                             OSEDIC, HSEDIM, self._HSUBG_AUCV_RC, OWARM,
                             KKA, KKU, KKL, self._KSPLITR,
                             timestep, KMI, 7 if hail else 6,
                             PDZZ, PRHODJ, PRHODREF, PEXN,
                             P, previous_state['ni'].reshape(shape), PCLDFR,
                             PTHT, PRV, PRC, PRR,
                             PRI, PRS, PRG,
                             PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,
                             previous_state['SIGS'].reshape(shape), PSEA, PTOWN,
                             *([] if not hail else [PRH, PRHS, True]))
        next_state = {}
        if self._ice_version in ['ICE3', 'OLD3']:
            (next_state['ni'], PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,
             PINPRC, PINPRR, PINPRR3D, PEVAP3D, PINPRS, PINPRG, PINDEP) = result
        else:
            (next_state['ni'], PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,
             PINPRC, PINPRR, PINPRR3D, PEVAP3D, PINPRS, PINPRG, PINDEP, PRHS, PINPRH) = result
        next_state['ni'] = next_state['ni'].reshape(shape_ori)
        next_state['T'] = (PTHS * timestep * PEXN).reshape(shape_ori)
        next_state['rv'] = PRVS.reshape(shape_ori) * timestep
        next_state['rc'] = PRCS.reshape(shape_ori) * timestep
        next_state['rr'] = PRRS.reshape(shape_ori) * timestep
        next_state['ri'] = PRIS.reshape(shape_ori) * timestep
        next_state['rs'] = PRSS.reshape(shape_ori) * timestep
        next_state['rg'] = PRGS.reshape(shape_ori) * timestep
        if hail:
            next_state['rh'] = PRHS.reshape(shape_ori) * timestep

        return next_state
