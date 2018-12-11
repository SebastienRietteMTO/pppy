#!/usr/bin/env python 3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

#Written by S. Riette

"""
This module contains the PPPY implementation for calling sedimentation
schemes of Meso-NH model.
"""

import numpy
import os
from pppy import ctypesForFortran
import pppy

class pppy_sedim_MNH54(pppy.PPPY):
    """
    PPPY implementation for calling sedimentation schemes of Meso-NH model.
    """

    def __init__(self, dt, method, name, tag, solib, hail, version, maxcfl=None):
        """
        In addition to dt, method, name and tag parameters
        defined in the PPPY class, this parameterization
        needs the following parameters:
        solib   : path to the shared lib
        hail    : True or False to activate or no hail
        version : string that can take the following values:
                   - 'STAT': statistical scheme
                   - 'SPLI': eulerian scheme (old version with constant time-splitting)
                   - 'SPLN': eulerian scheme (new version, future operational implementation)
                   - 'SPL2': eulerian scheme (new version with momentum transport)
        maxcfl  : maximum value of CFL for SPLN and SPL2
        """
        super().__init__(dt, method, name, tag,
                         solib=solib, hail=hail, version=version,
                         **(dict(maxcfl=maxcfl) if maxcfl is not None else {}))

        assert hail in [True, False], "hail must be set to True or False"
        assert version in ['STAT', 'SPLI', 'SPLN', 'SPL2'], "version is unknown"
        if version in ['STAT', 'SPLI']:
            assert maxcfl is None, "maxcfl must be None if version is STAT or SPLI"
        else:
            assert maxcfl is not None, "maxcfl must not be None if version is SPLN or SPL2"

        self._handle = None
        self._aroini_micro_py = None
        self._ini_cst_py = None
        self._sedim_py = None
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

        ctypesFF, self._handle = ctypesForFortran.ctypesForFortranFactory(self._options['solib'])

        @ctypesFF()
        def init_py(KULOUT,PTSTEP,LDWARM,CMICRO,CCSEDIM,LDCRIAUTI,
                    PCRIAUTI,PT0CRIAUTI,PCRIAUTC, PTSTEP_TS, CCSNOWRIMING, PMRSTEP, KMAXITER,
                    LDFEEDBACKT, LDEVLIMIT, LDNULLWETG, LDWETGPOST, LDNULLWETH, LDWETHPOST,
                    PFRACM90, LDCONVHG, CCSUBG_RC_RR_ACCR, CCSUBG_RR_EVAP, CCSUBG_PR_PDF,
                    LDCRFLIMIT, CCFRAC_ICE_ADJUST, PSPLIT_MAXCFL,
                    CCFRAC_ICE_SHALLOW_MF, LDSEDIM_AFTER):
            """
            This function calls the init_py fortran subroutine.
            """
            return ([KULOUT,PTSTEP,LDWARM,CMICRO,CCSEDIM,LDCRIAUTI,
                     PCRIAUTI,PT0CRIAUTI,PCRIAUTC, PTSTEP_TS, CCSNOWRIMING, PMRSTEP, KMAXITER,
                     LDFEEDBACKT, LDEVLIMIT, LDNULLWETG, LDWETGPOST, LDNULLWETH, LDWETHPOST,
                     PFRACM90, LDCONVHG, CCSUBG_RC_RR_ACCR, CCSUBG_RR_EVAP, CCSUBG_PR_PDF,
                     LDCRFLIMIT, CCFRAC_ICE_ADJUST, PSPLIT_MAXCFL,
                     CCFRAC_ICE_SHALLOW_MF, LDSEDIM_AFTER],
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
                     (numpy.float64, None, IN), #REAL, INTENT (IN) :: PTSTEP_TS
                     (numpy.str, (4, ), IN), #CHARACTER(4), INTENT (IN) :: CCSNOWRIMING
                     (numpy.float64, None, IN), #REAL, INTENT (IN) :: XMRSTEP
                     (numpy.int64, None, IN), #INTEGER, INTENT (IN) :: KMAXITER
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN) :: LDFEEDBACKT
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN) :: LDEVLIMIT
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN) :: LDNULLWETG
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN) :: LDWETGPOST
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN) :: LDNULLWETH
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN) :: LDWETHPOST
                     (numpy.float64, None, IN), #REAL, INTENT (IN) :: PFRACM90
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN) :: LDCONVHG
                     (numpy.str, (80, ), IN), #CHARACTER(LEN=80), INTENT(IN) :: CCSUBG_RC_RR_ACCR
                     (numpy.str, (80, ), IN), #CHARACTER(LEN=80), INTENT(IN) :: CCSUBG_RR_EVAP
                     (numpy.str, (80, ), IN), #CHARACTER(LEN=80), INTENT(IN) :: CCSUBG_RPR_PDF
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN) :: LDCRFLIMIT
                     (numpy.str, (1, ), IN), #CHARACTER(LEN=1), INTENT(IN) :: CCFRAC_ICE_ADJUST
                     (numpy.float64, None, IN), #REAL, INTENT (IN) :: PSPLIT_MAXCFL
                     (numpy.str, (1, ), IN), #CHARACTER(LEN=1), INTENT(IN) :: CCFRAC_ICE_SHALLOW_MF
                     (numpy.bool, None, IN), #LOGICAL, INTENT (IN) :: LDSEDIM_AFTER
                    ],
                    None)
        self._init_py = init_py

        #SEDIM_PY(HSEDIM, IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
        #                   &PTSTEP, KRR,KSPLITR, OMOMENTUM, &
        #                   &PSEA, PTOWN, PDZZ, &
        #                   &PRHODREF, PPABST, PTHT, PRHODJ, &
        #                   &PRCT, PRRT, PRIT, PRST, PRGT, PRHT, &
        #                   &PRCS, PRRS, PRIS, PRSS, PRGS, PRHS, &
        #                   &PINPRC, PINPRR, PINPRI, PINPRS, PINPRG, PINPRH)
        @ctypesFF()
        def sedim_py(HSEDIM, LHAIL, IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL,
                     PTSTEP, KRR,KSPLITR, OMOMENTUM,
                     PSEA, PTOWN, PDZZ,
                     PRHODREF, PPABST, PTHT, PRHODJ,
                     PRCT, PRRT, PRIT, PRST, PRGT, PRHT,
                     PRCS, PRRS, PRIS, PRSS, PRGS, PRHS):
            """
            This function calls the sedim_py fortran subroutine.
            """
            return ([HSEDIM, LHAIL, IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL,
                     PTSTEP, KRR,KSPLITR, OMOMENTUM,
                     PSEA, PTOWN, PDZZ,
                     PRHODREF, PPABST, PTHT, PRHODJ,
                     PRCT, PRRT, PRIT, PRST, PRGT, PRHT,
                     PRCS, PRRS, PRIS, PRSS, PRGS, PRHS],
                    [(numpy.str, (4, ), IN), #CHARACTER(LEN=4), INTENT(IN) :: HSEDIM
                     (numpy.bool, None, IN), #HAIL
                     (numpy.int64, None, IN), (numpy.int64, None, IN), (numpy.int64, None, IN), #INTEGER, INTENT(IN) :: IIB, IIE, IIT
                     (numpy.int64, None, IN), (numpy.int64, None, IN), (numpy.int64, None, IN), #INTEGER, INTENT(IN) :: IJB, IJE, IJT
                     (numpy.int64, None, IN), (numpy.int64, None, IN), (numpy.int64, None, IN), #INTEGER, INTENT(IN) :: IKB, IKE, IKTB
                     (numpy.int64, None, IN), (numpy.int64, None, IN), (numpy.int64, None, IN), #INTEGER, INTENT(IN) :: IKTE, IKT, KKL
                     (numpy.float64, None, IN), #REAL, INTENT(IN) :: PTSTEP
                     (numpy.int64, None, IN), (numpy.int64, None, IN), #INTEGER, INTENT(IN) :: KRR,KSPLITR
                     (numpy.bool, None, IN), #MOMENTUM
                     (numpy.float64, (IIT, IJT), IN), #REAL, DIMENSION(:,:), INTENT(IN)        :: PSEA ! Sea Mask
                     (numpy.float64, (IIT, IJT), IN), #REAL, DIMENSION(:,:), INTENT(IN)        :: PTOWN! Fraction that is town
                     (numpy.float64, (IIT, IJT, IKT), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
                     (numpy.float64, (IIT, IJT, IKT), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
                     (numpy.float64, (IIT, IJT, IKT), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! absolute pressure at t
                     (numpy.float64, (IIT, IJT, IKT), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
                     (numpy.float64, (IIT, IJT, IKT), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
                     (numpy.float64, (IIT, IJT, IKT), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
                     (numpy.float64, (IIT, IJT, IKT), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
                     (numpy.float64, (IIT, IJT, IKT), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
                     (numpy.float64, (IIT, IJT, IKT), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
                     (numpy.float64, (IIT, IJT, IKT), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel m.r. at t
                     (numpy.float64, (IIT, IJT, IKT), IN), #REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHT    ! Hail m.r. at t
                     (numpy.float64, (IIT, IJT, IKT), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
                     (numpy.float64, (IIT, IJT, IKT), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
                     (numpy.float64, (IIT, IJT, IKT), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
                     (numpy.float64, (IIT, IJT, IKT), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
                     (numpy.float64, (IIT, IJT, IKT), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source
                     (numpy.float64, (IIT, IJT, IKT), INOUT), #REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRHS    ! Hail m.r. source
                     (numpy.float64, (IIT, IJT), OUT), #REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRC! Cloud instant precip
                     (numpy.float64, (IIT, IJT), OUT), #REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRR! Rain instant precip
                     (numpy.float64, (IIT, IJT), OUT), #REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRI! Rain instant precip
                     (numpy.float64, (IIT, IJT), OUT), #REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRS! Snow instant precip
                     (numpy.float64, (IIT, IJT), OUT), #REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRG! Graupel instant precip
                     (numpy.float64, (IIT, IJT), OUT), #REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRH! Hail instant precip
                    ],
                    None)
        self._sedim_py = sedim_py

        #We call init_py with ICE3 or ICE4 microphysic scheme and
        #SPLI sedimentation scheme to force KSPLITR computation
        self._KSPLITR = self._init_py(1, self._dt, True,
                                      'ICE4' if self._options['hail'] else 'ICE3',
                                      'SPLI', False, 0., 0., 0.,
                                      0., 'OLD', 0.,
                                      1, True, True, True, True, True, True, 0.1, True,
                                      'NONE'.ljust(80), 'NONE'.ljust(80), 'SIGM'.ljust(80),
                                      True, 'T', self._options.get('maxcfl', 1.), 'T', True)
        #In init min(dz) is set to 20m, below is a correction to take into
        #account the true value of min(dz)
        Z = init_state['Z_half']
        self._KSPLITR = int(self._KSPLITR * 20. / (Z[..., 1:] - Z[..., :-1]).min())

    def finalize(self):
        """
        The method closes the shared library so we are able to load
        a new one with the same symbol names. Moreover, temporary file
        is removed.
        """
        super().finalize()
        del self._init_py, self._sedim_py
        ctypesForFortran.dlclose(self._handle)
        if os.path.exists('fort.1'):
            os.remove('fort.1')

    def build_init_state(self, state):
        """
        This method checks present variables, initializes hail
        and surface flux variables.
        """
        state = super().build_init_state(state)
        needed = ['Z_half', 'T', 'P', 'rc', 'rr', 'ri', 'rs', 'rg']
        if self._options['hail']:
            needed += ['rh']
        for var in needed:
            if var not in state:
                raise ValueError(var + " must be in state")
        if not self._options['hail']:
            if 'rh' in state:
                state['rg'] += state['rh']
                state['rh'] = state['rh'] * 0.
            else:
                state['rh'] = numpy.zeros(state['ri'].shape)
        for surf_flux in ['cum_rc', 'cum_rr', 'cum_ri', 'cum_rs', 'cum_rg', 'cum_rh']:
            if not surf_flux in state:
                state[surf_flux] = numpy.zeros(state['ri'].shape[:-1])
            else:
                assert state[surf_flux].shape == state['ri'].shape[:-1], \
                       "surface flux must have the same shape as " + \
                       "mixing ratio without vertical dimension"
        assert len(set([state[v].shape for v in (needed + ['rh']) if v != 'Z_half'])) == 1, \
               "all variables (except Z_half) must share the same shape"
        shape_z = tuple(list(state['rc'].shape[:-1]) + [state['rc'].shape[-1] + 1])
        assert state['Z_half'].shape == shape_z, \
               "Z_half must have the same shape as other variables " + \
               "except on vertical (must have one more point)"
        return state

    def execute(self, previous_state, timestep, timestep_number):
        """
        execute implementation to call the sedimentation schemes
        """
        super().execute(previous_state, timestep, timestep_number)

        #Dimensions
        PT = previous_state['T']
        shape = PT.shape
        shape_ori = shape
        if len(shape) > 3:
            raise ValueError("Maximum shape length is 3")
        if len(shape) == 1:
            shape = tuple([1, 1] + list(shape))
        elif len(shape) == 2:
            shape = tuple([1] + list(shape))
        assert len(shape) == 3, "Internal error"
        PT.reshape(shape)
        KKL = 1

        #Auxiliary fields (not really needed)
        PRHODJ = numpy.ones(shape)
        PRHODREF = numpy.ones(shape)
        PSEA = numpy.zeros(tuple(shape[:2]))
        PTOWN = numpy.zeros(tuple(shape[:2]))

        #Temperature
        P0 = 100000.
        Boltz = 1.380658E-23
        Avogadro = 6.0221367E+23
        Md = 28.9644E-3
        Rd = Avogadro * Boltz / Md
        Cpd = 7. * Rd / 2.
        PEXN = (previous_state['P'].reshape(shape) / P0) ** (Rd / Cpd)
        PTHT = PT / PEXN

        #Tendencies
        PRCS = previous_state['rc'].reshape(shape) / timestep
        PRRS = previous_state['rr'].reshape(shape) / timestep
        PRIS = previous_state['ri'].reshape(shape) / timestep
        PRSS = previous_state['rs'].reshape(shape) / timestep
        PRGS = previous_state['rg'].reshape(shape) / timestep
        PRHS = previous_state['rh'].reshape(shape) / timestep

        PDZZ = previous_state['Z_half'][..., 1:] - previous_state['Z_half'][..., :-1]
        PDZZ = PDZZ.reshape(shape)
        result = self._sedim_py(self._options['version'], self._options['hail'], 1,
                                shape[0], shape[0], 1, shape[1], shape[1], 1,
                                shape[2], 1, shape[2], shape[2], KKL,
                                timestep, 7 if self._options['hail'] else 6,
                                self._KSPLITR, self._options['version'] == 'SPL2',
                                PSEA, PTOWN, PDZZ,
                                PRHODREF, previous_state['P'].reshape(shape), PTHT, PRHODJ,
                                previous_state['rc'].reshape(shape),
                                previous_state['rr'].reshape(shape),
                                previous_state['ri'].reshape(shape),
                                previous_state['rs'].reshape(shape),
                                previous_state['rg'].reshape(shape),
                                previous_state['rh'].reshape(shape),
                                PRCS, PRRS, PRIS, PRSS, PRGS, PRHS)
        next_state = {}
        inst = {}
        (PRCS, PRRS, PRIS, PRSS, PRGS, PRHS,
         inst['rc'], inst['rr'], inst['ri'],
         inst['rs'], inst['rg'], inst['rh']) = result
        for var in ['rc', 'rr', 'ri', 'rs', 'rg', 'rh']:
            if self._method == 'step-by-step':
                next_state['cum_' + var] = previous_state['cum_' + var] + inst[var] * timestep
            else:
                next_state['cum_' + var] = inst[var] * timestep
        next_state['T'] = previous_state['T'].reshape(shape_ori)
        next_state['rc'] = PRCS.reshape(shape_ori) * timestep
        next_state['rr'] = PRRS.reshape(shape_ori) * timestep
        next_state['ri'] = PRIS.reshape(shape_ori) * timestep
        next_state['rs'] = PRSS.reshape(shape_ori) * timestep
        next_state['rg'] = PRGS.reshape(shape_ori) * timestep
        next_state['rh'] = PRHS.reshape(shape_ori) * timestep

        return next_state
