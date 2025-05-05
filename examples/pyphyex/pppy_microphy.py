#!/usr/bin/env python3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

"""
This module contains the PPPY implementation for the ICE microphysics of PHYEX
"""

import os
import tempfile
import json
import logging
import numpy
import f90nml
import pppy
from ctypesForFortran import MISSING

# pylint: disable=invalid-name #Disabled because many names are taken from the FORTRAN routines

class pppy_microphy(pppy.PPPY):
    """
    PPPY implementation for calling the ICE microphysical scheme of PHYEX.
    """

    def __init__(self, dt, method, name, tag, adjustment=True, **namel):
        """
        In addition to dt, method, name and tag parameters
        defined in the PPPY class, this parameterization
        needs the following parameters:
        namel: json string representing the namelists to use as a dict of dict
        """
        super().__init__(dt, method, name, tag, adjustment=adjustment, **namel)

        self._init = None
        self._close = None
        self._paramICE = None
        self._paramOLD = None
        self._paramADJU = None
        self._full_nml = None

    def setup(self, init_state, duration):
        """
        This method opens the shared library and do the init stuff.
        """
        super().setup(init_state, duration)

        #Opening of shared lib is done here to avoid loading it at initialization
        #Because if shared lib is loaded at init, all declared schemes will have
        #their shared lib loaded simultaneously which can be a source of error.
        #This way, shared lib is loaded in a sub-process (if class instance is used
        #through a PPPYComp instance).

        try:
            from pyphyex import PYRAIN_ICE, PYRAIN_ICE_OLD, PYICE_ADJUST, PYINI_PHYEX, close
        except ModuleNotFoundError:
            logging.error("The module pyphyex has not been found. " + \
                          "This module can be built with the PHYEX package " + \
                          "(available on github). Details on the compilation " + \
                          "process can be found in the PHYEX documentation.")
            raise

        self._init = PYINI_PHYEX
        self._close = close
        self._paramICE = PYRAIN_ICE
        self._paramOLD = PYRAIN_ICE_OLD
        self._paramADJU = PYICE_ADJUST

        #Setup
        with tempfile.NamedTemporaryFile() as namel:
            nml = f90nml.read(namel.name)
            for k, v in json.loads(self._options['namel']).items():
                nml[k] = v
            nml.write(namel.name, force=True)
            #First call to really initialize
            self._init('AROME', 33, namel.name, False, 20, 0, 1, self._dt, 20.,
                       'ICE3', 'EDKF', 'TKEL',
                       LDCHANGEMODEL=True, LDDEFAULTVAL=True, LDREADNAM=True, LDCHECK=True,
                       KPRINT=0, LDINIT=True)
            os.remove('fort.20')
            #Second call to only write the namelist
            numnml = 21
            self._init('AROME', 33, namel.name, False, numnml, 0, 1, self._dt, 20.,
                       'ICE3', 'EDKF', 'TKEL',
                       LDCHANGEMODEL=False, LDDEFAULTVAL=False, LDREADNAM=False, LDCHECK=False,
                       KPRINT=1, LDINIT=False)
            self._full_nml = f90nml.read(f'fort.{numnml}')
            os.remove(f'fort.{numnml}')

    def finalize(self):
        """
        We close the shared library so we are able to load a new one with the same symbol names.
        """
        super().finalize()
        self._close()

    def execute(self, previous_state, timestep, timestep_number):
        """
        This method does the actual work
        """
        super().execute(previous_state, timestep, timestep_number)
        ps = previous_state

        #Thermo
        XG = 9.80665
        XBOLTZ = 1.380658E-23
        XAVOGADRO = 6.0221367E+23
        XMD = 28.9644E-3
        XMV    = 18.0153E-3
        XRD = XAVOGADRO * XBOLTZ / XMD
        XRV    = XAVOGADRO * XBOLTZ / XMV
        XCPD = 7.* XRD / 2.

        #Dimensions
        NKT, NIJT = ps['Theta'].shape
        KRR = 7 if 'rh' in ps else 6

        #Derived arrays
        exner = (ps['P'] / 1.E5) ** (XRD / XCPD)
        PRHODREF = ps['P'] / ((XRD + ps['rv'] * XRV) * ps['Theta'] * exner)
        PMFCONV = numpy.zeros((NKT, NIJT))
        PWEIGHT_MF_CLOUD = numpy.ones((NKT, NIJT))
        PRVS = ps['rv'] / timestep
        PRCS = ps['rc'] / timestep
        PRRS = ps['rr'] / timestep
        PRIS = ps['ri'] / timestep
        PRSS = ps['rs'] / timestep
        PRGS = ps['rg'] / timestep
        PRHS = (ps['rh'] / timestep) if KRR == 7 else MISSING
        PTHS = ps['Theta'] / timestep
        rhodj = ps['dzz'] * PRHODREF
        CF = ps['CF']
        HLC_HRC, HLC_HCF = ps['HLC_HRC'], ps['HLC_HCF']
        HLI_HRI, HLI_HCF = ps['HLI_HRI'], ps['HLI_HCF']

        #Copy
        PTHT, PRVT, PRCT = ps['Theta'], ps['rv'], ps['rc']

        ns = {}
        if self._options['adjustment']:
            PSIGQSAT = numpy.ones((NIJT,)) * self._full_nml['NAM_NEBn']['VSIGQSAT']
            result = self._paramADJU(NIJT, NKT, 1, 0, 6, 'DEPO', timestep, PSIGQSAT, rhodj,
                                     exner, PRHODREF, ps['sigs'], False, PMFCONV,
                                     ps['P'], ps['Z_mass'], exner,
                                     ps['CF_MF'], ps['rc_MF'],
                                     ps['ri_MF'], PWEIGHT_MF_CLOUD, ps['rv'],
                                     ps['rc'], PRVS, PRCS,
                                     ps['Theta'], PTHS,
                                     True, ps['rr'], ps['ri'],
                                     PRIS, ps['rs'], ps['rg'], 0)
            (_, _, _, _, _, _, _, _, ns['src'], CF, _, PRVT, PRCT, PRIT, PTHT,
             HLC_HRC, HLC_HCF, HLI_HRI, HLI_HCF) = result
            ns['CF'] = CF
            ns['HLC_HRC'], ns['HLC_HCF'] = HLC_HRC, HLC_HCF
            ns['HLI_HRI'], ns['HLI_HCF'] = HLI_HRI, HLI_HCF
            PRVS = PRVT / timestep
            PRCS = PRCT / timestep
            PRIS = PRIT / timestep

        inst = {}
        if self._full_nml['NAM_PARAM_ICEn']['LRED']:
            PRT = numpy.ndarray((KRR, NKT, NIJT))
            PRT[0] = PRVT
            PRT[1] = PRCT
            PRT[2] = ps['rr']
            PRT[3] = PRIT
            PRT[4] = ps['rs']
            PRT[5] = ps['rg']
            if KRR == 7:
                PRT[6] = ps['rh']
            PRS = numpy.ndarray((KRR, NKT, NIJT))
            PRS[0] = PRVS
            PRS[1] = PRCS
            PRS[2] = PRRS
            PRS[3] = PRIS
            PRS[4] = PRSS
            PRS[5] = PRGS
            if KRR == 7:
                PRS[6] = PRHS
            result = self._paramICE(NIJT, NKT, 1, 0, False, False, 0., timestep, KRR, exner,
                                    ps['dzz'], rhodj, PRHODREF, exner, ps['P'], ps['ni'], CF,
                                    HLC_HRC, HLC_HCF, HLI_HRI, HLI_HCF,
                                    PTHT, PRT, PTHS, PRS, ps['sigs'], 0,
                                    PSEA=ps['sea'], PTOWN=ps['town'],
                                    missingOUT=['PFPR'] if KRR == 7 else ['PINPRH', 'PFPR'])
            (ns['ci'], _, _, _, _, PTHS, PRS,
             inst['c'], inst['r'], _, inst['s'], inst['g'], _, PRAINFR) = result[:14]
            if KRR == 7:
                inst['rh'] = result[14:]
            PRVS = PRS[0]
            PRCS = PRS[1]
            PRRS = PRS[2]
            PRIS = PRS[3]
            PRSS = PRS[4]
            PRGS = PRS[5]
            if KRR == 7:
                PRHS = PRS[6]
        else:
            GMICRO = numpy.ndarray((NKT, NIJT), dtype=bool)
            GMICRO[...] = True
            PICLDFR = PSSIO = PSSIU = PIFR = PCLDROP = PIFNNC = numpy.ndarray((NKT, NIJT))
            PICENU = PKGN_ACON = PKGN_SBGR = numpy.ndarray((NIJT,))
            result = self._paramOLD(NIJT, NKT, 1, 0, self._full_nml['NAM_PARAM_ICEn']['LSEDIC'],
                                    False, False, False, self._full_nml['NAM_PARAM_ICEn']['CSEDIM'],
                                    self._full_nml['NAM_PARAM_ICEn']['CSUBG_AUCV_RC'],
                                    True, 1, 1, 1,
                                    timestep, KRR, NIJT * NKT, GMICRO, ps['dzz'], rhodj,
                                    PRHODREF, exner, ps['P'], ps['ni'], ps['CF'],
                                    PICLDFR, PSSIO, PSSIU, PIFR, PTHT, PRVT, PRCT, ps['rr'], PRIT,
                                    ps['rs'], ps['rg'], PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,
                                    ps['sigs'], ps['sea'], ps['town'], False, False, PCLDROP, PIFNNC, 0,
                                    PICENU, PKGN_ACON, PKGN_SBGR,
                                    PRHT=ps['rh'] if KRR == 7 else MISSING, PRHS=PRHS,
                                    missingOUT=['PFPR'] if KRR == 7 else ['PINPRH', 'PFPR'])
            (ns['ci'], _, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,
             inst['c'], inst['r'], _, inst['s'], inst['g']) = result[:14]
            if KRR == 7:
                PRHS, inst['rh'] = result[14:]

        ns = {}
        ns['rv'] = PRVS * timestep
        ns['rc'] = PRCS * timestep
        ns['rr'] = PRRS * timestep
        ns['ri'] = PRIS * timestep
        ns['rs'] = PRSS * timestep
        ns['rg'] = PRGS * timestep
        if KRR == 7:
            ns['rh'] = PRHS * timestep
        ns['Theta'] = PTHS  * timestep

        for spe in ['c', 'r', 's', 'g']:
            if self._method == 'step-by-step':
                ns['cum_' + spe] = ps['cum_' + spe] + inst[spe] * timestep
            else:
                ns['cum_' + spe] = inst[spe] * timestep

        return ns
