#!/usr/bin/env python3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

"""
This module contains the PPPY implementation for ice_adjust
"""

import os
import tempfile
import json
import logging
import numpy
import f90nml
import pppy

# pylint: disable=invalid-name #Disabled because many names are taken from the FORTRAN routines

class pppy_ice_adjust(pppy.PPPY):
    """
    PPPY implementation for ice_adjust
    """

    def __init__(self, dt, method, name, tag, **namel):
        """
        In addition to dt, method, name and tag parameters
        defined in the PPPY class, this parameterization
        needs the following parameters:
        namel: json string representing the namelists to use as a dict of dict
               namel must contain ['NAM_NEBn']['VSIGQSAT']
        """
        super().__init__(dt, method, name, tag, **namel)

        self._init = None
        self._close = None
        self._param = None
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
            from pyphyex import PYICE_ADJUST, PYINI_PHYEX, close
        except ModuleNotFoundError:
            logging.error("The module pyphyex has not been found. " + \
                          "This module can be built with the PHYEX package " + \
                          "(available on github). Details on the compilation " + \
                          "process can be found in the PHYEX documentation.")
            raise
        self._init = PYINI_PHYEX
        self._close = close
        self._param = PYICE_ADJUST

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

        #Derived arrays
        #VSIGQSAT = json.loads(self._options['namel'])['NAM_NEBn']['VSIGQSAT']
        PSIGQSAT = numpy.ones((NIJT,)) * self._full_nml['NAM_NEBn']['VSIGQSAT']
        exner = (ps['P'] / 1.E5) ** (XRD / XCPD)
        PRHODREF = ps['P'] / ((XRD + ps['rv'] * XRV) * ps['Theta'] * exner)
        PMFCONV = numpy.zeros((NKT, NIJT))
        PRVS = ps['rv'] / timestep
        PRCS = ps['rc'] / timestep
        PRIS = ps['ri'] / timestep
        PTHS = ps['Theta'] / timestep
        rhodj = ps['dzz'] * PRHODREF

        result = self._param(NIJT, NKT, 1, 0, 6, 'DEPO', timestep, PSIGQSAT, rhodj,
                             exner, PRHODREF, ps['sigs'], False, PMFCONV,
                             ps['P'], ps['Z_mass'], exner,
                             ps['CF_MF'], ps['rc_MF'],
                             ps['ri_MF'], ps['rv'],
                             ps['rc'], PRVS, PRCS,
                             ps['Theta'], PTHS,
                             True, ps['rr'], ps['ri'],
                             PRIS, ps['rs'], ps['rg'], 0)

        ns = {}
        (_, _, _, _, _, _, _, _, ns['src'], ns['CF'], _, ns['rv'], ns['rc'], ns['ri'], ns['Theta'],
        ns['HLC_HRC'], ns['HLC_HCF'], ns['HLI_HRI'], ns['HLI_HCF']) = result

        return ns
