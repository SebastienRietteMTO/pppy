#!/usr/bin/env python3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

"""
This module contains the PPPY implementation for calling the shallow convection scheme
"""

import os
import tempfile
import json
import logging
import numpy
import f90nml
import pppy

# pylint: disable=invalid-name #Disabled because many names are taken from the FORTRAN routines

class pppy_shallow(pppy.PPPY):
    """
    PPPY implementation for calling the shallow convection scheme
    """

    def __init__(self, dt, method, name, tag, dx, dy, **namel):
        """
        In addition to dt, method, name and tag parameters
        defined in the PPPY class, this parameterization
        needs the following parameters:
        dx, dy: horizontal grid size
        namel: json string representing the namelists to use as a dict of dict
        """
        super().__init__(dt, method, name, tag, **namel)

        self._init = None
        self._close = None
        self._param = None
        self._dx, self._dy = dx, dy

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
            from pyphyex import PYSHALLOW_MF, PYINI_PHYEX, close
        except ModuleNotFoundError:
            logging.error("The module pyphyex has not been found. " + \
                          "This module can be built with the PHYEX package " + \
                          "(available on github). Details on the compilation " + \
                          "process can be found in the PHYEX documentation.")
            raise

        self._init = PYINI_PHYEX
        self._close = close
        self._param = PYSHALLOW_MF

        #Setup
        with tempfile.NamedTemporaryFile() as namel:
            nml = f90nml.read(namel.name)
            for k, v in json.loads(self._options['namel']).items():
                nml[k] = v
            nml.write(namel.name, force=True)
            self._init('AROME', 33, namel.name, False, 20, 0, 1, self._dt, 20.,
                       'ICE3', 'EDKF', 'TKEL',
                       LDCHANGEMODEL=True, LDDEFAULTVAL=True, LDREADNAM=True, LDCHECK=True,
                       KPRINT=0, LDINIT=True)
            os.remove('fort.20')

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
        XCPV   = 4.* XRV
        XCL    = 4.218E+3
        XCI    = 2.106E+3
        XTT    = 273.16
        XLVTT  = 2.5008E+6
        XLSTT  = 2.8345E+6

        #Dimensions
        NKT, NIJT = ps['Theta'].shape
        KRR = 6

        #Derived arrays
        exner = (ps['P'] / 1.E5) ** (XRD / XCPD)
        PRHODREF = ps['P'] / ((XRD + ps['rv'] * XRV) * ps['Theta'] * exner)
        rhodj = ps['dzz'] * PRHODREF

        PTHL_UP = numpy.zeros((NKT, NIJT))
        PRT_UP = numpy.zeros((NKT, NIJT))
        PRV_UP = numpy.zeros((NKT, NIJT))
        PRC_UP = numpy.zeros((NKT, NIJT))
        PRI_UP = numpy.zeros((NKT, NIJT))
        PU_UP = numpy.zeros((NKT, NIJT))
        PV_UP = numpy.zeros((NKT, NIJT))
        PTKE_UP = numpy.zeros((NKT, NIJT))
        PTHV_UP = numpy.zeros((NKT, NIJT))
        PW_UP = numpy.zeros((NKT, NIJT))
        PFRAC_UP = numpy.zeros((NKT, NIJT))
        PEMF = numpy.zeros((NKT, NIJT))
        PSVM = numpy.zeros((0, NKT, NIJT))
        PRM = numpy.ndarray((KRR, NKT, NIJT))
        for i, var in enumerate(['rv', 'rc', 'rr', 'ri', 'rs', 'rg']):
            PRM[i, :, :] = ps[var]

        result = self._param(NIJT, NKT, 1, 0, KRR, 2, 3, 0, False, 0, 0,
                             timestep, ps['dzz'], ps['Z_flux'], rhodj, PRHODREF,
                             ps['P'], exner, ps['sfth'], ps['sfrv'], ps['Theta'],
                             PRM, ps['u'], ps['v'], ps['tke'], PSVM,
                             PTHL_UP, PRT_UP, PRV_UP, PRC_UP, PRI_UP, PU_UP,
                             PV_UP, PTKE_UP, PTHV_UP, PW_UP, PFRAC_UP, PEMF,
                             self._dx, self._dy, KBUDGETS=0)

        ns = {}
        (PDUDT_MF, PDVDT_MF, PDTKEDT_MF, PDTHLDT_MF, PDRTDT_MF, PDSVDT_MF,
         ns['sigs_MF'], ns['rc_MF'], ns['ri_MF'], ns['CF_MF'],
         PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF, PWEIGHT_MF_CLOUD,
         PFLXZTHVMF, PFLXZTHMF, PFLXZRMF, PFLXZUMF, PFLXZVMF, PFLXZTKEMF,
         PTHL_UP, PRT_UP, PRV_UP, PRC_UP, PRI_UP, PU_UP, PV_UP, PTKE_UP, PTHV_UP, PW_UP, PFRAC_UP, PEMF,
         PDETR, PENTR, KKLCL, KKETL, KKCTL) = result
        ns['tke'] = ps['tke'] + PDTKEDT_MF *timestep
        ns['Theta'] = ps['Theta'] + PDTHLDT_MF * timestep
        ns['rv'] = ps['rv'] + PDRTDT_MF * timestep
        ns['u'] = ps['u'] + PDUDT_MF * timestep
        ns['v'] = ps['v'] + PDVDT_MF * timestep

        return ns
