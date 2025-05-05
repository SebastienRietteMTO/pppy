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

class pppy_turb(pppy.PPPY):
    """
    PPPY implementation for calling the shallow convection scheme
    """

    def __init__(self, dt, method, name, tag, **namel):
        """
        In addition to dt, method, name and tag parameters
        defined in the PPPY class, this parameterization
        needs the following parameters:
        namel: json string representing the namelists to use as a dict of dict
        """
        super().__init__(dt, method, name, tag, **namel)

        self._init = None
        self._close = None
        self._param = None

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
            from pyphyex import PYTURB, PYINI_PHYEX, close
        except ModuleNotFoundError:
            logging.error("The module pyphyex has not been found. " + \
                          "This module can be built with the PHYEX package " + \
                          "(available on github). Details on the compilation " + \
                          "process can be found in the PHYEX documentation.")
            raise

        self._init = PYINI_PHYEX
        self._close = close
        self._param = PYTURB

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
        NVEXT = 1 #for turb!!!
        NKT, NIJT = previous_state['Theta'].shape
        NKT += 2 * NVEXT
        KRR, KRRL, KRRI = 6, 2, 3
        KSV, KSV_LGBEG, KSV_LGEND = 0, 0, 0
        KSV_LIMA_NR, KSV_LIMA_NS, KSV_LIMA_NG, KSV_LIMA_NH = 0, 0, 0, 0

        #Add extra levels
        ps = {var:previous_state[var] for var in ['sfth', 'sfrv', 'sfu', 'sfv']} #2D var untouched
        for var in ['Theta', 'rv', 'rc', 'rr', 'ri', 'rs', 'rg',
                    'u', 'v', 'w', 'tke', 'src', 'flxz_thv_MF', 'flxz_u_MF', 'flxz_v_MF',
                    'Z_flux', 'dzz', 'P']:
            #3D array, copy value in extra levels
            newval = numpy.ndarray((NKT, NIJT))
            newval[1:-1, :] = previous_state[var]
            newval[0, :] = newval[1, :]
            newval[-1, :] = newval[-2, :]
            ps[var] = newval
        #Update Z_flux and P extra values by extrapolation
        for var in ['Z_flux', 'P']:
            ps[var][0, :] = 2 * ps[var][1, :] - ps[var][2, :]
            ps[var][-1, :] = 2 * ps[var][-2, :] - ps[var][-3, :]

        #Other values
        HLBCX = HLBCY = numpy.array(['CYCL', 'CYCL'], dtype=('S', 4))
        KGRADIENTSLEO, KGRADIENTSGOG, KHALO, KSPLIT = 0, 0, 1, 1
        OCLOUDMODIFLM = False
        O2D, ONOMIXLG, OFLAT, OCOUPLES = False, False, False, False
        OBLOWSNOW, OIBM, OFLYER, OCOMPUTE_SRC = False, False, True, True
        PRSNOW = 1.
        OOCEAN, ODEEPOC, ODIAG_IN_RUN = False, False, False
        HTURBLEN_CL, HCLOUD, HELEC = 'DELT', 'ICE3', 'NONE'
        PSVT = numpy.zeros((KSV, NKT, NIJT))
        PSFSV = numpy.zeros((KSV, NIJT))
        PDXX = PDYY = PDZX = PDZY = numpy.ndarray((NKT, NIJT))
        PDIRCOSXW = PDIRCOSYW = PDIRCOSZW = PCOSSLOPE = numpy.ones((NIJT,))
        PSINSLOPE = numpy.zeros((NIJT,))
        PHGRADLEO = numpy.ndarray((KGRADIENTSLEO, NKT, NIJT))
        PHGRADGOG = numpy.ndarray((KGRADIENTSGOG, NKT, NIJT))
        PLENGTHM, PLENGTHH = numpy.zeros((NKT, NIJT)), numpy.zeros((NKT, NIJT))
        MFMOIST = numpy.zeros((NKT, NIJT))
        PCEI = numpy.zeros((NKT, NIJT))
        PCEI_MIN, PCEI_MAX = 0.001E-06, 0.01E-06
        PCOEF_AMPL_SAT = 5.
        KBUDGETS = 12
        KUNIT = 999
        PBL_DEPTH, PSBL_DEPTH = numpy.zeros((NIJT,)), numpy.zeros((NIJT,))
        zs = ps['Z_flux'][0, :]


        exner = (ps['P'] / 1.E5) ** (XRD / XCPD)
        PRHODREF = ps['P'] / ((XRD + ps['rv'] * XRV) * ps['Theta'] * exner)
        rhodj = ps['dzz'] * PRHODREF

        PSVM = numpy.zeros((0, NKT, NIJT))
        PRT = numpy.ndarray((KRR, NKT, NIJT))
        for i, var in enumerate(['rv', 'rc', 'rr', 'ri', 'rs', 'rg']):
            PRT[i, :, :] = ps[var]
        PTHVREF = ps['Theta'] * (1 + ps['rv'] * XRV / XRD) / (1 + PRT.sum(axis=0))
        PRUS = ps['u'] / timestep * rhodj
        PRVS = ps['v'] / timestep * rhodj
        PRWS = ps['w'] / timestep * rhodj
        PRTHLS = ps['Theta'] / timestep * rhodj
        PRRS = PRT / timestep * rhodj
        PRSVS = PSVM / timestep * rhodj
        PRTKES = ps['tke'] / timestep * rhodj

        result = self._param(NIJT, NKT, 1, NVEXT, KRR, KRRL, KRRI, HLBCX, HLBCY, KGRADIENTSLEO, KGRADIENTSGOG,
                             KHALO, KSPLIT, OCLOUDMODIFLM, KSV, KSV_LGBEG, KSV_LGEND,
                             KSV_LIMA_NR, KSV_LIMA_NS, KSV_LIMA_NG, KSV_LIMA_NH,
                             O2D, ONOMIXLG, OFLAT, OCOUPLES, OBLOWSNOW, OIBM, OFLYER,
                             OCOMPUTE_SRC, PRSNOW, OOCEAN, ODEEPOC, ODIAG_IN_RUN,
                             HTURBLEN_CL, HCLOUD, HELEC, timestep, KUNIT,
                             PDXX, PDYY, ps['dzz'], PDZX, PDZY, ps['Z_flux'],
                             PDIRCOSXW, PDIRCOSYW, PDIRCOSZW, PCOSSLOPE, PSINSLOPE,
                             rhodj, PTHVREF, PHGRADLEO, PHGRADGOG, zs, ps['sfth'], ps['sfrv'],
                             PSFSV, ps['sfu'], ps['sfv'], ps['P'], ps['u'], ps['v'], ps['w'],
                             ps['tke'], PSVT, ps['src'], PLENGTHM, PLENGTHH, MFMOIST,
                             PBL_DEPTH, PSBL_DEPTH, PCEI, PCEI_MIN, PCEI_MAX, PCOEF_AMPL_SAT,
                             ps['Theta'], PRT, PRUS, PRVS, PRWS, PRTHLS, PRRS, PRSVS, PRTKES,
                             ps['flxz_thv_MF'], ps['flxz_u_MF'], ps['flxz_v_MF'],
                             KBUDGETS, missingOUT=['PEDR', 'PLEM', 'PDPMF',
                             'PTPMF', 'PDRUS_TURB', 'PDRVS_TURB', 'PDRTHLS_TURB', 'PDRRTS_TURB',
                             'PDRSVS_TURB', 'PTR', 'PDISS', 'PIBM_XMUT'])

        ns = {}
        (_, _, _, _, PRUS, PRVS, PRWS, PRTHLS, PRRS,
         _, PRTKES, PSIGS, _, _, _, _, _, _, _) = result

        ns['Theta'] = PRTHLS[1:-1, :] * timestep / rhodj[1:-1, :]
        for i, var in enumerate(['rv', 'rc', 'rr', 'ri', 'rs', 'rg']):
            ns[var] = PRRS[i, 1:-1, :] * timestep / rhodj[1:-1, :]
        ns['u'] = PRUS[1:-1, :] * timestep / rhodj[1:-1, :]
        ns['v'] = PRVS[1:-1, :] * timestep / rhodj[1:-1, :]
        ns['w'] = PRWS[1:-1, :] * timestep / rhodj[1:-1, :]
        ns['tke'] = PRTKES[1:-1, :] * timestep / rhodj[1:-1, :]
        ns['sigs'] = PSIGS[1:-1, :]

        return ns
