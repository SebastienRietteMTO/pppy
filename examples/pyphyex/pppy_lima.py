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

class pppy_lima(pppy.PPPY):
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
        self._param = None
        self._paramADJU = None
        self._full_nml = None

    def setup(self, init_state, duration):
        """
        This method opens the shared library and do the init stuff.
        """
        super().setup(init_state, duration)

        # Opening of shared lib is done here to avoid loading it at initialization
        # Because if shared lib is loaded at init, all declared schemes will have
        # their shared lib loaded simultaneously which can be a source of error.
        # This way, shared lib is loaded in a sub-process (if class instance is used
        # through a PPPYComp instance).

        try:
            from pyphyex import PYLIMA, PYLIMA_ADJUST_SPLIT, PYINI_PHYEX, close
        except ModuleNotFoundError:
            logging.error("The module pyphyex has not been found. " +
                          "This module can be built with the PHYEX package " +
                          "(available on github). Details on the compilation " +
                          "process can be found in the PHYEX documentation.")
            raise

        self._init = PYINI_PHYEX
        self._close = close
        self._param = PYLIMA
        self._paramADJU = PYLIMA_ADJUST_SPLIT

        # Setup
        with tempfile.NamedTemporaryFile() as namel:
            nml = f90nml.read(namel.name)
            for k, v in json.loads(self._options['namel']).items():
                nml[k] = v
            nml.write(namel.name, force=True)
            # First call to really initialize
            self._init('AROME', 33, namel.name, False, 20, 0, 1, self._dt, 20.,
                       'LIMA', 'EDKF', 'TKEL',
                       LDCHANGEMODEL=True, LDDEFAULTVAL=True, LDREADNAM=True, LDCHECK=True,
                       KPRINT=0, LDINIT=True)
            os.remove('fort.20')
            # Second call to only write the namelist
            numnml = 21
            self._init('AROME', 33, namel.name, False, numnml, 0, 1, self._dt, 20.,
                       'LIMA', 'EDKF', 'TKEL',
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

        # Thermo
        XG = 9.80665
        XBOLTZ = 1.380658E-23
        XAVOGADRO = 6.0221367E+23
        XMD = 28.9644E-3
        XMV = 18.0153E-3
        XRD = XAVOGADRO * XBOLTZ / XMD
        XRV = XAVOGADRO * XBOLTZ / XMV
        XCPD = 7. * XRD / 2.

        # Dimensions
        NKT, NJT, NIT = ps['Theta'].shape
        NSV = 7
        KRR = 7 if 'rh' in ps else 6

        # Derived arrays
        exner = (ps['P'] / 1.E5) ** (XRD / XCPD)
        PRHODREF = ps['P'] / ((XRD + ps['rv'] * XRV) * ps['Theta'] * exner)
        PMFCONV = numpy.zeros((NKT, NJT, NIT))
        PDTHRAD = numpy.zeros((NKT, NJT, NIT))
        PW_NU = numpy.zeros((NKT, NJT, NIT))
        PRT = numpy.zeros((KRR, NKT, NJT, NIT))
        PRT[0, :, :, :] = ps['rv']
        PRT[1, :, :, :] = ps['rc']
        PRT[2, :, :, :] = ps['rr']
        PRT[3, :, :, :] = ps['ri']
        PRT[4, :, :, :] = ps['rs']
        PRT[5, :, :, :] = ps['rg']
        if KRR == 7:
            PRT[6, :, :, :] = ps['rh']
        PRS = PRT / timestep
        PTHS = ps['Theta'] / timestep
        rhodj = ps['dzz'] * PRHODREF
        PSVT = numpy.zeros((NSV, NKT, NJT, NIT))
        isv = 0
        if self._full_nml['NAM_PARAM_LIMA']['NMOM_C'] >= 2:
            PSVT[isv, :, :, :] = ps.get('nc', numpy.zeros((NKT, NJT, NIT)))
            isv += 1
        if self._full_nml['NAM_PARAM_LIMA']['NMOM_R'] >= 2:
            PSVT[isv, :, :, :] = ps.get('nr', numpy.zeros((NKT, NJT, NIT)))
            isv += 1
        NCCN = self._full_nml['NAM_PARAM_LIMA']['NMOD_CCN']
        if NCCN > 0:
            for imodccn in range(NCCN):
                PSVT[isv, :, :, :] = ps.get(f'ccn{imodccn+1}ft', numpy.zeros((NKT, NJT, NIT)))
                isv += 1
            for imodccn in range(NCCN):
                PSVT[isv, :, :, :] = ps.get(f'ccn{imodccn+1}at', numpy.zeros((NKT, NJT, NIT)))
                isv += 1
        if self._full_nml['NAM_PARAM_LIMA']['LSCAV'] and \
           self._full_nml['NAM_PARAM_LIMA']['LAERO_MASS']:
            raise NotImplementedError('Scavenging')
            PSVT[isv, :, :, :] = X
            isv += 1
        if self._full_nml['NAM_PARAM_LIMA']['NMOM_I'] >= 2:
            PSVT[isv, :, :, :] = ps.get('ni', numpy.zeros((NKT, NJT, NIT)))
            isv += 1
        if self._full_nml['NAM_PARAM_LIMA']['NMOM_S'] >= 2:
            PSVT[isv, :, :, :] = ps.get('ns', numpy.zeros((NKT, NJT, NIT)))
            isv += 1
        if self._full_nml['NAM_PARAM_LIMA']['NMOM_G'] >= 2:
            PSVT[isv, :, :, :] = ps.get('ng', numpy.zeros((NKT, NJT, NIT)))
            isv += 1
        if self._full_nml['NAM_PARAM_LIMA']['NMOM_H'] >= 2:
            PSVT[isv, :, :, :] = ps.get('nh', numpy.zeros((NKT, NJT, NIT)))
            isv += 1
        NIFN = self._full_nml['NAM_PARAM_LIMA']['NMOD_IFN']
        if NIFN > 0:
            for imodccn in range(NIFN):
                PSVT[isv, :, :, :] = ps.get(f'ifn{imodccn+1}ft', numpy.zeros((NKT, NJT, NIT)))
                isv += 1
            for imodccn in range(NIFN):
                PSVT[isv, :, :, :] = ps.get(f'ifn{imodccn+1}at', numpy.zeros((NKT, NJT, NIT)))
                isv += 1
        NIMM = self._full_nml['NAM_PARAM_LIMA']['NMOD_IMM']
        if NIMM > 0:
            raise NotImplementedError('IMM')
            for imodimm in range(NIMM):
                PSVT[isv, :, :, :] = X
                isv += 1
        if self._full_nml['NAM_PARAM_LIMA']['LHHONI']:
            raise NotImplementedError('Homogeneous freezing of CCN')
            PSVT[isv, :, :, :] = X
            isv += 1
        if self._full_nml['NAM_PARAM_LIMA']['LSPRO']:
            raise NotImplementedError('Supersaturation')
            PSVT[isv, :, :, :] = X
            isv += 1
        PSVS = PSVT / timestep

        ns = {}
        if self._options['adjustment']:
            PSIGQSAT = self._full_nml['NAM_NEBn']['VSIGQSAT']
            HCONDENS = self._full_nml['NAM_NEBn']['CCONDENS']
            HLAMBDA3 = self._full_nml['NAM_NEBn']['CLAMBDA3']
            OSUBG_COND = self._full_nml['NAM_NEBn']['LSUBG_COND']
            OSIGMAS = self._full_nml['NAM_NEBn']['LSIGMAS']
            result = self._paramADJU(NIT * NJT, NKT, 1, 0, 0, KRR, 1,
                                     HCONDENS, HLAMBDA3,
                                     OSUBG_COND, OSIGMAS, timestep, PSIGQSAT,
                                     PRHODREF, rhodj, exner, ps['sigs'], True, PMFCONV,
                                     ps['P'], ps['P'], ps['Z_mass'], True, PDTHRAD, PW_NU,
                                     PRT, PRS, PSVT, PSVS, PTHS, True, ps['CFw'], ps['CFi'],
                                     ps['rc_MF'], ps['ri_MF'], ps['CF_MF'])
            PRS, _, PTHS, ns['src'], ns['CFw'], ns['CFi'] = result
            PTHT = PTHS * timestep
            PRT = PRS * timestep
            PCLDFR, PICEFR = ns['CFw'], ns['CFi']
        else:
            PTHT = ps['Theta']
            PCLDFR, PICEFR = ps['CFw'], ps['CFi']

        inst = {}
        PTHVREFZIKB = 300.
        PPRCFR = numpy.zeros((NKT, NJT, NIT))
        PFPR = numpy.zeros((KRR, NKT, NJT, NIT))
        result = self._param(NIT * NJT, NKT, 1, 0, 0, KRR, timestep, False,
                             PRHODREF, exner, ps['dzz'], PTHVREFZIKB, rhodj,
                             ps['P'], NCCN, NIFN, NIMM, True, PDTHRAD,
                             PTHT, PRT, PSVT, PW_NU, PTHS, PRS, PSVS,
                             PCLDFR, PICEFR, PPRCFR, PFPR)
        (PTHS, PRS, PSVS, inst['c'], _, inst['r'], inst['i'], inst['s'],
         inst['g'], inst['h'], _, ns['CFw'], ns['CFi'], _, _) = result

        ns['rv'] = PRS[0, :, :, :] * timestep
        ns['rc'] = PRS[1, :, :, :] * timestep
        ns['rr'] = PRS[2, :, :, :] * timestep
        ns['ri'] = PRS[3, :, :, :] * timestep
        ns['rs'] = PRS[4, :, :, :] * timestep
        ns['rg'] = PRS[5, :, :, :] * timestep
        if KRR == 7:
            ns['rh'] = PRS[6, :, :, :] * timestep
        ns['Theta'] = PTHS * timestep

        isv = 0
        if self._full_nml['NAM_PARAM_LIMA']['NMOM_C'] >= 2:
            ns['nc'] = PSVS[isv, :, :, :] * timestep
            isv += 1
        if self._full_nml['NAM_PARAM_LIMA']['NMOM_R'] >= 2:
            ns['nr'] = PSVS[isv, :, :, :] * timestep
            isv += 1
        if NCCN > 0:
            for imodccn in range(NCCN):
                ns[f'ccn{imodccn+1}ft'] = PSVS[isv, :, :, :] * timestep
                isv += 1
            for imodccn in range(NCCN):
                ns[f'ccn{imodccn+1}at'] = PSVS[isv, :, :, :] * timestep
                isv += 1
        if self._full_nml['NAM_PARAM_LIMA']['LSCAV'] and \
           self._full_nml['NAM_PARAM_LIMA']['LAERO_MASS']:
            raise NotImplementedError('Scavenging')
            X = PSVS[isv, :, :, :] * timestep
            isv += 1
        if self._full_nml['NAM_PARAM_LIMA']['NMOM_I'] >= 2:
            ns['ni'] = PSVS[isv, :, :, :] * timestep
            isv += 1
        if self._full_nml['NAM_PARAM_LIMA']['NMOM_S'] >= 2:
            ns['ns'] = PSVS[isv, :, :, :] * timestep
            isv += 1
        if self._full_nml['NAM_PARAM_LIMA']['NMOM_G'] >= 2:
            ns['ng'] = PSVS[isv, :, :, :] * timestep
            isv += 1
        if self._full_nml['NAM_PARAM_LIMA']['NMOM_H'] >= 2:
            ns['nh'] = PSVS[isv, :, :, :] * timestep
            isv += 1
        if NIFN > 0:
            for imodccn in range(NIFN):
                ns[f'ifn{imodccn+1}ft'] = PSVS[isv, :, :, :] * timestep
                isv += 1
            for imodccn in range(NIFN):
                ns[f'ifn{imodccn+1}at'] = PSVS[isv, :, :, :] * timestep
                isv += 1
        if NIMM > 0:
            raise NotImplementedError('IMM')
            for imodimm in range(NIMM):
                X = PSVS[isv, :, :, :] * timestep
                isv += 1
        if self._full_nml['NAM_PARAM_LIMA']['LHHONI']:
            raise NotImplementedError('Homogeneous freezing of CCN')
            X = PSVS[isv, :, :, :] * timestep
            isv += 1
        if self._full_nml['NAM_PARAM_LIMA']['LSPRO']:
            raise NotImplementedError('Supersaturation')
            X = PSVS[isv, :, :, :] * timestep
            isv += 1

        for spe in ['c', 'r', 's', 'g']:
            if self._method == 'step-by-step':
                ns['cum_' + spe] = ps['cum_' + spe] + inst[spe] * timestep
            else:
                ns['cum_' + spe] = inst[spe] * timestep

        return ns
