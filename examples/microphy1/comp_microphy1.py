#!/usr/bin/env python3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

#Written by S. Riette

"""
This an example of comparison involving microphysics parameterization
schemes from AROME and Meso-NH models.
"""

import numpy
import os
import matplotlib.pyplot as plt

from pppy import PPPYComp
from pppy_ice3ice4 import pppy_ice3ice4
from pppy_microPhyIce_MNH53 import pppy_microPhyIce_MNH53
from pppy_microPhyLima_MNH53 import pppy_microPhyLima_MNH53

#Configuration of comparison
cwd = os.path.dirname(os.path.abspath(__file__))
output_dir = cwd
solib_AROME = cwd + "/lib/AROME/libpppyAROME.so"
solib_MNH = cwd + "/lib/MNH/libMNH53.so"

#Configuration of schemes
iv = '4' #4 to use hail in ICE and LIMA, 3 otherwise

ice_1 = pppy_ice3ice4(solib=solib_AROME,
                      dt=1., method='step-by-step',
                      name="ICE" + iv + " scheme with dt=1.s",
                      tag="ICE" + iv + "_dt=1.",
                      ice_version="ICE" + iv, XTSTEP_TS=0., CSNOWRIMING='M90 ',
                      XMRSTEP=1.E-4, NMAXITER=500, HSUBG_AUCV_RC='PDF ',
                      HSUBG_RC_RR_ACCR='NONE'.ljust(80), HSUBG_RR_EVAP='NONE'.ljust(80),
                      HSUBG_PR_PDF='SIGM'.ljust(80), LCRFLIMIT=True)

ice_10 = pppy_ice3ice4(solib=solib_AROME,
                       dt=10., method='step-by-step',
                       name="ICE" + iv + " scheme with dt=10.s",
                       tag="ICE" + iv + "_dt=10.",
                       ice_version="ICE" + iv, XTSTEP_TS=0., CSNOWRIMING='M90 ',
                       XMRSTEP=1.E-4, NMAXITER=500, HSUBG_AUCV_RC='PDF ',
                       HSUBG_RC_RR_ACCR='NONE'.ljust(80), HSUBG_RR_EVAP='NONE'.ljust(80),
                       HSUBG_PR_PDF='SIGM'.ljust(80), LCRFLIMIT=True)

ice_20 = pppy_ice3ice4(solib=solib_AROME,
                       dt=20., method='step-by-step',
                       name="ICE" + iv + " scheme with dt=20.s",
                       tag="ICE" + iv + "_dt=20.",
                       ice_version="ICE" + iv, XTSTEP_TS=0., CSNOWRIMING='M90 ',
                       XMRSTEP=1.E-4, NMAXITER=500, HSUBG_AUCV_RC='PDF ',
                       HSUBG_RC_RR_ACCR='NONE'.ljust(80), HSUBG_RR_EVAP='NONE'.ljust(80),
                       HSUBG_PR_PDF='SIGM'.ljust(80), LCRFLIMIT=True)

ice_30 = pppy_ice3ice4(solib=solib_AROME,
                       dt=30., method='step-by-step',
                       name="ICE" + iv + " scheme with dt=30.s",
                       tag="ICE" + iv + "_dt=30.",
                       ice_version="ICE" + iv, XTSTEP_TS=0., CSNOWRIMING='M90 ',
                       XMRSTEP=1.E-4, NMAXITER=500, HSUBG_AUCV_RC='PDF ',
                       HSUBG_RC_RR_ACCR='NONE'.ljust(80), HSUBG_RR_EVAP='NONE'.ljust(80),
                       HSUBG_PR_PDF='SIGM'.ljust(80), LCRFLIMIT=True)

old_1 = pppy_microPhyIce_MNH53(solib=solib_MNH,
                               dt=1., method='step-by-step',
                               name="OLD" + iv + " scheme with dt=1.s",
                               tag="OLD" + iv + "_dt=1.",
                               ice_version="ICE" + iv, HSUBG_AUCV_RC='SIGM')

old_10 = pppy_microPhyIce_MNH53(solib=solib_MNH,
                                dt=10., method='step-by-step',
                                name="OLD" + iv + " scheme with dt=10.s",
                                tag="OLD" + iv + "_dt=10.",
                                ice_version="ICE" + iv, HSUBG_AUCV_RC='SIGM')

old_20 = pppy_microPhyIce_MNH53(solib=solib_MNH,
                                dt=20., method='step-by-step',
                                name="OLD" + iv + " scheme with dt=20.s",
                                tag="OLD" + iv + "_dt=20.",
                                ice_version="ICE" + iv, HSUBG_AUCV_RC='SIGM')

old_30 = pppy_microPhyIce_MNH53(solib=solib_MNH,
                                dt=30., method='step-by-step',
                                name="OLD" + iv + " scheme with dt=30.s",
                                tag="OLD" + iv + "_dt=30.",
                                ice_version="ICE" + iv, HSUBG_AUCV_RC='SIGM')

lima_1 = pppy_microPhyLima_MNH53(solib=solib_MNH,
                                dt=1., method='step-by-step',
                                name="LIMA" + iv + " scheme with dt=1.s",
                                tag="LIMA" + iv+ "_dt=1.")

lima_10 = pppy_microPhyLima_MNH53(solib=solib_MNH,
                                 dt=10., method='step-by-step',
                                 name="LIMA" + iv + " scheme with dt=10.s",
                                 tag="LIMA" + iv+ "_dt=10.")

lima_20 = pppy_microPhyLima_MNH53(solib=solib_MNH,
                                 dt=20., method='step-by-step',
                                 name="LIMA" + iv + " scheme with dt=20.s",
                                 tag="LIMA" + iv+ "_dt=20.")

lima_30 = pppy_microPhyLima_MNH53(solib=solib_MNH,
                                 dt=30., method='step-by-step',
                                 name="LIMA" + iv + " scheme with dt=30.s",
                                 tag="LIMA" + iv+ "_dt=30.")

#Comparison and plots
conf = {
        'schemes' : [ice_1, ice_10, old_1, old_10, lima_1, lima_10],
        'output_dir': output_dir,
        'duration': 180.,
        'init_state': dict(P=numpy.array([100000.]), T=numpy.array([270.]),
                           rv=numpy.array([1.E-2]), rc=numpy.array([1.E-3]),
                           rr=numpy.array([1.E-3]), ri=numpy.array([1.E-3]),
                           rs=numpy.array([1.E-3]), rg=numpy.array([1.E-3]),
                           rh=numpy.array([0.]),
                           ccn1ft=numpy.array([1.E12]), ccn1at=numpy.array([0.]),
                           ifn1ft=numpy.array([0.]), ifn1at=numpy.array([0.]),
                           nc=numpy.array([300.E6]), nr=numpy.array([3572.29]),
                           ni=numpy.array([0.1E6])),
        'experiment_name': "First test",
        'experiment_tag': "firstTest"
       }
comp = PPPYComp(**conf)
comp.run(force=False)
r_all = ['rv', 'rc', 'rr', 'ri', 'rs', 'rg'] + (['rh'] if iv == "4" else [])
figs = []
plot1 = 'evol', dict(var_names=r_all)
plot2 = 'evol', dict(var_names=['T'])
plot3 = 'evol', dict(var_names=['nc', 'nr', 'ni'])
fig, plots = comp.plot_multi((3, 1), [plot1, plot2, plot3])
figs.append(fig)
plots = [('evol', dict(var_names=[s])) for s in r_all]
fig, plots = comp.plot_multi((3, 3) if iv == "4" else (3, 2), plots)
figs.append(fig)
plt.show()
for fig in figs:
    plt.close(fig)

