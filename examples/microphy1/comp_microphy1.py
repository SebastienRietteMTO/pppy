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
from pppy_microphyAROME_ICE import pppy_microphyAROME_ICE
from pppy_microPhyMNH53_ICE import pppy_microPhyMNH53_ICE
from pppy_microPhyMNH53_LIMA import pppy_microPhyMNH53_LIMA

#Configuration of comparison
cwd = os.path.dirname(os.path.abspath(__file__))
output_dir = cwd
solib_AROME = cwd + "/lib/AROME/libpppyAROME.so"
solib_MNH = cwd + "/lib/MNH/libMNH53.so"

#Configuration of schemes
iv = '4' #4 to use hail in ICE and LIMA, 3 otherwise

dtList = [1, 10, 20, 30]
ice = {dt:pppy_microphyAROME_ICE(solib=solib_AROME,
                                 dt=float(dt), method='step-by-step',
                                 name="ICE" + iv + " scheme with dt=" + str(dt) + ".s",
                                 tag="ICE" + iv + "_dt=" + str(dt) + ".",
                                 ice_version="ICE" + iv, XTSTEP_TS=0., CSNOWRIMING='M90 ',
                                 XMRSTEP=1.E-4, NMAXITER=500, HSUBG_AUCV_RC='PDF ',
                                 HSUBG_RC_RR_ACCR='NONE'.ljust(80), HSUBG_RR_EVAP='NONE'.ljust(80),
                                 HSUBG_PR_PDF='SIGM'.ljust(80), LCRFLIMIT=True)
       for dt in dtList}

old = {dt:pppy_microPhyMNH53_ICE(solib=solib_MNH,
                                 dt=float(dt), method='step-by-step',
                                 name="OLD" + iv + " scheme with dt=" + str(dt) + ".s",
                                 tag="OLD" + iv + "_dt=" + str(dt) + ".",
                                 ice_version="ICE" + iv, HSUBG_AUCV_RC='SIGM')
       for dt in dtList}

lima = {dt:pppy_microPhyMNH53_LIMA(solib=solib_MNH,
                                   dt=float(dt), method='step-by-step',
                                   name="LIMA" + iv + " scheme with dt=" + str(dt) + ".s",
                                   tag="LIMA" + iv+ "_dt=" + str(dt) + ".")
        for dt in dtList}


#Comparison and plots
conf = {
        'schemes' : [ice[1], ice[10], old[1], old[10], lima[1], lima[10]],
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
        'name': "First test",
        'tag': "firstTest"
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

