#!/usr/bin/env python3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

"""
This is an example of comparison for ice_adjust
"""

import os
import json
import numpy
import matplotlib.pyplot as plt

from pppy import PPPYComp
from pppy_ice_adjust import pppy_ice_adjust

#Configuration of comparison
output_dir = os.path.dirname(os.path.abspath(__file__))

#Configuration of schemes
namel = {'NAM_NEBn':{'CCONDENS':'CB02', 'LSUBG_COND':True, 'LSIGMAS':True}}#, 'VSIGQSAT':0.02}}
ice_adjust_CB02 = pppy_ice_adjust(dt=60., #timestep to use with this parameterization
                                  method='step-by-step', #like a true simulation
                                  name="Ice adjust with CB02 option", #name to use for plots
                                  tag="ice_adjust_CB02", #tag to use for file names
                                  namel=json.dumps(namel))
namel = {'NAM_NEBn':{'CCONDENS':'GAUS', 'LSUBG_COND':True, 'LSIGMAS':True, 'VSIGQSAT':0.02}}
ice_adjust_GAUS = pppy_ice_adjust(dt=60., #timestep to use with this parameterization
                                  method='step-by-step', #like a true simulation
                                  name="Ice adjust with a Gaussian", #name to use for plots
                                  tag="ice_adjust_GAUS", #tag to use for file names
                                  namel=json.dumps(namel))
#Initial state
NIJT, NKT = 1, 1
PRV = numpy.linspace(0.003, 0.01, num=NKT)
PTH = numpy.linspace(280., 295., num=NIJT)
PTH, PRV = numpy.meshgrid(PTH, PRV)
init_state = {'Theta': PTH,
              'P': numpy.ones((NKT, NIJT)) * 101325.,
              'rv': PRV,
              'rc': numpy.zeros((NKT, NIJT)),
              'ri': numpy.zeros((NKT, NIJT)),
              'rr': numpy.zeros((NKT, NIJT)),
              'rs': numpy.zeros((NKT, NIJT)),
              'rg': numpy.zeros((NKT, NIJT)),
              'CF_MF': numpy.zeros((NKT, NIJT)),
              'rc_MF': numpy.zeros((NKT, NIJT)),
              'ri_MF': numpy.zeros((NKT, NIJT)),
              'dzz': numpy.ones((NKT, NIJT)),
              'sigs': numpy.ones((NKT, NIJT)) * 0.0005,
              'Z_mass': numpy.ones((NKT, NIJT)),
              'src': numpy.zeros((NKT, NIJT)),
              'CF': numpy.zeros((NKT, NIJT)),
              'HLC_HRC': numpy.zeros((NKT, NIJT)),
              'HLC_HCF': numpy.zeros((NKT, NIJT)),
              'HLI_HRI': numpy.zeros((NKT, NIJT)),
              'HLI_HCF': numpy.zeros((NKT, NIJT)),
             }

#Comparison and plots
comp = PPPYComp(schemes=[ice_adjust_CB02, ice_adjust_GAUS], #List of parameterizations to compare
                output_dir=output_dir, #directory in which results are stored (hdf5 files)
                duration=180., #duration of simulation
                init_state=init_state, #initial state
                name="CB02 and GUASS comparison", #name to use in plots
                tag="CB02_GAUS") #tag to use for file names
comp.run(force=True)
plot = 'evol', dict(var_names=['CF'])
fig, plots = comp.plot_multi((1, 1), [plot])
plt.show()
plt.close(fig)
