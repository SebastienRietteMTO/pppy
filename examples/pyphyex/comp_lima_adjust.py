#!/usr/bin/env python3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

"""
This is an example of comparison for lima_adjust
"""

import os
import json
import numpy
import matplotlib.pyplot as plt

from pppy import PPPYComp
from pppy_lima_adjust import pppy_lima_adjust

# Configuration of comparison
output_dir = os.path.dirname(os.path.abspath(__file__))

# Configuration of schemes
namel = {'NAM_NEBn': {'CCONDENS': 'CB02', 'LSUBG_COND': True, 'LSIGMAS': True}}  # , 'VSIGQSAT':0.02}}
lima_adjust = pppy_lima_adjust(dt=60.,  # timestep to use with this parameterization
                               method='step-by-step',  # like a true simulation
                               name="LIMA adjust",  # name to use for plots
                               tag="lima_adjust",  # tag to use for file names
                               namel=json.dumps(namel))

# Initial state
NIT, NJT, NKT = 1, 1, 1
PRV = numpy.linspace(0.003, 0.01, num=NKT)
PTH = numpy.linspace(280., 295., num=NIT * NJT)
PTH, PRV = numpy.meshgrid(PTH, PRV)
PRV = PRV.reshape((NKT, NJT, NIT))
PTH = PTH.reshape((NKT, NJT, NIT))
init_state = {'Theta': PTH,
              'P': numpy.ones((NKT, NJT, NIT)) * 101325.,
              'rv': PRV,
              'rc': numpy.zeros((NKT, NJT, NIT)),
              'ri': numpy.zeros((NKT, NJT, NIT)),
              'rr': numpy.zeros((NKT, NJT, NIT)),
              'rs': numpy.zeros((NKT, NJT, NIT)),
              'rg': numpy.zeros((NKT, NJT, NIT)),
              'nc': numpy.ones((NKT, NJT, NIT)) * 3.E8,
              'nr': numpy.zeros((NKT, NJT, NIT)),
              'ni': numpy.ones((NKT, NJT, NIT)),
              'ns': numpy.zeros((NKT, NJT, NIT)),
              'ng': numpy.zeros((NKT, NJT, NIT)),
              'ccn1ft': numpy.ones((NKT, NJT, NIT)) * 1.E8,
              'ccn1at': numpy.ones((NKT, NJT, NIT)) * 3.E8,
              'ifn1ft': numpy.ones((NKT, NJT, NIT)) * 1000.,
              'ifn1at': numpy.ones((NKT, NJT, NIT)) * 1000.,
              'CF_MF': numpy.zeros((NKT, NJT, NIT)),
              'rc_MF': numpy.zeros((NKT, NJT, NIT)),
              'ri_MF': numpy.zeros((NKT, NJT, NIT)),
              'dzz': numpy.ones((NKT, NJT, NIT)),
              'sigs': numpy.ones((NKT, NJT, NIT)) * 0.0005,
              'Z_mass': numpy.ones((NKT, NJT, NIT)),
              'src': numpy.zeros((NKT, NJT, NIT)),
              'CFw': numpy.zeros((NKT, NJT, NIT)),
              'CFi': numpy.zeros((NKT, NJT, NIT)),
             }

# Comparison and plots
comp = PPPYComp(schemes=[lima_adjust],  # List of parameterizations to compare
                output_dir=output_dir,  # directory in which results are stored (hdf5 files)
                duration=180.,  # duration of simulation
                init_state=init_state,  # initial state
                name="LIMA test",  # name to use in plots
                tag="LIMA")  # tag to use for file names
comp.run(force=True)
plot_CF = 'evol', dict(var_names=['CFw', 'CFi'])
plot_r = 'evol', dict(var_names=['rv', 'rc', 'ri'])
plot_n = 'evol', dict(var_names=['nc', 'ni'])
fig, plots = comp.plot_multi((1, 3), [plot_CF, plot_r, plot_n])
plt.show()
plt.close(fig)
