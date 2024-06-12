#!/usr/bin/env python3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

"""
This is an example of comparison involving fake parameterizations
"""

import numpy
import os
import matplotlib.pyplot as plt

from pppy import PPPYComp
from pppy_param1 import pppy_param1
from pppy_param2 import pppy_param2

#Configuration of comparison
cwd = os.path.dirname(os.path.abspath(__file__))
output_dir = cwd
solib = cwd + "/lib/param.so"

#Configuration of schemes
iconf = 0
param_1 = pppy_param1(dt=60., #timestep to use with this parameterization
                      method='step-by-step', #like a true simulation
                      name="Param #1", #name to use for plots
                      tag="param_1", #tag to use for file names
                      solib=solib, #first option to give to pppy_param1: filename of the shared lib
                      iconf=iconf) #second option to give to pppy_param1: configuration
param_2 = pppy_param2(dt=60.,
                      method='step-by-step',
                      name="Param #2",
                      tag="param_3",
                      solib=solib, iconf=iconf)

#Comparison and plots
comp = PPPYComp(schemes=[param_1, param_2], #List of parameterizations to compare
                output_dir=output_dir, #directory in which results are stored (hdf5 files)
                duration=180., #duration of simulation
                init_state=dict(x=numpy.array([0.], order='F')), #initial state
                name="First test", #name to use in plots
                tag="firstTest") #tag to use for file names
comp.run(force=True)
plot = 'evol', dict(var_names=['x'])
fig, plots = comp.plot_multi((1, 1), [plot])
plt.show()
plt.close(fig)