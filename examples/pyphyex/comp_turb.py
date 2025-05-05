#!/usr/bin/env python3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

"""
This modules is an example of comparison for the sedimentation schemes
"""

import os
import json
import numpy
import matplotlib.pyplot as plt

from pppy import PPPYComp
from pppy_turb import pppy_turb

#Configuration of comparison
output_dir = os.path.dirname(os.path.abspath(__file__))

#Configuration of schemes
turb = pppy_turb(dt=60., method='step-by-step',
                 name="TKEL scheme with dt=60.s",
                 tag="TKEL=60.",
                 namel=json.dumps({'NAM_NEBn': {'LSUBG_COND':True}}))

#Initial state
NIJT, NKT = 1, 100
zmax = 5000.
dzz = numpy.ones((NKT, NIJT)) * zmax / NKT #distance between two flux levels
Z_flux = dzz.cumsum(axis=0) - dzz[0,:]
Z_mass = numpy.ndarray((NKT, NIJT))
Z_mass[:-1, :] = .5*(Z_flux[1:, :]+Z_flux[:-1, :])
Z_mass[-1, :] = 2 * Z_mass[-2, :] - Z_mass[-3, :]
Theta = numpy.array(numpy.where(Z_mass < 2000.,
                                280. - Z_mass * 5. / 2000.,
                                275. + (Z_mass - 2000.) * 0.006))
#Approximative computation of pressure
XG = 9.80665
XBOLTZ = 1.380658E-23
XAVOGADRO = 6.0221367E+23
XMD = 28.9644E-3
XRD = XAVOGADRO * XBOLTZ / XMD
XCPD = 7.* XRD / 2.
Ps = 101325.
P = numpy.ndarray((NKT, NIJT))
P[0, :] = numpy.log(Ps)
for i in range(NKT - 1):
    T = Theta[i, :] * (numpy.exp(P[i, :]) / 1.E5) ** (XRD / XCPD)
    P[i + 1, :] = P[i, :] - XG / (XRD * T) * dzz[i, :]
P = numpy.exp(P)

init_state = dict(P=numpy.ones((NKT, NIJT)) * 101325.,
                  Theta=Theta,
                  rv=numpy.ones((NKT, NIJT)) * 0.001,
                  rc=numpy.zeros((NKT, NIJT)),
                  rr=numpy.zeros((NKT, NIJT)),
                  ri=numpy.zeros((NKT, NIJT)),
                  rs=numpy.zeros((NKT, NIJT)),
                  rg=numpy.zeros((NKT, NIJT)),
                  dzz=dzz,
                  sfth=numpy.ones((NIJT,)) * 100.,
                  sfrv=numpy.ones((NIJT,)) * 100.,
                  sfu=numpy.ones((NIJT,)),
                  sfv=numpy.ones((NIJT,)),
                  u=numpy.zeros((NKT, NIJT)),
                  v=numpy.zeros((NKT, NIJT)),
                  w=numpy.zeros((NKT, NIJT)),
                  tke=numpy.ones((NKT, NIJT)) * 0.00001,
                  src=numpy.zeros((NKT, NIJT)),
                  flxz_thv_MF=numpy.zeros((NKT, NIJT)),
                  flxz_u_MF=numpy.zeros((NKT, NIJT)),
                  flxz_v_MF=numpy.zeros((NKT, NIJT)),
                  sigs=numpy.zeros((NKT, NIJT)),
                  Z_flux=Z_flux,
                  Z_mass=Z_mass,
                 )

#Comparison and plots
conf = {
        'schemes': [turb],
        'output_dir': output_dir,
        'duration': 6000.,
        'init_state': init_state,
        'name': "First turb test",
        'tag': "firstTurbTest"
       }
comp = PPPYComp(**conf)
comp.run(force=True)
plot1 = 'evol', dict(var_names=['Theta'], ylim=(0, 500), y_var='Z_mass')
plot2 = 'evol', dict(var_names=['rv'], ylim=(0, 500), y_var='Z_mass', mfactor=1000.)
fig, plots = comp.plot_multi((1, 2), [plot1, plot2])

plt.show()
plt.close(fig)
