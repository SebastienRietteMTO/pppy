#!/usr/bin/env python3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

#Written by S. Riette

"""
This modules is an example of use of PPPY.
It compares sedimentation schemes of Meso-NH.
"""

import numpy
import os
import matplotlib.pyplot as plt

from pppy import PPPYComp
from pppy_sedim_MNH54 import pppy_sedim_MNH54

#Configuration of comparison
cwd = os.path.dirname(os.path.abspath(__file__))
output_dir = cwd
solib = cwd + "/lib/libMNH54.so"
comp = 1

#Configuration of schemes
sedim_STAT = pppy_sedim_MNH54(solib=solib,
                              dt=60., method='step-by-step',
                              name="Statistical scheme with dt=60.s",
                              tag="STAT_dt=60.",
                              hail=False, version='STAT')
sedim_SPLI = pppy_sedim_MNH54(solib=solib,
                              dt=60., method='step-by-step',
                              name="Old eulerian scheme with dt=60.s",
                              tag="SPLI_dt=60.",
                              hail=False, version='SPLI')
sedim_SPLN = pppy_sedim_MNH54(solib=solib,
                              dt=60., method='step-by-step',
                              name="New eulerian scheme with dt=60.s",
                              tag="SPLN_dt=60.",
                              hail=False, version='SPLN', maxcfl=.8)
sedim_SPL2 = pppy_sedim_MNH54(solib=solib,
                              dt=60., method='step-by-step',
                              name="Eulerian scheme with momentum transport with dt=60.s",
                              tag="SPL2_dt=60.",
                              hail=False, version='SPL2', maxcfl=.8)
sedim_SREF = pppy_sedim_MNH54(solib=solib,
                              dt=60., method='one-step',
                              name="Reference scheme (eulerian scheme with momentum " +
                                   "transport run in one-step mode) with dt=60.s",
                              tag="SREF_dt=60.",
                              hail=False, version='SPL2', maxcfl=.8)

dt_list = [1., 2., 3., 4., 5., 6., 8., 10., 12., 15., 20., 30., 45., 60.]
sedim_STAT2 = [pppy_sedim_MNH54(solib=solib,
                                dt=dt, method='step-by-step',
                                name="Statistical scheme with dt=" + str(dt) + "s",
                                tag="STAT_dt=" + str(dt),
                                hail=False, version='STAT')
               for dt in dt_list]

#Initial state
op_profile = numpy.zeros((20, )) #vertical profile of hydrometeor
op_profile[op_profile.shape[0]-3] = 1 #Only one point with not null content
Z_half = numpy.arange(0, 21 * 50., 50.) #Vertical grid, flux levels
Z_mass = 0.5 * (Z_half[..., 1:] + Z_half[..., :-1]) #Vertical grid, mass levels
init_state = dict(P=numpy.ones(op_profile.shape) * 100000.,
                  T=numpy.ones(op_profile.shape) * 280.,
                  rc=op_profile.copy() * 1.E-4,
                  rr=op_profile.copy() * 1.E-4,
                  ri=op_profile.copy() * 1.E-4,
                  rs=op_profile.copy() * 1.E-4,
                  rg=op_profile.copy() * 1.E-4,
                  Z_half=Z_half, Z_mass=Z_mass)

#Comparison and plots
if comp == 1:
    #Comparison of different schemes
    conf = {
            'schemes': [sedim_STAT, sedim_SPLI, sedim_SPLN, sedim_SPL2, sedim_SREF],
            'output_dir': output_dir,
            'duration': 3600.,
            'init_state': init_state,
            'name': "First sedimentation test",
            'tag': "firstSedimTest"
           }
    comp = PPPYComp(**conf)
    comp.run(force=False)
    figs = []
    levels = [1.E-4, 1.E-3, 5.E-3, 1.E-2, 2.E-2, 3.E-2, 4.E-2, 5.E-2, 6.E-2, 7.E-2, 8.E-2, 9.E-2]

    #rr for the different schemes
    plot1 = 'evol', dict(var_names=['rr'], only_param_tag=[sedim_STAT.tag],
                         y_var='Z_mass', mfactor=1000., contourlevels=levels)
    plot2 = 'evol', dict(var_names=['rr'], only_param_tag=[sedim_SPLI.tag],
                         y_var='Z_mass', mfactor=1000., contourlevels=levels)
    plot3 = 'evol', dict(var_names=['rr'], only_param_tag=[sedim_SPLN.tag],
                         y_var='Z_mass', mfactor=1000., contourlevels=levels)
    plot4 = 'evol', dict(var_names=['rr'], only_param_tag=[sedim_SPL2.tag],
                         y_var='Z_mass', mfactor=1000., contourlevels=levels)
    plot5 = 'evol', dict(var_names=['rr'], only_param_tag=[sedim_SREF.tag],
                         y_var='Z_mass', mfactor=1000., contourlevels=levels)
    fig, plots = comp.plot_multi((5, 1), [plot1, plot2, plot3, plot4, plot5])
    figs.append(fig)

    #cum_rr for the different schemes
    plot = 'evol', dict(var_names=['cum_rr'])
    fig, plots = comp.plot_multi((1, 1), [plot])
    figs.append(fig)

    plt.show()
    for fig in figs:
        plt.close(fig)

elif comp == 2:
    #Comparison of different time step for the same scheme
    conf = {
            'schemes': sedim_STAT2,
            'output_dir': output_dir,
            'duration': 1080.,
            'init_state': init_state,
            'name': "First sedimentation test",
            'tag': "firstSedimTest"
           }

    comp = PPPYComp(**conf)
    comp.run(force=False)
    figs = []
    levels = [1.E-4, 1.E-3, 5.E-3, 1.E-2, 2.E-2, 3.E-2, 4.E-2, 5.E-2, 6.E-2, 7.E-2, 8.E-2, 9.E-2]

    #cum_rr
    plot = 'comp', dict(var_names=['cum_rr'], only_times=[120., 180., 360., 720., 1080.],
                        x_var=dt_list)
    fig, plots = comp.plot_multi((1, 1), [plot])
    figs.append(fig)

    #rr at different times on same plot
    plot = 'comp', dict(var_names=['rr'], only_times=[120., 180., 360., 720.],
                        x_var=dt_list, mfactor=1000., contourlevels=levels)
    fig, plots = comp.plot_multi((1, 1), [plot])
    figs.append(fig)

    #rr at different times on different plots
    plot1 = 'comp', dict(var_names=['rr'], only_times=[120.], y_var='Z_mass',
                         x_var=dt_list, mfactor=1000., contourlevels=levels)
    plot2 = 'comp', dict(var_names=['rr'], only_times=[180.], y_var='Z_mass',
                         x_var=dt_list, mfactor=1000., contourlevels=levels)
    plot3 = 'comp', dict(var_names=['rr'], only_times=[360.], y_var='Z_mass',
                         x_var=dt_list, mfactor=1000., contourlevels=levels)
    plot4 = 'comp', dict(var_names=['rr'], only_times=[720.], y_var='Z_mass',
                         x_var=dt_list, mfactor=1000., contourlevels=levels)
    plot5 = 'comp', dict(var_names=['rr'], only_times=[1080.], y_var='Z_mass',
                         x_var=dt_list, mfactor=1000., contourlevels=levels)
    fig, plots = comp.plot_multi((5, 1), [plot1, plot2, plot3, plot4, plot5])
    figs.append(fig)

    #rr and rs
    plot = 'comp', dict(var_names=['rr', 'rs'], only_times=[120., 180.], y_var='Z_mass',
                        x_var=dt_list, mfactor=1000., contourlevels=levels)
    fig, plots = comp.plot_multi((1, 1), [plot])
    figs.append(fig)

    plt.show()
    for fig in figs:
        plt.close(fig)

