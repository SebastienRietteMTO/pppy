#!/usr/bin/env python3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

#Written by S. Riette

"""
This an example of comparison involving microphysics parameterization
schemes from Meso-NH and WRF models.
"""

import numpy
import os
import matplotlib.pyplot as plt

import logging
logging.basicConfig(level=getattr(logging, 'DEBUG', None))

from pppy import PPPYComp, VAR_NAME, VAR_UNIT
from pppy_microPhyIce_MNH54 import pppy_microPhyIce_MNH54
from pppy_microPhyLima_MNH54 import pppy_microPhyLima_MNH54
from pppy_microPhyWRF import pppy_microPhyWRF

#Configuration of comparison
cwd = os.path.dirname(os.path.abspath(__file__))
output_dir = cwd
solib_MNH = cwd + "/lib/MNH/libMNH54.so"
solib_WRF = cwd + "/lib/WRF/libWRF.so"

#Configuration of schemes
dtList = [0.001, 0.01, 0.1, 1., 5., 10., 15., 20.]
dtList = [1., 5., 10., 15.]
iv = '4' #4 to use hail in ICE and LIMA, 3 otherwise

old = {dt:pppy_microPhyIce_MNH54(solib=solib_MNH,
                                 dt=dt, method='step-by-step',
                                 name="old ICE" + iv + " scheme with dt=" + str(dt) +"s",
                                 tag="OLD" + iv + "_dt=" + str(dt),
                                 ice_version="OLD" + iv, HSUBG_AUCV_RC='PDF ',
                                 HSUBG_PR_PDF='SIGM'.ljust(80))
          for dt in dtList}

ice = {dt:pppy_microPhyIce_MNH54(solib=solib_MNH,
                                 dt=dt, method='step-by-step',
                                 name="new ICE" + iv + " scheme with dt=" + str(dt) +"s",
                                 tag="ICE" + iv + "_dt=" + str(dt),
                                 ice_version="ICE" + iv, HSUBG_AUCV_RC='PDF ',
                                 HSUBG_PR_PDF='SIGM'.ljust(80),
                                 maxiter=100, mrstep=5.E-5, tstep_ts=0., frac_ice_adjust='S')
          for dt in dtList}

lima = {dt:pppy_microPhyLima_MNH54(solib=solib_MNH,
                                   dt=dt, method='step-by-step',
                                   name="LIMA" + iv + " scheme with dt=" + str(dt) +"s",
                                   tag="LIMA_" + iv+ "_dt=" + str(dt),
                                   LHAIL=iv == '4')
           for dt in dtList}

WRF_Kessler = {dt:pppy_microPhyWRF(solib=solib_WRF,
                                   dt=dt, method='step-by-step',
                                   name="Kessler scheme with dt=" + str(dt) +"s",
                                   tag="Kessler_dt=" + str(dt),
                                   mp_physics=pppy_microPhyWRF.KESSLERSCHEME)
                  for dt in dtList}

WRF_Thompson_wo_aerosols = {dt:pppy_microPhyWRF(solib=solib_WRF,
                                                dt=dt, method='step-by-step',
                                                name="Thompson scheme without aerosols with dt=" +
                                                     str(dt) + "s",
                                                tag="Thompson_wo_aerosols_dt=" + str(dt),
                                                mp_physics=pppy_microPhyWRF.THOMPSON)
                               for dt in dtList}

WRF_Thompson_with_aerosols = {dt:pppy_microPhyWRF(solib=solib_WRF,
                                                  dt=dt, method='step-by-step',
                                                  name="Thompson scheme with aerosols with dt=" +
                                                       str(dt) + "s",
                                                  tag="Thompson_with_aerosols_dt=" + str(dt),
                                                  mp_physics=pppy_microPhyWRF.THOMPSONAERO,
                                                  dx=1000., dy=1000.,
                                                  enableCCNsource=False)
                                 for dt in dtList}

WRF_Full_SBM = {dt:pppy_microPhyWRF(solib=solib_WRF,
                                    dt=dt, method='step-by-step',
                                    name="Full SBM scheme with dt=" + str(dt) +"s",
                                    tag="FullSBM_dt=" + str(dt),
                                    mp_physics=pppy_microPhyWRF.FULL_KHAIN_LYNN,
                                    dx=1000., dy=1000.,
                                    enableCCNinit=False)
                   for dt in dtList}

scheme_names = {'old_ICE':old, 'new_ICE':ice, 'lima':lima,
                'WRF_Kessler':WRF_Kessler,
                'WRF_Thompson_wo_aerosols':WRF_Thompson_wo_aerosols,
                'WRF_Thompson_with_aerosols':WRF_Thompson_with_aerosols,
                'WRF_Full_SBM':WRF_Full_SBM}
long_names = {'old_ICE':'old ICE', 'new_ICE':'new ICE', 'lima':'lima',
              'WRF_Kessler':'WRF Kessler',
              'WRF_Thompson_wo_aerosols':'WRF Thompson without aerosols',
              'WRF_Thompson_with_aerosols':'WRF Thompson with aerosols',
              'WRF_Full_SBM':'WRF Full SBM'}
keys = sorted(scheme_names.keys(), key=str.lower)
scheme_list_all = [scheme[dt] for dt in dtList for scheme in scheme_names.values()]
scheme_list_dt = {dt:[scheme_names[key][dt] for key in keys] for dt in dtList}
scheme_list_sch = {scheme:[scheme_names[scheme][dt] for dt in dtList] for scheme in scheme_names.keys()}

#Comparison and plots
warm_state = dict(P=numpy.array([100000.]), T=numpy.array([290.]),
                  rv=numpy.array([1.E-2]), rc=numpy.array([1.E-2]),
                  rr=numpy.array([1.E-4]), ri=numpy.array([0.]),
                  ric =numpy.array([0.]), rid=numpy.array([0.]), rip=numpy.array([0.]),
                  rs=numpy.array([0.]), rg=numpy.array([0.]), rh=numpy.array([0.]),
                  nc=numpy.array([3.E8]), nr=numpy.array([2000.]),
                  ni=numpy.array([0.]), ns=numpy.array([0.]),
                  ng=numpy.array([0.]), nh=numpy.array([0.]),
                  nic=numpy.array([0.]), nid=numpy.array([0.]), nip=numpy.array([0.]),
                  ccn1ft=numpy.array([1.E8]), ccn1at=numpy.array([1.E8]), #CCN free/activated
                  ifn1ft=numpy.array([0.]), ifn1at=numpy.array([0.]), #IFN free/activated
                  u=numpy.array([0.]), v=numpy.array([0.]), w=numpy.array([0.]),
                  xland=numpy.array([1.]))
cold_state = dict(P=numpy.array([100000.]), T=numpy.array([270.]),
                  rv=numpy.array([0.007]), rc=numpy.array([1.E-3]),
                  rr=numpy.array([1.E-3]), ri=numpy.array([9.E-4]),
                  ric =numpy.array([3.E-4]), rid=numpy.array([3.E-4]), rip=numpy.array([3.E-4]),
                  rs=numpy.array([1.E-3]), rg=numpy.array([1.E-3]), rh=numpy.array([0.]),
                  nc=numpy.array([3.E8]), nr=numpy.array([2000.]),
                  ni=numpy.array([30000.]), ns=numpy.array([21000.]),
                  ng=numpy.array([8500.]), nh=numpy.array([0.]),
                  nic=numpy.array([10000.]), nid=numpy.array([10000.]), nip=numpy.array([10000.]),
                  ccn1ft=numpy.array([1.E8]), ccn1at=numpy.array([3.E8]), #CCN free/activated
                  ifn1ft=numpy.array([1000.]), ifn1at=numpy.array([1000.]), #IFN free/activated
                  u=numpy.array([0.]), v=numpy.array([0.]), w=numpy.array([0.]),
                  xland=numpy.array([1.]))
cold2_state = dict(P=numpy.array([100000.]), T=numpy.array([270.]),
                   rv=numpy.array([0.004]), rc=numpy.array([1.E-3]),
                   rr=numpy.array([1.E-3]), ri=numpy.array([9.E-4]),
                   ric =numpy.array([3.E-4]), rid=numpy.array([3.E-4]), rip=numpy.array([3.E-4]),
                   rs=numpy.array([1.E-3]), rg=numpy.array([1.E-3]), rh=numpy.array([0.]),
                   nc=numpy.array([3.E8]), nr=numpy.array([2000.]),
                   ni=numpy.array([30000.]), ns=numpy.array([21000.]),
                   ng=numpy.array([8500.]), nh=numpy.array([0.]),
                   nic=numpy.array([10000.]), nid=numpy.array([10000.]), nip=numpy.array([10000.]),
                   ccn1ft=numpy.array([1.E8]), ccn1at=numpy.array([3.E8]), #CCN free/activated
                   ifn1ft=numpy.array([1000.]), ifn1at=numpy.array([1000.]), #IFN free/activated
                   u=numpy.array([0.]), v=numpy.array([0.]), w=numpy.array([0.]),
                   xland=numpy.array([1.]))
state = 'cold'
comp_name = "Test with MNH and WRF, " + state + " state"
conf = {
        'schemes' : scheme_list_all,
        'output_dir': output_dir,
        'duration': 100.,
        'init_state': {'cold':cold_state, 'cold2':cold2_state, 'warm':warm_state}[state],
        'experiment_name': comp_name,
        'experiment_tag': "test_MNH_WRF_" + state
       }

#Run all simulations
comp = PPPYComp(**conf)
comp.run(force=False)

#Plots by timestep
for dt in dtList:
    conf['schemes'] = scheme_list_dt[dt]
    conf['experiment_name'] = comp_name + ", dt=" + str(dt)
    comp = PPPYComp(**conf)
    fig1, plots = comp.plot_multi((2, 3), [('evol', dict(var_names=[s])
                                           ) for s in ['rv','rc', 'rr', 'T', 'nc', 'nr']],
                                  figsize=(24, 14))
    fig2, plots = comp.plot_multi((2, 4), [('evol', dict(var_names=[s])
                                           ) for s in ['ri','rs', 'rg', 'rh',
                                                       'ni', 'ns', 'ng', 'nh']],
                                  figsize=(24, 14))
    fig1.savefig(conf['experiment_tag'] + "/fig1_dt=" + str(dt) + ".png")
    fig2.savefig(conf['experiment_tag'] + "/fig2_dt=" + str(dt) + ".png")
    plt.close(fig1)
    plt.close(fig2)

#Plots by scheme
for schemeTag in scheme_list_sch.keys():
    conf['schemes'] = scheme_list_sch[schemeTag]
    conf['experiment_name'] = comp_name + ", " + long_names[schemeTag]
    comp = PPPYComp(**conf)
    fig1, plots = comp.plot_multi((2, 3), [('evol', dict(var_names=[s])
                                           ) for s in ['rv','rc', 'rr', 'T', 'nc', 'nr']],
                                  figsize=(24, 14))
    fig2, plots = comp.plot_multi((2, 4), [('evol', dict(var_names=[s])
                                           ) for s in ['ri','rs', 'rg', 'rh',
                                                       'ni', 'ns', 'ng', 'nh']],
                                  figsize=(24, 14))
    fig1.savefig(conf['experiment_tag'] + "/fig1_sch=" + schemeTag + ".png")
    fig2.savefig(conf['experiment_tag'] + "/fig2_sch=" + schemeTag + ".png")
    plt.close(fig1)
    plt.close(fig2)

#Plots by parameter
for scheme in scheme_list_all:
    scheme.name = scheme.name.split()[-1]
conf['schemes'] = scheme_list_all
conf['experiment_name'] = comp_name
comp = PPPYComp(**conf)
for var in ['rv','rc', 'rr', 'ri','rs', 'rg', 'rh', 'T', 'nc', 'nr', 'ni', 'ns', 'ng', 'nh']:
    plots = [('evol', dict(var_names=[var],
                           title=long_names[schemeTag],
                           only_param_tag=[scheme.tag for scheme in scheme_list_sch[schemeTag]],
                           switch_cls=True,
                           linewidth=3,
                          )) for schemeTag in keys]
    fig, plots = comp.plot_multi((3, 3), plots,
                                 title=VAR_NAME[var] + " (" + VAR_UNIT[var] + ")",
                                 figsize=(24, 14), sharey=True)
    fig.tight_layout(pad=5)
    fig.savefig(conf['experiment_tag'] + "/var=" + var + ".png")

