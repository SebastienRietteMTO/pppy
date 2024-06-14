#!/usr/bin/env python3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

#Written by S. Riette

"""
This an example of comparison involving microphysics parameterization
"""

import os
import json
import numpy
import matplotlib.pyplot as plt

from pppy import PPPYComp, VAR_NAME, VAR_UNIT
from pppy_microphy import pppy_microphy

#Configuration of comparison
output_dir = os.path.dirname(os.path.abspath(__file__))

#Configuration of schemes
dtList = [0.001, 0.01, 0.1, 1., 5., 10., 15., 20.]
dtList = [1., 5., 10., 15.]
hail = True

icev = '4' if hail else '3'
old = {dt:pppy_microphy(dt=dt, method='step-by-step',
                        name="old ICE" + icev + " scheme with dt=" + str(dt) +"s",
                        tag="OLD" + icev + "_dt=" + str(dt),
                        namel=json.dumps({'NAM_NEBn': {'CFRAC_ICE_ADJUST': 'S',
                                                       'VSIGQSAT': 0.02, 'CCONDENS':'CB02',
                                                       'LSUBG_COND': True, 'LSIGMAS': True},
                                          'NAM_PARAM_ICEn': {'CSEDIM': 'NONE', 'LRED': False,
                                                             'LSEDIC': False,
                                                             'CSUBG_AUCV_RC': 'PDF',
                                                             'CSUBG_PR_PDF': 'SIGM',
                                                             'NMAXITER_MICRO': 100,
                                                             'XMRSTEP': 5.E-5, 'XTSTEP_TS': 0.}}))
          for dt in dtList}

ice = {dt:pppy_microphy(dt=dt, method='step-by-step',
                        name="new ICE" + icev + " scheme with dt=" + str(dt) +"s",
                        tag="ICE" + icev + "_dt=" + str(dt),
                        namel=json.dumps({'NAM_NEBn': {'CFRAC_ICE_ADJUST': 'S',
                                                       'VSIGQSAT': 0.02, 'CCONDENS':'CB02',
                                                       'LSUBG_COND': True, 'LSIGMAS': True},
                                          'NAM_PARAM_ICEn': {'CSEDIM': 'NONE', 'LRED': True,
                                                             'CSUBG_AUCV_RC': 'PDF',
                                                             'CSUBG_PR_PDF': 'SIGM',
                                                             'NMAXITER_MICRO': 100,
                                                             'XMRSTEP': 5.E-5, 'XTSTEP_TS': 0.}}))
          for dt in dtList}

scheme_names = {'old_ICE':old, 'new_ICE':ice}
long_names = {'old_ICE':'old ICE', 'new_ICE':'new ICE'}

keys = sorted(scheme_names.keys(), key=str.lower)
scheme_list_all = [scheme[dt] for dt in dtList for scheme in scheme_names.values()]
scheme_list_dt = {dt:[scheme_names[key][dt] for key in keys] for dt in dtList}
scheme_list_sch = {k: [v[dt] for dt in dtList] for (k, v) in scheme_names.items()}

#Comparison and plots
common = dict(dzz=numpy.array([[1.]]),
              Z_mass=numpy.array([[0.]]),
              sea=numpy.array([0.]),
              town=numpy.array([0.]),
              cum_c=numpy.array([[0.]]),
              cum_r=numpy.array([[0.]]),
              cum_s=numpy.array([[0.]]),
              cum_g=numpy.array([[0.]]),
              src=numpy.array([[0.]]),
              CF=numpy.array([[0.]]),
              HLC_HRC=numpy.array([[0.]]),
              HLC_HCF=numpy.array([[0.]]),
              HLI_HRI=numpy.array([[0.]]),
              HLI_HCF=numpy.array([[0.]]),
              sigs=numpy.array([[0.]]),
              CF_MF=numpy.array([[0.]]),
              rc_MF=numpy.array([[0.]]),
              ri_MF=numpy.array([[0.]]),
             )

warm_state = dict(P=numpy.array([[100000.]]), Theta=numpy.array([[290.]]),
                  rv=numpy.array([[1.E-2]]), rc=numpy.array([[1.E-2]]),
                  rr=numpy.array([[1.E-4]]), ri=numpy.array([[0.]]),
                  #ric =numpy.array([[0.]]), rid=numpy.array([[0.]]), rip=numpy.array([[0.]]),
                  rs=numpy.array([[0.]]), rg=numpy.array([[0.]]),
                  #nc=numpy.array([[3.E8]]), nr=numpy.array([[2000.]]),
                  ni=numpy.array([[0.]]), #ns=numpy.array([[0.]]),
                  #ng=numpy.array([[0.]]), nh=numpy.array([[0.]]),
                  #nic=numpy.array([[0.]]), nid=numpy.array([[0.]]), nip=numpy.array([[0.]]),
                  #ccn1ft=numpy.array([[1.E8]]), ccn1at=numpy.array([[1.E8]]), #CCN
                  #ifn1ft=numpy.array([[0.]]), ifn1at=numpy.array([[0.]]), #IFN
                  #u=numpy.array([[0.]]), v=numpy.array([[0.]]), w=numpy.array([[0.]]),
                  **common)
cold_state = dict(P=numpy.array([[100000.]]), Theta=numpy.array([[270.]]),
                  rv=numpy.array([[0.007]]), rc=numpy.array([[1.E-3]]),
                  rr=numpy.array([[1.E-3]]), ri=numpy.array([[9.E-4]]),
                  #ric =numpy.array([[3.E-4]]), rid=numpy.array([[3.E-4]]),
                  #rip=numpy.array([[3.E-4]]),
                  rs=numpy.array([[1.E-3]]), rg=numpy.array([[1.E-3]]),
                  #nc=numpy.array([[3.E8]]), nr=numpy.array([[2000.]]),
                  ni=numpy.array([[30000.]]), #ns=numpy.array([[21000.]]),
                  #ng=numpy.array([[8500.]]), nh=numpy.array([[0.]]),
                  #nic=numpy.array([[10000.]]), nid=numpy.array([[10000.]]),
                  #nip=numpy.array([[10000.]]),
                  #ccn1ft=numpy.array([[1.E8]]), ccn1at=numpy.array([[3.E8]]), #CCN
                  #ifn1ft=numpy.array([[1000.]]), ifn1at=numpy.array([[1000.]]), #IFN
                  #u=numpy.array([[0.]]), v=numpy.array([[0.]]), w=numpy.array([[0.]]),
                  **common)
cold2_state = dict(P=numpy.array([[100000.]]), Theta=numpy.array([[270.]]),
                   rv=numpy.array([[0.004]]), rc=numpy.array([[1.E-3]]),
                   rr=numpy.array([[1.E-3]]), ri=numpy.array([[9.E-4]]),
                   #ric =numpy.array([[3.E-4]]), rid=numpy.array([[3.E-4]]),
                   #rip=numpy.array([[3.E-4]]),
                   rs=numpy.array([[1.E-3]]), rg=numpy.array([[1.E-3]]),
                   nc=numpy.array([[3.E8]]), nr=numpy.array([[2000.]]),
                   ni=numpy.array([[30000.]]), #ns=numpy.array([[21000.]]),
                   #ng=numpy.array([[8500.]]), nh=numpy.array([[0.]]),
                   #nic=numpy.array([[10000.]]), nid=numpy.array([[10000.]]),
                   #nip=numpy.array([[10000.]]),
                   #ccn1ft=numpy.array([[1.E8]]), ccn1at=numpy.array([[3.E8]]), #CCN
                   #ifn1ft=numpy.array([[1000.]]), ifn1at=numpy.array([[1000.]]), #IFN
                   #u=numpy.array([[0.]]), v=numpy.array([[0.]]), w=numpy.array([[0.]]),
                   **common)
state_name = 'cold'
state = {'cold':cold_state, 'cold2':cold2_state, 'warm':warm_state}[state_name]
if hail:
    rh = ['rh']
    nh = ['nh']
    state['rh'] = numpy.array([[0.]])
else:
    rh = []
    nh = []
rwarm = ['rc', 'rr']
nwarm = ['nc', 'nr']
rcold = ['ri', 'rs', 'rg'] + rh
ncold = ['ni', 'ns', 'ng'] + nh
ncold = nwarm = [] #Preparation for LIMA
comp_name = f"Test with PHYEX, {state_name} state"
conf = {
        'schemes' : scheme_list_all,
        'output_dir': output_dir,
        'duration': 100.,
        'init_state': state,
        'name': comp_name,
        'tag': f"test_microphy_{state_name}"
       }

#Run all simulations
comp = PPPYComp(**conf)
comp.run(force=False)

#Plots by timestep
for dt in dtList:
    conf['schemes'] = scheme_list_dt[dt]
    conf['name'] = comp_name + ", dt=" + str(dt)
    comp = PPPYComp(**conf)
    fig1, plots = comp.plot_multi((2, 3), [('evol', dict(var_names=[s])
                                           ) for s in ['rv'] + rwarm + ['Theta'] + nwarm],
                                  figsize=(24, 14))
    fig2, plots = comp.plot_multi((2, 4), [('evol', dict(var_names=[s])
                                           ) for s in rcold + ncold],
                                  figsize=(24, 14))
    fig1.savefig(conf['tag'] + "/fig1_dt=" + str(dt) + ".png")
    fig2.savefig(conf['tag'] + "/fig2_dt=" + str(dt) + ".png")
    plt.show()
    plt.close(fig1)
    plt.close(fig2)

#Plots by scheme
for schemeTag in scheme_list_sch.keys():
    conf['schemes'] = scheme_list_sch[schemeTag]
    conf['name'] = comp_name + ", " + long_names[schemeTag]
    comp = PPPYComp(**conf)
    fig1, plots = comp.plot_multi((2, 3), [('evol', dict(var_names=[s])
                                           ) for s in ['rv'] + rwarm + ['Theta'] + nwarm],
                                  figsize=(24, 14))
    fig2, plots = comp.plot_multi((2, 4), [('evol', dict(var_names=[s])
                                           ) for s in rcold + ncold],
                                  figsize=(24, 14))
    fig1.savefig(conf['tag'] + "/fig1_sch=" + schemeTag + ".png")
    fig2.savefig(conf['tag'] + "/fig2_sch=" + schemeTag + ".png")
    plt.show()
    plt.close(fig1)
    plt.close(fig2)

#Plots by parameter
for scheme in scheme_list_all:
    scheme.name = scheme.name.split()[-1]
conf['schemes'] = scheme_list_all
conf['name'] = comp_name
comp = PPPYComp(**conf)
for var in ['rv'] + rwarm + rcold + ['Theta'] + nwarm + ncold:
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
    fig.savefig(conf['tag'] + "/var=" + var + ".png")
    plt.show()
    plt.close(fig)
