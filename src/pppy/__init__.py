#!/usr/bin/env python3
"""
This module aims at calling different implementations
of an individual parameterization to compare them.

``pppy`` stands for Physical Parameterizations with PYthon.
"""

__all__ = []

__version__ = '1.1.0'

__license__ = 'CeCILL-C'

__authors__ = ['Sébastien Riette']

__contributors__ = []

#List of the known variable names
#A variable not listed in VAR_NAME will have an empty legend
VAR_NAME = {
            'qv': "Specific content of vapor",
            'qc': "Specific content of cloud droptlets",
            'qr': "Specific content of rain",
            'qi': "Specific content of cloud ice",
            'qs': "Specific content of snow",
            'qg': "Specific content of graupel",
            'qh': "Specific content of hail",
            'rv': "Mixing-ratio of vapor",
            'rc': "Mixing-ratio of cloud droptlets",
            'rr': "Mixing-ratio of rain",
            'ri': "Mixing-ratio of cloud ice",
            'rs': "Mixing-ratio of snow",
            'rg': "Mixing-ratio of graupel",
            'rh': "Mixing-ratio of hail",
            'nc': "Number concentration of cloud droplets",
            'nr': "Number concentration of rain",
            'ni': "Number concentration of cloud ice",
            'ns': "Number concentration of snow",
            'ng': "Number concentration of graupel",
            'nh': "Number concentration of hail",
            'T': "Temperature",
            'Theta': "Potential temperature",
            'P': "Pressure",
            't': "Time",
            'Z_half': "Altitude of flux levels",
            'Z_mass': "Altitude of mass levels",
            'ccn1ft': "Water-friendly free aerosol number (CCN)",
            'ccn1at': "Water-friendly activated aerosol number (CCN)",
            'ifn1ft': "Ice-friendly free aerosol number (IFN)",
            'ifn1at': "Ice-friendly activated aerosol number (IFN)",
            'u': "U-component of wind",
            'v': "V-component of wind",
            'w': "Vertical velocity (positive for updraft)",
            "CF": "Cloud fraction",
            'Ps': 'Surface pressure',
            'Zs': 'Surface altitude above sea level',
            'slhf': 'Surface latent heat flux',
            'sshf': 'Surface sensible heat flux',
           }

#Unit for variables
#A variable not listed in VAR_UNIT will have an empty unit in legends
VAR_UNIT = {
            'qv': "kg/kg",
            'qc': "kg/kg",
            'qr': "kg/kg",
            'qi': "kg/kg",
            'qs': "kg/kg",
            'qg': "kg/kg",
            'qh': "kg/kg",
            'rv': "kg/kg",
            'rc': "kg/kg",
            'rr': "kg/kg",
            'ri': "kg/kg",
            'rs': "kg/kg",
            'rg': "kg/kg",
            'rh': "kg/kg",
            'nc': "#/kg",
            'nr': "#/kg",
            'ni': "#/kg",
            'ns': "#/kg",
            'ng': "#/kg",
            'nh': "#/kg",
            'T': "K",
            'Theta': 'K',
            'P': "Pa",
            't': 's',
            'Z_half': 'm',
            'Z_mass': 'm',
            'ccn1ft': "#/kg",
            'ccn1at': "#/kg",
            'ifn1ft': "#/kg",
            'ifn1at': "#/kg",
            'u': 'm/s',
            'v': 'm/s',
            'w': 'm/s',
            'CF': '0-1',
            'Ps': 'Pa',
            'Zs': 'm',
            'slhf': 'W/m²',
            'sshf': 'W/m²'
           }

COLORS = ['black', 'red', 'darksalmon', 'gold', 'olivedrab', 'silver',
          'chartreuse', 'skyblue', 'darkblue', 'purple', 'magenta']

STYLES = ['-', ':', '--', '-.', (20.0, 20.0), (30.0, 10.0),
          (20., 5., 1., 5., 1., 5., 1., 5.), (20., 5., 5., 5., 5., 5.)]

from .pppy import PPPY
from .pppyComp import PPPYComp
