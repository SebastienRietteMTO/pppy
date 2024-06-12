#!/usr/bin/env python3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

"""
This module contains the PPPY implementation for a test parameterization
"""

import numpy
import ctypesForFortran
import pppy

class pppy_param2(pppy.PPPY):
    """PPPY implementation for a test parameterization"""

    def __init__(self, dt, method, name, tag, solib, iconf):
        """
        In addition to dt, method, name and tag parameters
        defined in the PPPY class, this parameterization
        needs the following parameters:
        solib : path to the shared library to use
        iconf : useless parameter
        """
        super().__init__(dt, method, name, tag,
                         solib=solib, iconf=iconf)
        self._solib = solib

        #Check received options
        assert iconf in [0, 1], "iconf must be 0 or 1"
        self._iconf = iconf

        self._init = None
        self._handle = None
        self._param = None

    def setup(self, init_state, duration):
        """
        This method opens the shared library and do the init stuff.
        """
        super().setup(init_state, duration)

        print("A new simulation begins with a duration of " + str(duration) + "s " +
              "for the parameterization #2")

        #Opening of shared lib is done here to avoid loading it at initialization
        #Because if shared lib is loaded at init, all declared schemes will have
        #their shared lib loaded simultaneously which can be a source of error.
        #This way, shared lib is loaded in a sub-process (if class instance is used
        #through a PPPYComp instance).

        #Define fortran routine to use through ctypesForFortran
        IN = ctypesForFortran.IN
        OUT = ctypesForFortran.OUT
        INOUT = ctypesForFortran.INOUT

        ctypesFF, self._handle = ctypesForFortran.ctypesForFortranFactory(self._solib)

        @ctypesFF()
        def init(ICONF): #Name of the function must be the name of the actual fortran function
            "init function"
            return ([ICONF],
                    [(numpy.int64, None, IN), #INTEGER, INTENT (IN)
                    ],
                    None)
        self._init = init

        @ctypesFF()
        def param2_py(x):
            "Function that actually call the parameterization"
            return ([x, x.shape[0], x.shape[1]],
                    [(numpy.float64, x.shape, IN),
                     (numpy.float64, x.shape, OUT),
                     (numpy.int64, None, IN),
                     (numpy.int64, None, IN)
                     ], None)
        self._param = param2_py

        #Setup
        self._init(self._iconf)

    def finalize(self):
        """
        We close the shared library so we are able to load a new one with the same symbol names.
        """
        super().finalize()
        ctypesForFortran.dlclose(self._handle)

    def build_init_state(self, state):
        """
        This method can modify the initial state to include variables
        needed by the parameterization
        """
        state = super().build_init_state(state)
        #Nothing to do here, in this example
        return state

    def execute(self, previous_state, timestep, timestep_number):
        """
        This method do the actual work
        """
        super().execute(previous_state, timestep, timestep_number)

        #previous value for variable x
        previous = previous_state['x']

        #Dimensions
        #During the comparison, x is a 1D vector but param need a 2D matrix
        x = previous.reshape((len(previous), 1))

        #We call the parameterization
        y = self._param(x)

        #Dimensions
        result = y.reshape((len(previous), ))

        return {'x': result}



