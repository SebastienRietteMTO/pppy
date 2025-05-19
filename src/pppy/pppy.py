#!/usr/bin/env python3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

"""
The PPPY class is an abstract class to represent an individual parameterization.
The methods that must (or can) be defined are:
    - __init__ (optional) to list the options allowed for the parameterization
                          with default values. super().__init__(...) is then
                          called with all arguments.
    - setup (optional) this is the place where the parameterization can be
                       initialized. The setup method is called in the
                       subprocess reserved for the parameterization whereas
                       __init__ is called in the main process. It can be
                       important to delay some initialisations here when dealing
                       with shared library to prevent interferences. Moreover,
                       setup is called each time the parameterization is used
                       whereas __init__ is called only once.
    - finalize (optional) to remove temporary files for instance
    - build_init_state (optional) this is the place to add variables to the
                                  provided common state. For instance, if a
                                  parameterization needs to track a specific
                                  variable that is not in the variables tracked
                                  by the tool.
    - diag (optional) to compute and save diagnostics. This methods is called
                      at the end.
    - execute The method receives the state (dictionary holding the different
              variables) and must return the same dictionary containing the
              different variables after a one-time-step advance. If conversions
              (units, array shape, 32bit/64bit...) are needed to call the
              parameterization, this is the place to perform them. Be careful
              to do the reverse conversion after the actual call to return the
              same variables in the output state as the one received.
If these methods are defined, they must call the method of the abstract class
through the super() mechanism.

The instances of PPPY are intended to be used with an instance of the PPPYComp class
"""

import numpy
import os
import h5py
import copy

from pppy import VAR_NAME, VAR_UNIT

class PPPY():
    """
    Abstract class used to call a single version of a parameterization.

    To create a usable class, one needs to implement the following methods:
     - setup (optional): called first, may be used to initialize some values
     - build_init_state (optional): called after setup to complete initial state
     - execute (mandatory): called to compute evolution of state variables
       over one time step
     - finalize (optional): called at the end to enable operations
       like unloading library from memory

    To be as general as possible, individual variables followed by the
    parameterization are not defined. Only one variable is followed (named
    state in the different methods) and this variable is a dictionary
    whose keys are implementation dependent and can vary from one
    parameterization to the other.
    """

    def __init__(self, dt, method, name, tag, **options):
        """
        At instantiation, we define the "time integration method" and
        the timestep.
        Moreover, we set name and tag to identify the parameterization.

        :param dt:        time step to use (in seconds)
        :param method:    among ('step-by-step', 'one-step'). Two kind of
                          simulations are implemented. 'step-by-step' is
                          like a true simulation, output of a time step is
                          computed from output of previous time step.
                          For 'one-step', output (at all output times)
                          is computed by direct integration from the initial state
        :param name:      name of the parameterization, used in legends
        :param tag:       string to identify uniquely the parameterization, the
                          string must be usable to build filenames
        :param \*\*options: a specific implementation may add some args and
                          kwargs. Those transmitted to this __init__
                          method (through this \*\*options argument) will be
                          stored in the output file as metadata.

        This method can be extended by a specific implementation to add,
        in particular, specific options.

        Instances are intended to be used by a PPPYComp instance.
        """

        assert isinstance(dt, float), "dt must be a float"
        assert dt > 0., "dt must be > 0."
        assert method in ('step-by-step', 'one-step'), "method must be 'step-by-step' or 'one-step'"

        self._dt = dt
        self._method = method
        self.name= name
        self.tag = tag
        self._options = options

    def setup(self, init_state, duration):
        """
        This method is called each time the __call__ method is called.

        :param init_state: dictionary holding the 'state' variables
        :param duration:   the simulation total duration in seconds
        """

        pass

    def finalize(self):
        """
        This method is called at the end of each call to the __call__ method.
        This method can be implemented if something needs to be done after execution.
        """

        pass

    def execute(self, previous_state, timestep, timestep_number):
        """
        This method do the computational part of the time advance
        (setup and building of the initial state are done before).

        :param previous_state:  dictionary holding all the followed 'state' variables.
        :param timestep:        is the simulation duration in seconds (this is the time step
                                defined in __init__ if method is 'step-by-step')
        :param timestep_number: is the number of current time step (always 1 if method is 'one-step')
        :returns: a 'state' dictionary containing variable values after time integration

        The method must be implemented.
        """

        assert isinstance(previous_state, dict), "state must be a dictionary"
        assert isinstance(timestep, float), "duration must be a float"
        assert timestep > 0., "duration must be > 0."
        assert isinstance(timestep_number, int), "timestep_number must be an integer"
        assert timestep_number >= 1, "timestep_number must be at least 1"

    def _open_file(self, output_file, times, init_state):
        """
        This internal method opens an hdf5 file and creates all dataset.

        :param output_file: filename of the hdf5 file to open
        :param times: output times
        :param init_state: dictionary containing the initial 'state'

        :returns: The method returns the hdf5 file and a dictionary containing the datasets.
        """
        output = h5py.File(output_file, 'w')
        output.attrs['dt'] = self._dt
        output.attrs['method'] = self._method
        for key, value in self._options.items():
            output.attrs[key] = value
        dset = {}
        dset['t'] = output.create_dataset("t", shape=(len(times),), dtype=float, data=times)
        dset['t'].attrs['name'] = VAR_NAME.get('t', '')
        dset['t'].attrs['unit'] = VAR_UNIT.get('t', '')
        for key, value in init_state.items():
            shape = tuple([len(times)] + list(value.shape))
            dset[key] = output.create_dataset(key, shape=shape, dtype=value.dtype)
            dset[key].set_fill_value = numpy.nan
            dset[key][...] = numpy.nan
            dset[key].attrs['name'] = VAR_NAME.get(key, '')
            if VAR_UNIT.get(key, '') is not None:
                dset[key].attrs['unit'] = VAR_UNIT.get(key, '')
            dset[key][0] = value
        return output, dset

    def __call__(self, init_state, duration, output_file):
        """
        This method do the time advance.

        :param init_state:  dictionary holding the 'state' variables
        :param duration:    the simulation total duration in seconds
        :param output_file: the filename of the hdf5 file that will
                            contain the output variables
        """

        assert isinstance(init_state, dict), "state must be a dictionary"
        assert isinstance(duration, float), "duration must be a float"
        assert duration > 0., "duration must be > 0."

        assert all([isinstance(value, numpy.ndarray) for value in list(init_state.values())]), \
               "All init_state item must be ndarrays"

        if os.path.exists(output_file):
            raise IOError("output file already exists")

        #We prepare state and output file
        self.setup(init_state, duration)
        state = self.build_init_state(init_state)
        old_state = copy.deepcopy(state)
        times = numpy.array(list(numpy.arange(0, duration, self._dt)) + [duration])
        output, dset = self._open_file(output_file, times, state)

        try:
            for i in range(1, len(times)):
                if self._method == 'step-by-step':
                    state = self.execute(old_state, times[i] - times[i-1], i)
                else:
                    state = self.execute(copy.deepcopy(old_state), times[i], 1)
                for key, value in old_state.items():
                    if key not in state:
                        state[key] = value
                if self._method == 'step-by-step':
                    old_state = state
                for key, value in state.items():
                    dset[key][i] = value
            self.diag(output, dset)
            output.close()
            self.finalize()
        except:
            try:
                output.close()
            except KeyboardInterrupt:
                raise
            except:
                pass
            if os.path.exists(output_file):
                os.remove(output_file)
            self.finalize()
            import traceback
            traceback.print_exc()
            raise

    def build_init_state(self, state):
        """
        Method used to modify the initial state.
        Initial state can be incomplete. For instance, when comparing microphysical
        schemes, we can only know the initial mixing ratios and this method
        must compute initial concentration number if the scheme needs it.
        This method must add to the dictionary all the missing variables.

        :param state: dictionary holding the initial 'state' that can be incomplete
                      for the actual scheme
        :returns: complete 'state' dictionary
        """

        assert isinstance(state, dict), "state must be a dictionay"

        return state

    def diag(self, hdf5file, dsets):
        """
        Method used to acompute diagnostics to be stored in the output file
        :param hdf5file: the output file (in hdf5 format)
        :param dsets: a dictionary containing the existing datasets
        """
