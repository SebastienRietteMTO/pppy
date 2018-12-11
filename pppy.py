#!/usr/bin/env python3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

"""
The PPPY class represent an individual parameterization.
The PPPYComp class is used to run the different parameterizations,
save the results and plot some variables
"""

import numpy
import os
import h5py
import copy
import multiprocessing
import logging

#List of the known variable names
#A variable not listed in VAR_NAME will have an empty legend
VAR_NAME = {'rv': "Mixing-ratio of vapor",
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
            'w': "Vertical velocity (positive for updraft)"
           }

#Unit for variables
#A variable not listed in VAR_UNIT will have an empty unit in legends
VAR_UNIT = {'rv': "kg/kg",
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
            'w': 'm/s'
           }

COLORS = ['black', 'red', 'darksalmon', 'gold', 'olivedrab', 'silver',
          'chartreuse', 'skyblue', 'darkblue', 'purple', 'magenta']
STYLES = ['-', ':', '--', '-.', (20.0, 20.0), (30.0, 10.0),
          (20., 5., 1., 5., 1., 5., 1., 5.), (20., 5., 5., 5., 5., 5.)]

class PPPY():
    """
    Abstract class used to call a single version of a parameterization.

    To create a usable class, one needs to implement the following methods:
     - setup (optional): called first, may be used to initialize some values
     - build_init_state (optional): called after setup to complete initial state
     - execute (mandatory): called to compute evolution of state variables
       over one timestep
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

        :param dt:        timestep to use (in seconds)
        :param method:    among ('step-by-step', 'one-step'). Two kind of
                          simulations are implemented. 'step-by-step' is
                          like a true simulation, output of a timestep is
                          computed from output of previous timestep.
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
        :param timestep:        is the simulation duration in seconds (this is the timestep
                                defined in __init__ if method is 'step-by-step')
        :param timestep_number: is the number of current timestep (always 1 if method is 'one-step')
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

class PPPYComp():
    """
    This class is used to perform a comparison between several implementation
    of a parameterization.
    This class can be directly instantiated, there is no need to extend it,
    unless you need to write some diagnostic methods (e.g. plots, statistics...).
    """

    def __init__(self, schemes, output_dir, duration, init_state, experiment_name, experiment_tag):
        """
        At instantiation, we define the schemes to compare, the total duration of the
        time integration, the initial conditions and the name of the experiment.

        :param schemes:         a list of PPPY instances
        :param output_dir:      directory in which results are stored
        :param duration:        total duration of the time integration (in seconds)
        :param init_state:      dictionary holding the initial state
        :param experiment_name: name of the comparison experiment used in legends
        :param experiment_tag:  string to identify uniquely the experiment, the
                                string must be usable to build filenames

        After instantiation, the different parameterizations can be run with the
        run method. Results can be plotted (plot_evol, plot_comp and plot_multi
        methods) or time series can be accessed (get_series method).
        """

        assert isinstance(schemes, list), "schemes must be a list"
        for scheme in schemes:
            assert isinstance(scheme, PPPY), "third element of tuple must be a PPPY instance"
        if not len(schemes) == len(set([scheme.name for scheme in schemes])):
            logging.warning("Scheme names are not unique")
        assert len(schemes) == len(set([scheme.tag for scheme in schemes])), \
               "all scheme tags must be different"

        assert isinstance(duration, float) and duration > 0., \
               "duration must be a float and must be > 0."
        assert isinstance(init_state, dict), "init_state must be a dict"
        assert isinstance(experiment_name, str) and isinstance(experiment_tag, str), \
               "experiment_name and experiment_tag must be strings"

        values = list(init_state.values())
        assert all([isinstance(value, numpy.ndarray) for value in values]), \
               "All init_state item must be ndarrays"

        if not os.path.exists(output_dir):
            raise IOError("output dir does not exist")
        self._output_dir = os.path.join(output_dir, experiment_tag)
        if not os.path.exists(self._output_dir):
            os.mkdir(self._output_dir)

        self._schemes = schemes
        self._duration = duration
        self._init_state = init_state
        self._experiment_name = experiment_name
        self._experiment_tag = experiment_tag

    def run(self, only_param_tag=None, force=False):
        """
        This methods runs the different parameterizations.

        :param only_param_tag: (optional) list of parameterization tags to actually run.
                               By default (None), all parameterizations are run.
        :param force:          (optional) By default (False), parameterizations are run only
                               if needed (output file for same experiment_tag and
                               same parameterization tag not found). If force is True,
                               parameterizations are run anyway.
        """
        assert isinstance(only_param_tag, list) or only_param_tag is None, \
               "only_param_tag must be a list or None"
        if not only_param_tag is None:
            assert all([isinstance(item, str) for item in only_param_tag]), \
                   "if only_param_tag is a list, each item must be a string"

        for scheme in [scheme for scheme in self._schemes
                       if only_param_tag is None or scheme.tag in only_param_tag]:
            output_file = os.path.join(self._output_dir, scheme.tag + '.hdf5')
            if os.path.exists(output_file) and force:
                os.remove(output_file)
            if not os.path.exists(output_file):
                #To isolate calls, we use multiprocessing instead of calling directly
                #scheme(copy.deepcopy(self._initState), self._duration, output_file)
                #Without this trick we can have issues with different schemes using
                #same symbol names and eventually with other python libs (such as plt)
                child = multiprocessing.Process(target=scheme,
                                                args=(copy.deepcopy(self._init_state),
                                                      self._duration, output_file))
                child.start()
                child.join()

    @staticmethod
    def _var_legend(var, mfactor):
        """
        Internal method that return the string to use as a legend.

        :param var:     name of the variable
        :param mfactor: multiplying factor applied to data

        :returns: The method returns the string to use.
        """
        legend = VAR_NAME.get(var, var)
        unit = VAR_UNIT.get(var, "??")
        if mfactor == 1.:
            legend += " (" + unit + ")"
        else:
            legend += " (" + str(1./mfactor) + " " + unit + ")"
        return legend

    def plot_multi(self, plot_shape, plot_spec, title=None, figsize=None, sharey=False):
        """
        This method is a wrapper around matplotlib features to deal with a
        figure containing several plots.

        :param plot_shape: tuple (nb_rows, nb_cols) to ask for a figure containing
                           nb_rows * nb_cols subplots
        :param plot_spec:  list of plot definitions. Each plot definition is a tuple
                           (kind, args, kwargs). Method will plot the subfigure using
                           method plot_<kind> which must accept a matplotlib axis
                           as first parameter then args and kwargs (args and/or
                           kwargs can be omitted)
        :param title:      (optional) main title to use for the figure
        :param figsize:    (optional) figsize argument to build the matplotlib figure
        :param sharey:     (optional) must we share the y axis (this option is active
                           only if all plots have the same kind)

        :returns: The method returns the figure and the list of return values of each plot.
        """
        assert isinstance(plot_spec, list), "plot_spec must be a list"
        assert isinstance(plot_shape, tuple), "plot_shape must be a tuple"
        assert len(plot_shape) == 2, "plot_shape must have two elements"
        assert len(plot_spec) <= plot_shape[0] * plot_shape[1], \
               "number of plots must be inferior to plot_shape[0] * plot_shape[1]"

        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=figsize)
        if title is None:
            fig.suptitle(self._experiment_name)
        else:
            fig.suptitle(title)
        result = []
        for iplot, plot_def in enumerate(plot_spec):
            assert isinstance(plot_def, tuple), \
                   "Each plot definition must be a tuple"
            kind = plot_def[0]
            if len(plot_def) > 3:
                raise ValueError("A plot definition is a tuple with a length of 1, 2 or 3")
            elif len(plot_def) == 3:
                plot_args, plot_kwargs = plot_def[1:]
                assert isinstance(plot_args, list) and isinstance(plot_kwargs, dict), \
                       "Second item of a plot definition must be a list and third " + \
                       "item must be a dictionary"
            elif len(plot_def) == 2:
                assert isinstance(plot_def[1], list) or isinstance(plot_def[1], dict), \
                       "Second item of a plot definition must be a list or a dictionary"
                if isinstance(plot_def[1], list):
                    plot_args = plot_def[1]
                    plot_kwargs = {}
                else:
                    plot_kwargs = plot_def[1]
                    plot_args = []
            else:
                plot_args = []
                plot_kwargs = {}
            kwargs = {}
            if iplot != 0 and len(set([item[0] for item in plot_spec])) == 1:
                kwargs = {'sharex':ax}
                if sharey:
                    kwargs['sharey'] = ax
            ax = fig.add_subplot(plot_shape[0], plot_shape[1], iplot + 1, **kwargs)
            method = getattr(self, "plot_" + kind)
            result.append(method(ax, *plot_args, **plot_kwargs))
        return fig, result

    def common_times(self, only_param_tag=None, only_times=None, common_times=1):
        """
        This method extracts output times of the different parameterizations and
        looks for common times.

        :param only_param_tag: (optional) list of parameterization tags to actually consider.
                               By default (None), all parameterizations are used.
        :param only_times:     (optional) If set, result is limited to times listed in this list
        :param common_times:   (optional, defaults to 1) If 1, the time series is the list of all
                               encountered time values in all parameterizations; if 2, time values
                               are limited to values common to all parameterizations
        :returns: The method returns a numpy array containing the time values.
        """

        assert isinstance(only_param_tag, list) or only_param_tag is None, \
               "only_param_tag must be a list or None"
        if not only_param_tag is None:
            assert all([isinstance(item, str) for item in only_param_tag]), \
                   "if only_param_tag is a list, each item must be a string"
        assert isinstance(only_times, list) or only_times is None, \
            "only_times must be a list or None"
        if not only_times is None:
            assert all([isinstance(item, float) for item in only_times]), \
                   "if only_times is a list, each item must be a float"

        schemes = [scheme for scheme in self._schemes
                   if only_param_tag is None or scheme.tag in only_param_tag]
        files = {scheme.tag: h5py.File(os.path.join(self._output_dir, scheme.tag + '.hdf5'), 'r')
                 for scheme in schemes}
        if common_times == 1:
            times = []
            for scheme in self._schemes:
                times.extend(list(files[scheme.tag]['t'].value))
            times = sorted(list(set(times)))
            times = [t for t in times if only_times is None or t in only_times]
        elif common_times == 2:
            times = []
            for scheme in self._schemes:
                times = [t for t in list(files[scheme.tag]['t'].value) if t in times]
            times = sorted(times)
            times = [t for t in times if only_times is None or t in only_times]
        for hdf5_file in files.values():
            hdf5_file.close()

        return times

    def get_series(self, var_names, slicing=None,
                   mfactor=1., only_param_tag=None,
                   only_times=None, common_times=0):
        """
        This method read output files and return time series of variables.

        :param var_names:      list of variables to read
        :param slicing:        (optional) slicing to apply to the data series.
                               This can be used to run a parameterization on a
                               multidimensional initial state and extract result
                               for only one dimension.
        :param mfactor:        (optional) multiplying factor to apply to data
        :param only_param_tag: (optional) list of parameterization tags to actually consider.
                               By default (None), all parameterizations are used.
        :param only_times:     (optional) If set, result is limited to times listed
        :param common_times:   (optional, defaults to 0) If 0, each parameterization
                               has its own time values; If 1, the time series is the list
                               of all encountered time values in all parameterizations,
                               resulting values are masked where undefined; If 2, time values
                               are limited to values common to all parameterizations.
        :returns: The method returns a 3-tuple whose components are:

                  - var_names
                  - parameterization tags
                  - a dictionary with:
                      - keys of the form (<var_name>, <param_tag>)
                      - values being (time_serie, values)
        """
        assert isinstance(only_param_tag, list) or only_param_tag is None, \
               "only_param_tag must be a list or None"
        if not only_param_tag is None:
            assert all([isinstance(item, str) for item in only_param_tag]), \
                   "if only_param_tag is a list, each item must be a string"
        assert isinstance(only_times, list) or only_times is None, \
               "only_times must be a list or None"
        if not only_times is None:
            assert all([isinstance(item, float) for item in only_times]), \
                   "if only_times is a list, each item must be a float"
        assert isinstance(var_names, list), "var_names must be a list"
        assert all([isinstance(item, str) for item in var_names]), \
               "var_names must be a list of variable names"
        assert common_times in [0, 1, 2], "common_times must be 0, 1 or 2"

        #common_times
        if common_times in [1, 2]:
            times = self.common_times(only_param_tag, only_times, common_times)

        schemes = [scheme for scheme in self._schemes
                   if only_param_tag is None or scheme.tag in only_param_tag]
        files = {scheme.tag: h5py.File(os.path.join(self._output_dir, scheme.tag + '.hdf5'), 'r')
                 for scheme in schemes}

        result = {}
        for scheme in schemes:
            simul_time = files[scheme.tag]['t'].value
            if common_times == 0:
                if only_times is None:
                    time = simul_time
                else:
                    time = numpy.intersect1d(time, numpy.array(only_times))
            else:
                time = times
            used_time = numpy.intersect1d(time, simul_time, True)
            index_simul = numpy.searchsorted(simul_time, used_time)
            index_result = numpy.searchsorted(time, used_time)
            for var in var_names:
                simul_values = files[scheme.tag][var].value * mfactor
                if slicing is not None:
                    simul_values = simul_values[:, slicing]
                shape = tuple([len(time)] + list(simul_values.shape[1:]))
                serie = numpy.ma.masked_all(shape)
                serie[index_result] = simul_values[index_simul]
                result[(var, scheme.tag)] = (time, serie)

        for hdf5_file in files.values():
            hdf5_file.close()

        return var_names, [scheme.tag for scheme in schemes], result

    def plot_evol(self, ax, var_names, slicing=None, y_var=None,
                  mfactor=1., only_param_tag=None, enable_contourf=True,
                  contourlevels=None, title=None,
                  handlelength=4, framealpha=0.5, fontsize=8,
                  switch_cls=False, linewidth=None):
        """
        This method plots a time evolution (x-axis is time) of different variables
        for different schemes.

        If data is 0D, each couple variable/scheme is represented by a line.

        If data is 1D, each couple variable/scheme is represented by contour lines.

        :param ax:              matplolib axis to use for plotting
        :param var_names:       list of variables to plot
        :param slicing:         (optional) slicing to apply to the data series.
                                This can be used to run a parameterization on a
                                multidimensional initial state and extract result
                                for only one dimension.
        :param y_var:           (optional) If data is 1D, use this variable on the y-axis
        :param mfactor:         (optional) multiplying factor to apply to data
        :param only_param_tag:  (optional) list of parameterization tags to actually consider.
                                By default (None), all parameterizations are used.
        :param enable_contourf: (optional) If only one 1D variable is plotted, enables
                                use of contourf instead of contour
        :param contourlevels:   (optional) List of values to use with contour and contourf
        :param title:           (optional) Title to use on subplot
        :param handlelength:    (optional) The length of the legend handles
        :param framealpha:      (optional) Alpha value for legend box
        :param fontsize:        (optional) Font size for legends
        :param switch_cls:      (optional, default to False) By default, A line style is assigned
                                to each parameterization and a color is assign to each variable.
                                Setting this variable to True, reverse this default behavior.
        :param linewidth:       Control the line width used in plots

        :returns: The method returns the output of the plotting function (plot, contour or contourf)
        """
        assert isinstance(only_param_tag, list) or only_param_tag is None, \
               "only_param_tag must be a list or None"
        if not only_param_tag is None:
            assert all([isinstance(item, str) for item in only_param_tag]), \
                   "if only_param_tag is a list, each item must be a string"
        assert isinstance(var_names, list), "var_names must be a list"
        assert all([isinstance(item, str) for item in var_names]), \
               "var_names must be a list of variable names"

        plot_schemes = [scheme for scheme in self._schemes
                        if only_param_tag is None or scheme.tag in only_param_tag]

        #line style and color
        def filter_color(var, scheme):
            "returns the element (var or scheme) that controls the color"
            if switch_cls:
                return scheme
            else:
                return var
        def filter_ls(var, scheme):
            "returns the element (var or scheme) that controls the style"
            if switch_cls:
                return var
            else:
                return scheme
        if switch_cls:
            #color by parameterization
            if len(COLORS) >= len(self._schemes):
                #We keep the same association between
                #schemes and colors even if all schemes are not plotted
                colors = {k:COLORS[i] for i, k in enumerate([scheme.tag
                                                             for scheme in self._schemes])}
            else:
                colors = {k:COLORS[i] for i, k in enumerate([scheme.tag
                                                             for scheme in plot_schemes])}
            #color by variable
            styles = {k:STYLES[i] for i, k in enumerate(var_names)}
        else:
            #line style by parameterization
            if len(STYLES) >= len(self._schemes):
                #We keep the same association between
                #schemes and styles even if all schemes are not plotted
                styles = {k:STYLES[i] for i, k in enumerate([scheme.tag
                                                             for scheme in self._schemes])}
            else:
                styles = {k:STYLES[i] for i, k in enumerate([scheme.tag
                                                             for scheme in plot_schemes])}
            #color by variable
            colors = {k:COLORS[i] for i, k in enumerate(var_names)}

        #Plotting output
        _, _, series = self.get_series(var_names, slicing, mfactor, only_param_tag, common_times=0)
        if y_var is not None:
            _, _, series_y = self.get_series([y_var], slicing, 1., only_param_tag, common_times=0)
        if title is not None:
            ax.set_title(title)
        lines = []
        result = []
        for var in var_names:
            for scheme in plot_schemes:
                time, serie = series[(var, scheme.tag)]
                shape = serie.shape
                dim = sum([1 if l > 1 else 0 for l in shape])
                if dim not in [1, 2]:
                    raise NotImplementedError("We only plot 0D and 1D variables, please " +
                                              "use the slicing argument to reduce the problem size")
                if dim == 1:
                    if len(var_names) == 1:
                        if title is None:
                            ax.set_title(self._var_legend(var, mfactor))
                        label = scheme.name
                    elif len(plot_schemes) == 1:
                        if title is None:
                            ax.set_title(scheme.name)
                        label = self._var_legend(var, mfactor)
                    else:
                        label = self._var_legend(var, mfactor) + " - " + scheme.name
                    ls_is_str = isinstance(styles[filter_ls(var, scheme.tag)], str)
                    line, = ax.plot(time, serie.squeeze(),
                                    color=colors[filter_color(var, scheme.tag)],
                                    linestyle=styles[filter_ls(var, scheme.tag)] if ls_is_str
                                                                                 else '',
                                    label=label,
                                    **({'linewidth':linewidth} if linewidth is not None else {}))
                    if not ls_is_str:
                        line.set_dashes(styles[filter_ls(var, scheme.tag)])
                    lines.append(line)
                    result.append(line)
                else:
                    X = numpy.repeat(time, serie.squeeze().shape[-1]).reshape(serie.squeeze().shape)
                    if y_var is None:
                        Y = numpy.mgrid[0:len(time), 0:serie.squeeze().shape[-1]][1]
                    else:
                        _, Y = series_y[(y_var, scheme.tag)]
                        Y = Y.squeeze()
                    if len(plot_schemes) == 1 and len(var_names) == 1 and enable_contourf:
                        cs = ax.contourf(X, Y, serie.squeeze(), levels=contourlevels)
                        if title is None:
                            ax.set_title(self._var_legend(var, mfactor) + " - " + scheme.name)
                        cs.clabel(inline=1, fontsize=10, colors='black')
                    else:
                        ls_is_str = isinstance(styles[filter_ls(var, scheme.tag)], str)
                        cs = ax.contour(X, Y, serie.squeeze(),
                                        colors=colors[filter_color(var, scheme.tag)],
                                        linestyles=styles[filter_ls(var, scheme.tag)] if ls_is_str
                                                                                else None,
                                        label=self._var_legend(var, mfactor) + " - " + scheme.name,
                                        levels=contourlevels,
                                        **({'linewidths':linewidth} if linewidth is not None
                                                                    else {}))
                        cs.clabel(inline=1, fontsize=10)
                        if not ls_is_str:
                            raise NotImplementedError("plotting of contour lines with dashes " +
                                                      "is not yet implemented")
                            #code below does not work
                            for lc in cs.collections:
                                lc.set_dashes(styles[filter_ls(var, scheme.tag)])
                        ls_is_str = isinstance(styles[filter_ls(var, scheme.tag)], str)
                        line, = ax.plot([], [],
                                        color=colors[filter_color(var, scheme.tag)],
                                        linestyle=styles[filter_ls(var, scheme.tag)] if ls_is_str
                                                                                     else '',
                                        label=self._var_legend(var, mfactor) + " - " + scheme.name,
                                        **({'linewidth':linewidth} if linewidth is not None else{}))
                        if not ls_is_str:
                            line.set_dashes(styles[filter_ls(var, scheme.tag)])
                        lines.append(line)
                    result.append(cs)
        #Legend
        if len(var_names) == 1 or len(plot_schemes) == 1:
            if not (dim == 2 and len(plot_schemes) == 1 and
                    len(var_names) == 1 and enable_contourf):
                ax.legend(handlelength=handlelength, framealpha=framealpha, fontsize=fontsize)
        else:
            if switch_cls:
                lines = []
                for scheme in plot_schemes:
                    line, = ax.plot([], [], color=colors[filter_color('', scheme.tag)],
                                    linestyle='-', label=scheme.name,
                                     **({'linewidth':linewidth} if linewidth is not None else{}))
                    lines.append(line)
                ax.legend(handles=lines, title="Colors (schemes)", loc=1,
                          handlelength=handlelength, framealpha=framealpha, fontsize=fontsize)
                lines = []
                ax2 = ax.twinx()
                ax2.get_xaxis().set_visible(False)
                ax2.get_yaxis().set_visible(False)
                for var in var_names:
                    ls_is_str = isinstance(styles[filter_ls(var, '')], str)
                    line, = ax2.plot([], [], color='black', label=self._var_legend(var, mfactor),
                                     linestyle=styles[filter_ls(var, '')] if ls_is_str else '',
                                     **({'linewidth':linewidth} if linewidth is not None else{}))
                    if not ls_is_str:
                        line.set_dashes(styles[filter_ls(var, '')])
                    lines.append(line)
                ax2.legend(handles=lines, title="Line styles (parameters)", loc=2,
                           handlelength=handlelength, framealpha=framealpha, fontsize=fontsize)
            else:
                lines = []
                for var in var_names:
                    line, = ax.plot([], [], color=colors[filter_color(var, '')],
                                    linestyle='-', label=self._var_legend(var, mfactor),
                                    **({'linewidth':linewidth} if linewidth is not None else{}))
                    lines.append(line)
                ax.legend(handles=lines, title="Colors (parameters)", loc=1,
                          handlelength=handlelength, framealpha=framealpha, fontsize=fontsize)
                lines = []
                ax2 = ax.twinx()
                ax2.get_xaxis().set_visible(False)
                ax2.get_yaxis().set_visible(False)
                for scheme in plot_schemes:
                    ls_is_str = isinstance(styles[filter_ls('', scheme.tag)], str)
                    line, = ax2.plot([], [], color='black', label=scheme.name,
                                     linestyle=styles[filter_ls('', scheme.tag)] if ls_is_str
                                                                                 else '',
                                     **({'linewidth':linewidth} if linewidth is not None else{}))
                    if not ls_is_str:
                        line.set_dashes(styles[filter_ls('', scheme.tag)])
                    lines.append(line)
                ax2.legend(handles=lines, title="Line styles (schemes)", loc=2,
                           handlelength=handlelength, framealpha=framealpha, fontsize=fontsize)

        ax.set_xlabel("Time")
        if dim != 1:
            ax.set_ylabel("Index" if y_var is None else y_var)
        else:
            pass

        return result

    def plot_comp(self, ax, var_names, slicing=None, y_var=None,
                  mfactor=1., only_param_tag=None, enable_contourf=True,
                  contourlevels=None, title=None, x_var=None, only_times=None,
                  handlelength=4, framealpha=0.5, fontsize=8,
                  switch_cls=False, linewidth=None):
        """
        This method plots a scheme comparison (x-axis is for scheme) of different variables
        for different times.

        If data is 0D, each couple variable/scheme is represented by a line.

        If data is 1D, each couple variable/scheme is represented by contour lines.

        :param ax:              matplolib axis to use for plotting
        :param var_names:       list of variables to plot
        :param slicing:         (optional) slicing to apply to the data series.
                                This can be used to run a parameterization on a
                                multidimensional initial state and extract result
                                for only one dimension.
        :param y_var:           (optional) If data is 1D, use this variable on the y-axis
        :param mfactor:         (optional) multiplying factor to apply to data
        :param only_param_tag:  (optional) list of parameterization tags to actually consider.
                                By default (None), all parameterizations are used.
        :param enable_contourf: (optional) If only one 1D variable is plotted, enables
                                use of contourf instead of contour
        :param contourlevels:   (optional) List of values to use with contour and contourf
        :param title:           (optional) Title to use on subplot
        :param x_var:           (optional) By default, the schemes are evenly spaced on the x-axis.
                                If this variable is set, it must be a list of position to use
                                (specified in the same order as schemes).
        :param only_times:      (optional) By default, all output times are plotted.
                                If set, plotting is limited to times listed.
        :param handlelength:    (optional) The length of the legend handles
        :param framealpha:      (optional) Alpha value for legend box
        :param fontsize:        (optional) Font size for legends
        :param switch_cls:      (optional, default to False) By default, A line style is assigned
                                to each time and a color is assign to each variable.
                                Setting this variable to True, reverse this default behavior.
        :param linewidth:       Control the line width used in plots

        :returns: The method returns the output of the plotting function (plot, contour or contourf)
        """

        assert isinstance(only_param_tag, list) or only_param_tag is None, \
               "only_param_tag must be a list or None"
        if not only_param_tag is None:
            assert all([isinstance(item, str) for item in only_param_tag]), \
                   "if only_param_tag is a list, each item must be a string"
        assert isinstance(only_times, list) or only_times is None, \
               "only_times must be a list or None"
        if not only_times is None:
            assert all([isinstance(item, float) for item in only_times]), \
                   "if only_times is a list, each item must be a float"
        assert isinstance(var_names, list), "var_names must be a list"
        assert all([isinstance(item, str) for item in var_names]), \
               "var_names must be a list of variable names"
        assert x_var is None or isinstance(x_var, list) or isinstance(x_var, numpy.ndarray), \
               "x_var must be None, a list or a numpy array"

        plot_schemes = [scheme for scheme in self._schemes
                        if only_param_tag is None or scheme.tag in only_param_tag]
        if x_var is not None:
            assert len(x_var) == len(plot_schemes), \
                   "x_var must be None or must have the same " + \
                   "length as the number of schemes to plot"
        else:
            x_var = range(len(plot_schemes))

        #Time and values
        times = self.common_times(only_param_tag, None, common_times=1) #All possible times
        plot_times = self.common_times(only_param_tag,
                                       only_times,
                                       common_times=1) #Only selected times
        _, _, series = self.get_series(var_names, slicing, mfactor, only_param_tag, common_times=1)
        if y_var is not None:
            _, _, series_y = self.get_series([y_var], slicing, 1., only_param_tag, common_times=1)

        #line style and color
        def filter_color(var, time):
            "returns the element (var or time) that controls the color"
            if switch_cls:
                return time
            else:
                return var
        def filter_ls(var, time):
            "returns the element (var or time) that controls the style"
            if switch_cls:
                return var
            else:
                return time
        if switch_cls:
            #color by time
            if len(COLORS) >= len(times):
                #We keep the same association between
                #times and colors even if all times are not plotted
                colors = {k:COLORS[i] for i, k in enumerate(times)}
            else:
                colors = {k:COLORS[i] for i, k in enumerate(plot_times)}
            #color by variable
            styles = {k:STYLES[i] for i, k in enumerate(var_names)}
        else:
            #line style by time
            if len(STYLES) >= len(times):
                #We keep the same association between
                #times and styles even if all times are not plotted
                styles = {k:STYLES[i] for i, k in enumerate(times)}
            else:
                styles = {k:STYLES[i] for i, k in enumerate(plot_times)}
            #color by variable
            colors = {k:COLORS[i] for i, k in enumerate(var_names)}

        #Plotting output
        if title is not None:
            ax.set_title(title)
        lines = []
        result = []
        for var in var_names:
            for time in plot_times:
                #time serie of the first scheme for var (all time series are equal)
                time_values = series[(var, plot_schemes[0].tag)][0]
                index_time = numpy.nonzero(time_values==time)[0][0]
                serie = numpy.ma.array([series[(var, scheme.tag)][1][index_time]
                                        for scheme in plot_schemes]).squeeze()
                if y_var is not None:
                    Y = numpy.ma.array([series_y[(y_var, scheme.tag)][1][index_time]
                                        for scheme in plot_schemes])
                dim = len(serie.shape)
                if dim not in [1, 2]:
                    raise NotImplementedError("We only plot 0D and 1D variables, please " +
                                              "use the slicing argument to reduce the problem size")
                if dim == 1:
                    ls_is_str = isinstance(styles[filter_ls(var, time)], str)
                    line, = ax.plot(x_var, serie,
                                    color=colors[filter_color(var, time)],
                                    linestyle=styles[filter_ls(var, time)] if ls_is_str else '',
                                    label=self._var_legend(var, mfactor) + " - t="+ str(time),
                                    **({'linewidth':linewidth} if linewidth is not None else {}))
                    if not ls_is_str:
                        line.set_dashes(styles[filter_ls(var, time)])
                    lines.append(line)
                    result.append(line)
                else:
                    X = numpy.repeat(x_var, serie.shape[-1]).reshape(serie.shape)
                    if y_var is None:
                        Y = numpy.mgrid[0:len(x_var), 0:serie.shape[-1]][1]
                    else:
                        Y = Y.squeeze()
                    if len(plot_times) == 1 and len(var_names) == 1 and enable_contourf:
                        cs = ax.contourf(X, Y, serie, levels=contourlevels)
                        if title is None:
                            ax.set_title(self._var_legend(var, mfactor) + " - t="+ str(time))
                        cs.clabel(inline=1, fontsize=10, colors='black')
                    else:
                        ls_is_str = isinstance(styles[filter_ls(var, time)], str)
                        cs = ax.contour(X, Y, serie,
                                        colors=colors[filter_color(var, time)],
                                        linestyles=styles[filter_ls(var, time)] if ls_is_str
                                                                                else None,
                                        label=self._var_legend(var, mfactor) + " - t="+ str(time),
                                        levels=contourlevels,
                                        **({'linewidths':linewidth} if linewidth is not None
                                                                    else {}))
                        cs.clabel(inline=1, fontsize=10)
                        if not ls_is_str:
                            raise NotImplementedError("plotting of contour lines with dashes " +
                                                      "is not yet implemented")
                            #code below does not work
                            for lc in cs.collections:
                                lc.set_dashes(styles[filter_ls(var, time)])
                        ls_is_str = isinstance(styles[filter_ls(var, time)], str)
                        line, = ax.plot([], [],
                                        color=colors[filter_color(var, time)],
                                        linestyle=styles[filter_ls(var, time)] if ls_is_str else '',
                                        label=self._var_legend(var, mfactor) + " - t="+ str(time),
                                        **({'linewidth':linewidth} if linewidth is not None
                                                                   else {}))
                        if not ls_is_str:
                            line.set_dashes(styles[filter_ls(var, time)])
                        lines.append(line)
                    result.append(cs)
        #Legend
        if len(var_names) == 1 or len(plot_times) == 1:
            if not (dim == 2 and len(plot_times) == 1 and len(var_names) == 1 and enable_contourf):
                ax.legend(handlelength=handlelength, framealpha=framealpha, fontsize=fontsize)
        else:
            if switch_cls:
                lines = []
                for time in plot_times:
                    line, = ax.plot([], [], color=colors[filter_color('', time)],
                                    linestyle='-', label="t=" + str(time),
                                    **({'linewidth':linewidth} if linewidth is not None else {}))
                    lines.append(line)
                ax.legend(handles=lines, title="Colors (times)", loc=1,
                          handlelength=handlelength, framealpha=framealpha, fontsize=fontsize)
                lines = []
                ax2 = ax.twinx()
                ax2.get_xaxis().set_visible(False)
                ax2.get_yaxis().set_visible(False)
                for var in var_names:
                    ls_is_str = isinstance(styles[filter_ls(var, -1)], str)
                    line, = ax2.plot([], [], color='black', label=self._var_legend(var, mfactor),
                                     linestyle=styles[filter_ls(var, -1)] if ls_is_str else '',
                                        **({'linewidth':linewidth} if linewidth is not None
                                                                   else {}))
                    if not ls_is_str:
                        line.set_dashes(styles[filter_ls(var, -1)])
                    lines.append(line)
                ax2.legend(handles=lines, title="Line styles (parameters)", loc=2,
                           handlelength=handlelength, framealpha=framealpha, fontsize=fontsize)
            else:
                lines = []
                for var in var_names:
                    line, = ax.plot([], [], color=colors[filter_color(var, 0)],
                                    linestyle='-', label=self._var_legend(var, mfactor),
                                    **({'linewidth':linewidth} if linewidth is not None
                                                               else {}))
                    lines.append(line)
                ax.legend(handles=lines, title="Colors (parameters)", loc=1,
                          handlelength=handlelength, framealpha=framealpha, fontsize=fontsize)
                lines = []
                ax2 = ax.twinx()
                ax2.get_xaxis().set_visible(False)
                ax2.get_yaxis().set_visible(False)
                for time in plot_times:
                    ls_is_str = isinstance(styles[filter_ls('', time)], str)
                    line, = ax2.plot([], [], color='black', label="t=" + str(time),
                                     linestyle=styles[filter_ls('', time)] if ls_is_str else '',
                                        **({'linewidth':linewidth} if linewidth is not None
                                                                   else {}))
                    if not ls_is_str:
                        line.set_dashes(styles[filter_ls('', time)])
                    lines.append(line)
                ax2.legend(handles=lines, title="Line styles (times)", loc=2,
                           handlelength=handlelength, framealpha=framealpha, fontsize=fontsize)

        ax.set_xlabel("Scheme")
        if dim != 1:
            ax.set_ylabel("Index" if y_var is None else y_var)
        else:
            pass

        return result
