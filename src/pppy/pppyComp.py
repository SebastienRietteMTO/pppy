#!/usr/bin/env python3

# Copyright (c) Météo France (2018-)
# This software is governed by the CeCILL-C license under French law.
# http://www.cecill.info

"""
The PPPYComp class is used to run the different parameterizations (instances of PPPY),
save the results and plot some variables.
Variables are defined using a string (e.g. temperature is commonly associated
with the 'T' string which is used as a key in the sate dictionary) but the plot
methods try to build legends using the long name of the variables. To do this, the
variable name is searched in the VAR_NAME dictionary and, if found, the long name
is set to the associated key, otherwise, the long name is set to the name.
The same mechanism is used for the unit of the variable which is searched in the
VAR_UNIT dictionary.
Both dictionaries can be modified before plotting to configure the long name
and the unit to use in the plot legends.

>>> import pppy
>>> pppy.VAR_NAME['T'] = "Temperature"
>>> pppy.VAR_UNIT['T'] = "K"
"""

import numpy
import os
import h5py
import copy
import multiprocessing
import logging
from .pppy import PPPY
from pppy import VAR_NAME, VAR_UNIT, COLORS, STYLES

class PPPYComp():
    """
    This class is used to perform a comparison between several implementation
    of a parameterization.
    This class can be directly instantiated, there is no need to extend it,
    unless you need to write some diagnostic methods (e.g. plots, statistics...).
    In particular, the get_scheme_output method can be extended to compute
    diagnostics from the variable computed by a scheme and make these
    diagnostics available for plotting.
    """

    def __init__(self, schemes, output_dir, duration, init_state, name, tag):
        """
        At instantiation, we define the schemes to compare, the total duration of the
        time integration, the initial conditions and the name of the experiment.

        :param schemes:         a list of PPPY instances
        :param output_dir:      directory in which results are stored
        :param duration:        total duration of the time integration (in seconds)
        :param init_state:      dictionary holding the initial state
        :param name:            name of the comparison experiment used in legends
        :param tag:             string to identify uniquely the experiment, the
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
        assert isinstance(name, str) and isinstance(tag, str), \
               "name and tag must be strings"

        values = list(init_state.values())
        assert all([isinstance(value, numpy.ndarray) for value in values]), \
               "All init_state item must be ndarrays"

        if not os.path.exists(output_dir):
            raise IOError("output dir does not exist")
        self._output_dir = os.path.join(output_dir, tag)
        if not os.path.exists(self._output_dir):
            os.mkdir(self._output_dir)

        self._schemes = schemes
        self._duration = duration
        self._init_state = init_state
        self.name = name
        self.tag = tag
        self._files = {}

    def run(self, only_param_tag=None, force=False):
        """
        This methods runs the different parameterizations.

        :param only_param_tag: (optional) list of parameterization tags to actually run.
                               By default (None), all parameterizations are run.
        :param force:          (optional) By default (False), parameterizations are run only
                               if needed (output file for same tag and
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
        if unit is not None:
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
            fig.suptitle(self.name)
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
        if common_times == 1:
            times = []
            for scheme in schemes:
                times.extend(list(self.get_scheme_output(scheme, 't')))
            times = sorted(list(set(times)))
            times = [t for t in times if only_times is None or t in only_times]
        elif common_times == 2:
            times = []
            for scheme in schemes:
                times = [t for t in list(self.get_scheme_output(scheme, 't')) if t in times]
            times = sorted(times)
            times = [t for t in times if only_times is None or t in only_times]

        return times

    def get_scheme_output(self, scheme, var):
        """
        This method returns the output of a scheme
        :param scheme: scheme
        :param val: variable to read
        :returns: The corresponding array
        """
        filename = os.path.join(self._output_dir, scheme.tag + '.hdf5')
        if not os.path.exists(filename):
            raise IOError(f'Filename does not exist: {filename}')
        if scheme.tag not in self._files:
            try:
                f =  h5py.File(filename, 'r')
            except OSError:
                #This error can be due to too many files opened
                self.close_files()
                f =  h5py.File(filename, 'r')
            self._files[scheme.tag] =  f
        if var in self._files[scheme.tag]:
            return self._files[scheme.tag][var][...]
        else:
            print(self._files[scheme.tag].items())
            raise ValueError(("This variable (%s) is not present in " + \
                              "this output file (%s)") % (var, filename))

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

        result = {}
        for scheme in schemes:
            simul_time = self.get_scheme_output(scheme, 't')
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
                simul_values = self.get_scheme_output(scheme, var) * mfactor
                if slicing is not None:
                    simul_values = simul_values[:, slicing]
                shape = tuple([len(time)] + list(simul_values.shape[1:]))
                serie = numpy.ma.masked_all(shape)
                serie[index_result] = simul_values[index_simul]
                result[(var, scheme.tag)] = (time, serie)

        return var_names, [scheme.tag for scheme in schemes], result

    def plot_evol(self, ax, var_names, slicing=None, y_var=None,
                  mfactor=1., only_param_tag=None, enable_contourf=True,
                  contourlevels=None, title=None,
                  switch_cls=False, linewidth=None,
                  xlabel=None, ylabel=None, colorbar=True,
                  clabel=False, xlim=None, ylim=None,
                  legend='min', colors=None, styles=None,
                  kwargs_contourf=None, kwargs_contour=None, kwargs_plot=None,
                  kwargs_legend=None):
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
        :param title:           (optional) Title to use on subplot, if None title is automatic
        :param switch_cls:      (optional, default to False) By default, A line style is assigned
                                to each parameterization and a color is assign to each variable.
                                Setting this variable to True, reverse this default behavior.
        :param linewidth:       Control the line width used in plots
        :param xlabel:          To impose an x-axis label (None to generate it automatically)
        :param ylabel:          To impose an y-axis label (None to generate it automatically)
        :param colorbar:        In case of contourf use, adds a colorbar
        :param clabel:          To add labels to line contours
        :param xlim:            tuple (xmin, xmax) or None
        :param ylim:            tuple (ymin, ymax) or None
        :param legend:          'min' to label curves with the minimum needed to distinguish them
                                'full' to label curves with the entire information
                                'no' to not plot legend
        :param colors:          list of color names to use, None to use the default list
        :param styles:          list of line styles to use, None to use the default list
                                a line style is a string (recognized by matplotlib) or a
                                even length tuple of on and off ink in points
        :param kwargs_contourf, kwargs_contour and kwargs_plot: optional dictionaries of
                                arguments to be passed to the contourf, contour and plot methods
        :param kwargs_legend:   optional dictionaries of arguments to be passed to ax.legend

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
        assert legend in ['min', 'full', 'no'], "legend must be 'min', 'full' or 'no'"

        kwargs_contourf = {} if kwargs_contourf is None else kwargs_contourf
        kwargs_contour = {} if kwargs_contour is None else kwargs_contour
        kwargs_plot = {} if kwargs_plot is None else kwargs_plot
        kwargs_legend = {} if kwargs_legend is None else kwargs_legend

        plot_schemes = [scheme for scheme in self._schemes
                        if only_param_tag is None or scheme.tag in only_param_tag]

        #line style and color
        colors = COLORS if colors is None else colors
        styles = STYLES if styles is None else styles
        def get_color(var, scheme):
            "returns the color associated to var or scheme"
            if switch_cls:
                #color by parameterization
                if len(colors) >= len(self._schemes):
                    #We keep the same association between
                    #schemes and colors even if all schemes are not plotted
                    return {k:colors[i] for i, k in enumerate([sch.tag
                                                               for sch in self._schemes])}[scheme]
                else:
                    return {k:colors[i] for i, k in enumerate([sch.tag
                                                               for sch in plot_schemes])}[scheme]
            else:
                #color by variable
                return {k:colors[i] for i, k in enumerate(var_names)}[var]
        def get_ls(var, scheme):
            "returns the line style associated to var or scheme"
            if switch_cls:
                #line style by variable
                return {k:styles[i] for i, k in enumerate(var_names)}[var]
            else:
                #line style by parameterization
                if len(styles) >= len(self._schemes):
                    #We keep the same association between
                    #schemes and styles even if all schemes are not plotted
                    return {k:styles[i] for i, k in enumerate([sch.tag
                                                               for sch in self._schemes])}[scheme]
                else:
                    return {k:styles[i] for i, k in enumerate([sch.tag
                                                               for sch in plot_schemes])}[scheme]

        #Plotting output
        _, _, series = self.get_series(var_names, slicing, mfactor, only_param_tag, common_times=0)
        if y_var is not None:
            _, _, series_y = self.get_series([y_var], slicing, 1., only_param_tag, common_times=0)
        if title is not None:
            ax.set_title(title)
        result = []
        handles_var = {}
        handles_scheme = {}
        plot_nb = 0
        for var in var_names:
            for scheme in plot_schemes:
                plot_nb += 1
                time, serie = series[(var, scheme.tag)]
                shape = serie.shape
                dim = sum([1 if l > 1 else 0 for l in shape])
                if dim not in [1, 2]:
                    raise NotImplementedError("We only plot 0D and 1D variables, please " +
                                              "use the slicing argument to reduce the problem size")

                label = self._var_legend(var, mfactor) + " - " + scheme.name
                if len(var_names) == len(plot_schemes) == 1:
                    if title is None:
                        ax.set_title(self._var_legend(var, mfactor) + " - " + scheme.name)
                    if legend == 'min':
                        legend = 'no'
                elif len(var_names) == 1:
                    if title is None:
                        ax.set_title(self._var_legend(var, mfactor))
                    if legend != 'full':
                        label = scheme.name
                elif len(plot_schemes) == 1:
                    if title is None:
                        ax.set_title(scheme.name)
                    if legend != 'full':
                        label = self._var_legend(var, mfactor)

                if dim == 1:
                    ls_is_str = isinstance(get_ls(var, scheme.tag), str)
                    if linewidth is not None:
                        kwargs_plot['linewidth'] = linewidth
                    line, = ax.plot(time, serie.squeeze(),
                                    color=get_color(var, scheme.tag),
                                    linestyle=get_ls(var, scheme.tag) if ls_is_str else '',
                                    label=label,
                                    **kwargs_plot)
                    if not ls_is_str:
                        line.set_dashes(get_ls(var, scheme.tag))
                    handles_var[var] = line
                    handles_scheme[scheme.name] = line
                    result.append(line)
                else:
                    X = numpy.repeat(time, serie.squeeze().shape[-1]).reshape(serie.squeeze().shape)
                    if y_var is None:
                        Y = numpy.mgrid[0:len(time), 0:serie.squeeze().shape[-1]][1]
                    else:
                        _, Y = series_y[(y_var, scheme.tag)]
                        Y = Y.squeeze()
                    if len(plot_schemes) == 1 and len(var_names) == 1 and enable_contourf:
                        qcs = ax.contourf(X, Y, serie.squeeze(), levels=contourlevels,
                                          **kwargs_contourf)
                        if clabel:
                            qcs.clabel(inline=1, fontsize=10, colors='black')
                        if colorbar:
                            ax.figure.colorbar(qcs, ax=ax)
                    else:
                        ls_is_str = isinstance(get_ls(var, scheme.tag), str)
                        if linewidth is not None:
                            kwargs_contour['linewidths'] = linewidth
                        qcs = ax.contour(X, Y, serie.squeeze(),
                                         colors=get_color(var, scheme.tag),
                                         linestyles=get_ls(var, scheme.tag) if ls_is_str else None,
                                         levels=contourlevels,
                                         **kwargs_contour)
                        lc0 = qcs.collections[0]
                        lc0.set_label(label)
                        handles_var[var] = lc0
                        handles_scheme[scheme.name] = lc0
                        if clabel:
                            qcs.clabel(inline=1, fontsize=10)
                        if not ls_is_str:
                            for line in qcs.collections:
                                line.set_dashes([(0., get_ls(var, scheme.tag))])
                    result.append(qcs)
        #Legend
        if plot_nb > 0 and legend != 'no':
            if len(var_names) == 1 or len(plot_schemes) == 1:
                if not (dim == 2 and len(plot_schemes) == len(var_names) == 1 and enable_contourf):
                    ax.legend(**kwargs_legend)
            else:
                if switch_cls:
                    handles1 = [handles_scheme[scheme.name] for scheme in plot_schemes]
                    labels1 = [scheme.name for scheme in plot_schemes]
                    handles2 = [handles_var[var] for var in var_names]
                    labels2 = [self._var_legend(var, mfactor) for var in var_names]
                    title1 = "Schemes (colors)"
                    title2 = "Parameters (line styles)"
                else:
                    handles1 = [handles_var[var] for var in var_names]
                    labels1 = [self._var_legend(var, mfactor) for var in var_names]
                    handles2 = [handles_scheme[scheme.name] for scheme in plot_schemes]
                    labels2 = [scheme.name for scheme in plot_schemes]
                    title1 = "Parameters (colors)"
                    title2 = "Schemes (line styles)"
                legend1 = ax.legend(handles=handles1, labels=labels1,
                                    title=title1, loc=1, **kwargs_legend)
                for handle in legend1.legend_handles:
                    handle.set_linestyle('-')
                legend2 = ax.legend(handles=handles2, labels=labels2,
                                    title=title2, loc=2, **kwargs_legend)
                for handle in legend2.legend_handles:
                    handle.set_color('black')
                ax.add_artist(legend1)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        ax.set_xlabel("Time" if xlabel is None else xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        elif plot_nb > 0 and dim != 1:
            ax.set_ylabel("Index" if y_var is None else y_var)

        return result

    def plot_comp(self, ax, var_names, slicing=None, y_var=None,
                  mfactor=1., only_param_tag=None, enable_contourf=True,
                  contourlevels=None, title=None, x_var=None, only_times=None,
                  switch_cls=False, linewidth=None,
                  xlabel=None, ylabel=None, remove_ma=False, colorbar=True,
                  clabel=False, xlim=None, ylim=None,
                  legend='min', colors=None, styles=None,
                  kwargs_contourf=None, kwargs_contour=None, kwargs_plot=None,
                  kwargs_legend=None):
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
        :param title:           (optional) Title to use on subplot, if None title is automatic
        :param x_var:           (optional) By default, the schemes are evenly spaced on the x-axis.
                                If this variable is set, it must be a list of position to use
                                (specified in the same order as schemes).
        :param only_times:      (optional) By default, all output times are plotted.
                                If set, plotting is limited to times listed.
        :param switch_cls:      (optional, default to False) By default, A line style is assigned
                                to each time and a color is assign to each variable.
                                Setting this variable to True, reverse this default behavior.
        :param linewidth:       Control the line width used in plots
        :param xlabel:          To impose an x-axis label (None to generate it automatically)
        :param ylabel:          To impose an y-axis label (None to generate it automatically)
        :param remove_ma:       If True, masked data are removed before plotting. This enable
                                a smoother plot but it is then not guaranteed that all
                                superimposed plots are plotted for the same values (if one scheme
                                does not provide output for a given time whereas other schemes do)
        :param colorbar:        In case of contourf use, adds a colorbar
        :param clabel:          To add labels to line contours
        :param xlim:            tuple (xmin, xmax) or None
        :param ylim:            tuple (ymin, ymax) or None
        :param legend:          'min' to label curves with the minimum needed to distinguish them
                                'full' to label curves with the entire information
                                'no' to not plot legend
        :param colors:          list of color names to use, None to use the default list
        :param styles:          list of line styles to use, None to use the default list
                                a line style is a string (recognized by matplotlib) or a
                                even length tuple of on and off ink in points
        :param kwargs_contourf, kwargs_contour and kwargs_plot: optional dictionaries of
                                arguments to be passed to the contourf, contour and plot methods
        :param kwargs_legend:   optional dictionaries of arguments to be passed to ax.legend

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
        assert legend in ['min', 'full', 'no'], "legend must be 'min', 'full' or 'no'"

        kwargs_contourf = {} if kwargs_contourf is None else kwargs_contourf
        kwargs_contour = {} if kwargs_contour is None else kwargs_contour
        kwargs_plot = {} if kwargs_plot is None else kwargs_plot
        kwargs_legend = {} if kwargs_legend is None else kwargs_legend

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
        colors = COLORS if colors is None else colors
        styles = STYLES if styles is None else styles
        def get_color(var, time):
            "returns the color associated to var or time"
            if switch_cls:
                #color by time
                if len(colors) >= len(times):
                    #We keep the same association between
                    #times and colors even if all times are not plotted
                    return {k:colors[i] for i, k in enumerate(times)}[time]
                else:
                    return {k:colors[i] for i, k in enumerate(plot_times)}[time]
            else:
                #color by variable
                return {k:colors[i] for i, k in enumerate(var_names)}[var]
        def get_ls(var, time):
            "returns the line style associated to var or time"
            if switch_cls:
                #line style by variable
                return {k:styles[i] for i, k in enumerate(var_names)}[var]
            else:
                #line style by time
                if len(styles) >= len(times):
                    #We keep the same association between
                    #times and styles even if all times are not plotted
                    return {k:styles[i] for i, k in enumerate(times)}[time]
                else:
                    return {k:styles[i] for i, k in enumerate(plot_times)}[time]

        #Plotting output
        if title is not None:
            ax.set_title(title)
        result = []
        handles_var = {}
        handles_time = {}
        plot_nb = 0
        for var in var_names:
            for time in plot_times:
                plot_nb += 1
                #time serie of the first scheme for var (all time series are equal)
                time_values = series[(var, plot_schemes[0].tag)][0]
                index_time = numpy.nonzero(time_values==time)[0][0]
                serie = numpy.ma.array([series[(var, scheme.tag)][1][index_time]
                                        for scheme in plot_schemes]).squeeze()
                if y_var is not None:
                    Y = numpy.ma.array([series_y[(y_var, scheme.tag)][1][index_time]
                                        for scheme in plot_schemes])
                x_var_plt = numpy.array(x_var)
                if remove_ma:
                    mask = numpy.array([numpy.all(numpy.ma.getmaskarray(serie[i]))
                                        for i in range(len(plot_schemes))])
                    mask = numpy.logical_not(mask)
                    serie = serie[mask]
                    x_var_plt = x_var_plt[mask]
                    if y_var is not None:
                        Y = Y[mask]

                dim = len(serie.shape)
                if dim not in [1, 2]:
                    raise NotImplementedError("We only plot 0D and 1D variables, please " +
                                              "use the slicing argument to reduce the problem size")

                label = self._var_legend(var, mfactor) + " - t="+ str(time)
                if len(var_names) == len(plot_times) == 1:
                    if title is None:
                        ax.set_title(self._var_legend(var, mfactor) + " - t=" + str(time))
                    if legend == 'min':
                        legend = 'no'
                elif len(var_names) == 1:
                    if title is None:
                        ax.set_title(self._var_legend(var, mfactor))
                    if legend != 'full':
                        label = "t=" + str(time)
                elif len(plot_times) == 1:
                    if title is None:
                        ax.set_title("t=" + str(time))
                    if legend != 'full':
                        label = self._var_legend(var, mfactor)

                if dim == 1:
                    ls_is_str = isinstance(get_ls(var, time), str)
                    if linewidth is not None:
                        kwargs_plot['linewidth'] = linewidth
                    line, = ax.plot(x_var_plt, serie,
                                    color=get_color(var, time),
                                    linestyle=get_ls(var, time) if ls_is_str else '',
                                    label=label,
                                    **kwargs_plot)
                    if not ls_is_str:
                        line.set_dashes(get_ls(var, time))
                    handles_var[var] = line
                    handles_time[time] = line
                    result.append(line)
                else:
                    X = numpy.repeat(x_var_plt, serie.shape[-1]).reshape(serie.shape)
                    if y_var is None:
                        Y = numpy.mgrid[0:len(x_var_plt), 0:serie.shape[-1]][1]
                    else:
                        Y = Y.squeeze()
                    if len(plot_times) == 1 and len(var_names) == 1 and enable_contourf:
                        qcs = ax.contourf(X, Y, serie, levels=contourlevels,
                                          **kwargs_contourf)
                        if clabel:
                            qcs.clabel(inline=1, fontsize=10, colors='black')
                        if colorbar:
                            ax.figure.colorbar(qcs, ax=ax)
                    else:
                        ls_is_str = isinstance(get_ls(var, time), str)
                        if linewidth is not None:
                            kwargs_contour['linewidths'] = linewidth
                        qcs = ax.contour(X, Y, serie,
                                         colors=get_color(var, time),
                                         linestyles=get_ls(var, time) if ls_is_str else None,
                                         levels=contourlevels,
                                         **kwargs_contour)
                        lc0 = qcs.collections[0]
                        lc0.set_label(label)
                        handles_var[var] = lc0
                        handles_time[time] = lc0
                        if clabel:
                            qcs.clabel(inline=1, fontsize=10)
                        if not ls_is_str:
                            for line in qcs.collections:
                                line.set_dashes([(0., get_ls(var, time))])
                    result.append(qcs)
        #Legend
        if plot_nb > 0 and legend != 'no':
            if len(var_names) == 1 or len(plot_times) == 1:
                if not (dim == 2 and len(plot_times) == len(var_names) == 1 and enable_contourf):
                    ax.legend(**kwargs_legend)
            else:
                if switch_cls:
                    handles1 = [handles_time[time] for time in plot_times]
                    labels1 = ["t=" + str(time) for time in plot_times]
                    handles2 = [handles_var[var] for var in var_names]
                    labels2 = [self._var_legend(var, mfactor) for var in var_names]
                    title1 = "Times (colors)"
                    title2 = "Parameters (line styles)"
                else:
                    handles1 = [handles_var[var] for var in var_names]
                    labels1 = [self._var_legend(var, mfactor) for var in var_names]
                    handles2 = [handles_time[time] for time in plot_times]
                    labels2 = ["t=" + str(time) for time in plot_times]
                    title1 = "Parameters (colors)"
                    title2 = "Times (line styles)"
                legend1 = ax.legend(handles=handles1, labels=labels1,
                                    title=title1, loc=1, **kwargs_legend)
                for handle in legend1.legend_handles:
                    handle.set_linestyle('-')
                legend2 = ax.legend(handles=handles2, labels=labels2,
                                    title=title2, loc=2, **kwargs_legend)
                for handle in legend2.legend_handles:
                    handle.set_color('black')
                ax.add_artist(legend1)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        ax.set_xlabel("Scheme" if xlabel is None else xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        elif plot_nb > 0 and dim != 1:
            ax.set_ylabel("Index" if y_var is None else y_var)

        return result

    def close_files(self):
        """
        This methods close the different hdf5 opened files
        """
        for hdf5_file in list(self._files.keys()):
            self._files[hdf5_file].close()
            del self._files[hdf5_file]

    def __del__(self):
        try:
            self.close_files()
        except:
            pass
