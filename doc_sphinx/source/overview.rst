.. _overview:

Overview
========

General design
--------------

Each of the parameterizations to compare is defined by a class that must be an extension
of an abstract class (:ref:`PPPY <PPPY>`) provided by the package. Actual implementation
of the parameterization can be written in python or must be interfaced to be called from this object.
If the existing implementation is in fortran, a utility (:ref:`ctypesForFortran <ctypesForFortran>`)
is provided to help the interfacing task. It is possible to use other interfacing method
(e.g. f2py) but this is not covered by this documentation. A parameterization is defined
by a time step and the "time integration method" (see the class documentation).

These objects representing parameterizations can be used to perform several comparisons.
Each of these comparisons is described by an instance of the :ref:`PPPYComp <PPPYComp>`
class by specifying the schemes to compare, the length of the simulation and the
initial conditions.

To actually perform the comparison, the parameterizations are run and results are saved in
hdf5 files. Then, some plot functions are provided to visualize the differences between
the parameterizations.

The complete documentation of the library can be found :ref:`here <library>` and examples are
provided in the package (their description can be found :ref:`here <examples>`).


Compilation
-----------
To call fortran code from python using ctypesForFortran, we need a shared library (.so under
linux with object codes compiled using the -fPIC option).
To insure compatibility between fortran and python, fortran routines directly called by the python code
must not use "assumed shape" arrays (arrays declared with ':' as dimensions).
In case some of the routines contain such elements, one should write a wrapper around these
routines to declare explicit dimensions.

For example, if the parameterization is this fortran subroutine:

.. code-block:: fortran

    SUBROUTINE PARAM(X, Y)
    REAL, DIMENSION(:, :), INTENT(IN) :: X
    REAL, DIMENSION(:, :), INTENT(OUT) :: Y
    Y=X+1
    END SUBROUTINE PARAM

It is necessary to write, compile and use this wrapper:

.. code-block:: fortran

    SUBROUTINE PARAM_PY(X, Y, I1, I2)
    INTEGER, INTENT(IN) :: I1, I2
    REAL, DIMENSION(I1, I2), INTENT(IN) :: X
    REAL, DIMENSION(I1, I2), INTENT(OUT) :: Y
    CALL PARAM(X, Y)
    END SUBROUTINE PARAM_PY

In the "examples" directory, several comparison examples can be found.
In each of these examples, there is a lib directory which contain the explanations
and the extra source code needed to compile AROME, Meso-NH or WRF to use some of their
parameterizations.

The compilation stage is certainly the most tricky with some traps (big/little endian,
size of the variables, fortran/C compatibility...). Interfacing examples using the
ctypesForFortran utility can be found in the :ref:`module <ctypesForFortran>` (and a trivial
example is in the examples/test/lib directory).

Parameterization implementation
-------------------------------
A class (inherited from the PPPY class) is written for each of the parameterizations.
During the instantiation of this class, the time step, the integration method, its
identification and some possible options are given.
The class contains a __call__ method which calls the different methods in this order:
setup, build_init_state, execute (called in a loop) then finalize.

The methods that can or must be implemented to represent a parameterization are described below
(see examples/test/pppy_param1.py for a trivial example):

__init__
++++++++
This method can be implemented to deal with optional arguments which represent options for the
parameterization. In this case, optional arguments are defined in the argument list of the method
and a call to super().__init__ must be done.

setup
+++++
This method do the initialization part that cannot be done earlier (__init__) or need to be done
again before each of the execution.
In the provided examples, this is the place where signatures of fortran routines are defined,
where the shared library is opened and where the initialization of the fortran modules are done.
To suppress interferences between parameterizations, each computation is done in a separate
process, the setup method is called in this sub-process whereas the __init__ method is called
in the main process.

finalize
++++++++
This method can be useful to clean memory or disk after having run a simulation. 

build_init_state
++++++++++++++++
The PPPYComp instance that performs the comparison calls each of the parameterizations
with the same initial state (this is a dictionary whose keys are the variable names and values
are classically numpy arrays). This method allows the parameterization to adapt the content of the
initial state to its particular need. The parameterization must modify, during the simulation,
the variables that are in the initial state but can follow, in addition, other variables.
All these new variables must be added in this method so that the output file is dimensioned
accordingly.

A concrete example: we want to compare microphysical schemes in warm conditions. Therefore, we
set values for the vapor, the cloud water and the rain contents. Each parameterization computes
the time evolution of these variables and values are compared at the end. If one of the
parameterizations is, in fact, a mixed scheme which needs values for ice species, this method is
the place to create corresponding numpy arrays and fill them with zeros. 

execute
+++++++
This method receive three arguments: the current state (dictionary containing the different
state variables) the time step length to use and the number of the current time step. 
The method must return the new state after having applied the parameterization during the
given time step. All variables that are in the current state but not in the output state
are considered to be constant over the time step duration.

In this method, it can be necessary to modify the shape of the arrays to fulfill the
requirements of the parameterization but variables in the output state must keep the
same shape as those received in the current state. For example, one can need to transform
a 1D array of length 1 into a 2D or 3D array with a border to feed the parameterization with.
Other conversions may be needed, such as change of unit, of physical variables (temperature
versus potential temperature) or memory representation (32bit versus 64bit...).

Comparison
----------
The goal is to compare several parameterizations.
A parameterization is defined with a PPPY instance then the comparison itself is described
by a PPPYComp instance which controls the execution and the comparison.

Parameterization definition
+++++++++++++++++++++++++++
The PPPY instance is created by defining the time step length, the time integration method,
a name (for the legends on plots) and a tag (a string used for filenames).
If we want to compare a parameterization using several time step lengths and/or
time integration methods, we must define several instances of a same PPPY class.

Synthetic example:

.. code-block:: python

    class param1(PPPY): pass #class definition
    class param2(PPPY): pass #class definition
    param1_dt1 = param1(dt=1, method='step-by-step', name='param 1, dt=1', tag='param1_1')
    param1_dt5 = param1(dt=5, method='step-by-step', name='param 1, dt=5', tag='param1_5')
    param2_dt1 = param2(dt=1, method='step-by-step', name='param 2, dt=1', tag='param2_1')
    param2_dt5 = param2(dt=5, method='step-by-step', name='param 2, dt=5', tag='param2_5')

The time step length is expressed in seconds, an output is computed and stored
every time step. The time integration method can take two values:

- ‘step-by-step’: This is the classical time integration method;
  the state computed after one time step is used as the beginning state for the next iteration;
- ‘one-step’: each output is computed by an integration starting from the same initial state
  with varying time step length.

Comparison definition
+++++++++++++++++++++
The PPPYComp instance is created with a list of PPPY instances, a simulation duration, and the initial state.
A name (used for plot legends) and a tag (used for filenames) is also associated to the
comparison.

Synthetic example:

.. code-block:: python

    comp = PPPYComp(schemes=[param1_dt1, param1_dt5, param2_dt1, param2_dt5],
                    output_dir= <directory_for_output>,
                    duration= <simulation_length_in_seconds>,
                    init_state={'var1' :..., 'var2' :..., ...},
                    name=<name_of_the_experiment>,
                    tag=<identifying_string_used_in_filenames>)

The execution of the different parameterizations (and the storage of the results)
is performed by the run method.

Diagnostics
+++++++++++
It is possible to compute diagnostics from the simulated values in order to
allow the plotting of these diagnostics.
This is done by extending the get_scheme_output method of the PPPYComp class
as in this example:

..  code-block:: python

    class MyComp(PPPYComp):
        def get_scheme_output(self, scheme, var):
            if var == 'rt':
                return numpy.array([self.get_scheme_output(scheme, v)
                                    for v in ['rv', 'rc', 'ri', 'rr', 'rs', 'rg', 'rh']]
                                  ).sum(axis=0)
            else:
                return super().get_scheme_output(scheme, var)

Plotting
++++++++
It is possible to implement diagnostics and plotting functions by creating a new
PPPYComp class (inherited from PPPYComp) or outside of this class by
directly reading the hdf5 file.

Two plotting functions already exist in the PPPYComp class:

- plot_evol: plots a time evolution for one or several parameterizations and for
  one or several variables in 1D (curves) or 2D (contour plots).
- plot_comp: plots a "parameterization evolution" for one or several output times
  and for one or several variables in 1D (curves) or 2D (contour plots). This is
  like the plot_evol plots apart that time and parameterization dimensions are exchanged.
  This can be useful when a parameter of the parameterization can take several values; we
  then plot the output with respect to this parameter.

A reasonable number of options and arguments are available through both methods.
The names and units of the different variables, that are
printed in the legend of the plots, can be modified through the VAR_NAME and
VAR_UNIT constants. Example:

.. code-block:: python

	import pppy
	pppy.VAR_NAME['new_var'] = "My super variable"
	pppy.VAR_UNIT['new_var'] = "J/K/m/s/year"

Moreover, the plot_multi method is a utility simplifying the plotting of several plots
over a same figure (matplotlib subplots).







