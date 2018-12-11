Introduction
============

The ``pppy`` package facilitates the comparison of physical parameterizations in an "off-line" way.
The package controls the execution of each of the parameterizations and, then, compares them.

What can we do with this?
-------------------------
At the beginning, this tool was written to compare the output of a single microphysical
scheme when using different time steps and a single point (0D simulation). Then, the tool
was applied to compare different microphysical schemes, still on a single point.
The tool was also used to compare sedimentation (vertical advection of hydrometeors) in 1D simulations.

The tool is written without reference to a given set of variables (except two constants that are set in
the code to provide descriptions and units for well known microphysical variables, but this can
be easily extended to other variables, and the tool is totally usable with unknown variables).
It can be used on simulations with greater number of dimensions than the 0D and 1D cases listed above.

How can I use it?
-----------------
The tool comes with some examples (their description can be found :ref:`here <examples>`).
These examples contain instructions for compilation and can be used directly.
And you can extend the tool with other parameterizations (see the :ref:`overview <overview>`
to begin).
