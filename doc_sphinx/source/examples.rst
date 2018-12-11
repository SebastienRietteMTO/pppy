.. _examples:

Examples
========

Four examples are provided in the ``example`` directory in the package folder.

- ``test``: this example is a very simple example that demonstrates how the package can be used.
- ``microphy1``: this is an example of use of microphysical parameterizations coming from AROME and Meso-NH
  models.
- ``microphy2``: this is an example of use of microphysical parameterizations coming from WRF and Meso-NH
  models (both model can be downloaded for free)
- ``sedimentation_AROME``: this is an example of 1D comparisons using sedimentation parameterizations from the
  AROME model.
- ``sedimentation_MNH``: this is the same example as sedimentation_AROME but build upon the
  Meso-NH model.

All these examples involve parameterizations from meteorological models (AROME, WRF and Meso-NH);
you need a copy of the corresponding source code to be able to run the examples.

To run one of these examples, you must copy the entire directory somewhere, then go
reading the README file in the lib directory that describes how to compile and link
a shared library. Then, you can execute the comp_*.py script to perform the comparison.

test
----
This example uses two parameterizations (param1 and param2) written in fortran that do almost
nothing. In the lib directory, a README file describes how to compile with gfortran the source
codes provided.

This is a very simple example but a good starting point before trying to interface a
real fortran subroutine. It uses the ctypesForFortran utility.

microphy1
---------
This is an example kept here as a reference on how to interface the AROME microphysical
parameterization with this tool. To play with microphysical parameterizations, the following
example is more interesting.

microphy2
---------
This example can be run almost everywhere as the two models involved are freely available
(Meso-Nh and WRF).
The lib directory provides the source code modifications that must be applied to the
official release of these models and the compilation instructions.
Not all the microphysical parameterizations of the WRF model are interfaced but the
extension to another parameterization must be relatively easy.

sedimentation_AROME and sedimentation_MNH
-----------------------------------------
This are examples of 1D simulations using the AROME or the Meso-NH model.
They use the sedimentation algorithms available for ICE3/ICE4 microphysical scheme.
