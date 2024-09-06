.. _examples:

Examples
========

Several examples are provided in the ``example`` directory in the package folder.

- ``test``: this example is a very simple example that demonstrates how the package can be used.
- ``pyphyex``: examples with the parametrisations from PHYEX (https://github.com/UMR-CNRM/PHYEX) involving
  the ICE3/ICE4 and LIMA microphysics, the saturation adjustement used with these microphysics schemes,
  the sedimentation schemes used with ICE3/ICE4, the shallow convection and the turbulence schemes.
  To be used, you need to compile the PHYEX source code.
- ``old``: old example using the WRF microphysics schemes (no more maintained)

To run one of these examples, you must copy the entire directory somewhere, then go
reading the README file in the directory (or the lib subdirectory) that describes how to compile and link
a shared library. Then, you can execute the comp_*.py script to perform the comparison.

test
----
This example uses two parameterizations (param1 and param2) written in fortran that do almost
nothing. In the lib directory, a README file describes how to compile with gfortran the source
codes provided.

This is a very simple example but a good starting point before trying to interface a
real fortran subroutine. It uses the ctypesForFortran utility.

pyphyex
---------
The README file describes hot to compile the PHYEX package.
The same shared lib contain the entry points for the different tests:

- microphy: example using the ICE3/ICE4 microphysics scheme
- sedimentation: the sedimentation schemes of ICE3/ICE4
- ice_adjust: the saturation adjustment used with ICE3/ICE4
- shallow: the shallow convection scheme
- turb: the turbulance scheme
- lima_adjust: the saturation adjustment used with the LIMA scheme
- lima: the LIMA microphysics scheme
