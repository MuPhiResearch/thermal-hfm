Examples
========

The ``examples/`` directory contains runnable scripts that demonstrate the
module's capabilities.  Each script is self-contained and can be executed
directly in MATLAB.


2D Thermal EDFM
---------------

``example_2d_thermal_edfm.m``: two intersecting fractures in a 2D Cartesian
grid, with hot-water injection on the west and production on the east.  It
demonstrates the basic 2D workflow of grid generation, NNC setup, thermal
preprocessing, model construction, and simulation.


3D Thermal EDFM
---------------

``example_3d_thermal_edfm.m``: a single vertical fracture plane in a 3D domain,
using the 3D ``hfm`` pipeline (``preProcessingFractures``,
``assembleGlobalGrid``, ``computeEffectiveTrans``).


DFM vs EDFM comparison
----------------------

``example_dfm_vs_edfm_comparison.m``: a side-by-side comparison of a fully
resolved Discrete Fracture Model (unstructured Delaunay mesh) against the EDFM
approach on the same single-fracture problem.


EDFM vs pEDFM (2D and 3D)
-------------------------

``example_edfm_vs_pedfm_2d.m`` and ``example_edfm_vs_pedfm_3d.m``: compare
standard EDFM and pEDFM on a domain with a barrier fracture and a conduit
fracture.  They demonstrate that pEDFM correctly blocks flow across
low-permeability barriers where EDFM fails.


Brine EOS (EDFM and pEDFM)
--------------------------

``example_edfm_brine_eos.m`` and ``example_pedfm_brine_eos.m``: thermal
simulation with ``CompositionalBrineFluid`` (H2O + NaCl), the Spivey EOS for
temperature- and pressure-dependent density and viscosity, and dynamic fluid
thermal conductivity via ``EDFMDynamicHeatTransmissibility``.


Thermal doublet
---------------

``example_thermal_doublet.m``: a realistic geothermal doublet (cold injection
and hot production) in a heterogeneous 3D reservoir with log-normal
permeability, two barrier fractures, and two conduit fractures.  A 20-year
simulation that compares EDFM and pEDFM results.
