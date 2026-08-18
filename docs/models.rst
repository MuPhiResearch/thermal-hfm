Model Class
===========

The thermal-hfm module provides a single model class, ``GeothermalHFMModel``,
which extends MRST's ``GeothermalModel`` with embedded and projection-based
discrete fracture (EDFM and pEDFM) support for coupled flow and heat transfer.

EDFM and pEDFM are **not** separate classes.  They are two configurations of
the same model, differing only in how the flow operators are built.  That
choice is driven by the **grid**, specifically whether the projection
preprocessing (``pMatFracNNCs3D``) has been run, which adds the projected
connection field ``G.nnc.pMMneighs``.  It is selected through the
``fractureMethod`` property.


Class architecture
------------------

.. only:: html

   .. code-block:: text

      ReservoirModel                       (ad-core)
        |
        +-- GeothermalModel                (geothermal)
              |
              +-- GeothermalHFMModel       fractureMethod: auto | edfm | pedfm
                    |
                    |  setupOperators()  dispatches by fractureMethod:
                    |     edfm  -> setupEDFMOperatorsTPFA  (hfm)
                    |     pedfm -> setupPEDFMOpsTPFA        (shale)
                    |
                    +-- setupStateFunctionGroupings() installs:
                          * EDFMDynamicFlowTransmissibility
                          * EDFMDynamicHeatTransmissibility
                          * EDFMDynamicMolecularTransmissibility   (all <- StateFunction)

      Preprocessing (set G.nnc.transHr / transHf / CI_dyn):
        computeThermalNNCTransFracMatrix , computeThermalNNCTransFracFrac

   The dynamic state functions self-detect pEDFM via ``G.nnc.pMMneighs``, so
   the same grouping serves both methods.

.. raw:: latex

   \input{chart.tex}


GeothermalHFMModel
------------------

Unified geothermal model for embedded and projection-based discrete fractures.

Inherits from ``GeothermalModel`` and overrides ``setupOperators`` to:

- Resolve ``fractureMethod`` (``'auto'`` selects ``'pedfm'`` when
  ``G.nnc.pMMneighs`` is present, ``'edfm'`` otherwise) and call the matching
  operator builder:

  - ``'edfm'`` uses ``setupEDFMOperatorsTPFA``: the extended neighbour list,
    flow transmissibilities with ``transmultEDFM``, and NNC flow
    transmissibilities.
  - ``'pedfm'`` uses ``setupPEDFMOpsTPFA``, which applies ``transmultpEDFM``
    and the projected NNCs on top of the EDFM operators.

- Compute rock and fluid heat face transmissibilities with the **same**
  geometric area-reduction multiplier as flow (``transmultEDFM``, and
  ``transmultpEDFM`` for pEDFM), then append the NNC thermal transmissibilities
  (``G.nnc.transHr``, ``G.nnc.transHf``).  The fracture-covered face area
  carries its cross-fracture conduction through the fracture/matrix heat NNCs,
  so a flow **barrier still conducts heat**.  This follows the coupled
  mass/energy finite-volume formulation of HosseiniMehr et al. (2022).

- Ensure the ``AccDiv`` operator exists for the extended connectivity.

It also overrides ``setupStateFunctionGroupings`` to replace the standard
``DynamicTransmissibility`` state functions with EDFM/pEDFM-aware versions when
dynamic (function-handle) conductivities or permeabilities are used.  These
versions apply the same area-reduction multiplier to dynamic heat conduction
and molecular diffusion.

**Properties:**

- ``fractureMethod`` (default ``'auto'``): one of ``'auto'``, ``'edfm'``, or
  ``'pedfm'``.  ``'pedfm'`` raises an error if the grid lacks projection data.
- ``edfmTolerance`` (default ``1e-5``): geometric tolerance for NNC and
  multiplier computations.

**Required grid fields:**

- ``G.nnc.cells``: NNC cell pairs
- ``G.nnc.T``: NNC flow transmissibilities
- ``G.nnc.transHr`` / ``G.nnc.transHf``: NNC rock / fluid heat transmissibilities
- ``G.nnc.CI_dyn``: contact index for dynamic updates
- ``G.nnc.type``: NNC type strings (for example ``'fracmat boundary'`` or
  ``'fracmat interior'``)
- *(pEDFM only)* ``G.nnc.pMMneighs`` and ``G.nnc.normal``: projected
  matrix/matrix neighbours and fracture-plane normals, set by ``pMatFracNNCs3D``

**Usage:**

.. code-block:: matlab

   % Auto-detect from the grid (pEDFM if pMatFracNNCs3D was run, else EDFM)
   model = GeothermalHFMModel(G, rock, fluid);

   % Force a method explicitly
   model = GeothermalHFMModel(G, rock, fluid, 'fractureMethod', 'edfm');
   model = GeothermalHFMModel(G, rock, fluid, 'fractureMethod', 'pedfm');

   % With a compositional (brine) fluid for salinity transport
   model = GeothermalHFMModel(G, rock, fluid, compFluid, ...
                              'fractureMethod', 'pedfm');

.. note::

   ``'auto'`` is convenient but silent: if you intend pEDFM but forget to run
   ``pMatFracNNCs3D``, the model runs EDFM without warning.  Pass
   ``'fractureMethod', 'pedfm'`` to require the projection data and raise a
   clear error when it is missing.


.. rubric:: Note on two-phase water and steam

The class inherits the fluid systems of ``GeothermalModel``.  Single-phase
liquid water and H2O + NaCl brine (a ``CompositionalMixture`` passed as the
optional ``compFluid`` argument) are fully supported with wells and fractures.

Single-component **two-phase water and steam** (boiling) is **not currently
supported** by this module.  It would rely on the base geothermal flash, which
is not yet ready for fractured or well-driven two-phase problems:

*What the base geothermal module has now.*  A physically sound single-component
flash that assigns the phase state from specific enthalpy and obtains the
two-phase saturation from the lever rule, with IAPWS-IF97 water and steam
properties.  It is validated for one-dimensional, convection /
boundary-condition-driven problems (Weis et al., 2014), which by construction
contain no fractures and no wells.

*What it would need for a proper implementation.*  Before two-phase could be
exposed here and extended to fractures, the base flash needs:

- active phase-appearance / disappearance handling (enthalpy-jump limiting and
  flip-flop damping across the saturation line), which is currently inactive;
- robust primary-variable switching between the single-phase and two-phase
  regions, as mature geothermal simulators (e.g. TOUGH2) use;
- reliable property tables near the saturation and critical regions, without
  ad-hoc patches, and a trustworthy thermal-conductivity table;
- stable well source/sink terms under flashing.

With those in place the two-phase flash could be coupled to the EDFM/pEDFM
connectivity exactly as the single-phase energy balance is.  Until then, use
single-phase liquid or brine flow.
