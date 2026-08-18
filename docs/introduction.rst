Introduction
============

Overview
--------

The **thermal-hfm** module simulates geothermal flow and heat transport
(pressure, temperature, and optionally salinity) through fractured reservoirs
in **MRST 2026a**.  It supports single-phase liquid water and H2O + NaCl brine
on the fracture grid.  It couples two existing MRST modules:

- **geothermal**: the energy-balance equation, thermal rock and fluid
  properties, and the ``GeothermalModel`` base class.
- **hfm / shale**: embedded-fracture grid generation, fracture/matrix
  intersection geometry, and Non-Neighbouring Connection (NNC) computation,
  including the projection-based (pEDFM) connections.

It bridges the two by computing the thermal NNC transmissibilities, providing
a single EDFM/pEDFM-aware model class (``GeothermalHFMModel``), and supplying
state functions that handle the extended connectivity (regular faces plus
NNCs) that fractured grids require.


Embedded discrete fracture model (EDFM)
---------------------------------------

EDFM represents fractures as additional cells appended to the matrix grid.
Each fracture segment intersecting a matrix cell becomes a fracture cell with
its own volume (area × aperture), porosity, and permeability.  Fracture cells
do not share regular grid faces with the matrix.  Instead they communicate
through **Non-Neighbouring Connections (NNCs)**, extra transmissibility entries
that link each fracture cell to the matrix cell that hosts it, and to other
fracture cells at intersections.

The governing geometric quantity is the **Contact Index (CI)**:

.. math::

   \CI = \frac{A_{\mathrm{intersection}}}{d_{\mathrm{avg}}}

where :math:`A_{\mathrm{intersection}}` is the fracture/cell intersection area
and :math:`d_{\mathrm{avg}}` is the volume-weighted average normal distance
from the matrix cell to the fracture plane.  An NNC flow transmissibility is
:math:`T = \CI \, k`, and the thermal NNC transmissibility uses the effective
rock or fluid conductivity in place of permeability.

EDFM (Lee et al., 2001; Moinfar et al., 2014) is accurate for conductive
fractures, but it has a well-known limitation.  A fracture lying in the
interior of a matrix cell cannot block flow between that cell and its
neighbours, because the matrix/matrix faces stay fully open.  Low-permeability
features such as barriers and sealed faults are therefore effectively
invisible to the flow field.


Projection-based EDFM (pEDFM)
-----------------------------

The projection-based extension (Ţene et al., 2017) removes this limitation in
three steps:

1. identify the matrix faces that an interior fracture projects onto (the
   faces between the host cell and its neighbours);
2. **reduce** the matrix/matrix transmissibility of each such face by the
   projected fracture area; and
3. add a **projected fracture/matrix** NNC that carries the displaced flux.

The reduction is purely geometric.  Its multiplier is
:math:`(A_{\mathrm{face}} - A_{\mathrm{proj}})/A_{\mathrm{face}}`, with no
permeability term, so it is applied identically whether the fracture is a
conduit or a barrier.  It only *matters* for barriers.  For a conductive
fracture the displaced flux re-enters through the high-permeability fracture
and the two effects nearly cancel, so pEDFM reproduces EDFM.  For a barrier
there is no bypass, so the reduction correctly seals the matrix.

.. note::

   The projection adjusts the same connection geometry that both the mass and
   the energy equations use.  The area-reduction multiplier (``transmultEDFM``
   on boundary fractures, ``transmultpEDFM`` on projected interior fractures)
   is therefore applied to the **conductive heat** and **molecular-diffusion**
   transmissibilities exactly as it is to fluid flow, following the coupled
   mass/energy finite-volume formulation of HosseiniMehr et al. (2022).  The
   fracture-covered face area carries its cross-fracture conduction through the
   fracture/matrix heat NNCs, which avoids double-counting the conduction path.
   Because that heat NNC is weighted by thermal conductivity rather than
   permeability, a flow **barrier still conducts heat**, as it should
   physically.


Overall workflow
----------------

A typical simulation involves the following steps.

**1. Grid generation.**  Create a Cartesian matrix grid, define fracture
planes (corner points, aperture, porosity, permeability), and embed them into
the global grid.

- 2D: ``processFracture2D`` |rarr| ``CIcalculator2D`` |rarr| ``gridFracture2D``
- 3D: ``EDFMshalegrid`` or ``preProcessingFractures`` |rarr| ``assembleGlobalGrid``

**2. NNC computation (flow).**  Compute fracture/matrix and fracture/fracture
flow NNCs.

- 2D: ``defineNNCandTrans``
- 3D: ``fracturematrixShaleNNC3D`` |rarr| ``fracturefractureShaleNNCs3D``
- pEDFM: additionally ``pMatFracNNCs3D`` (this is what makes a grid "pEDFM")

**3. Thermal properties.**  Attach thermal conductivities and heat capacities
to rock and fluid with MRST's ``addThermalRockProps`` and
``addThermalFluidProps``.

**4. Thermal NNC computation.**  Run ``computeThermalNNCTransFracMatrix`` (and
``computeThermalNNCTransFracFrac`` for the 2D pipeline) to obtain the rock and
fluid heat NNC transmissibilities and store the CI for dynamic updates.

**5. Model construction.**  Create a ``GeothermalHFMModel``.  Its
``fractureMethod`` property (``'auto'``, ``'edfm'``, or ``'pedfm'``) selects
the operator builder, ``setupEDFMOperatorsTPFA`` or ``setupPEDFMOpsTPFA``,
which assembles the extended neighbour list, applies the flow and thermal
multipliers, and appends the NNC transmissibilities.  ``'auto'`` picks pEDFM
when the grid carries ``G.nnc.pMMneighs`` (that is, ``pMatFracNNCs3D`` was
run), and EDFM otherwise.

**6. Simulation.**  Define wells, schedule, and initial state, then run
``simulateScheduleAD``.

.. |rarr| unicode:: U+2192 .. right arrow


Design notes
------------

1. **One model class.**  EDFM and pEDFM are two configurations of a single
   ``GeothermalHFMModel``, not separate classes.  The choice is a property of
   the *grid*, namely whether the projection preprocessing (``pMatFracNNCs3D``)
   has run and left ``G.nnc.pMMneighs``.  It is selected through the
   ``fractureMethod`` property: ``'auto'`` detects it, while ``'edfm'`` and
   ``'pedfm'`` force a choice (``'pedfm'`` raises an error when the projection data
   is absent).  The dynamic state functions detect pEDFM the same way, so the
   static and dynamic transmissibility paths agree by construction.

2. **Conduction and diffusion follow the flow connectivity.**  The geometric
   area-reduction multiplier is applied to heat and molecular diffusion
   exactly as to flow (HosseiniMehr et al., 2022).  It is the identity on
   fracture-free and interior-fracture faces, so EDFM results are unchanged.
   For pEDFM it removes the matrix/matrix conduction double-count, while the
   fracture/matrix heat NNCs preserve cross-fracture conduction so barriers
   still conduct heat.

3. **CI back-computation (3D).**  The 3D pipeline does not store CI directly.
   It is back-computed from ``G.nnc.T`` via the inverse flow-transmissibility
   relation and stored as ``G.nnc.CI_dyn`` for reuse by the dynamic state
   functions.

4. **Consistent property interface.**  Thermal conductivities are supplied
   through ``rock.lambdaR`` and ``fluid.lambdaF`` (the ``addThermalRockProps``
   and ``addThermalFluidProps`` convention).  Scalars, per-cell vectors, and
   ``@(T)`` or ``@(p,T)`` function handles are all accepted; the
   function-handle forms activate the dynamic state functions.


References
----------

- Lee, S.H., Lough, M.F., Jensen, C.L. (2001).  Hierarchical modeling of flow
  in naturally fractured formations with multiple length scales.
  *Water Resources Research*, 37(3), 443-455.
- Moinfar, A., Varavei, A., Sepehrnoori, K., Johns, R.T. (2014).  Development
  of an efficient embedded discrete fracture model for 3D compositional
  reservoir simulation in fractured reservoirs.  *SPE Journal*, 19(2), 289-303.
- Ţene, M., Bosma, S.B.M., Al Kobaisi, M.S., Hajibeygi, H. (2017).
  Projection-based Embedded Discrete Fracture Model (pEDFM).
  *Advances in Water Resources*, 105, 205-216.
- HosseiniMehr, M., Piguave Tomala, J., Vuik, C., Al Kobaisi, M.,
  Hajibeygi, H. (2022).  Projection-based embedded discrete fracture model
  (pEDFM) for flow and heat transfer in real-field geological formations with
  corner-point grids.  *Advances in Water Resources*, 159, 104091.
