Dynamic State Functions
=======================

When thermal conductivity or permeability is provided as a function handle
(e.g., ``@(p,T)``), the transmissibilities must be recomputed at every Newton
iteration.  The standard MRST ``DynamicTransmissibility`` state function
produces vectors of length ``G.faces.num``, but EDFM operators expect
``G.faces.num + nNNC``.  The classes below fix this size mismatch.


EDFMDynamicHeatTransmissibility
-------------------------------

Dynamic heat transmissibility for EDFM/pEDFM grids.

Activated automatically by ``GeothermalHFMModel.setupStateFunctionGroupings``
when ``rock.lambdaR`` is a ``@(T)`` function handle or ``fluid.lambdaF``
is a ``@(p,T)`` function handle.

**Constructor:**

.. code-block:: matlab

   prop = EDFMDynamicHeatTransmissibility(model, 'hf');  % fluid
   prop = EDFMDynamicHeatTransmissibility(model, 'hr');  % rock

**Evaluation** (called each Newton iteration):

1. **Face part** (size ``G.faces.num``): per-cell effective conductivity
   is converted to face transmissibilities via harmonic averaging of
   two-point half-transmissibilities.
2. **NNC part** (size ``nNNC``): uses the stored ``G.nnc.CI_dyn`` with a
   pore-volume-weighted (fluid) or solid-volume-weighted (rock) harmonic
   average.
3. The two parts are concatenated to produce a vector of length
   ``G.faces.num + nNNC``.

The same geometric area-reduction multiplier as flow (``transmultEDFM``, and
``transmultpEDFM`` for pEDFM) is applied to the face part; the
fracture-covered area carries its cross-fracture conduction through the NNC
part, so a flow barrier still conducts heat.


EDFMDynamicFlowTransmissibility
-------------------------------

Dynamic flow transmissibility for EDFM/pEDFM grids.

Activated when ``rock.perm`` is a ``@(p,T)`` function handle.

Applies ``transmultEDFM`` (and ``transmultpEDFM`` for pEDFM models) to
the face part, matching the static flow operator setup.

**Constructor:**

.. code-block:: matlab

   prop = EDFMDynamicFlowTransmissibility(model);


EDFMDynamicMolecularTransmissibility
-------------------------------------

Dynamic molecular (diffusive) transmissibility for EDFM/pEDFM grids.

Produces a cell array of transmissibility vectors (one per component),
each of length ``G.faces.num + nNNC``.

The same geometric area-reduction multiplier as flow is applied to the face
part, by analogy with heat conduction (molecular diffusion is a Laplacian
flux through the pore space).

**Constructor:**

.. code-block:: matlab

   prop = EDFMDynamicMolecularTransmissibility(model);
