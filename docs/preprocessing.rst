Preprocessing Functions
=======================

These functions compute thermal NNC transmissibilities.  They must be called
**after** the standard ``hfm`` NNC preprocessing (which sets up the flow NNCs)
and **before** constructing the model.


computeThermalNNCTransFracMatrix
--------------------------------

.. code-block:: matlab

   G = computeThermalNNCTransFracMatrix(G, rock, fluid)
   G = computeThermalNNCTransFracMatrix(G, rock, fluid, 'TRef', 393.15)

Compute rock and fluid heat NNC transmissibilities for all fracture/matrix
connections.

It uses a harmonic-average formula analogous to flow, but with thermal
conductivity instead of permeability:

.. math::

   T^{\mathrm{Hr}}_{\mathrm{nnc}} = \CI \cdot
     \frac{s_1 + s_2}
          {s_1 / \tilde{\lambda}_{\mathrm{R},1}
         + s_2 / \tilde{\lambda}_{\mathrm{R},2}}

where :math:`s_i` is the solid volume (cell volume minus pore volume) and
:math:`\tilde{\lambda}_{\mathrm{R},i} = \lambda_{\mathrm{R}} \cdot s_i / V_i`
is the effective rock conductivity.  An analogous formula applies for the
fluid heat transmissibility using pore volume and :math:`\lambda_{\mathrm{F}}`.

**Outputs stored on G:**

- ``G.nnc.transHr``: rock heat NNC transmissibilities
- ``G.nnc.transHf``: fluid heat NNC transmissibilities
- ``G.nnc.CI_dyn``: contact index for the dynamic state functions

**Optional parameters** (name-value pairs):

- ``'pRef'``: reference pressure for ``@(p,T)`` conductivities
  (default ``100*barsa``)
- ``'TRef'``: reference temperature (default ``293.15`` K)


computeThermalNNCTransFracFrac
------------------------------

.. code-block:: matlab

   G = computeThermalNNCTransFracFrac(G, F, fracture, rock, fluid)

Compute rock and fluid heat NNC transmissibilities for fracture/fracture
intersection connections using a star/delta transformation.

This function is only needed for the **2D pipeline**, where fracture/fracture
intersections are handled by ``defineNNCandTrans``.  In the 3D pipeline the
fracture/fracture NNCs are already included in the global NNC list and handled
by ``computeThermalNNCTransFracMatrix``.

**Parameters:**

- ``G``: grid with existing ``G.nnc`` from ``defineNNCandTrans``
- ``F``: fracture grid structure from ``gridFracture2D``
- ``fracture``: fracture definition structure
- ``rock``: with ``rock.lambdaR``
- ``fluid``: with ``fluid.lambdaF``
