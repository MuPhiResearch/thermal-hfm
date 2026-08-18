function G = computeThermalNNCTransFracMatrix(G, rock, fluid, varargin)
% Compute thermal NNC transmissibilities for fracture-matrix connections.
%
% SYNOPSIS:
%   G = computeThermalNNCTransFracMatrix(G, rock, fluid)
%   G = computeThermalNNCTransFracMatrix(G, rock, fluid, 'pRef', pRef, 'TRef', TRef)
%
% DESCRIPTION:
%   Computes rock and fluid thermal transmissibilities for fracture-matrix
%   NNC connections. Must be called AFTER the stock hfm NNC preprocessing
%   (defineNNCandTrans or fracturematrixNNC3D + computeEffectiveTrans).
%
%   The function uses the same harmonic-averaging approach as the hydraulic
%   transmissibility but substitutes thermal conductivity for permeability:
%     - Rock:  lambda_R * (V_cell - V_pore) / V_cell  [solid fraction]
%     - Fluid: lambda_F * V_pore / V_cell              [pore fraction]
%
%   If G.nnc.CI is available (2D pipeline), it is used directly.
%   Otherwise, CI is back-computed from G.nnc.T using the flow
%   transmissibility relationship.
%
%   If fluid.lambdaF is a function handle @(p,T) (temperature/pressure
%   dependent conductivity), it is evaluated at the reference state
%   (pRef, TRef) to produce the static NNC transmissibilities. This is an
%   acceptable engineering approximation because the NNC geometric factor
%   dominates over lambda_F variation for typical geothermal temperature
%   ranges. Dynamic face transmissibilities are handled separately by
%   GeothermalModel's FluidThermalConductivity state function.
%
%   Similarly, if rock.lambdaR is a function handle @(T), it is evaluated
%   at TRef.
%
% PARAMETERS:
%   G     - Grid structure with G.nnc.cells and either G.nnc.CI or
%           G.nnc.T already computed by the hfm module.
%   rock  - Rock structure with fields:
%             rock.lambdaR - Rock thermal conductivity [W/(m*K)].
%                            Scalar, per-cell array, or @(T) handle.
%             rock.perm    - Permeability (needed if G.nnc.CI absent)
%             rock.poro    - Porosity
%   fluid - Fluid structure with field:
%             fluid.lambdaF - Fluid thermal conductivity [W/(m*K)].
%                             Scalar or @(p,T) function handle.
%
% OPTIONAL PARAMETERS (name-value):
%   pRef  - Reference pressure for evaluating function-handle lambdaF
%           [Pa]. Default: 100*barsa.
%   TRef  - Reference temperature for evaluating function-handle lambdaF/R
%           [K]. Default: 293.15 (= 20 degC).
%
% RETURNS:
%   G - Grid with added fields:
%         G.nnc.transHr - Rock heat NNC transmissibilities [W/K]
%         G.nnc.transHf - Fluid heat NNC transmissibilities [W/K]

%{
Thermal EDFM Module
Couples hfm (EDFM fractures) with geothermal (thermal transport).
%}

% --- Optional reference state for function-handle properties ---
opt = struct('pRef', 100*barsa, 'TRef', 293.15);
opt = merge_options(opt, varargin{:});

% --- Input validation ---
assert(isfield(G, 'nnc'), ...
    'G.nnc not found. Run hfm NNC preprocessing first.');
assert(isfield(G.nnc, 'cells') && ~isempty(G.nnc.cells), ...
    'G.nnc.cells not found or empty.');
assert(isfield(rock, 'lambdaR'), ...
    'rock.lambdaR (rock thermal conductivity) is required.');
assert(isfield(fluid, 'lambdaF'), ...
    'fluid.lambdaF (fluid thermal conductivity) is required.');

nc = G.cells.num;

% --- Pore volume and solid volume ---
if isfield(rock, 'poro')
    pv = poreVolume(G, rock);
else
    pv = G.cells.volumes;
end
vol   = G.cells.volumes;
solid = vol - pv;

% --- Determine which NNCs are frac-matrix (CI-based) ---
% G.nnc.CI only covers frac-matrix NNCs. Frac-frac NNCs (star-delta)
% are appended later by frac_frac_nnc and have no CI entry.
% We compute thermal trans only for the CI-covered NNCs here.
if isfield(G.nnc, 'CI') && ~isempty(G.nnc.CI)
    % 2D pipeline: CI is directly available for frac-matrix NNCs
    nfm = numel(G.nnc.CI);
    CI  = G.nnc.CI;
elseif isfield(G.nnc, 'T') && ~isempty(G.nnc.T)
    % 3D pipeline: back-compute CI from flow transmissibility
    % Use G.nnc.cells as the authoritative NNC count (pMatFracNNCs3D
    % may leave G.nnc.T slightly longer than G.nnc.cells).
    nfm = size(G.nnc.cells, 1);
    assert(isfield(rock, 'perm'), ...
        'rock.perm needed to back-compute CI from G.nnc.T');
    T_nnc = G.nnc.T(1:nfm);
    c1_all = G.nnc.cells(:, 1);
    c2_all = G.nnc.cells(:, 2);
    w1 = pv(c1_all) ./ rock.perm(c1_all);
    w2 = pv(c2_all) ./ rock.perm(c2_all);
    wt = pv(c1_all) + pv(c2_all);
    CI = T_nnc .* (w1 + w2) ./ wt;
else
    error('thermal_edfm:missingData', ...
        'Either G.nnc.CI or G.nnc.T must exist.');
end

% Cell indices for frac-matrix NNCs only
c1 = G.nnc.cells(1:nfm, 1);
c2 = G.nnc.cells(1:nfm, 2);

% --- Rock thermal conductivity (per-cell) ---
% rock.lambdaR may be: scalar, per-cell array, or @(T) function handle.
% Function handles are evaluated at the reference temperature TRef.
lambdaR = rock.lambdaR;
if isa(lambdaR, 'function_handle')
    lambdaR = lambdaR(opt.TRef);   % evaluate @(T) at reference T
end
if isscalar(lambdaR)
    lambdaR = repmat(lambdaR, nc, 1);
end
% Effective rock thermal conductivity: weighted by solid fraction
lambdaR_eff = lambdaR .* solid ./ vol;

% --- Fluid thermal conductivity (per-cell) ---
% fluid.lambdaF may be: scalar or @(p,T) function handle.
% Function handles are evaluated at the reference state (pRef, TRef).
% This is appropriate because the NNC geometric factor dominates;
% dynamic face transmissibilities are handled by FluidThermalConductivity.
lambdaF = fluid.lambdaF;
if isa(lambdaF, 'function_handle')
    lambdaF_val = lambdaF(opt.pRef, opt.TRef);
    if isscalar(lambdaF_val)
        lambdaF = repmat(lambdaF_val, nc, 1);
    else
        lambdaF = lambdaF_val;
    end
elseif isscalar(lambdaF)
    lambdaF = repmat(lambdaF, nc, 1);
end
% Effective fluid thermal conductivity: weighted by pore fraction
lambdaF_eff = lambdaF .* pv ./ vol;

% --- Rock thermal transmissibility ---
% HTr = CI * (solid1 + solid2) / (solid1/lambdaR1 + solid2/lambdaR2)
h1 = solid(c1) ./ lambdaR_eff(c1);
h2 = solid(c2) ./ lambdaR_eff(c2);
ht = solid(c1) + solid(c2);
G.nnc.transHr = CI .* (ht ./ (h1 + h2));

% --- Fluid thermal transmissibility ---
% HTf = CI * (pv1 + pv2) / (pv1/lambdaF1 + pv2/lambdaF2)
h1 = pv(c1) ./ lambdaF_eff(c1);
h2 = pv(c2) ./ lambdaF_eff(c2);
ht = pv(c1) + pv(c2);
G.nnc.transHf = CI .* (ht ./ (h1 + h2));


assert(all(isfinite(G.nnc.transHr) & G.nnc.transHr >= 0), ...
    'Non-finite or negative rock thermal NNC transmissibilities.');
assert(all(isfinite(G.nnc.transHf) & G.nnc.transHf >= 0), ...
    'Non-finite or negative fluid thermal NNC transmissibilities.');

% --- Store geometric CI for dynamic NNC recomputation ---
% EDFMDynamicHeatTransmissibility uses this to recompute NNC heat
% transmissibilities each Newton iteration when fluid.lambdaF is a
% @(p,T) handle (dynamicHeatTransFluid = true).
%
% The same geometric CI applies to both rock and fluid NNCs:
%   transHf = CI_dyn * (pv1+pv2) / (pv1/lambdaF_eff1 + pv2/lambdaF_eff2)
%   transHr = CI_dyn * (sv1+sv2)  / (sv1/lambdaR_eff1 + sv2/lambdaR_eff2)
% where sv = vol - pv (solid volume fraction).
G.nnc.CI_dyn = CI;

end
