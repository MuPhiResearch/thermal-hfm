classdef GeothermalHFMModel < GeothermalModel
%Unified geothermal model for embedded/projected discrete fractures (EDFM/pEDFM).
%
% SYNOPSIS:
%   model = GeothermalHFMModel(G, rock, fluid)                      % auto
%   model = GeothermalHFMModel(G, rock, fluid, compFluid)          % auto + EOS
%   model = GeothermalHFMModel(..., 'fractureMethod', 'edfm')      % force EDFM
%   model = GeothermalHFMModel(..., 'fractureMethod', 'pedfm')     % force pEDFM
%
% DESCRIPTION:
%   Single model class for geothermal flow and heat transfer through
%   fractured reservoirs discretized with the (projection-based) embedded
%   discrete fracture model.  EDFM and pEDFM are the SAME model differing only in
%   how the flow operators are built, and that choice is driven by the GRID
%   (whether the pEDFM projection preprocessing has been run), not by the
%   class.
%
%   fractureMethod selects the operator builder:
%     'auto'  (default) - pEDFM iff the grid carries projected connections
%                         (isfield(G.nnc,'pMMneighs')), EDFM otherwise.  This
%                         is the SAME gate EDFMDynamicFlowTransmissibility
%                         uses at runtime, so the static and dynamic
%                         transmissibility paths agree by construction.
%     'edfm'            - force standard EDFM (setupEDFMOperatorsTPFA).
%     'pedfm'           - force pEDFM (setupPEDFMOpsTPFA); raises an error if the
%                         grid lacks the projection data.
%
%   Heat conduction and molecular diffusion across the matrix-matrix
%   interfaces shadowed by a fracture receive the SAME geometric
%   area-reduction multiplier as fluid flow (transmultEDFM, and
%   transmultpEDFM for pEDFM), consistent with the coupled mass/energy
%   finite-volume formulation of HosseiniMehr et al. (2022). The
%   projected/boundary fracture-matrix heat NNCs carry the cross-fracture
%   conduction (weighted by thermal conductivity), so a flow barrier still
%   conducts heat.
%
%   FLUID SYSTEMS. Being a thin specialisation of GeothermalModel, the class
%   inherits its fluid systems and works through the fracture grid unchanged:
%     - single-phase liquid water (temperature formulation);
%     - H2O + NaCl brine via a CompositionalMixture passed as the optional
%       compFluid argument (salinity transport with the Spivey EOS).
%   Single-component two-phase water and steam (boiling) is NOT supported: it
%   would depend on the base geothermal flash, which is not yet ready for
%   fractured or well-driven two-phase problems.
%
% REQUIRED FIELDS (on G):
%   G.nnc.cells, G.nnc.T, G.nnc.transHr, G.nnc.transHf, G.nnc.type, G.nnc.area
%   pEDFM additionally requires: G.nnc.pMMneighs, G.nnc.normal
%
% SEE ALSO:
%   GeothermalModel, setupEDFMOperatorsTPFA, setupPEDFMOpsTPFA,
%   transmultEDFM, transmultpEDFM, pMatFracNNCs3D,
%   computeThermalNNCTransFracMatrix, EDFMDynamicFlowTransmissibility,
%   EDFMDynamicHeatTransmissibility, EDFMDynamicMolecularTransmissibility

    properties
        % Tolerance for EDFM/pEDFM transmissibility-multiplier computation
        edfmTolerance = 1e-5;
        % Fracture discretization method: 'auto' | 'edfm' | 'pedfm'
        fractureMethod = 'auto';
    end

    methods
        %-----------------------------------------------------------------%
        function model = GeothermalHFMModel(G, rock, fluid, varargin)
        % Class constructor.
        %
        % 'fractureMethod' is peeled off here BEFORE the parent constructor,
        % because GeothermalModel uses STRICT merge_options which errors on
        % unknown name/value pairs. The parent runs setupOperators() once
        % under the default fractureMethod='auto'; an explicit method then
        % triggers a single rebuild (MATLAB forbids assigning model.* before
        % the super-call, so this is the correct pattern).

            method = 'auto';
            keep   = {};
            i = 1;
            while i <= numel(varargin)
                if (ischar(varargin{i}) || isstring(varargin{i})) ...
                        && strcmpi(varargin{i}, 'fractureMethod') ...
                        && i < numel(varargin)
                    method = validatestring(lower(char(varargin{i+1})), ...
                        {'auto', 'edfm', 'pedfm'});
                    i = i + 2;
                else
                    keep{end+1} = varargin{i}; %#ok<AGROW>
                    i = i + 1;
                end
            end

            % Parent constructor (runs setupOperators() once under 'auto')
            model = model@GeothermalModel(G, rock, fluid, keep{:});

            % Honor an explicit method with a single operator rebuild
            if ~strcmpi(method, 'auto')
                model.fractureMethod = method;
                model = model.setupOperators();
            end
        end

        %-----------------------------------------------------------------%
        function method = resolveFractureMethod(model, G)
        % Resolve 'auto' to 'edfm'/'pedfm' from the grid contents.
            method = lower(model.fractureMethod);
            if strcmp(method, 'auto')
                if isfield(G, 'nnc') && isfield(G.nnc, 'pMMneighs')
                    method = 'pedfm';
                else
                    method = 'edfm';
                end
            end
        end

        %-----------------------------------------------------------------%
        function model = setupOperators(model, G, rock, varargin)
        % Set up flow + thermal operators with EDFM/pEDFM support.

            opt = struct('transHr', [], ...
                         'transHf', [], ...
                         'vol'    , []);
            [opt, ~] = merge_options(opt, varargin{:});
            if nargin < 3, rock = model.rock; end
            if nargin < 2, G = model.G;       end

            % --- Validate NNC fields ---
            assert(isfield(G, 'nnc'), ...
                'G.nnc required for GeothermalHFMModel');
            assert(isfield(G.nnc, 'transHr') && isfield(G.nnc, 'transHf'), ...
                'G.nnc must have transHr and transHf fields');

            method = model.resolveFractureMethod(G);

            % --- Dummy perm for dynamic flow transmissibility setup ---
            drock = rock;
            if model.dynamicFlowTrans()
                drock.perm = rock.perm(1*barsa, 273.15 + 20*Kelvin);
            end

            % --- Flow operators ---
            switch method
                case 'pedfm'
                    assert(isfield(G.nnc, 'pMMneighs') && isfield(G.nnc, 'normal'), ...
                        ['fractureMethod=''pedfm'' requires G.nnc.pMMneighs and ', ...
                         'G.nnc.normal. Run pMatFracNNCs3D in the shale-3D pipeline.']);
                    % setupPEDFMOpsTPFA re-reads G.rock (not its rock arg),
                    % so make G.rock consistent with the dummy perm.
                    G.rock = drock;
                    model.operators = setupPEDFMOpsTPFA(G, drock, model.edfmTolerance);
                case 'edfm'
                    model.operators = setupEDFMOperatorsTPFA(G, drock, model.edfmTolerance);
            end
            model.rock = rock;

            % Ensure AccDiv exists (setupPEDFMOpsTPFA provides it; the hfm
            % EDFM setup does not). Both equal acc + operators.C'*flux.
            if ~isfield(model.operators, 'AccDiv')
                C = model.operators.C;
                model.operators.AccDiv = @(acc, flux) acc + C'*flux;
            end

            % --- Cell volumes ---
            pv  = model.operators.pv;
            vol = opt.vol;
            if isempty(vol), vol = G.cells.volumes; end
            model.operators.vol = vol;

            % --- Geometric area-reduction multiplier for M-M conduction ---
            % Heat conduction across a matrix-matrix face shadowed by a
            % fracture is reduced by the SAME geometric factor as fluid flow
            % (the fracture occupies part of the face; that area's
            % cross-fracture conduction is carried by the fracture-matrix
            % heat NNCs in G.nnc.transHr/transHf). This mirrors the coupled
            % mass/energy finite-volume formulation of HosseiniMehr et al.
            % (2022): the projected/reduced connection areas modify the shared
            % mass+energy connection list.
            % tmult = 1 on every fracture-free face, so interior fractures
            % and regular faces are unaffected.
            tmult = ones(G.faces.num, 1);
            if ~model.dynamicHeatTransRock() || ~model.dynamicHeatTransFluid()
                tmult = transmultEDFM(G, model.edfmTolerance);
                if strcmp(method, 'pedfm')
                    tmult = tmult .* transmultpEDFM(G, model.edfmTolerance);
                end
            end

            % --- Rock heat transmissibility ---
            if ~model.dynamicHeatTransRock()
                Thr = opt.transHr;
                if isempty(Thr)
                    % Reference temperature for @(T) lambdaR evaluation
                    TRef = 273.15 + 20;
                    lambdaR = rock.lambdaR;
                    if isa(lambdaR, 'function_handle')
                        lambdaR = lambdaR(TRef);
                    end
                    if isscalar(lambdaR)
                        lambdaR = repmat(lambdaR, G.cells.num, 1);
                    end
                    lambdaR_eff = lambdaR .* (vol - pv) ./ vol;
                    r   = struct('perm', lambdaR_eff);
                    % Reduce M-M conduction by the fracture area fraction
                    Thr = getFaceTransmissibility(G, r) .* tmult;
                    % Append NNC thermal transmissibilities (carry the
                    % cross-fracture conduction, incl. barriers)
                    Thr = [Thr; G.nnc.transHr];
                end
                model.operators.Thr     = Thr(model.operators.internalConn);
                model.operators.Thr_all = Thr;
            end

            % --- Fluid heat transmissibility ---
            if ~model.dynamicHeatTransFluid()
                Thf = opt.transHf;
                if isempty(Thf)
                    % Reference state for @(p,T) lambdaF evaluation
                    pRef = 1*barsa; TRef = 273.15 + 20;
                    lambdaF = model.fluid.lambdaF;
                    if isa(lambdaF, 'function_handle')
                        lambdaF_val = lambdaF(pRef, TRef);
                        if isscalar(lambdaF_val)
                            lambdaF = repmat(lambdaF_val, G.cells.num, 1);
                        else
                            lambdaF = lambdaF_val;
                        end
                    elseif isscalar(lambdaF)
                        lambdaF = repmat(lambdaF, G.cells.num, 1);
                    end
                    lambdaF_eff = lambdaF .* pv ./ vol;
                    r   = struct('perm', lambdaF_eff);
                    % Reduce M-M conduction by the fracture area fraction
                    Thf = getFaceTransmissibility(G, r) .* tmult;
                    % Append NNC thermal transmissibilities (carry the
                    % cross-fracture conduction, incl. barriers)
                    Thf = [Thf; G.nnc.transHf];
                end
                model.operators.Thf     = Thf(model.operators.internalConn);
                model.operators.Thf_all = Thf;
            end

            % Compute hash
            model.operators.hashRock = obj2hash(rock);

        end
        %-----------------------------------------------------------------%
        function model = setupStateFunctionGroupings(model, varargin)
        % Override to replace DynamicTransmissibility with EDFM/pEDFM-aware
        % versions for heat, molecular diffusion, and flow transmissibility.
        %
        % DynamicTransmissibility builds its harmonic-average operator from
        % G.cells.faces (regular faces only) and produces G.faces.num
        % values, but operators.internalConn has length G.faces.num + nNNC.
        % The EDFM-aware replacements handle NNCs correctly. They self-detect
        % pEDFM (G.nnc.pMMneighs) internally, so this grouping is identical
        % for EDFM and pEDFM.

            % Pre-build the geothermal flux discretization so the parent grouping
            % does not emit "Assuming default flux discretization": GeothermalModel
            % reads model.FlowDiscretization at entry and would otherwise build
            % this very same default itself. Physics is unchanged.
            if isempty(model.FlowDiscretization) || ...
                    ~isa(model.FlowDiscretization, 'GeothermalFlowDiscretization')
                model.FlowDiscretization = GeothermalFlowDiscretization(model);
            end

            % Call parent to set up all standard groupings first
            model = setupStateFunctionGroupings@GeothermalModel(model, varargin{:});

            disc = model.FlowDiscretization;

            % --- Dynamic heat transmissibilities ---
            if model.dynamicHeatTransFluid()
                disc = disc.setStateFunction('FluidHeatTransmissibility', ...
                    EDFMDynamicHeatTransmissibility(model, 'hf'));
            end
            if model.dynamicHeatTransRock()
                disc = disc.setStateFunction('RockHeatTransmissibility', ...
                    EDFMDynamicHeatTransmissibility(model, 'hr'));
            end

            % --- Dynamic molecular transmissibility ---
            hasDiffusion = any(cellfun(@(c) c.molecularDiffusivity, ...
                model.Components) > 0);
            if hasDiffusion
                disc = disc.setStateFunction('MolecularTransmissibility', ...
                    EDFMDynamicMolecularTransmissibility(model));
            end

            % --- Dynamic flow transmissibility ---
            if model.dynamicFlowTrans()
                disc = disc.setStateFunction('Transmissibility', ...
                    EDFMDynamicFlowTransmissibility(model));
            end

            model.FlowDiscretization = disc;
        end

        %-----------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
        % Override to allocate flux arrays with correct EDFM/pEDFM sizing
        % (G.faces.num + nNNC).

            [state, report] = updateAfterConvergence@ReservoirModel(model, state0, state, dt, drivingForces);
            if model.extraStateOutput
                rho = model.getProps(state, 'Density');
                state.rho = horzcat(rho{:});
            end
            if model.outputFluxes
                nNNC = size(model.G.nnc.cells, 1);
                nf   = model.G.faces.num + nNNC;

                state_flow = model.FlowDiscretization.buildFlowState(model, state, state0, dt);
                f = model.getProp(state_flow, 'PhaseFlux');
                nph = numel(f);
                state.flux = zeros(nf, nph);
                state.flux(model.operators.internalConn, :) = [f{:}];

                [heatFluxAdv, heatFluxCond] = model.getProps(state, ...
                    'AdvectiveHeatFlux', 'ConductiveHeatFlux');
                state.heatFluxAdv  = zeros(nf, nph);
                state.heatFluxCond = zeros(nf, 1);
                state.heatFluxAdv(model.operators.internalConn,:) = horzcat(heatFluxAdv{:});
                state.heatFluxCond(model.operators.internalConn)  = heatFluxCond;

                if ~isempty(drivingForces.bc)
                    [p, s, mob, r, b] = model.getProps(state, ...
                        'PhasePressures', 's', 'Mobility', 'Density', 'ShrinkageFactors');
                    sat = expandMatrixToCell(s);
                    rho = expandMatrixToCell(r);
                    [~, ~, ~, fRes] = getBoundaryConditionFluxesAD(model, p, sat, mob, rho, b, drivingForces.bc);
                    idx = model.getActivePhases();
                    fWOG = cell(3, 1);
                    fWOG(idx) = fRes;
                    state = model.storeBoundaryFluxes(state, fWOG{1}, fWOG{2}, fWOG{3}, drivingForces);

                    if ~model.thermal, return; end
                    drivingForces.bc = getForceProperties(drivingForces.bc, model, state_flow);
                    src = struct();
                    phaseMass = cellfun(@(rho, q) rho.*q, ...
                        drivingForces.bc.propsRes.rho, fRes', 'UniformOutput', false);
                    src.bc.phaseMass = phaseMass;
                    src.bc.sourceCells = sum(model.G.faces.neighbors(drivingForces.bc.face,:), 2);
                    src.src.sourceCells = [];
                    src = getHeatFluxFromSources(model, src, drivingForces);

                    faces = drivingForces.bc.face;
                    sgn = 1 - 2*(model.G.faces.neighbors(faces, 2) == 0);
                    for ph = 1:model.getNumberOfPhases()
                        state.heatFluxAdv(faces, ph) = src.bc.advHeatFlux{ph}.*sgn;
                    end
                    state.heatFluxCond(faces) = src.bc.condHeatFlux.*sgn;
                end
            end
        end
        %-----------------------------------------------------------------%
    end

end
