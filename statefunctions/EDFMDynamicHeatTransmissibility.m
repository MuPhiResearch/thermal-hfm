classdef EDFMDynamicHeatTransmissibility < StateFunction
%Dynamic fluid or rock heat transmissibility for EDFM/pEDFM grids.
%
% SYNOPSIS:
%   prop = EDFMDynamicHeatTransmissibility(model, 'hf')   % fluid heat
%   prop = EDFMDynamicHeatTransmissibility(model, 'hr')   % rock heat
%
% DESCRIPTION:
%   Standard DynamicTransmissibility fails in EDFM/pEDFM because it builds
%   its harmonic-average operator over G.cells.faces (regular faces only),
%   then indexes the result with model.operators.internalConn which in EDFM
%   is extended to length G.faces.num + nNNC.  The result is a size mismatch.
%
%   This class produces the correct (nf_internal + nNNC) heat transmissibility
%   vector needed by ConductiveHeatFlux by computing two contributions
%   separately and concatenating:
%
%     1. FACE part (nf_internal): per-cell effective conductivity ->
%        two-point half-transmissibility -> harmonic average over faces ->
%        geometric area-reduction multiplier (transmultEDFM, and
%        transmultpEDFM for pEDFM, the SAME factor applied to fluid flow) ->
%        filter to internal faces.  The fracture-covered area carries its
%        cross-fracture conduction through the F-M heat NNCs (part 2).
%
%     2. NNC part (nNNC): harmonic average of per-cell effective conductivity
%        weighted by pore volume (fluid) or solid volume (rock), using the
%        geometric contact indicator G.nnc.CI_dyn stored by
%        computeThermalNNCTransFracMatrix.
%
%   The NNC terms carry full AD derivatives with respect to pressure and
%   temperature, so the Jacobian correctly accounts for conductivity changes
%   across Newton iterations.
%
% PREREQUISITES:
%   G.nnc.CI_dyn must exist.  Call computeThermalNNCTransFracMatrix before
%   constructing the model; it stores CI_dyn automatically.
%
% PARAMETERS:
%   model   - GeothermalHFMModel instance.
%   postfix - 'hf' (fluid heat) or 'hr' (rock heat).
%
% SEE ALSO:
%   GeothermalHFMModel,
%   computeThermalNNCTransFracMatrix, DynamicTransmissibility

    properties
        postfix           % 'hf' or 'hr'
        conductivity_name % 'FluidThermalConductivity' or 'RockThermalConductivity'
        use_pore_volume   % true = fluid (uses pv), false = rock (uses solid)
        % Pre-built geometric operators (constructed once from G.cells.faces)
        twoPointOperator
        harmonicAvgOperator
        % Geometric area-reduction multiplier for matrix-matrix conduction
        % (transmultEDFM, and transmultpEDFM for pEDFM) - the SAME factor
        % applied to fluid flow.
        transmult
    end

    methods
        %-----------------------------------------------------------------%
        function prop = EDFMDynamicHeatTransmissibility(model, postfix)
            prop@StateFunction(model);
            assert(ismember(postfix, {'hf', 'hr'}), ...
                'postfix must be ''hf'' (fluid) or ''hr'' (rock)');
            assert(isfield(model.G, 'nnc') && isfield(model.G.nnc, 'CI_dyn'), ...
                ['G.nnc.CI_dyn not found. ', ...
                 'Run computeThermalNNCTransFracMatrix before creating the model.']);

            prop.postfix = postfix;
            if strcmp(postfix, 'hf')
                prop.conductivity_name = 'FluidThermalConductivity';
                prop.use_pore_volume   = true;
            else
                prop.conductivity_name = 'RockThermalConductivity';
                prop.use_pore_volume   = false;
            end

            % Register dependency on the per-cell effective conductivity.
            % Use the local-dependency form (no grouping argument) so that
            % getEvaluatedDependencies can read it from state.FluxDisc,
            % matching the pattern used by DynamicTransmissibility.
            prop = prop.dependsOn(prop.conductivity_name);

            % Pre-build geometry-only operators from G.cells.faces
            prop.twoPointOperator    = buildTwoPointOp(model.G);
            prop.harmonicAvgOperator = buildHarmonicAvgOp(model.G);

            % Geometric area-reduction multiplier for matrix-matrix
            % conduction - the SAME factor applied to fluid flow. The
            % fracture-covered face area carries its cross-fracture
            % conduction through the F-M heat NNCs instead (consistent with
            % the coupled mass/energy FV formulation, HosseiniMehr 2022).
            tol = 1e-5;
            if isprop(model, 'edfmTolerance'), tol = model.edfmTolerance; end
            prop.transmult = transmultEDFM(model.G, tol);
            if isfield(model.G.nnc, 'pMMneighs')
                prop.transmult = prop.transmult .* transmultpEDFM(model.G, tol);
            end

            prop.label = ['T_{', postfix, '}'];
        end

        %-----------------------------------------------------------------%
        function Th = evaluateOnDomain(prop, model, state, allFaces)
        % Compute dynamic heat transmissibility.
        %
        % allFaces = false (default):
        %   Returns a vector of length nf_internal + nNNC suitable for use
        %   with model.operators.Grad and model.operators.AccDiv in the
        %   energy equation.
        %
        % allFaces = true:
        %   Returns a vector of length G.faces.num (all faces, no NNCs),
        %   indexed by absolute face number.  Required by getForceProperties
        %   which calls evaluateOnDomain(model, state, true) and then indexes
        %   the result by BC face indices.

            if nargin < 4, allFaces = false; end

            G  = model.G;
            pv = model.operators.pv;  % pore volume per cell

            % Per-cell effective conductivity from the state function:
            %   fluid: lambdaF(p,T) * pv/vol   [W/(m·K)]
            %   rock:  lambdaR(T)   * sv/vol   [W/(m·K)]
            lambda_eff = prop.getEvaluatedDependencies(state, ...
                prop.conductivity_name);

            % ---- Regular face contributions (size G.faces.num) ---------
            % Two-point half-transmissibilities at each (cell, face) entry
            T_half  = prop.twoPointOperator(lambda_eff);
            % Harmonic average: one value per face (size G.faces.num), with
            % the geometric area-reduction multiplier (same as flow).
            T_faces = prop.harmonicAvgOperator(T_half) .* prop.transmult;

            if allFaces
                % Return all G.faces.num values so getForceProperties can
                % index by absolute BC face number (same contract as
                % DynamicTransmissibility when called with allFaces=true).
                Th = T_faces;
                return;
            end

            % ---- Filter to internal faces (nf_internal) ----------------
            % model.operators.internalConn is logical [nf_all + nNNC].
            % The first nf_all entries are the standard internal-face mask.
            internalFacesMask = model.operators.internalConn(1:G.faces.num);
            Th_faces = T_faces(internalFacesMask);   % nf_internal values

            % ---- NNC contributions (nNNC) -------------------------------
            c1 = G.nnc.cells(:, 1);
            c2 = G.nnc.cells(:, 2);
            CI = G.nnc.CI_dyn;          % geometric contact indicator

            if prop.use_pore_volume
                % Fluid heat: harmonic average weighted by pore volume
                %   T_hf = CI * (pv1+pv2) / (pv1/lambdaF_eff1 + pv2/lambdaF_eff2)
                w1 = pv(c1) ./ lambda_eff(c1);
                w2 = pv(c2) ./ lambda_eff(c2);
                wt = pv(c1) + pv(c2);
            else
                % Rock heat: harmonic average weighted by solid volume
                %   T_hr = CI * (sv1+sv2) / (sv1/lambdaR_eff1 + sv2/lambdaR_eff2)
                sv = G.cells.volumes - pv;
                w1 = sv(c1) ./ lambda_eff(c1);
                w2 = sv(c2) ./ lambda_eff(c2);
                wt = sv(c1) + sv(c2);
            end
            Th_nnc = CI .* wt ./ (w1 + w2);

            % ---- Concatenate: nf_internal + nNNC -----------------------
            Th = [Th_faces; Th_nnc];
        end
    end
end

%--------------------------------------------------------------------------
function tp = buildTwoPointOp(G)
% Build half-transmissibility operator over all cell-face adjacencies.
% Returns a function handle:  T_half = tp(lambda_eff)
% where lambda_eff is a per-cell vector and T_half has numel(G.cells.faces)
% entries (one per cell-face adjacency).
    cells = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
    faces = G.cells.faces(:, 1);
    C     = G.faces.centroids(faces, :) - G.cells.centroids(cells, :);
    sgn   = 2*(cells == G.faces.neighbors(faces, 1)) - 1;
    N     = bsxfun(@times, sgn, G.faces.normals(faces, :));
    cn    = sum(C .* N, 2) ./ sum(C .* C, 2);
    tp    = @(lambda) cn .* lambda(cells);
end

%--------------------------------------------------------------------------
function ha = buildHarmonicAvgOp(G)
% Build harmonic-average operator: face value = 1 / sum(1/half-values).
% Returns a function handle:  T_face = ha(T_half)
% where T_half has numel(G.cells.faces) entries and T_face is G.faces.num.
    faces = G.cells.faces(:, 1);
    M     = sparse(faces, 1:numel(faces), 1, G.faces.num, numel(faces));
    ha    = @(T) 1 ./ (M * (1 ./ T));
end
