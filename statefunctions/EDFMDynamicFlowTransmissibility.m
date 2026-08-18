classdef EDFMDynamicFlowTransmissibility < StateFunction
%Dynamic flow transmissibility for EDFM/pEDFM grids.
%
% SYNOPSIS:
%   prop = EDFMDynamicFlowTransmissibility(model)
%
% DESCRIPTION:
%   Standard DynamicTransmissibility fails in EDFM/pEDFM because it builds
%   its harmonic-average operator over G.cells.faces (regular faces only),
%   then indexes the result with model.operators.internalConn which in EDFM
%   is extended to length G.faces.num + nNNC.  The result is a size mismatch.
%
%   This class produces the correct (nf_internal + nNNC) flow
%   transmissibility vector by computing two contributions separately and
%   concatenating:
%
%     1. FACE part (nf_internal): per-cell permeability ->
%        two-point half-transmissibility -> harmonic average over faces ->
%        apply EDFM (and pEDFM) area-reduction multiplier ->
%        filter to internal faces.
%
%     2. NNC part (nNNC): pore-volume-weighted harmonic average of per-cell
%        permeability, using the geometric contact indicator G.nnc.CI_dyn
%        stored by computeThermalNNCTransFracMatrix:
%          T_nnc = CI * (pv1+pv2) / (pv1/k1 + pv2/k2)
%
%   The NNC terms carry full AD derivatives with respect to pressure and
%   temperature, so the Jacobian correctly accounts for permeability changes
%   across Newton iterations.
%
% PREREQUISITES:
%   G.nnc.CI_dyn must exist.  Call computeThermalNNCTransFracMatrix before
%   constructing the model; it stores CI_dyn automatically.
%
% PARAMETERS:
%   model - GeothermalHFMModel instance.
%
% SEE ALSO:
%   EDFMDynamicHeatTransmissibility, EDFMDynamicMolecularTransmissibility,
%   DynamicTransmissibility, HeatPermeability

    properties
        % Pre-built geometric operators (constructed once from G.cells.faces)
        twoPointOperator
        harmonicAvgOperator
        % Precomputed face transmissibility multiplier (G.faces.num vector).
        % For EDFM: transmultEDFM.  For pEDFM: transmultEDFM * transmultpEDFM.
        transmult
    end

    methods
        %-----------------------------------------------------------------%
        function prop = EDFMDynamicFlowTransmissibility(model)
            prop@StateFunction(model);

            G = model.G;
            assert(isfield(G, 'nnc') && isfield(G.nnc, 'CI_dyn'), ...
                ['G.nnc.CI_dyn not found. ', ...
                 'Run computeThermalNNCTransFracMatrix before creating the model.']);

            % Depend on the per-cell Permeability state function
            % (HeatPermeability, installed by GeothermalFlowDiscretization)
            prop = prop.dependsOn('Permeability');

            % Pre-build geometry-only operators from G.cells.faces
            prop.twoPointOperator    = buildTwoPointOp(G);
            prop.harmonicAvgOperator = buildHarmonicAvgOp(G);

            % Precompute EDFM face-area multiplier (purely geometric)
            tol = 1e-5;
            if isprop(model, 'edfmTolerance')
                tol = model.edfmTolerance;
            end
            prop.transmult = transmultEDFM(G, tol);
            % For pEDFM: also apply projected-connection multiplier
            if isfield(G.nnc, 'pMMneighs')
                prop.transmult = prop.transmult .* transmultpEDFM(G, tol);
            end

            prop.label = 'T';
        end

        %-----------------------------------------------------------------%
        function T = evaluateOnDomain(prop, model, state, allFaces)
        % Compute dynamic flow transmissibility.
        %
        % allFaces = false (default):
        %   Returns a vector of length nf_internal + nNNC suitable for use
        %   with model.operators.Grad and model.operators.AccDiv in the
        %   flow equations.
        %
        % allFaces = true:
        %   Returns a vector of length G.faces.num (all faces, no NNCs),
        %   indexed by absolute face number.  Required by getForceProperties
        %   and by GeothermalModel.getFlowEquations which stores the result
        %   as model.operators.T_all for BC flux computation.

            if nargin < 4, allFaces = false; end

            G  = model.G;
            pv = model.operators.pv;  % pore volume per cell

            % Per-cell permeability from HeatPermeability state function
            k = prop.getEvaluatedDependencies(state, 'Permeability');

            % ---- Regular face contributions (size G.faces.num) ---------
            % Two-point half-transmissibilities at each (cell, face) entry
            T_half  = prop.twoPointOperator(k);
            % Harmonic average: one value per face (size G.faces.num)
            T_faces = prop.harmonicAvgOperator(T_half);
            % Fix negative transmissibilities (same as DynamicTransmissibility)
            fix = T_faces < 0;
            T_faces(fix) = -T_faces(fix);
            % Apply EDFM (and pEDFM) area-reduction multiplier
            T_faces = T_faces .* prop.transmult;

            if allFaces
                % Return all G.faces.num values so getForceProperties can
                % index by absolute BC face number and getFlowEquations
                % can store as model.operators.T_all.
                T = T_faces;
                return;
            end

            % ---- Filter to internal faces (nf_internal) ----------------
            % model.operators.internalConn is logical [nf_all + nNNC].
            % The first nf_all entries are the standard internal-face mask.
            internalFacesMask = model.operators.internalConn(1:G.faces.num);
            T_faces_int = T_faces(internalFacesMask);   % nf_internal values

            % ---- NNC contributions (nNNC) -------------------------------
            % Pore-volume-weighted harmonic average of permeability:
            %   T_nnc = CI * (pv1+pv2) / (pv1/k1 + pv2/k2)
            % Same formula as fracturematrixNNC3D uses for flow.
            % CI_dyn encodes the geometric contact indicator including the
            % adjustT correction (back-computed from G.nnc.T).
            c1 = G.nnc.cells(:, 1);
            c2 = G.nnc.cells(:, 2);
            CI = G.nnc.CI_dyn;

            w1 = pv(c1) ./ k(c1);
            w2 = pv(c2) ./ k(c2);
            wt = pv(c1) + pv(c2);
            T_nnc = CI .* wt ./ (w1 + w2);

            % ---- Concatenate: nf_internal + nNNC -----------------------
            T = [T_faces_int; T_nnc];
        end
    end
end

%--------------------------------------------------------------------------
function tp = buildTwoPointOp(G)
% Build half-transmissibility operator over all cell-face adjacencies.
% Returns a function handle:  T_half = tp(perm)
% where perm is a per-cell vector and T_half has numel(G.cells.faces)
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
