classdef EDFMDynamicMolecularTransmissibility < StateFunction
%Dynamic molecular transmissibility for EDFM/pEDFM grids.
%
% SYNOPSIS:
%   prop = EDFMDynamicMolecularTransmissibility(model)
%
% DESCRIPTION:
%   Standard DynamicTransmissibility fails in EDFM/pEDFM because it builds
%   its harmonic-average operator over G.cells.faces (regular faces only),
%   producing G.faces.num values.  The downstream ComponentTotalDiffusiveFlux
%   then indexes with model.operators.internalConn which in EDFM has length
%   G.faces.num + nNNC - a size mismatch.
%
%   This class produces transmissibility vectors of size G.faces.num + nNNC
%   for each component with non-zero molecular diffusivity, so that the
%   existing internalConn filtering in ComponentTotalDiffusiveFlux works
%   correctly.
%
%   For each component:
%     1. FACE part (G.faces.num): per-cell effective diffusivity ->
%        two-point half-transmissibility -> harmonic average over faces.
%        Same formula as DynamicTransmissibility.
%
%     2. NNC part (nNNC): pore-volume-weighted harmonic average of per-cell
%        effective diffusivity, scaled by G.nnc.CI_dyn.  Molecular
%        diffusion occurs through the pore space, hence pore-volume
%        weighting (same as fluid heat transmissibility).
%
%   All terms carry full AD derivatives so the Jacobian is correct.
%
% PREREQUISITES:
%   G.nnc.CI_dyn must exist.  Call computeThermalNNCTransFracMatrix before
%   constructing the model.
%
% PARAMETERS:
%   model - GeothermalHFMModel instance.
%
% SEE ALSO:
%   EDFMDynamicHeatTransmissibility, DynamicTransmissibility,
%   MolecularDiffusivity, ComponentTotalDiffusiveFlux

    properties
        % Pre-built geometric operators (constructed once from G.cells.faces)
        twoPointOperator
        harmonicAvgOperator
        % Geometric area-reduction multiplier for matrix-matrix diffusion
        % (transmultEDFM, and transmultpEDFM for pEDFM) - the SAME factor
        % applied to fluid flow, by analogy with heat conduction.
        transmult
    end

    methods
        %-----------------------------------------------------------------%
        function prop = EDFMDynamicMolecularTransmissibility(model)
            prop@StateFunction(model);
            assert(isfield(model.G, 'nnc') && isfield(model.G.nnc, 'CI_dyn'), ...
                ['G.nnc.CI_dyn not found. ', ...
                 'Run computeThermalNNCTransFracMatrix before creating the model.']);

            prop = prop.dependsOn('MolecularDiffusivity');

            % Pre-build geometry-only operators from G.cells.faces
            prop.twoPointOperator    = buildTwoPointOp(model.G);
            prop.harmonicAvgOperator = buildHarmonicAvgOp(model.G);

            % Geometric area-reduction multiplier for matrix-matrix
            % diffusion - the SAME factor applied to fluid flow. Molecular
            % diffusion is a Laplacian (Fickian) flux through the pore space,
            % so the fracture-shadowed face area is reduced as for heat
            % conduction (by analogy; HosseiniMehr 2022 covers flow + heat).
            tol = 1e-5;
            if isprop(model, 'edfmTolerance'), tol = model.edfmTolerance; end
            prop.transmult = transmultEDFM(model.G, tol);
            if isfield(model.G.nnc, 'pMMneighs')
                prop.transmult = prop.transmult .* transmultpEDFM(model.G, tol);
            end

            prop.label = 'T_{mol}';
        end

        %-----------------------------------------------------------------%
        function T = evaluateOnDomain(prop, model, state)
        % Compute dynamic molecular transmissibility.
        %
        % Returns a cell array of transmissibility vectors, one per
        % component.  Each non-empty entry has length G.faces.num + nNNC,
        % ready for internalConn filtering by ComponentTotalDiffusiveFlux.

            G  = model.G;
            pv = model.operators.pv;  % pore volume per cell

            % Per-cell effective diffusivity from the state function:
            %   d{i} = D_i * tau * porosity   [m^2/s]
            lambda = prop.getEvaluatedDependencies(state, ...
                'MolecularDiffusivity');

            % NNC geometry
            c1 = G.nnc.cells(:, 1);
            c2 = G.nnc.cells(:, 2);
            CI = G.nnc.CI_dyn;

            T = cell(numel(lambda), 1);
            for i = 1:numel(lambda)
                if ~isempty(lambda{i})
                    % ---- Face part (G.faces.num) --------------------------
                    T_half  = prop.twoPointOperator(lambda{i});
                    T_faces = prop.harmonicAvgOperator(T_half);
                    % Handle negative transmissibility
                    fix = T_faces < 0;
                    T_faces(fix) = -T_faces(fix);
                    % Geometric area-reduction multiplier (same as flow)
                    T_faces = T_faces .* prop.transmult;

                    % ---- NNC part (nNNC) ----------------------------------
                    % Pore-volume-weighted harmonic average:
                    %   T_nnc = CI * (pv1+pv2) / (pv1/d1 + pv2/d2)
                    w1 = pv(c1) ./ lambda{i}(c1);
                    w2 = pv(c2) ./ lambda{i}(c2);
                    wt = pv(c1) + pv(c2);
                    T_nnc = CI .* wt ./ (w1 + w2);

                    % ---- Concatenate: G.faces.num + nNNC ------------------
                    T{i} = [T_faces; T_nnc];
                end
            end
        end
    end
end

%--------------------------------------------------------------------------
function tp = buildTwoPointOp(G)
% Build half-transmissibility operator over all cell-face adjacencies.
% Returns a function handle:  T_half = tp(lambda)
% where lambda is a per-cell vector and T_half has numel(G.cells.faces)
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
