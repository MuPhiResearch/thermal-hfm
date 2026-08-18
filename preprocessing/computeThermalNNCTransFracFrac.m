function G = computeThermalNNCTransFracFrac(G, F, fracture, rock, fluid)
% Compute thermal NNC transmissibilities for fracture-fracture connections.
%
% SYNOPSIS:
%   G = computeThermalNNCTransFracFrac(G, F, fracture, rock, fluid)
%
% DESCRIPTION:
%   Computes rock and fluid thermal transmissibilities for
%   fracture-fracture NNC connections using the star-delta transformation.
%   Must be called AFTER defineNNCandTrans (which calls frac_frac_nnc).
%
%   The function re-runs the intersection logic from frac_frac_nnc but
%   computes half-transmissibilities using thermal conductivity instead of
%   permeability, then applies the same star-delta formula.
%
%   If G.nnc.transHr/transHf already exist (from computeThermalNNCTrans
%   FracMatrix), the frac-frac values are APPENDED in the same order as
%   the flow NNCs were appended by frac_frac_nnc.
%
% PARAMETERS:
%   G        - Grid with G.nnc already populated by defineNNCandTrans.
%   F        - Fracture structure from gridFracture2D.
%   fracture - Fracture structure from gridFracture2D.
%   rock     - Rock with rock.lambdaR [W/(m*K)].
%   fluid    - Fluid with fluid.lambdaF [W/(m*K)].
%
% RETURNS:
%   G - Grid with G.nnc.transHr and G.nnc.transHf updated to include
%       frac-frac thermal NNC transmissibilities.

%{
Thermal EDFM Module
%}

    % Check if there are fracture intersections
    if ~isfield(fracture, 'intersections')
        % No intersections: nothing to do, but ensure fields exist
        if ~isfield(G.nnc, 'transHr')
            G.nnc.transHr = [];
        end
        if ~isfield(G.nnc, 'transHf')
            G.nnc.transHf = [];
        end
        return;
    end

    % Initialize transHr/transHf if not already present
    if ~isfield(G.nnc, 'transHr')
        G.nnc.transHr = [];
    end
    if ~isfield(G.nnc, 'transHf')
        G.nnc.transHf = [];
    end

    % Get thermal conductivities (scalar -> per-cell)
    nc = G.cells.num;
    lambdaR = rock.lambdaR;
    if isscalar(lambdaR), lambdaR = repmat(lambdaR, nc, 1); end
    lambdaF = fluid.lambdaF;
    if isscalar(lambdaF), lambdaF = repmat(lambdaF, nc, 1); end

    % --- Re-run intersection logic (mirrors frac_frac_nnc.m) ---
    for i = 1:size(fracture.intersections.lines, 1)
        lines  = fracture.intersections.lines(i, :);
        coords = fracture.intersections.coords(i, :);

        [~, Gface(1)] = ismember(roundsd(coords, 5), ...
            roundsd(F(lines(1)).nodes.coords, 5), 'rows');
        [~, Gface(2)] = ismember(roundsd(coords, 5), ...
            roundsd(F(lines(2)).nodes.coords, 5), 'rows');

        diff1 = diff([F(lines(1)).nodes.coords(1,:); ...
                      F(lines(1)).nodes.coords(end,:)]);
        diff2 = diff([F(lines(2)).nodes.coords(1,:); ...
                      F(lines(2)).nodes.coords(end,:)]);

        % Identify cells at the intersection
        cells_l = cell(2, 1);
        for j = 1:numel(Gface)
            if Gface(j) == 1
                cells_l{j, 1} = Gface(j);
            elseif Gface(j) == size(F(lines(j)).nodes.coords, 1)
                cells_l{j, 1} = Gface(j) - 1;
            else
                cells_l{j, 1} = [Gface(j)-1, Gface(j)];
            end
        end

        Gf1 = G.FracGrid.(['Frac', num2str(lines(1))]);
        Gf2 = G.FracGrid.(['Frac', num2str(lines(2))]);

        cells_l{1, 1} = Gf1.cells.start - 1 + cells_l{1, 1};
        cells_l{2, 1} = Gf2.cells.start - 1 + cells_l{2, 1};

        % --- Compute thermal half-transmissibilities ---
        % Rock thermal conductivity for fracture cells
        [Thr1, Thf1] = computeFracHalfTransThermal(Gf1, lambdaR, lambdaF);
        [Thr2, Thf2] = computeFracHalfTransThermal(Gf2, lambdaR, lambdaF);

        % Select faces at intersection (same logic as frac_frac_nnc)
        if diff1(1) == 0
            Gface(1) = Gf1.faces.num - ...
                (size(F(lines(1)).nodes.coords, 1) - Gface(1));
        end
        if diff2(1) == 0
            Gface(2) = Gf2.faces.num - ...
                (size(F(lines(2)).nodes.coords, 1) - Gface(2));
        end

        Thr1 = Thr1(Gf1.cells.faces(:, 1) == Gface(1));
        Thf1 = Thf1(Gf1.cells.faces(:, 1) == Gface(1));
        Thr2 = Thr2(Gf2.cells.faces(:, 1) == Gface(2));
        Thf2 = Thf2(Gf2.cells.faces(:, 1) == Gface(2));

        % --- Star-delta transformation ---
        [~, ~] = meshgrid(cells_l{1, 1}, cells_l{2, 1});
        [tr1, tr2] = meshgrid(Thr1, Thr2);
        [tf1, tf2] = meshgrid(Thf1, Thf2);
        tpairs_r = [tr1(:), tr2(:)];
        tpairs_f = [tf1(:), tf2(:)];

        if numel(cells_l{2, 1}) > 1
            tpairs_r = [Thr2'; tpairs_r]; %#ok
            tpairs_f = [Thf2'; tpairs_f]; %#ok
        end
        if numel(cells_l{1, 1}) > 1
            tpairs_r = [Thr1'; tpairs_r]; %#ok
            tpairs_f = [Thf1'; tpairs_f]; %#ok
        end

        % Star-delta: T_pair = prod(T_half_pair) / sum(all_T_halves)
        n_new = size(tpairs_r, 1);
        transHr_new = prod(tpairs_r, 2) / sum([sum(Thr1), sum(Thr2)]);
        transHf_new = prod(tpairs_f, 2) / sum([sum(Thf1), sum(Thf2)]);

        G.nnc.transHr = [G.nnc.transHr; transHr_new];
        G.nnc.transHf = [G.nnc.transHf; transHf_new];
    end

    % --- Sanity checks ---
    assert(numel(G.nnc.transHr) == size(G.nnc.cells, 1), ...
        'transHr size (%d) does not match G.nnc.cells (%d)', ...
        numel(G.nnc.transHr), size(G.nnc.cells, 1));
    assert(numel(G.nnc.transHf) == size(G.nnc.cells, 1), ...
        'transHf size (%d) does not match G.nnc.cells (%d)', ...
        numel(G.nnc.transHf), size(G.nnc.cells, 1));

end


function [Thr, Thf] = computeFracHalfTransThermal(Gf, lambdaR, lambdaF)
% Compute thermal half-transmissibilities for a fracture grid.
%
% Uses the same geometry as computeTrans but with thermal conductivity
% instead of permeability. The effective conductivity is volume-weighted:
%   Rock:  lambdaR * (vol - pv) / vol
%   Fluid: lambdaF * pv / vol

    vol = Gf.cells.volumes;
    if isfield(Gf.rock, 'poro')
        pv = poreVolume(Gf, Gf.rock);
    else
        pv = vol;
    end

    % Use global lambdaR value (scalar applied to all frac cells)
    % The scalar has already been expanded to per-cell before calling
    lambdaR_val = lambdaR(1); % Fractures use a single value
    lambdaF_val = lambdaF(1);

    % Rock thermal conductivity: weighted by solid fraction
    lambdaR_eff = lambdaR_val .* (vol - pv) ./ vol;
    rock_r = struct('perm', lambdaR_eff);
    Thr = computeTrans(Gf, rock_r);

    % Fluid thermal conductivity: weighted by pore fraction
    lambdaF_eff = lambdaF_val .* pv ./ vol;
    rock_f = struct('perm', lambdaF_eff);
    Thf = computeTrans(Gf, rock_f);
end
