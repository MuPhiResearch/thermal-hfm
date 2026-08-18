%% Example: DFM vs EDFM for Thermal Simulation
% Compares Discrete Fracture Model (DFM) and Embedded Discrete Fracture
% Model (EDFM) for thermal simulation in a 2D fractured reservoir using
% the SAME fracture geometry, physical properties, and wells.
%
%   - Single diagonal conductive fracture (half the domain diagonal)
%   - Injector at bottom-left, producer at top-right (along fracture)
%   - DFM: Unstructured Delaunay mesh with gradual refinement + hybrid cells
%   - EDFM: Regular Cartesian grid + Non-Neighboring Connections (NNCs)
%
% Sections:
%   Part 1: Shared problem definition
%   Part 2: DFM simulation
%   Part 3: EDFM simulation
%   Part 4: Physics comparison

clc; clear; close all;

%% 0. Load modules
mrstModule add dfm hfm geothermal thermal-hfm
mrstModule add ad-core ad-props incomp
checkLineSegmentIntersect;

%% =====================================================================
%  Part 1: Shared problem definition
%  =====================================================================
physdim = [200 200];           % domain size [m]
K0      = 273.15;              % 0 degC in Kelvin
pRef    = 50*barsa;            % reference pressure (shallow reservoir)
T0      = K0 + 80;            % initial reservoir temperature (80 degC)
T_inj   = K0 + 20;            % injection temperature (20 degC)

% --- Single diagonal conductive fracture (half the diagonal length) ---
% Fracture along the injector-producer diagonal, centered in the domain
% Diagonal length = physdim * sqrt(2) ~ 283 m; half = ~141 m
fl = [50, 50, 150, 150];      % from (50,50) to (150,150)

nFrac    = size(fl, 1);
fracAp   = 0.05;              % fracture aperture [m] (wider → better DFM conditioning)
K_frac   = 1000;              % fracture permeability [darcy]

% --- Rock properties ---
K_mat   = 10*milli*darcy;     % matrix permeability
phi_mat = 0.2;                % matrix porosity
lambdaR = 2.0;                % rock thermal conductivity [W/(m*K)]
rhoR    = 2700;               % rock density [kg/m3]
CpR     = 880;                % rock heat capacity [J/(kg*K)]

% --- Fluid properties ---
lambdaF = 0.6;                % fluid thermal conductivity [W/(m*K)]
Cp_f    = 4200;               % fluid heat capacity [J/(kg*K)]

% --- Well parameters ---
well_rate  = 5e-5;            % injection rate [m3/s] (~4.3 m3/day)
well_bhp   = pRef;            % production BHP
well_rad   = 0.1;             % wellbore radius [m]

% --- Simulation parameters ---
sim_time  = 10*year;
dt_target = 30*day;

fprintf('=== DFM vs EDFM Comparison ===\n');
fprintf('Domain: %.0f x %.0f m, 1 diagonal fracture (%.0f m)\n', ...
    physdim, norm(fl(3:4) - fl(1:2)));
fprintf('Reservoir T = %.0f degC, Injection T = %.0f degC\n', T0-K0, T_inj-K0);

%% =====================================================================
%  Part 2: DFM Simulation
%  =====================================================================
fprintf('\n========== DFM SIMULATION ==========\n');
tic_dfm = tic;

%% 2.1 Create Delaunay mesh with gradual refinement around fracture
h_coarse = 5;   % coarse cell size [m] - away from fracture
h_fine   = 1.5; % fine cell size [m] - near fracture

frac_pts = [fl(1:2); fl(3:4)];

% Generate point cloud with gradual density based on distance to fracture
p = [];
for hh = h_fine : 0.5 : h_coarse
    [X, Y] = meshgrid(0:hh:physdim(1), 0:hh:physdim(2));
    pts = [X(:), Y(:)];

    % Distance to fracture line segment
    v = frac_pts(2,:) - frac_pts(1,:);
    w = pts - frac_pts(1,:);
    t = max(0, min(1, (w * v') / (v * v')));
    proj = frac_pts(1,:) + t * v;
    d_frac = sqrt(sum((pts - proj).^2, 2));

    % Keep points where the local desired spacing matches this grid level
    h_desired = h_fine + (h_coarse - h_fine) * min(d_frac / 30, 1);
    keep = abs(hh - h_desired) < 0.5;
    p = [p; pts(keep, :)]; %#ok<AGROW>
end

% Remove duplicate/close points
p = unique(round(p / (h_fine * 0.3)) * (h_fine * 0.3), 'rows');

% Fracture vertices and constraints
vertices = [fl(:, 1:2); fl(:, 3:4)];
constraints = [reshape(1:2*nFrac, 2, [])' (1:nFrac)'];
box = [0 0; physdim];

% Remove points too close to fractures
perimiter = h_fine * 0.8;
p = remove_closepoints(vertices, constraints, p, perimiter);

% Handle fracture intersections (very small precision for straight lines)
args = struct('precision', 0.001);
[vertices, constraints] = removeFractureIntersections( ...
    vertices, constraints, box, args);

numOrdPt = size(p, 1);
p = [p; vertices];
constraints(:, 1:2) = constraints(:, 1:2) + numOrdPt;

% Add points along fracture edges
[p, constraints, ~] = partition_edges(p, constraints, h_fine * 0.4, box, args);

tags = constraints(:, 3);
constraints_c = constraints(:, 1:2);

% Constrained Delaunay triangulation
try
    delTri = DelaunayTri(p, constraints_c); %#ok<DDELTRI>
catch
    delTri = delaunayTriangulation(p, constraints_c);
end
G_dfm = triangleGrid(delTri.X, delTri.Triangulation);
G_dfm = computeGeometry(G_dfm);

%% 2.2 Tag fracture faces and add hybrid cells
G_dfm.faces.tags = zeros(G_dfm.faces.num, 1);
faceNodes = sort(reshape(G_dfm.faces.nodes, 2, [])', 2);
try
    constraints_sorted = sort(delTri.Constraints, 2);
catch
    constraints_sorted = sort(constraints_c, 2);
end

for iter = 1:size(constraints_sorted, 1)
    fracFace = find(ismember(faceNodes, constraints_sorted(iter, :), 'rows'));
    G_dfm.faces.tags(fracFace) = tags(iter);
end

aperture_dfm = zeros(G_dfm.faces.num, 1);
aperture_dfm(G_dfm.faces.tags > 0) = fracAp;

G_dfm = addhybrid(G_dfm, G_dfm.faces.tags > 0, aperture_dfm);

%% 2.3 Rock and fluid for DFM
hybridInd   = find(G_dfm.cells.hybrid);
nCells_dfm  = G_dfm.cells.num;
nMatrix_dfm = nCells_dfm - numel(hybridInd);
nFrac_dfm   = numel(hybridInd);

rock_dfm = makeRock(G_dfm, K_mat, phi_mat);
rock_dfm.perm(hybridInd, :) = K_frac * darcy;
rock_dfm.poro(hybridInd)    = 0.5;

rock_dfm = addThermalRockProps(rock_dfm, ...
    'lambdaR', lambdaR, 'rhoR', rhoR, 'CpR', CpR);

fluid_dfm = initSimpleADIFluid('phases', 'W', ...
    'mu', 1e-3, 'rho', 1000, 'pRef', pRef, 'c', 0, 'cR', 0);
fluid_dfm = addThermalFluidProps(fluid_dfm, ...
    'Cp', Cp_f, 'lambdaF', lambdaF, 'useEOS', false);

%% 2.4 Create GeothermalModel
gravity reset off;
model_dfm = GeothermalModel(G_dfm, rock_dfm, fluid_dfm);
% Pre-build the geothermal flux discretization so validateModel does not emit
% "Assuming default flux discretization". This is the same default it builds.
model_dfm.FlowDiscretization = GeothermalFlowDiscretization(model_dfm);
model_dfm.extraStateOutput     = true;
model_dfm.outputFluxes         = true;
model_dfm.minimumTemperature   = K0;
model_dfm.maximumTemperature   = K0 + 200*Kelvin;
model_dfm.maximumPressure      = 200e6*Pascal;

%% 2.5 Wells
cc_dfm = G_dfm.cells.centroids;
mat_cells_dfm = find(~G_dfm.cells.hybrid);

% Injector at bottom-left
[~, inj_idx] = min(cc_dfm(mat_cells_dfm,1).^2 + cc_dfm(mat_cells_dfm,2).^2);
inj_cell_dfm = mat_cells_dfm(inj_idx);

% Producer at top-right
[~, prod_idx] = min((cc_dfm(mat_cells_dfm,1)-physdim(1)).^2 + ...
                    (cc_dfm(mat_cells_dfm,2)-physdim(2)).^2);
prod_cell_dfm = mat_cells_dfm(prod_idx);

W_dfm = addWell([], G_dfm, rock_dfm, inj_cell_dfm, ...
    'type', 'rate', 'val', well_rate, 'comp_i', 1, ...
    'Name', 'Inj', 'sign', 1, 'Radius', well_rad, ...
    'InnerProduct', 'ip_simple');
W_dfm = addWell(W_dfm, G_dfm, rock_dfm, prod_cell_dfm, ...
    'type', 'bhp', 'val', well_bhp, 'comp_i', 1, ...
    'Name', 'Prod', 'sign', -1, 'Radius', well_rad, ...
    'InnerProduct', 'ip_simple');
W_dfm = addThermalWellProps(W_dfm, G_dfm, rock_dfm, fluid_dfm, 'T', T_inj);

%% 2.6 Initial state and run
state0_dfm   = initResSol(G_dfm, pRef, 1);
state0_dfm.T = T0 * ones(G_dfm.cells.num, 1);

dt_dfm       = rampupTimesteps(sim_time, dt_target);
schedule_dfm = simpleSchedule(dt_dfm, 'W', W_dfm);

nls_dfm = NonLinearSolver('maxIterations', 25, 'maxTimestepCuts', 12);

fprintf('DFM grid: %d cells (%d matrix + %d fracture)\n', ...
    nCells_dfm, nMatrix_dfm, nFrac_dfm);
fprintf('Running DFM simulation (%d timesteps)...\n', numel(dt_dfm));

[~, states_dfm, report_dfm] = simulateScheduleAD(state0_dfm, model_dfm, ...
    schedule_dfm, 'NonLinearSolver', nls_dfm);
time_dfm = toc(tic_dfm);
fprintf('DFM simulation complete (%.1f s).\n', time_dfm);

%% =====================================================================
%  Part 3: EDFM Simulation - SAME fracture and wells
%  =====================================================================
fprintf('\n========== EDFM SIMULATION ==========\n');
tic_edfm = tic;

%% 3.1 Cartesian grid
celldim_edfm = [40 40];
G_edfm = cartGrid(celldim_edfm, physdim);
G_edfm = computeGeometry(G_edfm);

%% 3.2 EDFM preprocessing - same fracture line fl
[G_edfm, fracture] = processFracture2D(G_edfm, fl);
fracture.aperture = fracAp;
G_edfm = CIcalculator2D(G_edfm, fracture);
[G_edfm, F, fracture] = gridFracture2D(G_edfm, fracture, ...
    'min_size', 0.05, 'cell_size', 0.1);

%% 3.3 Rock
G_edfm.rock.perm = ones(G_edfm.cells.num, 1) * K_mat;
G_edfm.rock.poro = phi_mat * ones(G_edfm.cells.num, 1);
G_edfm = makeRockFrac(G_edfm, K_frac, 'porosity', 0.5);

%% 3.4 Flow NNCs
[G_edfm, ~] = defineNNCandTrans(G_edfm, F, fracture);

nMatrix_edfm = celldim_edfm(1) * celldim_edfm(2);
nFrac_edfm   = G_edfm.cells.num - nMatrix_edfm;
nNNC_edfm    = size(G_edfm.nnc.cells, 1);

%% 3.5 Thermal properties
G_edfm.rock = addThermalRockProps(G_edfm.rock, ...
    'lambdaR', lambdaR, 'rhoR', rhoR, 'CpR', CpR);
rock_edfm = G_edfm.rock;

fluid_edfm = initSimpleADIFluid('phases', 'W', ...
    'mu', 1e-3, 'rho', 1000, 'pRef', pRef, 'c', 0, 'cR', 0);
fluid_edfm = addThermalFluidProps(fluid_edfm, ...
    'Cp', Cp_f, 'lambdaF', lambdaF, 'useEOS', false);

%% 3.6 Thermal NNCs
G_edfm = computeThermalNNCTransFracMatrix(G_edfm, rock_edfm, fluid_edfm);
G_edfm = computeThermalNNCTransFracFrac(G_edfm, F, fracture, rock_edfm, fluid_edfm);

%% 3.7 Create GeothermalHFMModel (EDFM)
model_edfm = GeothermalHFMModel(G_edfm, rock_edfm, fluid_edfm, 'fractureMethod', 'edfm');
model_edfm.extraStateOutput     = true;
model_edfm.outputFluxes         = true;
model_edfm.minimumTemperature   = K0;
model_edfm.maximumTemperature   = K0 + 200*Kelvin;
model_edfm.maximumPressure      = 200e6*Pascal;

%% 3.8 Wells - same locations as DFM
cc_edfm = G_edfm.cells.centroids;
inj_cell_edfm  = 1;                % bottom-left = cell 1
prod_cell_edfm = nMatrix_edfm;     % top-right = cell Nx*Ny

W_edfm = addWell([], G_edfm, rock_edfm, inj_cell_edfm, ...
    'type', 'rate', 'val', well_rate, 'comp_i', 1, ...
    'Name', 'Inj', 'sign', 1, 'Radius', well_rad);
W_edfm = addWell(W_edfm, G_edfm, rock_edfm, prod_cell_edfm, ...
    'type', 'bhp', 'val', well_bhp, 'comp_i', 1, ...
    'Name', 'Prod', 'sign', -1, 'Radius', well_rad);
W_edfm = addThermalWellProps(W_edfm, G_edfm, rock_edfm, fluid_edfm, 'T', T_inj);

%% 3.9 Run EDFM
state0_edfm   = initResSol(G_edfm, pRef, 1);
state0_edfm.T = T0 * ones(G_edfm.cells.num, 1);

dt_edfm       = rampupTimesteps(sim_time, dt_target);
schedule_edfm = simpleSchedule(dt_edfm, 'W', W_edfm);

fprintf('EDFM grid: %d cells (%d matrix + %d fracture), %d NNCs\n', ...
    G_edfm.cells.num, nMatrix_edfm, nFrac_edfm, nNNC_edfm);
fprintf('Running EDFM simulation (%d timesteps)...\n', numel(dt_edfm));

[~, states_edfm, report_edfm] = simulateScheduleAD(state0_edfm, model_edfm, schedule_edfm);
time_edfm = toc(tic_edfm);
fprintf('EDFM simulation complete (%.1f s).\n', time_edfm);

%% =====================================================================
%  Part 4: Physics Comparison
%  =====================================================================
times_days = cumsum(dt_edfm) / day;

%% Figure 1: Grid comparison
figure('Name', 'Grid Comparison', 'Position', [50, 50, 1200, 500]);

subplot(1, 2, 1);
plotGrid_DFM(G_dfm, 'FaceColor', [0.9 0.9 1], 'EdgeAlpha', 0.3);
plotFractures(G_dfm);
hold on;
plot(cc_dfm(inj_cell_dfm,1), cc_dfm(inj_cell_dfm,2), ...
    'rv', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'LineWidth', 2);
plot(cc_dfm(prod_cell_dfm,1), cc_dfm(prod_cell_dfm,2), ...
    'b^', 'MarkerSize', 14, 'MarkerFaceColor', 'b', 'LineWidth', 2);
hold off;
title(sprintf('DFM: %d cells (%d matrix + %d hybrid)', ...
    nCells_dfm, nMatrix_dfm, nFrac_dfm));
axis equal tight; set(gca, 'FontSize', 11);

subplot(1, 2, 2);
plotCellData(G_edfm, log10(G_edfm.rock.perm / darcy), 'EdgeAlpha', 0.1);
hold on;
line(fl(:,[1 3])', fl(:,[2 4])', 'Color', 'k', 'LineWidth', 2);
plot(cc_edfm(inj_cell_edfm,1), cc_edfm(inj_cell_edfm,2), ...
    'rv', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'LineWidth', 2);
plot(cc_edfm(prod_cell_edfm,1), cc_edfm(prod_cell_edfm,2), ...
    'b^', 'MarkerSize', 14, 'MarkerFaceColor', 'b', 'LineWidth', 2);
hold off;
colormap(gca, jet(25)); cbar = colorbar; title(cbar, 'log_{10}(K/D)');
title(sprintf('EDFM: %d cells (Cartesian + %d NNCs)', ...
    G_edfm.cells.num, nNNC_edfm));
axis equal tight; view(0, 90); set(gca, 'FontSize', 11);

%% Figure 2: Temperature snapshots
snap_targets = [365, 1825, 3650];
snap_idx = zeros(size(snap_targets));
for k = 1:numel(snap_targets)
    [~, snap_idx(k)] = min(abs(times_days - snap_targets(k)));
end

figure('Name', 'DFM vs EDFM Temperature', 'Position', [50, 50, 1200, 700]);
for k = 1:numel(snap_idx)
    si = snap_idx(k);

    % DFM (top row)
    subplot(2, numel(snap_idx), k);
    plotCellData_DFM(G_dfm, states_dfm{si}.T - K0);
    plotFractures(G_dfm, hybridInd, states_dfm{si}.T - K0);
    colormap(gca, jet(25)); caxis([20 80]);
    axis equal tight;
    if k == 1, ylabel('DFM', 'FontWeight', 'bold'); end
    title(sprintf('t = %d d', round(times_days(si))));
    if k == numel(snap_idx)
        cbar = colorbar; title(cbar, 'T [degC]');
    end

    % EDFM (bottom row)
    subplot(2, numel(snap_idx), numel(snap_idx) + k);
    plotCellData(G_edfm, states_edfm{si}.T - K0, 'EdgeColor', 'none');
    hold on;
    line(fl(:,[1 3])', fl(:,[2 4])', 'Color', 'w', 'LineWidth', 1);
    hold off;
    colormap(gca, jet(25)); caxis([20 80]);
    axis equal tight; view(0, 90);
    if k == 1, ylabel('EDFM', 'FontWeight', 'bold'); end
    if k == numel(snap_idx)
        cbar = colorbar; title(cbar, 'T [degC]');
    end
end

%% Figure 3: Temperature cross-section along diagonal (y = x)
x_dfm = G_dfm.cells.centroids(:, 1);
y_dfm = G_dfm.cells.centroids(:, 2);
x_edfm_c = G_edfm.cells.centroids(:, 1);
y_edfm_c = G_edfm.cells.centroids(:, 2);

mat_mask_dfm  = ~G_dfm.cells.hybrid;
mat_mask_edfm = (1:G_edfm.cells.num)' <= nMatrix_edfm;

% Select cells along diagonal (|y - x| < tolerance)
diag_tol_dfm  = h_coarse;
diag_tol_edfm = physdim(1) / celldim_edfm(1);
sel_dfm  = abs(y_dfm - x_dfm) < diag_tol_dfm & mat_mask_dfm;
sel_edfm = abs(y_edfm_c - x_edfm_c) < diag_tol_edfm & mat_mask_edfm;

% Distance along diagonal
d_dfm  = sqrt(x_dfm(sel_dfm).^2 + y_dfm(sel_dfm).^2);
d_edfm = sqrt(x_edfm_c(sel_edfm).^2 + y_edfm_c(sel_edfm).^2);
T_cs_dfm  = states_dfm{end}.T(sel_dfm) - K0;
T_cs_edfm = states_edfm{end}.T(sel_edfm) - K0;
[d_dfm, idx] = sort(d_dfm); T_cs_dfm = T_cs_dfm(idx);
[d_edfm, idx] = sort(d_edfm); T_cs_edfm = T_cs_edfm(idx);

figure('Name', 'Temperature Cross-Section', 'Position', [100, 100, 700, 450]);
plot(d_dfm, T_cs_dfm, 'b-o', 'LineWidth', 2, 'MarkerSize', 3, ...
    'DisplayName', 'DFM');
hold on;
plot(d_edfm, T_cs_edfm, 'r-s', 'LineWidth', 2, 'MarkerSize', 3, ...
    'DisplayName', 'EDFM');
hold off;
xlabel('Distance along diagonal [m]'); ylabel('Temperature [degC]');
title('Final Temperature Along Diagonal (Inj -> Prod)');
legend('Location', 'northeast'); grid on; set(gca, 'FontSize', 12);

%% Figure 4: Producer cell temperature over time
T_prod_dfm  = zeros(numel(states_dfm), 1);
T_prod_edfm = zeros(numel(states_edfm), 1);
for j = 1:numel(states_dfm)
    if ~isempty(states_dfm{j})
        T_prod_dfm(j) = states_dfm{j}.T(prod_cell_dfm) - K0;
    end
end
for j = 1:numel(states_edfm)
    if ~isempty(states_edfm{j})
        T_prod_edfm(j) = states_edfm{j}.T(prod_cell_edfm) - K0;
    end
end

figure('Name', 'Thermal Breakthrough', 'Position', [100, 100, 700, 450]);
plot(times_days/365, T_prod_dfm, 'b-', 'LineWidth', 2, 'DisplayName', 'DFM');
hold on;
plot(times_days/365, T_prod_edfm, 'r--', 'LineWidth', 2, 'DisplayName', 'EDFM');
yline(T_inj - K0, 'k:', 'T_{inj}', 'LineWidth', 1);
yline(T0 - K0, 'k:', 'T_0', 'LineWidth', 1);
hold off;
xlabel('Time [years]'); ylabel('Producer cell temperature [degC]');
title('Thermal Breakthrough: DFM vs EDFM');
legend('Location', 'southwest'); grid on; set(gca, 'FontSize', 12);

fprintf('\n=== Physics Comparison ===\n');
fprintf('DFM  final producer cell T: %.1f degC\n', T_prod_dfm(end));
fprintf('EDFM final producer cell T: %.1f degC\n', T_prod_edfm(end));
fprintf('Difference: %.1f degC\n', abs(T_prod_dfm(end) - T_prod_edfm(end)));

%% Physics comparison summary
fprintf('\n=== Summary ===\n');
fprintf('DFM and EDFM use the same fracture geometry, properties, and wells.\n');
fprintf('EDFM keeps a regular Cartesian matrix grid and represents the\n');
fprintf('fracture through Non-Neighbouring Connections (NNCs), whereas DFM\n');
fprintf('resolves the fracture with an unstructured mesh and hybrid cells.\n');
fprintf('The producer-temperature histories agree to within %.1f degC.\n', ...
    abs(T_prod_dfm(end) - T_prod_edfm(end)));
