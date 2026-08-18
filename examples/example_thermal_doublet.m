%% Example: Thermal Doublet - EDFM vs pEDFM in a Heterogeneous Fractured Reservoir
% Simulates a geothermal doublet (cold injection + hot production) in a
% 3D heterogeneous fractured reservoir and compares standard EDFM against
% projection-based EDFM (pEDFM).
%
% Physical setup:
%   - Domain: 600 x 400 x 60 m,  40 x 28 x 6 Cartesian grid (6720 cells).
%   - Reservoir at 120 degC (hot),  cold water injected at 20 degC.
%   - Permeability: 3D spatially correlated log-normal (Gaussian random field)
%     with K_mean = 80 mD, sigma_lnK = 2.5, correlation length 120 m (horiz).
%   - Fracture network (4 planes):
%       Frac 1  BARRIER  (1 nD)   - near injector (~45 m), blocks early sweeping.
%       Frac 2  CONDUIT  (5000 D) - main diagonal SW→NE bypass.
%       Frac 3  CONDUIT  (2000 D) - secondary NW→SE diagonal bypass.
%       Frac 4  BARRIER  (1 nD)   - near producer (~35 m), blocks direct approach.
%   - Injector: rate-controlled, perforates all 6 z-layers at x≈75 m.
%   - Producer: BHP-controlled, perforates all 6 z-layers at x≈525 m.
%   - Simulation: 20 years, dynamic lambda_F via EDFMDynamicHeatTransmissibility.
%
% Figures produced:
%   1. Setup:  Gaussian perm slice (mid-layer) + 3D fracture/well overview.
%   2. Cold-plume snapshots (pEDFM): 4 selected timesteps (with grid).
%   3. Final temperature comparison:  EDFM vs pEDFM (matrix, top view).
%   4. Thermal analysis: (a) producer breakthrough, (b) near-barrier zone.
%   5. Temperature difference pEDFM - EDFM (final state).
%   6. Final pressure comparison:  EDFM vs pEDFM (matrix, top view).
%   7. INTERACTIVE viewer: MRST plotToolbar with time-step player.

clc; clear; close all;

%% 1. Load modules
mrstModule add hfm shale geothermal thermal-hfm
mrstModule add ad-core ad-props incomp mrst-gui

%% 2. Create 3D matrix grid
rng(42);          % for reproducible permeability field
tol     = 1e-3;
celldim = [40 28 6];
physdim = [600 400 60];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

Nx = celldim(1);  Ny = celldim(2);  Nz = celldim(3);

%% 3. Heterogeneous rock properties - 3D Gaussian random field (log-normal)
% Spatially correlated log-normal permeability via Gaussian filtering of a
% white-noise field.  Parameters:
%   K_mean = 80 mD,  sigma_lnK = 2.5  (high contrast)
%   horizontal correlation length Lxy = 120 m,  vertical Lz = 15 m
K_mean_mD  = 80;      % geometric mean permeability [mD]
sigma_lnK  = 2.5;     % std-dev of ln(K)
corrLen_xy = 120;     % horizontal correlation length [m]
corrLen_z  = 15;      % vertical correlation length [m]

dx = physdim(1)/Nx;  dy = physdim(2)/Ny;  dz = physdim(3)/Nz;

% Generate 3D correlated field: white noise → Gaussian smooth → log-normal
Z  = randn(Nx, Ny, Nz);
gx = gauss_kernel(corrLen_xy / dx);
gy = gauss_kernel(corrLen_xy / dy);
gz = gauss_kernel(corrLen_z  / dz);
Z  = convn(Z, reshape(gx, [], 1, 1), 'same');
Z  = convn(Z, reshape(gy, 1, [], 1), 'same');
Z  = convn(Z, reshape(gz, 1, 1, []), 'same');
Z  = (Z - mean(Z(:))) / std(Z(:));        % standardise to N(0,1)
perm_matrix = exp(log(K_mean_mD * milli*darcy) + sigma_lnK * Z(:));

G.rock = makeShaleRock(G, perm_matrix, 0.20);

%% 4. Define fracture planes
fracplanes = struct;

% Frac 1 - BARRIER 1: near-injector (x ≈ 120-155 m, ~45 m from injector).
% pEDFM cuts M-M transmissibility immediately downstream of injection.
fracplanes(1).points   = [120 10 3; 120 10 57; 155 390 57; 155 390 3];
fracplanes(1).aperture = 1/25;
fracplanes(1).poro     = 0.2;
fracplanes(1).perm     = 1*nano*darcy;    % near-impermeable barrier

% Frac 2 - CONDUIT 1: main diagonal SW→NE bypass (high perm).
fracplanes(2).points   = [100 30 3; 100 30 57; 490 370 57; 490 370 3];
fracplanes(2).aperture = 1/25;
fracplanes(2).poro     = 0.5;
fracplanes(2).perm     = 5000*darcy;      % high-perm fast channel

% Frac 3 - CONDUIT 2: secondary NW→SE diagonal bypass (moderate-high perm).
fracplanes(3).points   = [80 360 3; 80 360 57; 430 40 57; 430 40 3];
fracplanes(3).aperture = 1/25;
fracplanes(3).poro     = 0.4;
fracplanes(3).perm     = 2000*darcy;      % secondary fast channel

% Frac 4 - BARRIER 2: near-producer, extended to touch the north boundary
% (y=400 m).  This seals off matrix flow from the north side, forcing all
% flow to pass either through or around the barrier's southern tip.
fracplanes(4).points   = [455 10 3; 455 10 57; 490 400 57; 490 400 3];
fracplanes(4).aperture = 1/25;
fracplanes(4).poro     = 0.2;
fracplanes(4).perm     = 1*nano*darcy;    % near-impermeable barrier

fprintf('Fracture setup:\n');
fprintf('  Frac 1 (BARRIER 1, near inj): perm = %.1g darcy\n',  fracplanes(1).perm/darcy);
fprintf('  Frac 2 (CONDUIT 1):           perm = %g darcy\n',    fracplanes(2).perm/darcy);
fprintf('  Frac 3 (CONDUIT 2):           perm = %g darcy\n',    fracplanes(3).perm/darcy);
fprintf('  Frac 4 (BARRIER 2, near prd): perm = %.1g darcy\n',  fracplanes(4).perm/darcy);

%% 5. EDFM preprocessing (shale pipeline)
fprintf('\nRunning EDFM preprocessing...\n');
[G, fracplanes] = EDFMshalegrid(G, fracplanes, ...
    'Tolerance', tol, 'plotgrid', false, ...
    'fracturelist', 1:numel(fracplanes));

G = fracturematrixShaleNNC3D(G, tol);
[G, fracplanes] = fracturefractureShaleNNCs3D(G, fracplanes, tol);

nMatrix = G.Matrix.cells.num;
nFrac   = G.cells.num - nMatrix;
nFrac1  = G.FracGrid.Frac1.cells.num;
nFrac2  = G.FracGrid.Frac2.cells.num;
nFrac3  = G.FracGrid.Frac3.cells.num;
nFrac4  = G.FracGrid.Frac4.cells.num;
fprintf('EDFM grid: %d cells (%d matrix + %d fracture)\n', G.cells.num, nMatrix, nFrac);
fprintf('  Frac1=%d  Frac2=%d  Frac3=%d  Frac4=%d cells\n', nFrac1, nFrac2, nFrac3, nFrac4);
fprintf('EDFM NNCs: %d\n', size(G.nnc.cells, 1));

G_edfm = G;   % save EDFM state

%% 6. pEDFM projected M-M connections
fprintf('\nComputing pEDFM projected NNCs (pMatFracNNCs3D)...\n');
G_pedfm = pMatFracNNCs3D(G, tol);
fprintf('pEDFM NNCs: %d  (%d added by pMatFracNNCs3D)\n', ...
    size(G_pedfm.nnc.cells, 1), ...
    size(G_pedfm.nnc.cells, 1) - size(G_edfm.nnc.cells, 1));

%% 7. Fluid properties (single-phase water, Spivey EOS, dynamic lambda_F)
pRef  = 150*barsa;
K0    = 273.15;
T0    = K0 + 120;    % initial reservoir temperature [K]
T_inj = K0 + 20;     % injection temperature [K]
qInj  = 5e-3;        % injection volumetric rate [m3/s]
pProd = pRef - 15*barsa;  % producer BHP [Pa]

% IAPWS 2008 lambda_F function handle - activates EDFMDynamicHeatTransmissibility
lambdaF_func = @(p, T) max(0.4, ...
    0.5563 + 2.31e-3.*(T - 273.15) - 8.7e-6.*(T - 273.15).^2);

fluid = initSimpleADIFluid('phases', 'W', ...
    'n'   , 1      , ...
    'mu'  , 1e-3   , ...
    'rho' , 1000   , ...
    'pRef', pRef   , ...
    'c'   , 4.5e-10, ...
    'cR'  , 0);
% addThermalFluidProps enforces double type for lambdaF via merge_options,
% so pass a scalar and then override with the function handle.
lambdaF_ref = lambdaF_func(pRef, T0);   % reference value at reservoir T
fluid = addThermalFluidProps(fluid, ...
    'Cp'     , 4200        , ...
    'lambdaF', lambdaF_ref , ...
    'useEOS' , true);
fluid.lambdaF = lambdaF_func;   % override: activates EDFMDynamicHeatTransmissibility

%% 8. Wells (vertical doublet, perforating all 6 z-layers)
xc = G_edfm.cells.centroids(:, 1);
yc = G_edfm.cells.centroids(:, 2);
mat_mask = (1:G_edfm.cells.num)' <= nMatrix;

[~, ref_inj] = min((xc - 75).^2  + (yc - 200).^2 + 1e12*(~mat_mask));
[~, ref_prd] = min((xc - 525).^2 + (yc - 200).^2 + 1e12*(~mat_mask));
x_inj = xc(ref_inj);  y_inj = yc(ref_inj);
x_prd = xc(ref_prd);  y_prd = yc(ref_prd);

inj_cells  = find(abs(xc - x_inj) < 1 & abs(yc - y_inj) < 1 & mat_mask);
prod_cells = find(abs(xc - x_prd) < 1 & abs(yc - y_prd) < 1 & mat_mask);
inj_cells  = sort(inj_cells);
prod_cells = sort(prod_cells);

fprintf('\nWells:\n');
fprintf('  Injector:  %d perforations  (x=%.0f m, y=%.0f m)\n', numel(inj_cells),  x_inj, y_inj);
fprintf('  Producer:  %d perforations  (x=%.0f m, y=%.0f m)\n', numel(prod_cells), x_prd, y_prd);

W = addWell([], G_edfm, G_edfm.rock, inj_cells, ...
    'type'  , 'rate', 'val', qInj,  'radius', 0.15, 'comp_i', 1, 'name', 'INJ');
W = addWell(W,  G_edfm, G_edfm.rock, prod_cells, ...
    'type'  , 'bhp',  'val', pProd, 'radius', 0.15, 'comp_i', 1, 'name', 'PRD');

lambdaR  = 2.5;
rock_tmp = G_edfm.rock;
rock_tmp = addThermalRockProps(rock_tmp, 'lambdaR', lambdaR, 'rhoR', 2650, 'CpR', 900);
% addThermalWellProps does fluid.lambdaF.*rock.poro - needs scalar lambdaF
fluid_well        = fluid;
fluid_well.lambdaF = lambdaF_ref;
W = addThermalWellProps(W, G_edfm, rock_tmp, fluid_well, 'T', T_inj);

%% 9. Timestep schedule (20-year simulation)
dt = rampupTimesteps(20*year, 120*day);
fprintf('\nSchedule: %d timesteps, total = %.0f years\n', numel(dt), sum(dt)/year);

%% 10. EDFM simulation
fprintf('\n--- EDFM Simulation ---\n');
rock_edfm = G_edfm.rock;
rock_edfm = addThermalRockProps(rock_edfm, 'lambdaR', lambdaR, 'rhoR', 2650, 'CpR', 900);
G_edfm.rock = rock_edfm;
G_edfm = computeThermalNNCTransFracMatrix(G_edfm, rock_edfm, fluid, 'pRef', pRef, 'TRef', T0);

gravity reset off;
model_edfm = GeothermalHFMModel(G_edfm, rock_edfm, fluid, 'fractureMethod', 'edfm');
model_edfm.extraStateOutput   = true;
model_edfm.outputFluxes       = true;
model_edfm.minimumTemperature = K0;
model_edfm.maximumTemperature = K0 + 400;
model_edfm.maximumPressure    = 600e6;
model_edfm = model_edfm.validateModel();

state0_edfm   = initResSol(G_edfm, pRef, 1);
state0_edfm.T = T0 * ones(G_edfm.cells.num, 1);
schedule_edfm = simpleSchedule(dt, 'W', W);

% Row-scaled direct solver: left diagonal scaling equilibrates the badly-scaled
% coupled flow+heat Jacobian (extreme 5000 D / 1 nD permeability contrast) so the
% backslash is not flagged as ill-conditioned. The scaling is exact, so the
% converged solution is unchanged.
lsolver = BackslashSolverAD('applyLeftDiagonalScaling', true);
t0 = tic;
[wellSols_edfm, states_edfm] = simulateScheduleAD(state0_edfm, model_edfm, schedule_edfm, ...
    'LinearSolver', lsolver);
fprintf('EDFM done (%.1f s).\n', toc(t0));

%% 11. pEDFM simulation
fprintf('\n--- pEDFM Simulation ---\n');
rock_pedfm = G_pedfm.rock;
rock_pedfm = addThermalRockProps(rock_pedfm, 'lambdaR', lambdaR, 'rhoR', 2650, 'CpR', 900);
G_pedfm.rock = rock_pedfm;
G_pedfm = computeThermalNNCTransFracMatrix(G_pedfm, rock_pedfm, fluid, 'pRef', pRef, 'TRef', T0);

model_pedfm = GeothermalHFMModel(G_pedfm, rock_pedfm, fluid, 'fractureMethod', 'pedfm');
model_pedfm.extraStateOutput   = true;
model_pedfm.outputFluxes       = true;
model_pedfm.minimumTemperature = K0;
model_pedfm.maximumTemperature = K0 + 400;
model_pedfm.maximumPressure    = 600e6;
model_pedfm = model_pedfm.validateModel();

state0_pedfm   = initResSol(G_pedfm, pRef, 1);
state0_pedfm.T = T0 * ones(G_pedfm.cells.num, 1);
schedule_pedfm = simpleSchedule(dt, 'W', W);

t0 = tic;
[wellSols_pedfm, states_pedfm] = simulateScheduleAD(state0_pedfm, model_pedfm, schedule_pedfm, ...
    'LinearSolver', lsolver);
fprintf('pEDFM done (%.1f s).\n', toc(t0));

%% ========== Visualization ==========
times_years = cumsum(dt) / year;
nM          = nMatrix;
dT_thresh   = 5;    % cooling threshold [degC] for "activated" cell
nSteps      = numel(states_pedfm);

% Output directory for report-ready PDF figures.
% Resolve from the module path (robust regardless of current directory).
modDir = mrstPath('thermal-hfm');   % .../mrst-2026a/solvers/thermal-hfm
figDir = fullfile(modDir, '..', '..', '..', 'report', 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
fprintf('\nSaving figures to: %s\n', figDir);

% Fracture cell index ranges (relative to global cell numbering)
frac1_cells = nM + (1:nFrac1);
frac2_cells = nM + nFrac1 + (1:nFrac2);
frac3_cells = nM + nFrac1 + nFrac2 + (1:nFrac3);
frac4_cells = nM + nFrac1 + nFrac2 + nFrac3 + (1:nFrac4);

%% Figure 1: Setup - Gaussian perm slice + 3D fracture/well overview
fig1 = figure('Name', 'Doublet Setup', 'Position', [50, 50, 1150, 500], 'Color', 'w', 'Visible', 'off');

subplot(1, 2, 1);
kz_mid     = ceil(Nz / 2);
midlay_idx = (kz_mid-1)*Nx*Ny + (1:Nx*Ny);
logK = log10(G_edfm.rock.perm(1:nM, 1) / darcy);
plotCellData(G_edfm, logK(midlay_idx), midlay_idx, 'EdgeAlpha', 0.15);
colormap(gca, parula(32));
cbar = colorbar; title(cbar, 'log_{10}(K [D])');
view(0, 90); axis equal tight;
title(sprintf('Mid-Layer Permeability  (K_{mean}=80 mD, \\sigma_{lnK}=2.5, L_{xy}=120 m)'));
xlabel('x [m]'); ylabel('y [m]');
hold on;
rectangle('Position', [x_inj-15, y_inj-15, 30, 30], ...
    'EdgeColor', 'b', 'LineWidth', 2.5, 'Curvature', [1,1]);
rectangle('Position', [x_prd-15, y_prd-15, 30, 30], ...
    'EdgeColor', 'r', 'LineWidth', 2.5, 'Curvature', [1,1]);
hold off;
set(gca, 'FontSize', 11);

subplot(1, 2, 2);
h_m  = plotGrid(G_edfm, 1:nM, 'FaceColor', [0.88 0.93 1], 'FaceAlpha', 0.07, 'EdgeAlpha', 0.03);
hold on;
h_b1 = plotCellData(G_edfm, zeros(nFrac1,1), frac1_cells, ...
    'FaceColor', [0.1 0.2 0.9], 'EdgeAlpha', 0.3, 'FaceAlpha', 0.75);   % barrier 1 blue
h_c1 = plotCellData(G_edfm, zeros(nFrac2,1), frac2_cells, ...
    'FaceColor', [0.9 0.1 0.1], 'EdgeAlpha', 0.3, 'FaceAlpha', 0.75);   % conduit 1 red
h_c2 = plotCellData(G_edfm, zeros(nFrac3,1), frac3_cells, ...
    'FaceColor', [0.95 0.55 0.0], 'EdgeAlpha', 0.3, 'FaceAlpha', 0.75); % conduit 2 orange
h_b2 = plotCellData(G_edfm, zeros(nFrac4,1), frac4_cells, ...
    'FaceColor', [0.0 0.80 0.80], 'EdgeAlpha', 0.3, 'FaceAlpha', 0.75); % barrier 2 cyan
zc  = G_edfm.cells.centroids(:, 3);
inj_top = inj_cells(G_edfm.cells.centroids(inj_cells,3) == max(zc(inj_cells)));
prd_top = prod_cells(G_edfm.cells.centroids(prod_cells,3) == max(zc(prod_cells)));
h_inj = plot3(G_edfm.cells.centroids(inj_top,1), G_edfm.cells.centroids(inj_top,2), physdim(3)+3, ...
    'b^', 'MarkerSize', 14, 'MarkerFaceColor', 'b');
h_prd = plot3(G_edfm.cells.centroids(prd_top,1), G_edfm.cells.centroids(prd_top,2), physdim(3)+3, ...
    'rv', 'MarkerSize', 14, 'MarkerFaceColor', 'r');
hold off;
legend([h_m(1), h_b1(1), h_c1(1), h_c2(1), h_b2(1), h_inj, h_prd], ...
    {'Matrix','Barrier 1 (near Inj)','Conduit 1 (5000 D)', ...
    'Conduit 2 (2000 D)','Barrier 2 (near Prd)','Injector','Producer'}, ...
    'Location','northeast','FontSize',8);
view([-35, 22]); axis equal tight;
title('3D Fracture Network + Well Locations  (barriers near wells)');
set(gca, 'FontSize', 10);

exportgraphics(fig1, fullfile(figDir, 'doublet_setup.png'), 'Resolution', 300);
fprintf('  Saved doublet_setup.png\n');

%% Figure 2: Cold-plume snapshots (pEDFM) - 4 selected timesteps
steps_show = unique(round(linspace(max(1, round(nSteps*0.15)), nSteps, 4)));
steps_show = steps_show(1:4);

fig2 = figure('Name', 'Cold-Plume Snapshots (pEDFM)', 'Position', [50, 120, 1200, 620], 'Color', 'w', 'Visible', 'off');
cmap_cool = flipud(cool(32));
for sp = 1:4
    st   = steps_show(sp);
    T_st = states_pedfm{st}.T(1:nM) - K0;
    dT   = (T0 - K0) - T_st;
    active = find(dT > dT_thresh);

    subplot(2, 2, sp);
    % Background grid in light grey so domain outline is always visible
    plotGrid(G_pedfm, 1:nM, 'FaceColor', [0.94 0.94 0.94], ...
        'EdgeColor', [0.75 0.75 0.75], 'EdgeAlpha', 0.35);
    hold on;
    if ~isempty(active)
        plotCellData(G_pedfm, dT(active), active, 'EdgeAlpha', 0.04);
        clim([dT_thresh, max(dT_thresh + 0.1, max(dT(active)))]);
        colormap(gca, cmap_cool); colorbar;
    else
        text(physdim(1)/2, physdim(2)/2, 30, 'No activated cells', ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    plot3(x_inj, y_inj, physdim(3)/2, 'b^', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    plot3(x_prd, y_prd, physdim(3)/2, 'rv', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    hold off;
    view(0, 90); axis equal tight;
    title(sprintf('t = %.1f yr  |  %d cells  (\\DeltaT > %d°C)', ...
        times_years(st), numel(active), dT_thresh));
    xlabel('x [m]'); ylabel('y [m]');
    set(gca, 'FontSize', 10);
end
sgtitle(sprintf('Cold-Water Plume - pEDFM  (\\DeltaT > %d°C threshold)', dT_thresh), ...
    'FontSize', 13, 'FontWeight', 'bold');

exportgraphics(fig2, fullfile(figDir, 'doublet_cold_plume.png'), 'Resolution', 300);
fprintf('  Saved doublet_cold_plume.png\n');

%% Figure 3: Final temperature - EDFM vs pEDFM side-by-side (top view)
T_edfm_fin  = states_edfm{end}.T(1:nM)  - K0;
T_pedfm_fin = states_pedfm{end}.T(1:nM) - K0;
Tmin_plot = floor(min([T_edfm_fin; T_pedfm_fin]));
Tmax_plot = T0 - K0;

fig3 = figure('Name', 'Final Temperature: EDFM vs pEDFM', 'Position', [50, 50, 1300, 560], 'Color', 'w', 'Visible', 'off');
subplot(1, 2, 1);
plotCellData(G_edfm, T_edfm_fin, 1:nM, 'EdgeAlpha', 0.06);
colormap(gca, jet(32)); clim([Tmin_plot, Tmax_plot]);
cbar = colorbar; title(cbar, 'T [°C]');
view(0, 90); axis equal tight;
title(sprintf('EDFM - Final Temperature  (t = %.0f yr)', times_years(end)));
xlabel('x [m]'); ylabel('y [m]'); set(gca, 'FontSize', 11);

subplot(1, 2, 2);
plotCellData(G_pedfm, T_pedfm_fin, 1:nM, 'EdgeAlpha', 0.06);
colormap(gca, jet(32)); clim([Tmin_plot, Tmax_plot]);
cbar = colorbar; title(cbar, 'T [°C]');
view(0, 90); axis equal tight;
title(sprintf('pEDFM - Final Temperature  (t = %.0f yr)', times_years(end)));
xlabel('x [m]'); ylabel('y [m]'); set(gca, 'FontSize', 11);

exportgraphics(fig3, fullfile(figDir, 'doublet_final_temperature.png'), 'Resolution', 300);
fprintf('  Saved doublet_final_temperature.png\n');

%% Figure 4: Thermal analysis - producer + near-barrier monitoring zone
prod_idx     = 2;
T_prod_edfm  = zeros(nSteps, 1);
T_prod_pedfm = zeros(nSteps, 1);
for k = 1:nSteps
    if ~isempty(wellSols_edfm{k}),  T_prod_edfm(k)  = wellSols_edfm{k}(prod_idx).T  - K0; end
    if ~isempty(wellSols_pedfm{k}), T_prod_pedfm(k) = wellSols_pedfm{k}(prod_idx).T - K0; end
end

% Define monitoring zone: matrix cells just downstream of Barrier 1
% Barrier 1 runs x ≈ 120-155 m;  monitor the zone x ∈ [160, 250] m
% which is immediately downstream of the barrier on the injector side.
cc_mat = G_edfm.cells.centroids(1:nM, :);
monitor_mask = cc_mat(:,1) >= 160 & cc_mat(:,1) <= 250;
n_monitor = sum(monitor_mask);

T_zone_edfm  = zeros(nSteps, 1);
T_zone_pedfm = zeros(nSteps, 1);
for k = 1:nSteps
    T_zone_edfm(k)  = mean(states_edfm{k}.T(monitor_mask))  - K0;
    T_zone_pedfm(k) = mean(states_pedfm{k}.T(monitor_mask)) - K0;
end

fig4 = figure('Name', 'Thermal Analysis', 'Position', [100, 100, 1200, 460], 'Color', 'w', 'Visible', 'off');

% Panel (a): Producer temperature breakthrough
subplot(1, 2, 1);
plot(times_years, T_prod_edfm,  'b-',  'LineWidth', 2.5, 'DisplayName', 'EDFM');
hold on;
plot(times_years, T_prod_pedfm, 'r--', 'LineWidth', 2.5, 'DisplayName', 'pEDFM');
yline(T0   - K0, 'k:', 'T_{res} = 120°C', 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.2);
yline(T_inj - K0, 'c:', 'T_{inj} = 20°C',  'LabelHorizontalAlignment', 'left', 'LineWidth', 1.2);
hold off;
xlabel('Time [years]'); ylabel('Producer temperature [°C]');
title('(a) Thermal Breakthrough at Producer');
legend('Location', 'southwest', 'FontSize', 10);
grid on; ylim([max(0, Tmin_plot-5), Tmax_plot+5]); set(gca, 'FontSize', 11);

% Panel (b): Mean temperature in monitoring zone downstream of Barrier 1
subplot(1, 2, 2);
plot(times_years, T_zone_edfm,  'b-',  'LineWidth', 2.5, 'DisplayName', 'EDFM');
hold on;
plot(times_years, T_zone_pedfm, 'r--', 'LineWidth', 2.5, 'DisplayName', 'pEDFM');
yline(T0 - K0, 'k:', 'T_{res}', 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.2);
hold off;
xlabel('Time [years]'); ylabel('Mean temperature [°C]');
title(sprintf('(b) Monitoring Zone Downstream of Barrier 1\n(x \\in [160, 250] m,  %d cells)', n_monitor));
legend('Location', 'southwest', 'FontSize', 10);
grid on; ylim([max(0, Tmin_plot-5), Tmax_plot+5]); set(gca, 'FontSize', 11);

exportgraphics(fig4, fullfile(figDir, 'doublet_breakthrough.png'), 'Resolution', 300);
fprintf('  Saved doublet_breakthrough.png\n');

%% Figure 5: Temperature difference pEDFM - EDFM (final state)
dT_models  = T_pedfm_fin - T_edfm_fin;
maxAbsDiff = max(abs(dT_models));

fig5 = figure('Name', 'Temp Difference: pEDFM - EDFM', 'Position', [100, 100, 820, 500], 'Color', 'w', 'Visible', 'off');
plotCellData(G_pedfm, dT_models, 1:nM, 'EdgeAlpha', 0.06);
t_cmap   = linspace(0, 1, 128)';
cmap_div = [bsxfun(@times, [0,0,1], (1-t_cmap)) + bsxfun(@times, ones(128,3), t_cmap); ...
            bsxfun(@times, ones(128,3), (1-t_cmap)) + bsxfun(@times, [1,0,0], t_cmap)];
colormap(gca, cmap_div); clim([-maxAbsDiff, maxAbsDiff]);
cbar = colorbar; title(cbar, '\DeltaT [°C]');
view(0, 90); axis equal tight;
title({'Temperature Difference: pEDFM - EDFM  (final state)', ...
       'Blue = pEDFM cooler (barriers block more),  Red = pEDFM warmer'});
xlabel('x [m]'); ylabel('y [m]'); set(gca, 'FontSize', 12);

exportgraphics(fig5, fullfile(figDir, 'doublet_temp_difference.png'), 'Resolution', 300);
fprintf('  Saved doublet_temp_difference.png\n');

%% Figure 6: Final pressure - EDFM vs pEDFM side-by-side (top view)
p_edfm_fin  = states_edfm{end}.pressure(1:nM)  / barsa;
p_pedfm_fin = states_pedfm{end}.pressure(1:nM) / barsa;
pmin_plot = floor(min([p_edfm_fin; p_pedfm_fin]));
pmax_plot = ceil(max([p_edfm_fin; p_pedfm_fin]));

fig6 = figure('Name', 'Final Pressure: EDFM vs pEDFM', 'Position', [50, 50, 1300, 560], 'Color', 'w', 'Visible', 'off');
subplot(1, 2, 1);
plotCellData(G_edfm, p_edfm_fin, 1:nM, 'EdgeAlpha', 0.06);
colormap(gca, parula(32)); clim([pmin_plot, pmax_plot]);
cbar = colorbar; title(cbar, 'p [bar]');
view(0, 90); axis equal tight;
title(sprintf('EDFM - Final Pressure  (t = %.0f yr)', times_years(end)));
xlabel('x [m]'); ylabel('y [m]'); set(gca, 'FontSize', 11);

subplot(1, 2, 2);
plotCellData(G_pedfm, p_pedfm_fin, 1:nM, 'EdgeAlpha', 0.06);
colormap(gca, parula(32)); clim([pmin_plot, pmax_plot]);
cbar = colorbar; title(cbar, 'p [bar]');
view(0, 90); axis equal tight;
title(sprintf('pEDFM - Final Pressure  (t = %.0f yr)', times_years(end)));
xlabel('x [m]'); ylabel('y [m]'); set(gca, 'FontSize', 11);

exportgraphics(fig6, fullfile(figDir, 'doublet_final_pressure.png'), 'Resolution', 300);
fprintf('  Saved doublet_final_pressure.png\n');

%% Figure 7: Interactive cold-plume viewer - MRST plotToolbar (not exported)
% plotToolbar (from mrst-gui) creates a figure with:
%   • A field-selector dropdown to choose between dT_cooling / T_celsius / pressure_bar
%   • A time-step slider at the bottom for interactive scrubbing
%
% Pre-build display states containing the fields we want to show.
fprintf('\nBuilding interactive plotToolbar viewer (%d states)...\n', nSteps);

states_tb = cell(nSteps, 1);
for k = 1:nSteps
    T_k = states_pedfm{k}.T;     % [K], all cells (matrix + fracture)
    states_tb{k} = struct( ...
        'dT_cooling_C',  max(0, T0 - T_k), ...          % cooling [K = degC]
        'T_celsius',     T_k - K0,         ...          % temperature [degC]
        'pressure_bar',  states_pedfm{k}.pressure / barsa);
end

fig7 = figure('Name', 'Interactive: pEDFM Cold Plume  [plotToolbar]', ...
    'Position', [100, 50, 1060, 720], 'Color', 'w', 'Visible', 'off');
plotToolbar(G_pedfm, states_tb, 'EdgeAlpha', 0.04);
colormap(flipud(cool(32)));
view(0, 90);  axis equal tight;
xlabel('x [m]');  ylabel('y [m]');
title({'pEDFM Interactive Viewer  (MRST plotToolbar)', ...
    'Dropdown: select field  |  Slider (bottom): step through time'}, ...
    'FontSize', 11);

% Export just the plot axes (not the whole figure) so the plotToolbar UI
% widgets are not in scope; this avoids the exporter's "UI components" warning.
ax_tb = findobj(fig7, 'type', 'axes');
exportgraphics(ax_tb(1), fullfile(figDir, 'doublet_plottoolbar.png'), 'Resolution', 300);
fprintf('  Saved doublet_plottoolbar.png\n');

%% Summary
fprintf('\n=== Thermal Doublet Results Summary ===\n');
fprintf('Grid: %d x %d x %d  (%d matrix cells)\n', Nx, Ny, Nz, nMatrix);
fprintf('Perm: log-normal Gaussian field  (K_mean = %d mD,  sigma_lnK = %.1f)\n', ...
    K_mean_mD, sigma_lnK);
fprintf('Fractures: 4  (2 barriers near wells + 2 conduits)\n');
fprintf('EDFM  producer T at t = %.0f yr: %.1f degC\n', times_years(end), T_prod_edfm(end));
fprintf('pEDFM producer T at t = %.0f yr: %.1f degC\n', times_years(end), T_prod_pedfm(end));
fprintf('Temperature difference (EDFM - pEDFM): %.1f degC\n', ...
    T_prod_edfm(end) - T_prod_pedfm(end));
fprintf('Max matrix dT (pEDFM-EDFM) = %.2f degC  at final state\n', maxAbsDiff);
fprintf('\nPhysical interpretation:\n');
fprintf('  Near-well barriers (Frac 1 at inj, Frac 4 at prd) force cold water\n');
fprintf('  through the diagonal conduit fractures.  pEDFM correctly cuts M-M\n');
fprintf('  transmissibility across the barriers; EDFM over-connects the matrix,\n');
fprintf('  causing earlier and more uniform thermal breakthrough.\n');

%% -----------------------------------------------------------------------
%  Local helper function
%% -----------------------------------------------------------------------
function k = gauss_kernel(sigma_cells)
%GAUSS_KERNEL  Normalised 1-D Gaussian smoothing kernel.
%   k = gauss_kernel(sigma_cells) returns a column vector with Gaussian
%   weights evaluated on an integer grid with std-dev sigma_cells (in
%   grid-cell units).  k sums to 1.  Used for separable 3-D smoothing of
%   random fields via convn.
    if sigma_cells < 0.5
        k = 1;   % correlation length smaller than one cell - no smoothing
        return;
    end
    n = ceil(3.5 * sigma_cells);
    x = (-n : n)';
    k = exp(-0.5 * (x / sigma_cells).^2);
    k = k / sum(k);
end
