%% Example: Thermal pEDFM with Brine EOS and NaCl Transport
% Demonstrates the thermal-hfm module using pEDFM (projection-based EDFM)
% with brine physics:
%
%   1. Two fractures: a low-perm BARRIER and a high-perm CONDUIT.
%      pEDFM correctly represents the barrier (blocks flow/heat transport)
%      by adding projected M-M connections that reduce transmissibility near
%      the barrier fracture.  Standard EDFM cannot capture this effect.
%   2. CompositionalBrineFluid (H2O + NaCl) - advective salinity transport.
%   3. useEOS = true  - Spivey (2004) T-p-dependent density and viscosity.
%   4. lambdaF_func   - @(p,T) function handle passed to addThermalFluidProps,
%      activating EDFMDynamicHeatTransmissibility for dynamic NNC thermal
%      transmissibilities (full AD derivatives) at every Newton iteration.
%   5. plotToolbar with a cell array of states for interactive time-stepping.
%
% pEDFM preprocessing pipeline:
%   EDFMshalegrid → fracturematrixShaleNNC3D → fracturefractureShaleNNCs3D
%   → pMatFracNNCs3D → computeThermalNNCTransFracMatrix
%
% Setup:
%   - 3D domain (120 x 80 x 40 m).
%   - Fracture 1 (BARRIER):  diagonal, 1e-6 darcy - blocks west-to-east flow.
%   - Fracture 2 (CONDUIT):  diagonal, 5000 darcy - accelerates transport.
%   - Hot saline brine (150 degC, 5 wt% NaCl) injected at west face.
%   - Production at east face (pressure BC).

clc; clear; close all;

%% 1. Load modules
mrstModule add hfm shale geothermal thermal-hfm
mrstModule add ad-core ad-props incomp mrst-gui

%% 2. Create 3D matrix grid
tol     = 1e-3;
celldim = [12 8 4];
physdim = [120 80 40];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);
G.rock = makeShaleRock(G, 50*milli*darcy, 0.15);

%% 3. Define fracture planes
fracplanes = struct;

% Fracture 1: BARRIER (low perm) - blocks direct matrix flow paths
fracplanes(1).points   = [35 5 2; 35 5 38; 35 75 38; 35 75 2];
fracplanes(1).aperture = 1/25;
fracplanes(1).poro     = 0.2;
fracplanes(1).perm     = 1e-6*darcy;   % near-impermeable barrier

% Fracture 2: CONDUIT (high perm) - fast transport channel
fracplanes(2).points   = [40 10 4; 40 10 36; 110 70 36; 110 70 4];
fracplanes(2).aperture = 1/25;
fracplanes(2).poro     = 0.5;
fracplanes(2).perm     = 5000*darcy;   % high-perm conduit

fprintf('Fracture 1 (barrier): perm = %.2g darcy\n', fracplanes(1).perm/darcy);
fprintf('Fracture 2 (conduit): perm = %.0f darcy\n', fracplanes(2).perm/darcy);

%% 4. pEDFM preprocessing pipeline
fprintf('\nProcessing fractures with EDFMshalegrid...\n');
[G, fracplanes] = EDFMshalegrid(G, fracplanes, ...
    'Tolerance', tol, 'plotgrid', false, ...
    'fracturelist', 1:numel(fracplanes));

% Fracture-matrix NNCs
G = fracturematrixShaleNNC3D(G, tol);

% Fracture-fracture NNCs (needed when multiple fractures may intersect)
[G, fracplanes] = fracturefractureShaleNNCs3D(G, fracplanes, tol);

nMatrix = G.Matrix.cells.num;
nFrac   = G.cells.num - nMatrix;
fprintf('EDFM grid: %d cells (%d matrix + %d fracture)\n', ...
    G.cells.num, nMatrix, nFrac);

% pEDFM projected M-M connections (key step for flow barrier representation)
fprintf('Computing pEDFM projected NNCs...\n');
G = pMatFracNNCs3D(G, tol);
nNNC = size(G.nnc.cells, 1);
fprintf('Total NNCs after pMatFracNNCs3D: %d\n', nNNC);

%% 5. Rock thermal properties
G.rock = addThermalRockProps(G.rock, ...
    'lambdaR', 2.5, 'rhoR', 2650, 'CpR', 900);
rock = G.rock;

%% 6. Fluid: brine with Spivey EOS + dynamic lambdaF handle
pRef = 200*barsa;

% T-dependent lambdaF (IAPWS 2008 polynomial, valid 20-200 degC).
% Passed as a function handle to activate EDFMDynamicHeatTransmissibility.
lambdaF_func   = @(p, T) max(0.4, ...
    0.5563 + 2.31e-3.*(T - 273.15) - 8.7e-6.*(T - 273.15).^2);
lambdaF_scalar = lambdaF_func(pRef, 273.15 + 115);   % ~0.707 W/(m·K), for Figure 4

% Base fluid
fluid = initSimpleADIFluid('phases', 'W', ...
    'n'   , 1      , ...
    'mu'  , 1e-3   , ...
    'rho' , 1000   , ...
    'pRef', pRef   , ...
    'c'   , 4.5e-10, ...
    'cR'  , 0);

% Spivey (2004) EOS.  Pass scalar lambdaF to satisfy merge_options type check,
% then override with function handle to activate EDFMDynamicHeatTransmissibility.
fluid = addThermalFluidProps(fluid, ...
    'Cp'     , 4200          , ...
    'lambdaF', lambdaF_scalar, ...
    'useEOS' , true);
fluid.lambdaF = lambdaF_func;   % override: activate dynamic NNC heat transmissibility

% Compositional brine: diffusivities = 0 (advection-only demonstration).
% Non-zero diffusivities are now supported via
% EDFMDynamicMolecularTransmissibility (see test_step7_brine_diffusion).
compFluid = CompositionalBrineFluid( ...
    {'H2O'              , 'NaCl'            }, ...
    [18.015281*gram/mol , 58.4428*gram/mol  ], ...
    [0                  , 0                 ]);

%% 7. Compute static NNC thermal transmissibilities
G = computeThermalNNCTransFracMatrix(G, rock, fluid, ...
    'pRef', pRef, 'TRef', 273.15 + 20);

%% 8. Create pEDFM model
gravity reset off;
model = GeothermalHFMModel(G, rock, fluid, compFluid, 'fractureMethod', 'pedfm');
model.extraStateOutput = true;
model.outputFluxes     = true;

K0 = 273.15;
model.minimumTemperature = K0;
model.maximumTemperature = K0 + 300;
model.maximumPressure    = 500e6;
% validateModel populates model.Components (needed by getMassFraction)
model = model.validateModel();

%% 9. Initial state: hot saline formation water
T0      = K0 + 80;     % 80 degC initial temperature
X0_NaCl = 0.15;        % 15 wt% NaCl

state0   = initResSol(G, pRef, 1);
state0.T = T0 * ones(G.cells.num, 1);
X0_mass  = repmat([1 - X0_NaCl, X0_NaCl], G.cells.num, 1);
state0.components = model.getMoleFraction(X0_mass);

%% 10. Boundary conditions
% All BC faces must receive non-nan component values (see example_edfm_brine_eos.m
% for detailed note). Outlet faces receive the initial reservoir composition.
faces = boundaryFaces(G);
fc    = G.faces.centroids(faces, :);
west  = abs(fc(:, 1) - 0)           < 1e-6;
east  = abs(fc(:, 1) - physdim(1))  < 1e-6;
wf    = faces(west);
ef    = faces(east);

T_inj      = K0 + 150;   % 150 degC
X_inj_NaCl = 0.05;       % 5 wt% NaCl

bc = addBC([], wf, 'flux'    , 3e-4);
bc = addBC(bc, ef, 'pressure', pRef);
Tbc   = [repmat(T_inj, numel(wf), 1); nan(numel(ef), 1)];
Hflux = [nan(numel(wf), 1)          ; zeros(numel(ef), 1)];
bc    = addThermalBCProps(bc, 'T', Tbc, 'Hflux', Hflux);

X_inj_mass = repmat([1 - X_inj_NaCl, X_inj_NaCl], numel(wf), 1);
X_out_mass = repmat([1 - X0_NaCl,    X0_NaCl   ], numel(ef), 1);
bc.components = [model.getMoleFraction(X_inj_mass); ...
                 model.getMoleFraction(X_out_mass)];

%% 11. Schedule: 2-year simulation with ramp-up
dt       = rampupTimesteps(2*year, 30*day);
schedule = simpleSchedule(dt, 'bc', bc);
fprintf('\nSchedule: %d timesteps, total %.0f days\n\n', numel(dt), sum(dt)/day);

%% 12. Run simulation
fprintf('Running brine/EOS pEDFM simulation...\n');
[~, states] = simulateScheduleAD(state0, model, schedule);
fprintf('Simulation complete.\n\n');

%% ========== Visualization ==========
times_days  = cumsum(schedule.step.val) / day;
frac_cells  = (nMatrix+1):G.cells.num;
x_all       = G.cells.centroids(:, 1);
east_mask   = x_all > (physdim(1) - physdim(1)/celldim(1));
east_matrix = east_mask & ((1:G.cells.num)' <= nMatrix);

%% Figure 1: Interactive plotToolbar with time slider
% Slider steps through timesteps; play button animates.
fprintf('Building interactive states for plotToolbar...\n');
states_data = cell(numel(states), 1);
for k = 1:numel(states)
    Xk = model.getMassFraction(states{k}.components);
    states_data{k} = struct( ...
        'Temperature_degC', states{k}.T - K0         , ...
        'Pressure_bar'    , states{k}.pressure / barsa, ...
        'NaCl_wt_frac'    , Xk(:, 2)                 );
end
mrstFigure();
plotToolbar(G, states_data, 'EdgeAlpha', 0.05);
view([-40, 25]); colormap(jet(32)); axis tight equal;
title('Brine/EOS pEDFM - Interactive (plotToolbar, use slider or play button)');
set(gca, 'FontSize', 11);

%% Figure 2: Final temperature - barrier and conduit highlighted
fig2 = mrstFigure();
ax2  = axes(fig2);

% Fracture cell ranges (fracture grids stored in G.FracGrid)
nFrac1 = G.FracGrid.Frac1.cells.num;
nFrac2 = G.FracGrid.Frac2.cells.num;
frac1_cells = (nMatrix+1):(nMatrix+nFrac1);
frac2_cells = (nMatrix+nFrac1+1):(nMatrix+nFrac1+nFrac2);

T_fin = states{end}.T - K0;
plotCellData(G, T_fin, 1:nMatrix, 'EdgeAlpha', 0.03, 'FaceAlpha', 0.5);
hold(ax2, 'on');
plotCellData(G, T_fin, frac1_cells, 'EdgeAlpha', 0.3, 'FaceColor', 'b'); % barrier
plotCellData(G, T_fin, frac2_cells, 'EdgeAlpha', 0.3);                   % conduit
hold(ax2, 'off');
colormap(ax2, hot(32)); clim(ax2, [T0 - K0, T_inj - K0]);
cbar2 = colorbar(ax2); title(cbar2, 'T [°C]');
view(ax2, [-40, 25]); axis(ax2, 'tight'); axis(ax2, 'equal');
title(ax2, sprintf('pEDFM Final Temperature  (t = %.0f d)  - barrier=blue, conduit=red', ...
    times_days(end)));
set(ax2, 'FontSize', 11);

%% Figure 3: Thermal breakthrough + EOS density + NaCl dilution
mean_T_out    = zeros(numel(states), 1);
mean_rho_out  = zeros(numel(states), 1);
mean_NaCl_out = zeros(numel(states), 1);
for j = 1:numel(states)
    if ~isempty(states{j})
        p_out  = mean(states{j}.pressure(east_matrix));
        T_out  = mean(states{j}.T(east_matrix));
        Xm_out = model.getMassFraction(states{j}.components);
        mean_T_out(j)    = T_out - K0;
        mean_rho_out(j)  = double(fluid.rhoW(p_out, T_out));
        mean_NaCl_out(j) = mean(Xm_out(east_matrix, 2));
    end
end

figure('Name', 'pEDFM Breakthrough', 'Position', [100, 100, 820, 380]);
yyaxis left;
plot(times_days, mean_T_out,   'r-',  'LineWidth', 2.5); hold on;
plot(times_days, mean_rho_out, 'b--', 'LineWidth', 2.5);
ylabel('T [°C]  /  \rho [kg/m^3]  (Spivey EOS)');
ylim([75, 1000]);
yyaxis right;
plot(times_days, mean_NaCl_out, 'm-.', 'LineWidth', 2.5);
ylabel('Outlet NaCl mass fraction [-]');
ylim([0.04, 0.16]);
xlabel('Time [days]');
title('Thermal Breakthrough, EOS Density, and Salinity Dilution  (Brine EOS pEDFM)');
legend({'Temperature [°C]','Density [kg/m^3]','NaCl fraction'}, 'Location', 'east');
grid on; set(gca, 'FontSize', 12);

%% Figure 4: T-dependent lambda_F profile
T_range  = linspace(20, 200, 200);
lam_vals = lambdaF_func(pRef, T_range + 273.15);
figure('Name', 'lambdaF(T)', 'Position', [150, 150, 560, 340]);
plot(T_range, lam_vals, 'k-', 'LineWidth', 2);
xline(T0 - K0,    'b--', 'T_0 = 80°C',       'LabelVerticalAlignment', 'bottom');
xline(T_inj - K0, 'r--', 'T_{inj} = 150°C',  'LabelVerticalAlignment', 'bottom');
xline(115,        'g--', 'T_{mean} (used)',   'LabelVerticalAlignment', 'top');
yline(lambdaF_scalar, 'g:', sprintf('%.4f W/(m·K)', lambdaF_scalar));
xlabel('Temperature [°C]');
ylabel('\lambda_F  [W m^{-1} K^{-1}]');
title({'T-Dependent Fluid Conductivity - IAPWS 2008', ...
       'Simulation uses @(p,T) handle (EDFMDynamicHeatTransmissibility)'});
grid on; set(gca, 'FontSize', 12);

%% Summary
fprintf('=== Results Summary (pEDFM Brine EOS) ===\n');
fprintf('Initial outlet T:      %.1f degC  (formation: %.1f degC)\n', ...
    mean_T_out(1),   T0 - K0);
fprintf('Final  outlet T:       %.1f degC  (injected:  %.1f degC)\n', ...
    mean_T_out(end), T_inj - K0);
fprintf('Initial outlet rho:    %.2f kg/m^3\n', mean_rho_out(1));
fprintf('Final  outlet rho:     %.2f kg/m^3  (Spivey EOS)\n', mean_rho_out(end));
fprintf('Initial outlet NaCl:   %.4f wt  (formation: %.4f wt)\n', ...
    mean_NaCl_out(1),   X0_NaCl);
fprintf('Final  outlet NaCl:    %.4f wt  (injected:  %.4f wt)\n', ...
    mean_NaCl_out(end), X_inj_NaCl);
fprintf('lambda_F at T_mean (115 degC): %.4f W/(m*K)  [reference]\n', ...
    lambdaF_scalar);
