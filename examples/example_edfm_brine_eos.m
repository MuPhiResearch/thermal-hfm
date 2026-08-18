%% Example: Thermal EDFM with Brine EOS and NaCl Transport
% Demonstrates the thermal-hfm module with brine physics:
%
%   1. CompositionalBrineFluid (H2O + NaCl) - advective salinity transport.
%      Molecular diffusivities set to 0 here (advection-only demonstration).
%      Non-zero diffusivities are now supported via
%      EDFMDynamicMolecularTransmissibility (see test_step7_brine_diffusion).
%   2. useEOS = true  - Spivey (2004) T-p-dependent density and viscosity.
%   3. lambdaF_func   - @(p,T) temperature-dependent fluid conductivity
%      (IAPWS 2008) passed directly to addThermalFluidProps, activating
%      EDFMDynamicHeatTransmissibility for full dynamic NNC thermal
%      transmissibilities at every Newton iteration.
%   4. plotToolbar with a cell array of states for interactive time-stepping
%      (use the slider or play/pause button in the MRST GUI toolbar).
%
% Setup:
%   - 3D domain (120 x 80 x 40 m), single diagonal high-perm fracture.
%   - Hot saline brine injected at the west face (flux BC).
%   - Production at east face (pressure BC).
%   - Initial: 80 degC, 15 wt% NaCl.  Injection: 150 degC, 5 wt% NaCl.
%
% NOTE on bc.components:
%   All BC faces (inlet AND outlet) must have non-nan component values.
%   For pressure/outflow BC faces, we specify the initial reservoir
%   composition - this value is only active for inflow and has no effect
%   on outflow (the equations use the upwind cell composition for outflow).
%
% NOTE on dynamic lambdaF:
%   Passing a @(p,T) function handle to addThermalFluidProps sets
%   fluid.lambdaF to the handle, which activates dynamicHeatTransFluid()
%   in GeothermalHFMModel.  The model then installs
%   EDFMDynamicHeatTransmissibility (instead of the standard
%   DynamicTransmissibility that does not support EDFM NNCs) to recompute
%   thermal transmissibilities including NNC contributions at every Newton
%   iteration.  G.nnc.CI_dyn (stored by computeThermalNNCTransFracMatrix)
%   provides the geometric NNC factor.

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

%% 3. Define fracture plane (diagonal, high-perm conduit)
fracplanes = struct;
fracplanes(1).points   = [20 10 4; 20 10 36; 100 70 36; 100 70 4];
fracplanes(1).aperture = 1/25;
fracplanes(1).poro     = 0.5;
fracplanes(1).perm     = 5000*darcy;

%% 4. EDFM preprocessing (shale pipeline)
fprintf('Processing fracture with EDFMshalegrid...\n');
[G, fracplanes] = EDFMshalegrid(G, fracplanes, ...
    'Tolerance', tol, 'plotgrid', false, ...
    'fracturelist', 1:numel(fracplanes));
G = fracturematrixShaleNNC3D(G, tol);

nMatrix = G.Matrix.cells.num;
nFrac   = G.cells.num - nMatrix;
nNNC    = size(G.nnc.cells, 1);
fprintf('Grid: %d cells (%d matrix + %d fracture), %d NNCs\n', ...
    G.cells.num, nMatrix, nFrac, nNNC);

%% 5. Rock thermal properties
G.rock = addThermalRockProps(G.rock, ...
    'lambdaR', 2.5, 'rhoR', 2650, 'CpR', 900);
rock = G.rock;

%% 6. Fluid: brine with Spivey EOS + dynamic lambdaF handle
pRef = 200*barsa;

% T-dependent lambdaF (IAPWS 2008 polynomial, valid 20-200 degC).
% Passed directly to addThermalFluidProps as a function handle, enabling
% EDFMDynamicHeatTransmissibility to recompute NNC thermal transmissibilities
% at each Newton iteration.
lambdaF_func   = @(p, T) max(0.4, ...
    0.5563 + 2.31e-3.*(T - 273.15) - 8.7e-6.*(T - 273.15).^2);
lambdaF_scalar = lambdaF_func(pRef, 273.15 + 115);   % ~0.707 W/(m·K), for Figure 4

% Base single-phase water fluid
fluid = initSimpleADIFluid('phases', 'W', ...
    'n'   , 1      , ...
    'mu'  , 1e-3   , ...
    'rho' , 1000   , ...
    'pRef', pRef   , ...
    'c'   , 4.5e-10, ...
    'cR'  , 0);

% Spivey (2004) EOS density/viscosity.
% addThermalFluidProps uses merge_options which enforces type matching on 'lambdaF',
% so we pass the scalar value (any valid double) and then override fluid.lambdaF
% with the function handle.  GeothermalHFMModel detects the handle via
% dynamicHeatTransFluid() and installs EDFMDynamicHeatTransmissibility.
fluid = addThermalFluidProps(fluid, ...
    'Cp'     , 4200          , ...
    'lambdaF', lambdaF_scalar, ...
    'useEOS' , true);
fluid.lambdaF = lambdaF_func;   % override: activate dynamic NNC heat transmissibility

% Compositional brine: H2O + NaCl, diffusivities = 0 (advection only).
% To enable NaCl diffusion, set the third argument to e.g. [0, 1.5e-9]
% and ensure rock.tau (tortuosity) is defined.
compFluid = CompositionalBrineFluid( ...
    {'H2O'              , 'NaCl'            }, ...
    [18.015281*gram/mol , 58.4428*gram/mol  ], ...
    [0                  , 0                 ]);

%% 7. Compute static NNC thermal transmissibilities
G = computeThermalNNCTransFracMatrix(G, rock, fluid, ...
    'pRef', pRef, 'TRef', 273.15 + 20);

%% 8. Create model
gravity reset off;
model = GeothermalHFMModel(G, rock, fluid, compFluid, 'fractureMethod', 'edfm');
model.extraStateOutput = true;
model.outputFluxes     = true;

K0 = 273.15;
model.minimumTemperature = K0;
model.maximumTemperature = K0 + 300;
model.maximumPressure    = 500e6;
% validateModel populates model.Components (needed by getMassFraction)
model = model.validateModel();

%% 9. Initial state: hot saline formation water
T0      = K0 + 80;     % 80 degC
X0_NaCl = 0.15;        % 15 wt% NaCl

state0   = initResSol(G, pRef, 1);
state0.T = T0 * ones(G.cells.num, 1);
X0_mass  = repmat([1 - X0_NaCl, X0_NaCl], G.cells.num, 1);
state0.components = model.getMoleFraction(X0_mass);

%% 10. Boundary conditions
% IMPORTANT: bc.components must cover ALL BC faces (no nan).
% Outlet (pressure BC) receives the initial reservoir composition;
% this value only activates for inflow - outflow uses the upwind cell.
faces = boundaryFaces(G);
fc    = G.faces.centroids(faces, :);
west  = abs(fc(:, 1) - 0)           < 1e-6;
east  = abs(fc(:, 1) - physdim(1))  < 1e-6;
wf    = faces(west);
ef    = faces(east);

T_inj      = K0 + 150;   % 150 degC injected brine
X_inj_NaCl = 0.05;       % 5 wt% NaCl (fresh brine)

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
fprintf('Schedule: %d timesteps, total %.0f days\n\n', numel(dt), sum(dt)/day);

%% 12. Run simulation
fprintf('Running brine/EOS EDFM simulation...\n');
[~, states] = simulateScheduleAD(state0, model, schedule);
fprintf('Simulation complete.\n\n');

%% ========== Visualization ==========
K0          = 273.15;
times_days  = cumsum(schedule.step.val) / day;
frac_cells  = (nMatrix+1):G.cells.num;
x_all       = G.cells.centroids(:, 1);
east_mask   = x_all > (physdim(1) - physdim(1)/celldim(1));
east_matrix = east_mask & ((1:G.cells.num)' <= nMatrix);

%% Figure 1: Interactive plotToolbar with time slider
% Pass the full states cell array so the slider lets you step through time.
% Use the play button (triangle) in the MRST figure toolbar to animate.
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
title('Brine/EOS EDFM - Interactive (plotToolbar, use slider or play button)');
set(gca, 'FontSize', 11);

%% Figure 2: Final state - fracture highlighted
fig2 = mrstFigure();
ax2  = axes(fig2);
T_fin = states{end}.T - K0;
plotCellData(G, T_fin, 1:nMatrix, 'EdgeAlpha', 0.03, 'FaceAlpha', 0.5);
hold(ax2, 'on');
plotCellData(G, T_fin, frac_cells, 'EdgeAlpha', 0.3);
hold(ax2, 'off');
colormap(ax2, hot(32)); clim(ax2, [T0 - K0, T_inj - K0]);
cbar2 = colorbar(ax2); title(cbar2, 'T [°C]');
view(ax2, [-40, 25]); axis(ax2, 'tight'); axis(ax2, 'equal');
title(ax2, sprintf('Final Temperature  (t = %.0f d)', times_days(end)));
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

figure('Name', 'Breakthrough', 'Position', [100, 100, 820, 380]);
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
title('Thermal Breakthrough, EOS Density, and Salinity Dilution  (Brine EOS EDFM)');
legend({'Temperature [°C]','Density [kg/m^3]','NaCl fraction'}, 'Location', 'east');
grid on; set(gca, 'FontSize', 12);

%% Figure 4: T-dependent lambda_F profile (IAPWS 2008)
T_range  = linspace(20, 200, 200);
lam_vals = lambdaF_func(pRef, T_range + 273.15);
figure('Name', 'lambdaF(T)', 'Position', [150, 150, 560, 340]);
plot(T_range, lam_vals, 'k-', 'LineWidth', 2);
xline(T0 - K0,          'b--', 'T_0 = 80°C',       'LabelVerticalAlignment', 'bottom');
xline(T_inj - K0,       'r--', 'T_{inj} = 150°C',  'LabelVerticalAlignment', 'bottom');
xline(115,              'g--', 'T_{mean} (used)',   'LabelVerticalAlignment', 'top');
yline(lambdaF_scalar,   'g:',  sprintf('%.4f W/(m·K)', lambdaF_scalar));
xlabel('Temperature [°C]');
ylabel('\lambda_F  [W m^{-1} K^{-1}]');
title({'T-Dependent Fluid Conductivity - IAPWS 2008', ...
       'Simulation uses @(p,T) handle (EDFMDynamicHeatTransmissibility)'});
grid on; set(gca, 'FontSize', 12);

%% Summary
fprintf('=== Results Summary ===\n');
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
