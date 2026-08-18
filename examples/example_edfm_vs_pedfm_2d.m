%% Example: Quasi-2D EDFM vs pEDFM Thermal Comparison
% Compares standard EDFM and projection-based EDFM (pEDFM) for thermal
% transport in a quasi-2D (single-layer 3D) fractured reservoir with both
% a conductive fracture and a non-conductive barrier at 45 degrees.
%
% Note: pEDFM is inherently 3D (projects onto 6 hex faces). A quasi-2D
% setup (nz=1) is used to produce 2D-like visualizations while using the
% full 3D pEDFM machinery.
%
% Key insight: EDFM treats low-permeability fractures as if they don't
% exist (the matrix cells on either side remain fully connected). pEDFM
% reduces the matrix transmissibility near barriers, correctly blocking
% flow and heat transport.
%
% Setup:
%   - Fracture 1 (barrier): 45-degree diagonal, very low perm
%   - Fracture 2 (conductive): 45-degree diagonal, very high perm
%   - Hot water injection from west, production from east

clc; clear; close all;

%% Load modules
mrstModule add hfm shale geothermal thermal-hfm
mrstModule add ad-core ad-props incomp

%% 1. Create quasi-2D matrix grid (single-layer 3D)
tol     = 1e-3;
celldim = [20 20 1];
physdim = [100 100 10];   % 10 m thick single layer
G = cartGrid(celldim, physdim);
G = computeGeometry(G);
G.rock = makeShaleRock(G, 100*milli*darcy, 0.2);

%% 2. Define fracture planes (45-degree diagonals)
fracplanes = struct;

% Barrier: 45-degree diagonal from lower-left to upper-right
fracplanes(1).points   = [20 15 0; 20 15 10; 75 70 10; 75 70 0];
fracplanes(1).aperture = 1/25;
fracplanes(1).poro     = 0.2;
fracplanes(1).perm     = 1*nano*darcy;   % barrier

% Conductive: 45-degree diagonal, offset
fracplanes(2).points   = [35 25 0; 35 25 10; 85 75 10; 85 75 0];
fracplanes(2).aperture = 1/25;
fracplanes(2).poro     = 0.5;
fracplanes(2).perm     = 10000*darcy;    % conduit

fprintf('Fracture 1 (barrier): perm = %.0e darcy\n', ...
    fracplanes(1).perm / darcy);
fprintf('Fracture 2 (conduit): perm = %.0e darcy\n', ...
    fracplanes(2).perm / darcy);

%% 3. EDFM preprocessing (shale pipeline)
fprintf('\nProcessing fractures...\n');
[G, fracplanes] = EDFMshalegrid(G, fracplanes, ...
    'Tolerance', tol, 'plotgrid', false, ...
    'fracturelist', 1:numel(fracplanes));
G = fracturematrixShaleNNC3D(G, tol);
[G, fracplanes] = fracturefractureShaleNNCs3D(G, fracplanes, tol);

nMatrix = G.Matrix.cells.num;
nFrac   = G.cells.num - nMatrix;
fprintf('Grid: %d cells (%d matrix + %d fracture)\n', ...
    G.cells.num, nMatrix, nFrac);

% Save EDFM state
G_edfm = G;

%% 4. Add pEDFM projected connections
fprintf('Computing pEDFM projected NNCs...\n');
G = pMatFracNNCs3D(G, tol);
G_pedfm = G;

fprintf('EDFM NNCs: %d,  pEDFM NNCs: %d\n', ...
    size(G_edfm.nnc.cells, 1), size(G_pedfm.nnc.cells, 1));

%% 5. Common thermal/fluid properties
pRef  = 100*barsa;
fluid = initSimpleADIFluid('phases', 'W', ...
    'mu', 1e-3, 'rho', 1000, 'pRef', pRef, 'c', 0, 'cR', 0);
fluid = addThermalFluidProps(fluid, ...
    'Cp', 4200, 'lambdaF', 0.6, 'useEOS', true);
K0 = 273.15;
T0    = K0 + 20;
T_inj = K0 + 100;

%% 6. Build and run EDFM model
fprintf('\n--- EDFM Simulation ---\n');
G = G_edfm;
G.rock = addThermalRockProps(G.rock, 'lambdaR', 2.0, 'rhoR', 2700, 'CpR', 880);
rock_edfm = G.rock;
G = computeThermalNNCTransFracMatrix(G, rock_edfm, fluid);

gravity reset off;
model_edfm = GeothermalHFMModel(G, rock_edfm, fluid, 'fractureMethod', 'edfm');
model_edfm.extraStateOutput     = true;
model_edfm.outputFluxes         = true;
model_edfm.stepFunctionIsLinear = true;
model_edfm.minimumTemperature   = K0;
model_edfm.maximumTemperature   = K0 + 200*Kelvin;
model_edfm.maximumPressure      = 200e6*Pascal;

state0_edfm   = initResSol(G, pRef, 1);
state0_edfm.T = T0 * ones(G.cells.num, 1);

% BCs
faces = boundaryFaces(G);
fc    = G.faces.centroids(faces, :);
west  = abs(fc(:, 1) - 0) < 1e-6;
east  = abs(fc(:, 1) - physdim(1)) < 1e-6;
wf    = faces(west);
ef    = faces(east);

bc = addBC([], wf, 'flux', 1e-4);
bc = addBC(bc, ef, 'pressure', pRef);
Tbc   = [repmat(T_inj, numel(wf), 1); nan(numel(ef), 1)];
Hflux = [nan(numel(wf), 1); zeros(numel(ef), 1)];
bc = addThermalBCProps(bc, 'T', Tbc, 'Hflux', Hflux);

dt       = rampupTimesteps(1*year, 30*day);
schedule = simpleSchedule(dt, 'bc', bc);

fprintf('Running EDFM simulation (%d timesteps)...\n', numel(dt));
[~, states_edfm] = simulateScheduleAD(state0_edfm, model_edfm, schedule);
fprintf('EDFM simulation complete.\n');

%% 7. Build and run pEDFM model
fprintf('\n--- pEDFM Simulation ---\n');
G = G_pedfm;
G.rock = addThermalRockProps(G.rock, 'lambdaR', 2.0, 'rhoR', 2700, 'CpR', 880);
rock_pedfm = G.rock;
G = computeThermalNNCTransFracMatrix(G, rock_pedfm, fluid);

model_pedfm = GeothermalHFMModel(G, rock_pedfm, fluid, 'fractureMethod', 'pedfm');
model_pedfm.extraStateOutput     = true;
model_pedfm.outputFluxes         = true;
model_pedfm.stepFunctionIsLinear = true;
model_pedfm.minimumTemperature   = K0;
model_pedfm.maximumTemperature   = K0 + 200*Kelvin;
model_pedfm.maximumPressure      = 200e6*Pascal;

state0_pedfm   = initResSol(G, pRef, 1);
state0_pedfm.T = T0 * ones(G.cells.num, 1);

faces = boundaryFaces(G);
fc    = G.faces.centroids(faces, :);
west  = abs(fc(:, 1) - 0) < 1e-6;
east  = abs(fc(:, 1) - physdim(1)) < 1e-6;
wf    = faces(west);
ef    = faces(east);

bc_p = addBC([], wf, 'flux', 1e-4);
bc_p = addBC(bc_p, ef, 'pressure', pRef);
Tbc   = [repmat(T_inj, numel(wf), 1); nan(numel(ef), 1)];
Hflux = [nan(numel(wf), 1); zeros(numel(ef), 1)];
bc_p = addThermalBCProps(bc_p, 'T', Tbc, 'Hflux', Hflux);

schedule_p = simpleSchedule(dt, 'bc', bc_p);

fprintf('Running pEDFM simulation (%d timesteps)...\n', numel(dt));
[~, states_pedfm] = simulateScheduleAD(state0_pedfm, model_pedfm, schedule_p);
fprintf('pEDFM simulation complete.\n');

%% ========== Comparison Visualization ==========
times_days = cumsum(schedule.step.val) / day;
nM_edfm  = G_edfm.Matrix.cells.num;
nM_pedfm = G_pedfm.Matrix.cells.num;

%% Figure 1: Fracture setup (2D top view)
figure('Name', 'Fracture Setup', 'Position', [50, 50, 600, 500]);
plotCellData(G_edfm, log10(G_edfm.rock.perm/darcy), 'EdgeAlpha', 0.1);
colormap(gca, jet(25)); cbar = colorbar; title(cbar, 'log_{10}(K/D)');
view(0, 90); axis equal tight;
title('Fracture Setup: Barrier + Conduit (45 degrees)');
set(gca, 'FontSize', 11);

%% Figure 2: Temperature evolution comparison (snapshots)
snap_targets = [30, 120, 365];
snap_idx = zeros(size(snap_targets));
for k = 1:numel(snap_targets)
    [~, snap_idx(k)] = min(abs(times_days - snap_targets(k)));
end

figure('Name', 'EDFM vs pEDFM Temperature', ...
    'Position', [50, 50, 1200, 700]);
for k = 1:numel(snap_idx)
    si = snap_idx(k);

    % EDFM
    subplot(2, numel(snap_idx), k);
    plotCellData(G_edfm, states_edfm{si}.T - K0, 'EdgeColor', 'none');
    colormap(gca, jet(25)); caxis([20 100]);
    view(0, 90); axis equal tight;
    if k == 1, ylabel('EDFM', 'FontWeight', 'bold'); end
    title(sprintf('t = %d d', round(times_days(si))));
    if k == numel(snap_idx)
        cbar = colorbar; title(cbar, 'T [degC]');
    end

    % pEDFM
    subplot(2, numel(snap_idx), numel(snap_idx) + k);
    plotCellData(G_pedfm, states_pedfm{si}.T - K0, 'EdgeColor', 'none');
    colormap(gca, jet(25)); caxis([20 100]);
    view(0, 90); axis equal tight;
    if k == 1, ylabel('pEDFM', 'FontWeight', 'bold'); end
    if k == numel(snap_idx)
        cbar = colorbar; title(cbar, 'T [degC]');
    end
end

%% Figure 3: Temperature cross-section along x (at y=50)
x_edfm  = G_edfm.cells.centroids(:, 1);
x_pedfm = G_pedfm.cells.centroids(:, 1);
y_edfm  = G_edfm.cells.centroids(:, 2);
y_pedfm = G_pedfm.cells.centroids(:, 2);

dy = physdim(2)/celldim(2);
mid_y_edfm  = abs(y_edfm - physdim(2)/2) < dy/2 + 1e-6;
mid_y_pedfm = abs(y_pedfm - physdim(2)/2) < dy/2 + 1e-6;
sel_edfm  = mid_y_edfm & (1:G_edfm.cells.num)' <= nM_edfm;
sel_pedfm = mid_y_pedfm & (1:G_pedfm.cells.num)' <= nM_pedfm;

figure('Name', 'Temperature Cross-Section', ...
    'Position', [100, 100, 700, 450]);

T_e = states_edfm{end}.T(sel_edfm) - K0;
T_p = states_pedfm{end}.T(sel_pedfm) - K0;
xe  = x_edfm(sel_edfm);
xp  = x_pedfm(sel_pedfm);
[xe, idx] = sort(xe); T_e = T_e(idx);
[xp, idx] = sort(xp); T_p = T_p(idx);

plot(xe, T_e, 'b-o', 'LineWidth', 2, 'MarkerSize', 5, ...
    'DisplayName', 'EDFM');
hold on;
plot(xp, T_p, 'r-s', 'LineWidth', 2, 'MarkerSize', 5, ...
    'DisplayName', 'pEDFM');
hold off;
xlabel('x [m]');
ylabel('Temperature [degC]');
title('Temperature Profile Along Flow (y = 50 m)');
legend('Location', 'northwest');
grid on;
set(gca, 'FontSize', 12);

%% Figure 4: Thermal breakthrough comparison
east_edfm  = x_edfm > (physdim(1) - physdim(1)/celldim(1));
east_pedfm = x_pedfm > (physdim(1) - physdim(1)/celldim(1));

mean_T_edfm  = zeros(numel(states_edfm), 1);
mean_T_pedfm = zeros(numel(states_pedfm), 1);
for j = 1:numel(states_edfm)
    if ~isempty(states_edfm{j})
        mean_T_edfm(j) = mean(states_edfm{j}.T(east_edfm)) - K0;
    end
    if ~isempty(states_pedfm{j})
        mean_T_pedfm(j) = mean(states_pedfm{j}.T(east_pedfm)) - K0;
    end
end

figure('Name', 'Breakthrough Comparison', 'Position', [100, 100, 700, 450]);
plot(times_days, mean_T_edfm, 'b-', 'LineWidth', 2, ...
    'DisplayName', 'EDFM');
hold on;
plot(times_days, mean_T_pedfm, 'r--', 'LineWidth', 2, ...
    'DisplayName', 'pEDFM');
hold off;
xlabel('Time [days]');
ylabel('Mean outlet temperature [degC]');
title('Thermal Breakthrough: EDFM vs pEDFM (Quasi-2D)');
legend('Location', 'southeast');
grid on;
ylim([15 105]);
set(gca, 'FontSize', 12);

fprintf('\n=== Results ===\n');
fprintf('EDFM  final outlet T: %.1f degC\n', mean_T_edfm(end));
fprintf('pEDFM final outlet T: %.1f degC\n', mean_T_pedfm(end));
fprintf('Difference: %.1f degC\n', abs(mean_T_edfm(end) - mean_T_pedfm(end)));
