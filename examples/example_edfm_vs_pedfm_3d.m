%% Example: 3D EDFM vs pEDFM Thermal Comparison
% Compares standard EDFM and projection-based EDFM (pEDFM) for thermal
% transport through a 3D fractured reservoir with both a conductive
% fracture and a non-conductive barrier at 45 degrees.
%
% Key insight: standard EDFM cannot accurately represent flow barriers
% because it does not modify the matrix-matrix transmissibility near
% fractures. pEDFM adds projected connections that reduce matrix
% transmissibility near low-permeability fractures, correctly capturing
% the barrier effect.
%
% Setup:
%   - Fracture 1 (barrier): 45-degree diagonal, very low perm
%   - Fracture 2 (conductive): 45-degree diagonal, very high perm
%   - Hot water injection from west, production from east
%
% Uses the shale pipeline: EDFMshalegrid -> fracturematrixShaleNNC3D
% -> fracturefractureShaleNNCs3D -> pMatFracNNCs3D

clc; clear; close all;

%% Load modules
mrstModule add hfm shale geothermal thermal-hfm
mrstModule add ad-core ad-props incomp

%% 1. Create 3D matrix grid
tol     = 1e-3;
celldim = [10 10 5];
physdim = [100 100 50];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);
G.rock = makeShaleRock(G, 100*milli*darcy, 0.2);

%% 2. Define fracture planes (45-degree diagonals)
fracplanes = struct;

% Barrier: 45-degree diagonal from lower-left to upper-right
fracplanes(1).points   = [20 15 5; 20 15 45; 75 70 45; 75 70 5];
fracplanes(1).aperture = 1/25;
fracplanes(1).poro     = 0.2;
fracplanes(1).perm     = 1*nano*darcy;   % barrier

% Conductive: 45-degree diagonal, offset
fracplanes(2).points   = [35 25 5; 35 25 45; 85 75 45; 85 75 5];
fracplanes(2).aperture = 1/25;
fracplanes(2).poro     = 0.5;
fracplanes(2).perm     = 10000*darcy;    % conduit

fprintf('Fracture 1 (barrier): perm = %.0e darcy\n', ...
    fracplanes(1).perm / darcy);
fprintf('Fracture 2 (conduit): perm = %.0e darcy\n', ...
    fracplanes(2).perm / darcy);

%% 3. EDFM preprocessing (shale pipeline)
fprintf('\nProcessing fractures with EDFMshalegrid...\n');
[G, fracplanes] = EDFMshalegrid(G, fracplanes, ...
    'Tolerance', tol, 'plotgrid', false, ...
    'fracturelist', 1:numel(fracplanes));

% Fracture-matrix NNCs
G = fracturematrixShaleNNC3D(G, tol);

% Fracture-fracture NNCs
[G, fracplanes] = fracturefractureShaleNNCs3D(G, fracplanes, tol);

nMatrix = G.Matrix.cells.num;
nFrac   = G.cells.num - nMatrix;
fprintf('Grid: %d cells (%d matrix + %d fracture)\n', ...
    G.cells.num, nMatrix, nFrac);

% Save EDFM state before adding pEDFM connections
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

% Boundary conditions
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

% Same BCs (recompute for pEDFM grid which has same faces)
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
x_edfm   = G_edfm.cells.centroids(:, 1);
x_pedfm  = G_pedfm.cells.centroids(:, 1);

%% Figure 1: Fracture network overview
figure('Name', 'Fracture Setup', 'Position', [50, 50, 900, 500]);
plotGrid(G_edfm, 1:nM_edfm, 'FaceColor', [0.9 0.9 0.9], ...
    'FaceAlpha', 0.1, 'EdgeAlpha', 0.05);
hold on;
% Color-code fractures by type
frac1_cells = (nM_edfm+1):(nM_edfm + G_edfm.FracGrid.Frac1.cells.num);
frac2_start = nM_edfm + G_edfm.FracGrid.Frac1.cells.num + 1;
frac2_cells = frac2_start:G_edfm.cells.num;
plotCellData(G_edfm, log10(G_edfm.rock.perm/darcy), frac1_cells, ...
    'EdgeAlpha', 0.3);
plotCellData(G_edfm, log10(G_edfm.rock.perm/darcy), frac2_cells, ...
    'EdgeAlpha', 0.3);
hold off;
colormap(gca, jet(25)); cbar = colorbar; title(cbar, 'log_{10}(K/D)');
view([-45, 25]); axis tight equal;
title('Fracture Setup: Barrier + Conduit (45 degrees)');
set(gca, 'FontSize', 11);

%% Figure 2: Final temperature comparison (side by side)
figure('Name', 'EDFM vs pEDFM: Final Temperature', ...
    'Position', [50, 50, 1200, 500]);

subplot(1,2,1);
T_edfm = states_edfm{end}.T - K0;
plotCellData(G_edfm, T_edfm, 1:nM_edfm, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.6);
colormap(gca, jet(25)); caxis([20 100]);
cbar = colorbar; title(cbar, 'T [degC]');
view([-45, 25]); axis tight equal;
title('EDFM - Final Temperature');
set(gca, 'FontSize', 11);

subplot(1,2,2);
T_pedfm = states_pedfm{end}.T - K0;
plotCellData(G_pedfm, T_pedfm, 1:nM_pedfm, 'EdgeAlpha', 0.1, 'FaceAlpha', 0.6);
colormap(gca, jet(25)); caxis([20 100]);
cbar = colorbar; title(cbar, 'T [degC]');
view([-45, 25]); axis tight equal;
title('pEDFM - Final Temperature');
set(gca, 'FontSize', 11);

%% Figure 3: Temperature cross-section along x (at y=50, z=25)
y_edfm  = G_edfm.cells.centroids(:, 2);
z_edfm  = G_edfm.cells.centroids(:, 3);
y_pedfm = G_pedfm.cells.centroids(:, 2);
z_pedfm = G_pedfm.cells.centroids(:, 3);

dy = physdim(2)/celldim(2);
dz = physdim(3)/celldim(3);
mid_y_edfm  = abs(y_edfm - physdim(2)/2) < dy/2 + 1e-6;
mid_z_edfm  = abs(z_edfm - physdim(3)/2) < dz/2 + 1e-6;
mid_y_pedfm = abs(y_pedfm - physdim(2)/2) < dy/2 + 1e-6;
mid_z_pedfm = abs(z_pedfm - physdim(3)/2) < dz/2 + 1e-6;

sel_edfm  = mid_y_edfm & mid_z_edfm & (1:G_edfm.cells.num)' <= nM_edfm;
sel_pedfm = mid_y_pedfm & mid_z_pedfm & (1:G_pedfm.cells.num)' <= nM_pedfm;

figure('Name', 'Temperature Cross-Section', ...
    'Position', [100, 100, 700, 450]);

T_line_edfm  = states_edfm{end}.T(sel_edfm) - K0;
T_line_pedfm = states_pedfm{end}.T(sel_pedfm) - K0;
x_line_edfm  = x_edfm(sel_edfm);
x_line_pedfm = x_pedfm(sel_pedfm);

[x_line_edfm, idx]  = sort(x_line_edfm);
T_line_edfm = T_line_edfm(idx);
[x_line_pedfm, idx] = sort(x_line_pedfm);
T_line_pedfm = T_line_pedfm(idx);

plot(x_line_edfm, T_line_edfm, 'b-o', 'LineWidth', 2, ...
    'MarkerSize', 6, 'DisplayName', 'EDFM');
hold on;
plot(x_line_pedfm, T_line_pedfm, 'r-s', 'LineWidth', 2, ...
    'MarkerSize', 6, 'DisplayName', 'pEDFM');
hold off;
xlabel('x [m]');
ylabel('Temperature [degC]');
title('Temperature Profile Along Flow (y=50, z=25)');
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
title('Thermal Breakthrough: EDFM vs pEDFM');
legend('Location', 'southeast');
grid on;
ylim([15 105]);
set(gca, 'FontSize', 12);

fprintf('\n=== Results ===\n');
fprintf('EDFM  final outlet T: %.1f degC\n', mean_T_edfm(end));
fprintf('pEDFM final outlet T: %.1f degC\n', mean_T_pedfm(end));
fprintf('Difference: %.1f degC\n', abs(mean_T_edfm(end) - mean_T_pedfm(end)));
