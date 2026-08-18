%% Example: 3D Thermal EDFM Simulation
% Demonstrates the thermal-hfm module on a 3D grid with a single vertical
% fracture plane. Hot water is injected from the west boundary and produced
% at the east boundary. The fracture acts as a preferential flow path,
% channeling heat transport.
%
% Uses the hfm 3D pipeline: preProcessingFractures -> assembleGlobalGrid
% -> computeEffectiveTrans, then thermal-hfm for thermal NNC coupling.

clc; clear; close all;

%% Load modules
mrstModule add hfm geothermal thermal-hfm
mrstModule add ad-core ad-props incomp
% The hfm-native 3D pipeline (preProcessingFractures -> markcells) needs
% getEnclosingCellsByFace from msrsb, and coarse-grid utilities from
% coarsegrid (same dependencies as the stock hfm 3D example
% simple2phHorizontalWell3D.m).
mrstModule add msrsb coarsegrid
checkLineSegmentIntersect;

%% 1. Create 3D matrix grid
celldim = [10 10 5];
physdim = [100 100 50];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

%% 2. Define fracture plane (vertical, aligned with flow direction)
fracplanes = struct;
% Offset off the matrix cell faces (y=45 sits inside a cell row, x=21/79 inside
% cell columns) so every fracture point locates cleanly in the matrix grid.
fracplanes(1).points = [21 45 5;    % corner 1
                        21 45 45;   % corner 2
                        79 45 45;   % corner 3
                        79 45 5];   % corner 4
fracplanes(1).aperture = 1/25;      % 0.04 m aperture
checkIfCoplanar(fracplanes);

%% 3. Process fractures (grid + NNC identification)
fprintf('Processing fractures...\n');
[G, fracplanes] = preProcessingFractures(G, fracplanes, ...
    'fractureCellSize', 0.3);

%% 4. Rock properties
G.rock.perm = ones(G.cells.num, 1) * 100 * milli * darcy;
G.rock.poro = 0.2 * ones(G.cells.num, 1);
K_frac = 10000;  % darcy
G = makeRockFrac(G, K_frac, 'porosity', 0.5);

%% 5. Assemble global grid and compute flow transmissibilities
G = assembleGlobalGrid(G);
G = computeEffectiveTrans(G);

% Add NNC type field (required by transmultEDFM in setupEDFMOperatorsTPFA)
nNNC = size(G.nnc.cells, 1);
G.nnc.type = repmat({'frac-matrix'}, nNNC, 1);

nMatrix = G.Matrix.cells.num;
nFrac   = G.cells.num - nMatrix;
fprintf('Grid: %d cells (%d matrix + %d fracture), %d NNCs\n', ...
    G.cells.num, nMatrix, nFrac, nNNC);

%% 6. Thermal properties
G.rock = addThermalRockProps(G.rock, ...
    'lambdaR', 2.0, 'rhoR', 2700, 'CpR', 880);
rock = G.rock;

pRef = 100*barsa;
fluid = initSimpleADIFluid('phases', 'W', ...
    'mu', 1e-3, 'rho', 1000, 'pRef', pRef, 'c', 0, 'cR', 0);
fluid = addThermalFluidProps(fluid, ...
    'Cp', 4200, 'lambdaF', 0.6, 'useEOS', true);

%% 7. Compute thermal NNC transmissibilities
G = computeThermalNNCTransFracMatrix(G, rock, fluid);

%% 8. Create model
gravity reset off;
model = GeothermalHFMModel(G, rock, fluid, 'fractureMethod', 'edfm');
model.extraStateOutput    = true;
model.outputFluxes        = true;
model.stepFunctionIsLinear = true;

K0 = 273.15;
model.minimumTemperature = K0;
model.maximumTemperature = K0 + 200*Kelvin;
model.maximumPressure    = 200e6*Pascal;

%% 9. Initial state
T0    = K0 + 20;   % 20 degC
T_inj = K0 + 100;  % 100 degC
state0   = initResSol(G, pRef, 1);
state0.T = T0 * ones(G.cells.num, 1);

%% 10. Boundary conditions
faces = boundaryFaces(G);
fc    = G.faces.centroids(faces, :);
west  = abs(fc(:, 1) - 0) < 1e-6;
east  = abs(fc(:, 1) - physdim(1)) < 1e-6;
wf    = faces(west);
ef    = faces(east);

bc = addBC([], wf, 'flux', 5e-4);       % injection per face [m3/s]
bc = addBC(bc, ef, 'pressure', pRef);    % production pressure
Tbc   = [repmat(T_inj, numel(wf), 1); nan(numel(ef), 1)];
Hflux = [nan(numel(wf), 1); zeros(numel(ef), 1)];
bc = addThermalBCProps(bc, 'T', Tbc, 'Hflux', Hflux);

%% 11. Schedule (1 year with ramp-up timesteps)
dt       = rampupTimesteps(1*year, 30*day);
schedule = simpleSchedule(dt, 'bc', bc);
fprintf('Schedule: %d timesteps, total %.0f days\n', numel(dt), sum(dt)/day);

%% 12. Run simulation
fprintf('Running simulation...\n');
[~, states] = simulateScheduleAD(state0, model, schedule);
fprintf('Simulation complete.\n');

%% ========== Visualization ==========
times_days = cumsum(schedule.step.val) / day;
frac_cells = (nMatrix+1):G.cells.num;
x = G.cells.centroids(:, 1);
east_cells = x > (physdim(1) - physdim(1)/celldim(1));

%% Figure 1: Fracture network and grid overview
figure('Name', 'Fracture Network', 'Position', [50, 50, 900, 500]);
% Semi-transparent matrix grid
plotGrid(G, 1:nMatrix, 'FaceColor', [0.9 0.9 0.9], 'FaceAlpha', 0.15, ...
    'EdgeAlpha', 0.08);
hold on;
% Fracture plane as a filled patch
pts = fracplanes(1).points;
patch(pts([1:end,1],1), pts([1:end,1],2), pts([1:end,1],3), ...
    'r', 'FaceAlpha', 0.5, 'EdgeColor', 'r', 'LineWidth', 2);
% Fracture cells colored by permeability
plotCellData(G, log10(G.rock.perm/darcy), frac_cells, ...
    'EdgeAlpha', 0.3);
hold off;
colormap(gca, jet(25));
cbar = colorbar; title(cbar, 'log_{10}(K/D)');
view([-45, 25]); axis tight equal;
title('3D Grid with Fracture Plane');
set(gca, 'FontSize', 11);

%% Figure 2: Animated 3D matrix temperature
fig2 = figure('Name', '3D Temperature Animation', ...
    'Position', [50, 50, 800, 600]);
ax2 = axes(fig2);
for k = 1:numel(states)
    if isempty(states{k}), continue; end
    cla(ax2);
    T_degC = states{k}.T - K0;
    plotCellData(G, T_degC, 1:nMatrix, 'EdgeAlpha', 0.05, ...
        'FaceAlpha', 0.6);
    hold(ax2, 'on');
    patch(ax2, pts([1:end,1],1), pts([1:end,1],2), pts([1:end,1],3), ...
        'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 2, ...
        'LineStyle', '--');
    hold(ax2, 'off');
    colormap(ax2, jet(25)); caxis(ax2, [20 100]);
    if k == 1, cbar2 = colorbar(ax2); title(cbar2, 'T [degC]'); end
    view(ax2, [-45, 25]); axis(ax2, 'tight'); axis(ax2, 'equal');
    title(ax2, sprintf('Matrix Temperature  |  t = %.0f d', ...
        times_days(k)));
    set(ax2, 'FontSize', 11);
    drawnow;
    pause(0.15);
end

%% Figure 3: Animated fracture cell temperature (3D view)
fig3 = figure('Name', 'Fracture Temperature Animation', ...
    'Position', [100, 100, 800, 600]);
ax3 = axes(fig3);
for k = 1:numel(states)
    if isempty(states{k}), continue; end
    cla(ax3);
    T_degC = states{k}.T - K0;
    plotGrid(G, 1:nMatrix, 'FaceColor', [0.9 0.9 0.9], ...
        'FaceAlpha', 0.08, 'EdgeAlpha', 0.03);
    hold(ax3, 'on');
    plotCellData(G, T_degC, frac_cells, 'EdgeAlpha', 0.3);
    hold(ax3, 'off');
    colormap(ax3, jet(25)); caxis(ax3, [20 100]);
    if k == 1, cbar3 = colorbar(ax3); title(cbar3, 'T [degC]'); end
    view(ax3, [-45, 25]); axis(ax3, 'tight'); axis(ax3, 'equal');
    title(ax3, sprintf('Fracture Temperature  |  t = %.0f d', ...
        times_days(k)));
    set(ax3, 'FontSize', 11);
    drawnow;
    pause(0.15);
end

%% Figure 4: Final pressure and temperature
figure('Name', 'Final State (3D)', 'Position', [100, 100, 1000, 400]);

subplot(1,2,1);
plotCellData(G, states{end}.pressure / barsa, 1:nMatrix, ...
    'EdgeAlpha', 0.1, 'FaceAlpha', 0.6);
colormap(gca, jet(25));
cbar = colorbar; title(cbar, 'P [bar]');
view([-45, 25]); axis tight equal;
title('Final Pressure');

subplot(1,2,2);
plotCellData(G, states{end}.T - K0, 1:nMatrix, ...
    'EdgeAlpha', 0.1, 'FaceAlpha', 0.6);
hold on;
plotCellData(G, states{end}.T - K0, frac_cells, 'EdgeAlpha', 0.3);
hold off;
colormap(gca, jet(25));
cbar = colorbar; title(cbar, 'T [degC]');
view([-45, 25]); axis tight equal;
title('Final Temperature (matrix + fracture)');

%% Figure 5: Thermal breakthrough curve
mean_T_east = zeros(numel(states), 1);
mean_T_frac = zeros(numel(states), 1);
for j = 1:numel(states)
    if ~isempty(states{j})
        mean_T_east(j) = mean(states{j}.T(east_cells)) - K0;
        mean_T_frac(j) = mean(states{j}.T(frac_cells)) - K0;
    end
end

figure('Name', 'Thermal Breakthrough', 'Position', [100, 100, 600, 400]);
plot(times_days, mean_T_east, 'b-', 'LineWidth', 2, ...
    'DisplayName', 'East boundary (avg)');
hold on;
plot(times_days, mean_T_frac, 'r--', 'LineWidth', 2, ...
    'DisplayName', 'Fracture cells (avg)');
hold off;
xlabel('Time [days]');
ylabel('Mean temperature [degC]');
title('Thermal Breakthrough Curves');
legend('Location', 'southeast');
grid on;
ylim([15 105]);
set(gca, 'FontSize', 12);

fprintf('\nFinal outlet temperature: %.1f degC\n', mean_T_east(end));
fprintf('Final fracture temperature: %.1f degC\n', mean_T_frac(end));
