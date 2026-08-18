%% Example: 2D Thermal EDFM Simulation
% Demonstrates the thermal-hfm module on a 2D grid with two intersecting
% fractures. Hot water is injected from the west boundary and produced at
% the east boundary, showing thermal breakthrough through the fracture
% network.
%
% It uses a modular approach: stock hfm for the EDFM grid and NNCs, stock
% geothermal for the thermal properties, and GeothermalHFMModel for the
% coupling.

clc; clear; close all;

%% Load modules
mrstModule add hfm geothermal thermal-hfm
mrstModule add ad-core ad-props incomp
checkLineSegmentIntersect;

%% 1. Create grid and fracture network
celldim = [40 40];
physdim = [10 10];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

% Two intersecting fractures
fl = [2, 5, 8, 5;   % horizontal fracture
      5, 2, 5, 8];  % vertical fracture

[G, fracture] = processFracture2D(G, fl);
fracture.aperture = 0.01;  % 1 cm aperture
G = CIcalculator2D(G, fracture);
[G, F, fracture] = gridFracture2D(G, fracture, ...
    'min_size', 0.05, 'cell_size', 0.1);

%% 2. Rock properties
% Matrix: moderate permeability
G.rock.perm = ones(G.cells.num, 1) * milli * darcy;
G.rock.poro = 0.2 * ones(G.cells.num, 1);

% Fractures: high permeability
K_frac = 1000; % darcy
G = makeRockFrac(G, K_frac, 'porosity', 0.5);

%% 3. Stock hfm NNC setup (flow transmissibilities)
[G, T] = defineNNCandTrans(G, F, fracture);

fprintf('Grid: %d cells (%d matrix + %d fracture), %d NNCs\n', ...
    G.cells.num, celldim(1)*celldim(2), ...
    G.cells.num - celldim(1)*celldim(2), size(G.nnc.cells, 1));

%% 4. Add thermal properties (using stock geothermal functions)
% Rock thermal properties
G.rock = addThermalRockProps(G.rock, ...
    'lambdaR', 2.0, ...    % W/(m*K) rock thermal conductivity
    'rhoR',    2700, ...    % kg/m3 rock density
    'CpR',     880);        % J/(kg*K) rock heat capacity
rock = G.rock;

% Fluid with EOS (p/T-dependent density and viscosity)
pRef = 10*barsa;
fluid = initSimpleADIFluid('phases', 'W', ...
    'mu', 1e-3, 'rho', 1000, 'pRef', pRef, 'c', 0, 'cR', 0);
fluid = addThermalFluidProps(fluid, ...
    'Cp',      4200, ...    % J/(kg*K) fluid heat capacity
    'lambdaF', 0.6, ...     % W/(m*K) fluid thermal conductivity
    'useEOS',  true);       % Spivey et al. (2004) EOS

%% 5. Compute thermal NNC transmissibilities (thermal-hfm functions)
G = computeThermalNNCTransFracMatrix(G, rock, fluid);
G = computeThermalNNCTransFracFrac(G, F, fracture, rock, fluid);

%% 6. Create the coupled thermal-EDFM model
gravity reset off;
model = GeothermalHFMModel(G, rock, fluid, 'fractureMethod', 'edfm');
model.extraStateOutput    = true;
model.outputFluxes        = true;
model.stepFunctionIsLinear = true;

% EOS validity bounds
K0 = 273.15;
model.minimumTemperature = K0;
model.maximumTemperature = K0 + 200*Kelvin;
model.maximumPressure    = 200e6*Pascal;

%% 7. Initial state (20 degC, hydrostatic pressure)
T0    = K0 + 20;   % 20 degC
T_inj = K0 + 100;  % 100 degC injection
state0   = initResSol(G, pRef, 1);
state0.T = T0 * ones(G.cells.num, 1);

%% 8. Boundary conditions
% West: flux BC with hot water injection
% East: pressure BC (outflow)
faces = boundaryFaces(G);
fc    = G.faces.centroids(faces, :);
west  = fc(:, 1) == min(G.faces.centroids(:, 1));
east  = fc(:, 1) == max(G.faces.centroids(:, 1));
wf    = faces(west);
ef    = faces(east);

bc = addBC([], wf, 'flux', 5e-9);       % injection rate per face [m3/s]
bc = addBC(bc, ef, 'pressure', pRef);    % production pressure
Tbc   = [repmat(T_inj, numel(wf), 1); nan(numel(ef), 1)];
Hflux = [nan(numel(wf), 1); zeros(numel(ef), 1)];
bc = addThermalBCProps(bc, 'T', Tbc, 'Hflux', Hflux);

%% 9. Schedule (5 years with ramp-up timesteps)
dt       = rampupTimesteps(5*year, 60*day);
schedule = simpleSchedule(dt, 'bc', bc);

fprintf('Schedule: %d timesteps, total %.0f days\n', numel(dt), sum(dt)/day);

%% 10. Run simulation
fprintf('Running simulation...\n');
[~, states] = simulateScheduleAD(state0, model, schedule);
fprintf('Simulation complete.\n');

%% ========== Visualization ==========

%% Figure 1: Fracture network and permeability
figure('Name', 'Grid and Properties', 'Position', [100, 100, 900, 400]);

subplot(1,2,1);
plotCellData(G, log10(G.rock.perm/darcy), 'EdgeColor', 'none');
line(fl(:,1:2:3)', fl(:,2:2:4)', 'Color', 'k', 'LineWidth', 2);
colormap(gca, jet(25));
cbar = colorbar; title(cbar, 'log_{10}(K/D)');
axis equal tight; view(0, 90);
title('Permeability field');

subplot(1,2,2);
plotCellData(G, G.rock.poro, 'EdgeColor', 'none');
line(fl(:,1:2:3)', fl(:,2:2:4)', 'Color', 'k', 'LineWidth', 2);
colormap(gca, jet(25));
cbar = colorbar; title(cbar, '\phi');
axis equal tight; view(0, 90);
title('Porosity field');

%% Figure 2: Temperature snapshots at selected times
times_days = cumsum(schedule.step.val) / day;
snap_targets = [30, 180, 365, 730, 1825]; % days
snap_idx = zeros(size(snap_targets));
for k = 1:numel(snap_targets)
    [~, snap_idx(k)] = min(abs(times_days - snap_targets(k)));
end

figure('Name', 'Temperature Evolution', 'Position', [100, 100, 1200, 500]);
for k = 1:numel(snap_idx)
    subplot(1, numel(snap_idx), k);
    T_degC = states{snap_idx(k)}.T - K0;
    plotCellData(G, T_degC, 'EdgeColor', 'none');
    line(fl(:,1:2:3)', fl(:,2:2:4)', 'Color', 'w', 'LineWidth', 1);
    colormap(gca, jet(25));
    caxis([20 100]);
    axis equal tight; view(0, 90);
    title(sprintf('t = %.0f d', times_days(snap_idx(k))));
    if k == numel(snap_idx)
        cbar = colorbar; title(cbar, 'T [degC]');
    end
end

%% Figure 3: Pressure and temperature at final time
figure('Name', 'Final State', 'Position', [100, 100, 900, 400]);

subplot(1,2,1);
plotCellData(G, states{end}.pressure / barsa, 'EdgeColor', 'none');
line(fl(:,1:2:3)', fl(:,2:2:4)', 'Color', 'k', 'LineWidth', 1.5);
colormap(gca, jet(25));
cbar = colorbar; title(cbar, 'P [bar]');
axis equal tight; view(0, 90);
title('Final Pressure');

subplot(1,2,2);
plotCellData(G, states{end}.T - K0, 'EdgeColor', 'none');
line(fl(:,1:2:3)', fl(:,2:2:4)', 'Color', 'w', 'LineWidth', 1.5);
colormap(gca, jet(25));
cbar = colorbar; title(cbar, 'T [degC]');
axis equal tight; view(0, 90);
title('Final Temperature');

%% Figure 4: Thermal breakthrough curve (east boundary avg temperature)
x = G.cells.centroids(:, 1);
east_cells = x > (physdim(1) - physdim(1)/celldim(1));

mean_T_east = zeros(numel(states), 1);
for j = 1:numel(states)
    if ~isempty(states{j})
        mean_T_east(j) = mean(states{j}.T(east_cells)) - K0;
    end
end

figure('Name', 'Thermal Breakthrough', 'Position', [100, 100, 600, 400]);
plot(times_days/365, mean_T_east, 'b-', 'LineWidth', 2);
xlabel('Time [years]');
ylabel('Mean outlet temperature [degC]');
title('Thermal Breakthrough Curve (East Boundary)');
grid on;
ylim([15 105]);
set(gca, 'FontSize', 12);

fprintf('\nFinal outlet temperature: %.1f degC\n', mean_T_east(end));
fprintf('Temperature range: [%.1f, %.1f] degC\n', ...
    min(states{end}.T - K0), max(states{end}.T - K0));
