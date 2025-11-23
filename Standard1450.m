%% Final Dual-Wavelength Simulation (Presentation Ready)
% Description:
% Runs sequential simulations for 1450nm (The Wall) and 1650nm (The Window).
% - Uses FLUSH geometry (Source/Detector at Z=0).
% - Calculates signals at 3mm and 7mm.
% - Forces Log10 visualization for valid 3D comparison.

clear; clc;
MCmatlab.closeMCmatlabFigures();

% Create Results Folder
folderName = 'Simulation_Results';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

%% Simulation Settings
wavelengths = [1450, 1650];
mua_values  = [23.0, 6.0];  % 23 for 1450nm (Water Peak), 6 for 1650nm (Valley)
titles      = {'1450 nm (The Wall)', '1650 nm (The Window)'};

fprintf('Starting Batch Simulation...\n');

%% Loop Through Wavelengths
for i = 1:length(wavelengths)
    wl = wavelengths(i);
    mua = mua_values(i);
    simTitle = titles{i};
    
    fprintf('\n--- Running %s ---\n', simTitle);
    
    % Initialize Model
    model = MCmatlab.model;
    
    % 1. Geometry (Flush Contact at Z=0)
    model.G.nx = 101; model.G.ny = 101; model.G.nz = 100;
    model.G.Lx = 2.0; model.G.Ly = 2.0; model.G.Lz = 1.0;
    
    % Pass the specific mu_a for this run into the parameter cell array
    model.G.mediaPropParams = {mua}; 
    model.G.mediaPropertiesFunc = @mediaPropertiesFunc; 
    model.G.geomFunc = @geometryDefinition;
    
    % 2. Monte Carlo Setup
    model.MC.simulationTimeRequested = 30.0; % [min]
    model.MC.matchedInterfaces = false;
    model.MC.boundaryType = 1;
    model.MC.wavelength = wl;
    
    % Source: Pencil Beam at Z=0 (Flush)
    model.MC.lightSource.sourceType = 0; 
    model.MC.lightSource.zFocus = 0.001; % Just inside tissue (10 microns deep)
    
    % Run MC
    model = runMonteCarlo(model);
    
    % 3. Extract Results
    % Get surface fluence (1st voxel layer)
    FluenceVol = model.MC.NFR;        % Full 3D Data
    FluenceMap = model.MC.NFR(:,:,1); % Surface Data
    
    % Calculate Signals
    x_axis = linspace(-model.G.Lx/2, model.G.Lx/2, model.G.nx);
    y_axis = linspace(-model.G.Ly/2, model.G.Ly/2, model.G.ny);
    z_axis = linspace(0, model.G.Lz, model.G.nz);
    [X, Y] = meshgrid(x_axis, y_axis);
    R = sqrt(X.^2 + Y.^2);
    
    mask_close = (R > 0.28) & (R < 0.32); % 3mm ring
    mask_far   = (R > 0.68) & (R < 0.72); % 7mm ring
    
    Signal_Close = sum(FluenceMap(mask_close), 'all');
    Signal_Far   = sum(FluenceMap(mask_far), 'all');
    Ratio = Signal_Far / Signal_Close;
    
    % Print Text Results (UPDATED to include Close Signal)
    fprintf('  Signal Close (3mm): %e\n', Signal_Close);
    fprintf('  Signal Far   (7mm): %e\n', Signal_Far);
    fprintf('  Ratio: %.6f\n', Ratio);
    
    %% 4. Visualization & Saving
    
    % --- A. Custom 3D Slice Plot (FORCED LOG SCALE) ---
    f3D = figure('Visible', 'off'); 
    [X3, Y3, Z3] = ndgrid(x_axis, y_axis, z_axis);
    
    % Calculate Log Volume (avoiding log(0))
    LogVol = log10(FluenceVol + 1e-15);
    
    % Create Slices at X=0 and Y=0
    h = slice(Y3, X3, Z3, LogVol, 0, 0, []); 
    set(h, 'EdgeColor', 'none', 'FaceColor', 'interp');
    
    colormap('jet');
    caxis([-10 0]); % Fixed color scale for comparison
    colorbar;
    xlabel('Y [cm]'); ylabel('X [cm]'); zlabel('Depth [cm]');
    set(gca, 'ZDir', 'reverse'); % Put depth going down
    title({simTitle, 'Log10 Fluence Cross-Section'});
    view(45, 30); 
    axis equal; axis tight;
    
    % Save 3D Image
    saveas(f3D, fullfile(folderName, sprintf('%dnm_3D_Slices_Log10.png', wl)));
    close(f3D);
    
    % --- B. Interactive Figure (Standard MCmatlab) ---
    fInteractive = figure('Visible', 'off');
    model = plot(model, 'MC');
    savefig(gcf, fullfile(folderName, sprintf('%dnm_Interactive.fig', wl)));
    close(gcf);

    % --- C. 2D Map (Top Down) ---
    f2D = figure('Visible', 'off');
    imagesc(x_axis, y_axis, log10(FluenceMap + 1e-15));
    colorbar; caxis([-8 0]); colormap('jet');
    title(sprintf('%s\nSurface Log10 Fluence', simTitle));
    xlabel('x [cm]'); ylabel('y [cm]');
    axis equal; hold on;
    
    % Draw Rings
    theta = linspace(0, 2*pi, 100);
    plot(0.3 * cos(theta), 0.3 * sin(theta), 'r-', 'LineWidth', 2); % 3mm
    plot(0.7 * cos(theta), 0.7 * sin(theta), 'w-', 'LineWidth', 2); % 7mm
    legend('Fluence', '3mm (Close)', '7mm (Far)');
    
    % Save 2D Image
    saveas(f2D, fullfile(folderName, sprintf('%dnm_2D_Map.png', wl)));
    close(f2D);
    
    fprintf('  Plots saved to %s\n', folderName);
end

fprintf('\nBatch Simulation Complete.\n');

%% Geometry Function
function M = geometryDefinition(X,Y,Z,parameters)
  % Gel starts exactly at Z=0. No Air layer on top.
  M = 2 * ones(size(X)); 
end

%% Media Properties Function
function mediaProperties = mediaPropertiesFunc(parameters)
  % Retrieve the dynamic mu_a
  current_mua = parameters{1}; 
  
  mediaProperties = MCmatlab.mediumProperties;
  
  % Medium 1: AIR (Used for outside boundary definition)
  j=1;
  mediaProperties(j).name = 'air';
  mediaProperties(j).mua = 1e-8; mediaProperties(j).mus = 1e-8; mediaProperties(j).g = 1; mediaProperties(j).n = 1.0;
  
  % Medium 2: GEL PHANTOM (Dynamic)
  j=2;
  mediaProperties(j).name = 'Gel Phantom';
  mediaProperties(j).mua = current_mua; 
  mediaProperties(j).mus = 100.0; 
  mediaProperties(j).g   = 0.9;   
  mediaProperties(j).n   = 1.33;  
end
