%% Description
% Double Loop Simulation:
% Loop 1: Wavelength Scenarios (1450nm vs 1650nm)
% Loop 2: Scattering Coefficients (10 to 10000)
% Output: Comparison plot of Absorbed Fraction vs Scattering.

%% MCmatlab abbreviations
% G: Geometry, MC: Monte Carlo, FMC: Fluorescence Monte Carlo, HS: Heat
% simulation, M: Media array, FR: Fluence rate, FD: Fractional damage.

%% Geometry definition
clc;
MCmatlab.closeMCmatlabFigures();
model = MCmatlab.model;

%% 1. Define Experimental Coefficients
ConcentrationWater = 0.8; 
ConcentrationGel   = 0.1;
ConcentrationIL    = 0.1;

ExtinctionWater1450    = 28.8; 
ExtinctionWater1650    = 5.7; 
ExtinctionGel      = 0; 
ExtinctionIL       = 0;

% Calculate Absorption (Mua) for both wavelengths
mua1450 = ConcentrationWater*ExtinctionWater1450 + ConcentrationGel*ExtinctionGel + ConcentrationIL*ExtinctionIL;
mua1650 = ConcentrationWater*ExtinctionWater1650 + ConcentrationGel*ExtinctionGel + ConcentrationIL*ExtinctionIL;

% Create lists for the Outer Loop
muaList = [mua1450, mua1650];
wavelengthNames = {'1450 nm', '1650 nm'};

% Create list for the Inner Loop
musValues = [10, 100, 1000, 10000]; 

% Initialize a matrix to store results for the final plot
% Rows = Wavelengths, Columns = Scattering values
resultsAbsorbed = zeros(2, length(musValues));

%% 2. Setup Base Model
model.G.nx                = 101; 
model.G.ny                = 101; 
model.G.nz                = 150; 
model.G.Lx                = 7.5; 
model.G.Ly                = 5.6; 
model.G.Lz                = 1.5; 
model.G.mediaPropertiesFunc = @mediaPropertiesFunc; 
model.G.geomFunc            = @geometryDefinition; 

% Plot geometry once
model = plot(model,'G');

% Simulation Parameters
model.MC.simulationTimeRequested  = .1; 
model.MC.matchedInterfaces        = true; 
model.MC.boundaryType             = 1; 
model.MC.lightSource.sourceType   = 0; 
model.MC.lightSource.xFocus       = 0; 
model.MC.lightSource.yFocus       = 0; 
model.MC.lightSource.zFocus       = 0; 
model.MC.lightSource.theta        = 0; 
model.MC.lightSource.phi          = 0; 

%% 3. The Double Simulation Loop (CORRECTED)

% Calculate Voxel Volume (needed for integration)
dx = model.G.Lx / model.G.nx;
dy = model.G.Ly / model.G.ny;
dz = model.G.Lz / model.G.nz;
voxelVolume = dx * dy * dz;

for w = 1:length(muaList) % --- OUTER LOOP (Wavelengths) ---
    currentMua = muaList(w);
    currentWavelengthName = wavelengthNames{w};
    
    for s = 1:length(musValues) % --- INNER LOOP (Scattering) ---
        currentMus = musValues(s);
        
        fprintf('------------------------------------------------\n');
        fprintf('Scenario: %s | Mus: %d | Mua: %.2f\n', currentWavelengthName, currentMus, currentMua);
        
        % Pass parameters
        model.MC.wavelength = str2double(currentWavelengthName(1:4)); 
        model.G.mediaPropParams = {currentMus, currentMua};
        
        % Run Simulation
        model = runMonteCarlo(model);
        
        % Plot (Optional - can comment out to speed up)
        model = plot(model,'MC');
        drawnow;
        
        % --- EXTRACT DATA (FIXED METHOD) ---
        % We calculate total absorption by summing: NFR * mua * Volume
        % Note: This assumes the whole grid is tissue (valid since zSurface=0)
        
        totalFluenceSum = sum(model.MC.NFR(:));
        totalAbsorbedFraction = totalFluenceSum * currentMua * voxelVolume;
        
        % Store result
        resultsAbsorbed(w, s) = totalAbsorbedFraction;
        
        fprintf('Result: %.1f%% Absorbed\n', totalAbsorbedFraction * 100);
    end
end

%% 4. Final Comparison Plot
figure('Name', 'Absorption vs Scattering Comparison', 'NumberTitle', 'off');
hold on;
grid on;

% Plot 1450nm Data (Row 1)
semilogx(musValues, resultsAbsorbed(1,:)*100, '-o', 'LineWidth', 2, 'DisplayName', '1450 nm (High Absorption)');

% Plot 1650nm Data (Row 2)
semilogx(musValues, resultsAbsorbed(2,:)*100, '-s', 'LineWidth', 2, 'DisplayName', '1650 nm (Low Absorption)');

xlabel('Scattering Coefficient \mu_s [cm^{-1}]');
ylabel('Total Absorbed Fraction [%]');
title('Impact of Scattering on Absorption Depth');
legend('show');
ylim([0 100]);

% Add annotation
annotation('textbox', [0.15, 0.2, 0.3, 0.1], 'String', 'High scattering causes reflectance (lower absorption)', 'FitBoxToText', 'on', 'BackgroundColor', 'white');

hold off;

%% Helper Functions

function M = geometryDefinition(X,Y,Z,parameters)
  zSurface = 0.0;
  M = ones(size(X)); % Air
  M(Z > zSurface) = 2; % "Standard" tissue
end

function mediaProperties = mediaPropertiesFunc(parameters)
  % Unpack parameters
  if isempty(parameters)
      dynamicMus = 100;
      dynamicMua = 1;
  else
      dynamicMus = parameters{1}; 
      dynamicMua = parameters{2}; 
  end

  mediaProperties = MCmatlab.mediumProperties;
  
  j=1;
  mediaProperties(j).name  = 'air';
  mediaProperties(j).mua   = 1e-8; 
  mediaProperties(j).mus   = 1e-8; 
  mediaProperties(j).g     = 1; 
  
  j=2;
  mediaProperties(j).name  = 'standard tissue';
  mediaProperties(j).mua   = dynamicMua; 
  mediaProperties(j).mus   = dynamicMus; 
  mediaProperties(j).g     = 0.9; 
end
