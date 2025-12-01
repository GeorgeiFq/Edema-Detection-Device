%% Description
% Simulation running loop with different scattering coefficients
% AND a custom calculated absorption coefficient.
% Run 1: mus = 10
% Run 2: mus = 100
% Run 3: mus = 1000
% Run 4: mus = 10000

%% MCmatlab abbreviations
% G: Geometry, MC: Monte Carlo, FMC: Fluorescence Monte Carlo, HS: Heat
% simulation, M: Media array, FR: Fluence rate, FD: Fractional damage.

%% Geometry definition
clc;
MCmatlab.closeMCmatlabFigures();
model = MCmatlab.model;

%% Estimation for absorption coefficient experimentally
ConcentrationWater = 0.8; 
ConcentrationGel   = 0.1;
ConcentrationIL    = 0.1;

% Note: 28.8 cm^-1 is very high absorption (typical for water at ~1450nm).
% This means light will not penetrate very deep!
ExtinctionWater    = 28.8; 
ExtinctionGel      = 0; 
ExtinctionIL       = 0;

% This calculates to roughly 23.04 cm^-1
absorptionCoeffExperimental = ConcentrationWater*ExtinctionWater + ConcentrationGel*ExtinctionGel + ConcentrationIL*ExtinctionIL;

% Define the scattering values we want to loop over
musValues = [10, 100, 1000, 10000]; 

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

%% Monte Carlo simulation Setup
model.MC.simulationTimeRequested  = .1; 
model.MC.matchedInterfaces        = true; 
model.MC.boundaryType             = 1; 
model.MC.wavelength               = 532; 
model.MC.lightSource.sourceType   = 0; 

% Source position
model.MC.lightSource.xFocus       = 0; 
model.MC.lightSource.yFocus       = 0; 
model.MC.lightSource.zFocus       = 0; 
model.MC.lightSource.theta        = 0; 
model.MC.lightSource.phi          = 0; 

%% Simulation Loop 
for i = 1:length(musValues)
    currentMus = musValues(i);
    
    fprintf('------------------------------------------------\n');
    fprintf('Run %d/%d: Mus = %d, Mua = %.2f\n', i, length(musValues), currentMus, absorptionCoeffExperimental);
    fprintf('------------------------------------------------\n');
    
    % --- UPDATE HERE ---
    % We now pass TWO parameters in the cell array:
    % 1. The Dynamic Scattering (currentMus)
    % 2. The Calculated Absorption (absorptionCoeffExperimental)
    model.G.mediaPropParams = {currentMus, absorptionCoeffExperimental};
    
    % Run
    model = runMonteCarlo(model);
    
    % Plot
    model = plot(model,'MC');
    drawnow; 
end

%% Geometry function
function M = geometryDefinition(X,Y,Z,parameters)
  zSurface = 0.0;
  M = ones(size(X)); % Air
  M(Z > zSurface) = 2; % "Standard" tissue
end

%% Media Properties function 
function mediaProperties = mediaPropertiesFunc(parameters)
  % --- UNPACK PARAMETERS HERE ---
  
  % Check if parameters are empty (initialization check)
  if isempty(parameters)
      dynamicMus = 100; % Default default
      dynamicMua = 1;   % Default default
  else
      % We unpack them in the same order we packed them in the loop
      dynamicMus = parameters{1}; 
      dynamicMua = parameters{2}; % <--- This is your absorptionCoeffExperimental
  end

  mediaProperties = MCmatlab.mediumProperties;
  
  % Air
  j=1;
  mediaProperties(j).name  = 'air';
  mediaProperties(j).mua   = 1e-8; 
  mediaProperties(j).mus   = 1e-8; 
  mediaProperties(j).g     = 1; 
  
  % Tissue
  j=2;
  mediaProperties(j).name  = 'standard tissue';
  mediaProperties(j).mua   = dynamicMua; % <--- NOW USES YOUR CALCULATION
  mediaProperties(j).mus   = dynamicMus; % <--- NOW USES YOUR LOOP VALUE
  mediaProperties(j).g     = 0.9; 
end
