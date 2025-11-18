function LUT = GelPhantom_1450_1650
%% GelPhantom_1450_1650
% Monte Carlo simulation of your gel phantom sweeping μs' values.
%
% For each μs' in mu_sp_list, it computes:
%   - λ = 1450 nm: air + gel geometry
%   - λ = 1650 nm: SAME air + gel geometry (for now, for stability)
% and measures:
%   - R_3mm = P_PD / P_incident at SDS = 3 mm
%   - R_7mm = P_PD / P_incident at SDS = 7 mm
%   - R_far/R_close at each wavelength
%
% Output:
%   LUT is an N×3 matrix:
%       [ μs' , (R_far/R_close)_1450 , (R_far/R_close)_1650 ]
%
% Main probe-physics features:
%   - LED emitting area ~0.3 mm × 0.3 mm (top-hat in x,y)
%   - LED angular spread from datasheets: 120° viewing angle → ±60° half-angle
%   - PD = APX-NG011SMD: 110 µm dia active area, ±65° half sensitivity
%       → NA ≈ sin(65°) ≈ 0.91
%   - Air (n=1.0) above gel (n=1.33) with a ~1 mm air layer.

    clc;
    MCmatlab.closeMCmatlabFigures();

    %% ===== User-tunable μs' sweep =====
    % Reduced scattering coefficients μs′ [cm^-1] to sweep over
    mu_sp_list = [5 10 15 20 25];   % EDIT this list as needed

    nMu = numel(mu_sp_list);

    % Lookup table: [mu_sp , Rratio1450 , Rratio1650]
    LUT = zeros(nMu, 3);

    %% ===== Common geometry sizes =====
    nx = 101;      % x bins
    ny = 101;      % y bins
    nz = 150;      % z bins

    Lx = 5.0;      % [cm] x size  (phantom ~5 cm)
    Ly = 7.0;      % [cm] y size  (phantom ~7 cm)
    Lz = 3.0;      % [cm] z size  (phantom ~3 cm thick)

    % Source–detector separations (cm)
    SDS_list_cm = [0.3, 0.7];      % 3 mm and 7 mm

    % Monte Carlo settings (same for all runs)
    simTime_min = 1.0;             % [min] MC simulation time

    % For debugging you can turn plots on/off
    doPlotGeom = false;
    doPlotMC   = false;

    %% ===== Main μs' sweep loop =====
    for iMu = 1:nMu
        mu_sp = mu_sp_list(iMu);        % current μs' [cm^-1]

        fprintf('\n==============================\n');
        fprintf(' Sweeping μs'' = %.2f cm^-1\n', mu_sp);
        fprintf('==============================\n');

        %% ----- 1450 nm model ----- 
        model1450 = MCmatlab.model;

        % Geometry
        model1450.G.nx = nx;
        model1450.G.ny = ny;
        model1450.G.nz = nz;
        model1450.G.Lx = Lx;
        model1450.G.Ly = Ly;
        model1450.G.Lz = Lz;

        % Air layer (0–0.10 cm) + gel below that
        model1450.G.geomFunc            = @geometryDefinition_air_gel;
        model1450.G.mediaPropertiesFunc = @mediaPropertiesFunc_1450;
        model1450.G.mediaPropParams     = {mu_sp};

        if doPlotGeom
            model1450 = plot(model1450,'G');
            title(sprintf('Geometry @ 1450 nm, μs'' = %.2f', mu_sp));
        end

        % MC settings
        model1450.MC.simulationTimeRequested = simTime_min;
        % Use refractive indices in mediaProperties (air vs gel)
        model1450.MC.matchedInterfaces       = false;
        model1450.MC.boundaryType            = 1;     % all faces escaping
        model1450.MC.wavelength              = 1450;  % [nm]
        % model1450.MC.calcNFR               = false;  % optional

        % Light source & collector (common config)
        model1450 = configureLightSource(model1450);
        model1450 = configureLightCollector(model1450);

        % Run for SDS = 3 mm and 7 mm
        R1450 = zeros(size(SDS_list_cm));

        fprintf('--- λ = 1450 nm ---\n');
        for k = 1:numel(SDS_list_cm)
            SDS = SDS_list_cm(k);

            % Set PD x-position (SDS along +x from source at x=0)
            model1450.MC.lightCollector.x = SDS;

            % Run MC
            model1450 = runMonteCarlo(model1450);

            % Store normalized collected power
            R1450(k) = model1450.MC.lightCollector.image;

            fprintf('  SDS = %.1f mm: R = %.3e W/W_incident\n', ...
                    SDS*10, R1450(k));
        end

        if doPlotMC
            model1450 = plot(model1450,'MC');
        end

        Rratio1450 = R1450(2) / R1450(1);   % R_far/R_close

        %% ----- 1650 nm model (same geometry, different μa/μs') -----
        model1650 = MCmatlab.model;

        % Geometry
        model1650.G.nx = nx;
        model1650.G.ny = ny;
        model1650.G.nz = nz;
        model1650.G.Lx = Lx;
        model1650.G.Ly = Ly;
        model1650.G.Lz = Lz;

        % Use the same air+gel geometry for stability/debugging
        model1650.G.geomFunc            = @geometryDefinition_air_gel;
        model1650.G.mediaPropertiesFunc = @mediaPropertiesFunc_1650;
        model1650.G.mediaPropParams     = {mu_sp};

        if doPlotGeom
            model1650 = plot(model1650,'G');
            title(sprintf('Geometry @ 1650 nm, μs'' = %.2f', mu_sp));
        end

        % MC settings
        model1650.MC.simulationTimeRequested = simTime_min;
        model1650.MC.matchedInterfaces       = false;
        model1650.MC.boundaryType            = 1;
        model1650.MC.wavelength              = 1650;  % [nm]
        % model1650.MC.calcNFR               = false;

        % Light source & collector (same probe geometry)
        model1650 = configureLightSource(model1650);
        model1650 = configureLightCollector(model1650);

        % Run for SDS = 3 mm and 7 mm
        R1650 = zeros(size(SDS_list_cm));

        fprintf('--- λ = 1650 nm ---\n');
        for k = 1:numel(SDS_list_cm)
            SDS = SDS_list_cm(k);

            % Set PD x-position
            model1650.MC.lightCollector.x = SDS;

            % Run MC
            model1650 = runMonteCarlo(model1650);

            % Store normalized collected power
            R1650(k) = model1650.MC.lightCollector.image;

            fprintf('  SDS = %.1f mm: R = %.3e W/W_incident\n', ...
                    SDS*10, R1650(k));
        end

        if doPlotMC
            model1650 = plot(model1650,'MC');
        end

        % Guard against exact zeros (if extremely low signal) to avoid NaNs
        if R1650(1) > 0
            Rratio1650 = R1650(2) / R1650(1);
        else
            Rratio1650 = NaN;
        end

        %% ----- Store in LUT -----
        LUT(iMu, :) = [mu_sp, Rratio1450, Rratio1650];

        fprintf(['=> μs'' = %.2f:  R1450_far/close = %.3f,  ', ...
                 'R1650_far/close = %.3f\n'], ...
                mu_sp, Rratio1450, Rratio1650);
    end

    %% ===== Print final LUT nicely =====
    fprintf('\n=========== Lookup Table (μs'' sweep) ===========\n');
    fprintf('   μs'' [cm^-1]   (R_far/R_close)_1450   (R_far/R_close)_1650\n');
    for iMu = 1:nMu
        fprintf('   %8.3f         %8.3f              %8.3f\n', ...
            LUT(iMu,1), LUT(iMu,2), LUT(iMu,3));
    end
    fprintf('=================================================\n');
end


%% ===== Helper: configure LED-like light source =====
function model = configureLightSource(model)
    % Rectangular LED-like emitter with datasheet angular spread.
    %
    %   - sourceType = 5: X/Y factorizable rectangular
    %   - 0.3 mm × 0.3 mm emitting area (top-hat in X,Y)
    %   - Angular: cosine (Lambertian-ish) truncated at ±60°
    %       (datasheets: 120° viewing angle)

    model.MC.lightSource.sourceType = 5;     % X/Y factorizable rectangular

    % Beam axis centered at (0,0), pointing along +z into the phantom
    model.MC.lightSource.xFocus = 0;
    model.MC.lightSource.yFocus = 0;
    model.MC.lightSource.zFocus = model.G.Lz/2;  % arbitrary focus depth
    model.MC.lightSource.theta  = 0;
    model.MC.lightSource.phi    = 0;

    % ---- Focal-plane intensity distribution (FPID): emitting area ----
    % Approximate LED die as 0.3 mm × 0.3 mm → 0.03 cm × 0.03 cm top-hat
    emitX = 0.03;   % [cm]
    emitY = 0.03;   % [cm]
    model.MC.lightSource.focalPlaneIntensityDistribution.XDistr = 1;  % top-hat
    model.MC.lightSource.focalPlaneIntensityDistribution.YDistr = 1;
    model.MC.lightSource.focalPlaneIntensityDistribution.XWidth = emitX;
    model.MC.lightSource.focalPlaneIntensityDistribution.YWidth = emitY;

    % ---- Angular intensity distribution (AID): 120° viewing angle ----
    % Viewing angle 120° → half-angle = 60°
    halfAngle_deg = 60;
    halfAngle_rad = halfAngle_deg * pi/180;

    % Use cosine (Lambertian-ish) in X and Y, truncated at ±60°
    model.MC.lightSource.angularIntensityDistribution.XDistr = 2;  % cosine
    model.MC.lightSource.angularIntensityDistribution.YDistr = 2;
    model.MC.lightSource.angularIntensityDistribution.XWidth = halfAngle_rad;
    model.MC.lightSource.angularIntensityDistribution.YWidth = halfAngle_rad;
end


%% ===== Helper: configure light collector (PD) =====
function model = configureLightCollector(model)
    % Configures a PD centered at (x,0) at the top boundary.
    % x (SDS) is set separately before each run.
    %
    % APX-NG011SMD (from datasheet):
    %   - Active area ≈ 0.0095 mm² → dia ≈ 110 µm → 0.011 cm
    %   - Angle of half sensitivity: ±65°
    %     → NA ≈ sin(65°) ≈ 0.91

    model.MC.useLightCollector   = true;
    model.MC.lightCollector.f    = Inf;        % PD / single-pixel collector
    model.MC.lightCollector.y    = 0.0;        % [cm], centered in y

    % PD placed at the top boundary (z ≈ 0+).
    % MCmatlab tallies photons escaping through this plane.
    model.MC.lightCollector.z    = -1e-4;      % [cm] just above z=0

    % PD diameter = 110 µm → 0.011 cm
    model.MC.lightCollector.diam = 0.011;      % [cm]

    % PD axis normal to surface, with finite numerical aperture
    model.MC.lightCollector.theta = 0;         % facing into +z
    model.MC.lightCollector.phi   = 0;

    % Half-sensitivity angle ≈ 65° → NA ≈ sin(65°) ≈ 0.91
    model.MC.lightCollector.NA    = 0.91;

    model.MC.lightCollector.res   = 1;         % single-pixel collector
end


%% ===== Geometry: air (0–0.10 cm) + gel below =====
function M = geometryDefinition_air_gel(X,Y,Z,parameters) %#ok<INUSD>
    % 0 <= z <= 0.10 cm → medium 1 (air)
    % z > 0.10 cm      → medium 2 (gel phantom)
    zSurface = 0.10;         % [cm] 1 mm air layer
    M = ones(size(X));       % start as air (1)
    M(Z > zSurface) = 2;     % below zSurface → gel (2)
end


%% ===== Media properties: 1450 nm (μs' passed via parameters) =====
function mediaProperties = mediaPropertiesFunc_1450(parameters)
    % Optical properties at 1450 nm
    % parameters{1} = μs' [cm^-1] for the gel phantom.

    mediaProperties = MCmatlab.mediumProperties;

    % ---- Air (medium 1) ----
    j = 1;
    mediaProperties(j).name = 'air';
    mediaProperties(j).mua  = 1e-8;  % [cm^-1]
    mediaProperties(j).mus  = 1e-8;  % [cm^-1]
    mediaProperties(j).g    = 1;
    mediaProperties(j).n    = 1.0;   % refractive index of air

    % ---- Gel phantom @ 1450 nm (medium 2) ----
    j = 2;
    g_gel = 0.9;

    if ~isempty(parameters)
        mu_sp = parameters{1};      % reduced scattering μs' [cm^-1]
    else
        mu_sp = 15.7;               % default value
    end
    mu_s    = mu_sp / (1 - g_gel);  % convert μs' -> μs

    % Your previous absorption estimate at 1450 nm
    mua_mix = 22.9;                 % [cm^-1]

    mediaProperties(j).name = 'gel phantom @1450nm';
    mediaProperties(j).mua  = mua_mix;
    mediaProperties(j).mus  = mu_s;
    mediaProperties(j).g    = g_gel;
    mediaProperties(j).n    = 1.33; % approximate refractive index of gel/water
end


%% ===== Media properties: 1650 nm (μs' passed via parameters) =====
function mediaProperties = mediaPropertiesFunc_1650(parameters)
    % Optical properties at 1650 nm
    % parameters{1} = μs' [cm^-1] for the gel phantom.

    mediaProperties = MCmatlab.mediumProperties;

    % ---- Air (medium 1) ----
    j = 1;
    mediaProperties(j).name = 'air';
    mediaProperties(j).mua  = 1e-8;
    mediaProperties(j).mus  = 1e-8;
    mediaProperties(j).g    = 1;
    mediaProperties(j).n    = 1.0;

    % ---- Gel phantom @ 1650 nm (medium 2) ----
    j = 2;
    g_gel = 0.9;

    if ~isempty(parameters)
        mu_sp = parameters{1};      % reduced scattering μs' [cm^-1]
    else
        mu_sp = 13.0;               % default value
    end
    mu_s    = mu_sp / (1 - g_gel);  % convert μs' -> μs

    % Your previous absorption estimate at 1650 nm
    mua_mix = 3.8;                  % [cm^-1]

    mediaProperties(j).name = 'gel phantom @1650nm';
    mediaProperties(j).mua  = mua_mix;
    mediaProperties(j).mus  = mu_s;
    mediaProperties(j).g    = g_gel;
    mediaProperties(j).n    = 1.33; % same gel index
end
