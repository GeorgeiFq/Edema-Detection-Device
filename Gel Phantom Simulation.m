function GelPhantom_1450_1650
%% GelPhantom_1450_1650
% Monte Carlo simulation of your gel phantom for:
%   - λ = 1450 nm: LED + PD ~1 mm above the gel (1 mm air layer)
%   - λ = 1650 nm: LED flush with the gel (no explicit air layer)
%
% For each wavelength, we compute:
%   - R_3mm = P_PD / P_incident at SDS = 3 mm
%   - R_7mm = P_PD / P_incident at SDS = 7 mm
%   - R_far / R_close ratio

    clc;
    MCmatlab.closeMCmatlabFigures();

    %% ===== Common geometry sizes =====
    nx = 101;      % x bins
    ny = 101;      % y bins
    nz = 150;      % z bins

    Lx = 5.0;      % [cm] x size  (phantom ~5 cm)
    Ly = 7.0;      % [cm] y size  (phantom ~7 cm)
    Lz = 3.0;      % [cm] z size  (phantom ~3 cm thick)

    %% ===== 1450 nm model: source + PD 1 mm above gel =====
    model1450 = MCmatlab.model;

    % Geometry
    model1450.G.nx = nx;
    model1450.G.ny = ny;
    model1450.G.nz = nz;
    model1450.G.Lx = Lx;
    model1450.G.Ly = Ly;
    model1450.G.Lz = Lz;

    % Air layer (0–0.10 cm) + gel below that
    model1450.G.geomFunc            = @geometryDefinition_1450;
    model1450.G.mediaPropertiesFunc = @mediaPropertiesFunc_1450;

    % Quick geometry plot (optional)
    model1450 = plot(model1450,'G');
    title('Geometry @ 1450 nm (air + gel)');

    % MC settings
    model1450.MC.simulationTimeRequested = 1.0;   % [min]
    model1450.MC.matchedInterfaces       = true;  % all n = 1
    model1450.MC.boundaryType            = 1;     % all faces escaping
    model1450.MC.wavelength              = 1450;  % [nm]

    % Light source: rectangular LED-like emitter
    model1450 = configureLightSource(model1450);

    % Light collector: PD ~1 mm above gel
    model1450 = configureLightCollector(model1450);

    % Run for SDS = 3 mm and 7 mm
    SDS_list_cm = [0.3, 0.7];          % [cm] 3 mm and 7 mm
    R1450       = zeros(size(SDS_list_cm));

    fprintf('\n=== Wavelength = 1450 nm (LED + PD ~1 mm above gel) ===\n');
    for k = 1:numel(SDS_list_cm)
        SDS = SDS_list_cm(k);

        % Set PD x-position (SDS along +x from source at x=0)
        model1450.MC.lightCollector.x = SDS;

        % Run MC
        model1450 = runMonteCarlo(model1450);

        % Store normalized collected power
        R1450(k) = model1450.MC.lightCollector.image;

        fprintf('SDS = %.1f mm: collected_norm = %.3e W/W_incident\n', ...
                SDS*10, R1450(k));
    end

    % Optional: MC plot for last 1450 run
    model1450 = plot(model1450,'MC');

    %% ===== 1650 nm model: source flush with gel (no air layer) =====
    model1650 = MCmatlab.model;

    % Geometry
    model1650.G.nx = nx;
    model1650.G.ny = ny;
    model1650.G.nz = nz;
    model1650.G.Lx = Lx;
    model1650.G.Ly = Ly;
    model1650.G.Lz = Lz;

    % Gel extends all the way to the top (z = 0) → source flush with gel
    model1650.G.geomFunc            = @geometryDefinition_1650;
    model1650.G.mediaPropertiesFunc = @mediaPropertiesFunc_1650;

    % Quick geometry plot (optional)
    model1650 = plot(model1650,'G');
    title('Geometry @ 1650 nm (gel to top)');

    % MC settings
    model1650.MC.simulationTimeRequested = 1.0;   % [min]
    model1650.MC.matchedInterfaces       = true;
    model1650.MC.boundaryType            = 1;
    model1650.MC.wavelength              = 1650;  % [nm]

    % Light source: same LED-like emitter configuration,
    % but now it is effectively at the gel surface.
    model1650 = configureLightSource(model1650);

    % Light collector: keep just above top surface (now gel surface)
    model1650 = configureLightCollector(model1650);

    % Run for SDS = 3 mm and 7 mm
    R1650 = zeros(size(SDS_list_cm));

    fprintf('\n=== Wavelength = 1650 nm (LED flush with gel) ===\n');
    for k = 1:numel(SDS_list_cm)
        SDS = SDS_list_cm(k);

        % Set PD x-position
        model1650.MC.lightCollector.x = SDS;

        % Run MC
        model1650 = runMonteCarlo(model1650);

        % Store normalized collected power
        R1650(k) = model1650.MC.lightCollector.image;

        fprintf('SDS = %.1f mm: collected_norm = %.3e W/W_incident\n', ...
                SDS*10, R1650(k));
    end

    % Optional: MC plot for last 1650 run
    model1650 = plot(model1650,'MC');

    %% ===== Summary =====
    fprintf('\n===== Summary (P_{PD} / P_{incident}) =====\n');
    fprintf('λ = 1450 nm:  R_3mm = %.3e,  R_7mm = %.3e,  R_far/R_close = %.3f\n', ...
        R1450(1), R1450(2), R1450(2)/R1450(1));
    fprintf('λ = 1650 nm:  R_3mm = %.3e,  R_7mm = %.3e,  R_far/R_close = %.3f\n', ...
        R1650(1), R1650(2), R1650(2)/R1650(1));
end


%% ===== Helper: configure LED-like light source =====
function model = configureLightSource(model)
    % Rectangular LED-like emitter with Lambertian angular distribution.

    model.MC.lightSource.sourceType = 5;     % X/Y factorizable rectangular

    % Beam axis at (0,0), pointing along +z into the phantom
    model.MC.lightSource.xFocus = 0;
    model.MC.lightSource.yFocus = 0;
    model.MC.lightSource.zFocus = model.G.Lz/2;  % arbitrary focus depth
    model.MC.lightSource.theta  = 0;
    model.MC.lightSource.phi    = 0;

    % Focal-plane intensity distribution (FPID):
    % uniform 0.3 mm × 0.3 mm emitting area
    emitX = 0.03;   % [cm]
    emitY = 0.03;   % [cm]
    model.MC.lightSource.focalPlaneIntensityDistribution.XDistr = 1;  % top-hat
    model.MC.lightSource.focalPlaneIntensityDistribution.YDistr = 1;
    model.MC.lightSource.focalPlaneIntensityDistribution.XWidth = emitX;
    model.MC.lightSource.focalPlaneIntensityDistribution.YWidth = emitY;

    % Angular intensity distribution (AID): Lambertian (cosine) in X and Y
    model.MC.lightSource.angularIntensityDistribution.XDistr = 2;  % cosine
    model.MC.lightSource.angularIntensityDistribution.YDistr = 2;
end


%% ===== Helper: configure light collector (PD) =====
function model = configureLightCollector(model)
    % Configures a 110 µm PD centered at (x,0) on top, with z just above surface.
    % x (SDS) is set separately before each run.

    model.MC.useLightCollector   = true;
    model.MC.lightCollector.f    = Inf;        % PD / fiber tip mode
    model.MC.lightCollector.y    = 0.0;        % [cm]

    % PD just above the top surface z=0
    % For 1450 nm: top is air, gel starts at z=0.10 → PD ~1 mm above gel.
    % For 1650 nm: top is gel → PD effectively just above gel surface.
    model.MC.lightCollector.z    = -1e-4;      % [cm] just above z=0

    model.MC.lightCollector.diam = 0.05;      % [cm] 110 µm -> 0.011
    model.MC.lightCollector.theta = 0;         % facing into +z
    model.MC.lightCollector.phi   = 0;
    model.MC.lightCollector.NA    = 1.0;
    model.MC.lightCollector.res   = 1;         % single-pixel collector
end


%% ===== Geometry: 1450 nm (air + gel) =====
function M = geometryDefinition_1450(X,Y,Z,parameters)
    % 0 <= z <= 0.10 cm → medium 1 (air)
    % z > 0.10 cm      → medium 2 (gel phantom)
    zSurface = 0.10;         % [cm] 1 mm air layer
    M = ones(size(X));       % start as air (1)
    M(Z > zSurface) = 2;     % below zSurface → gel (2)
end


%% ===== Geometry: 1650 nm (gel to top) =====
function M = geometryDefinition_1650(X,Y,Z,parameters)
    % No explicit air layer: gel occupies almost entire cuboid.
    % For all z > 0, set medium = gel.
    M = ones(size(X));       % start as "air"
    M(Z > 0) = 2;            % everything inside cuboid → gel (2)
end


%% ===== Media properties: 1450 nm =====
function mediaProperties = mediaPropertiesFunc_1450(parameters)
    % Optical properties at 1450 nm
    % Medium 1: air (μa ≈ μs ≈ 0)
    % Medium 2: 80% water, 10% gelatin, 10% IL phantom

    mediaProperties = MCmatlab.mediumProperties;

    % Air
    j = 1;
    mediaProperties(j).name = 'air';
    mediaProperties(j).mua  = 1e-8;  % [cm^-1]
    mediaProperties(j).mus  = 1e-8;  % [cm^-1]
    mediaProperties(j).g    = 1;

    % Gel phantom @ 1450 nm
    j = 2;
    g_gel   = 0.9;
    mu_sp   = 15.7;      % [cm^-1] reduced scattering μs'
    mu_s    = mu_sp / (1 - g_gel);
    mua_mix = 22.9;      % [cm^-1] absorption μa

    mediaProperties(j).name = 'gel phantom @1450nm';
    mediaProperties(j).mua  = mua_mix;
    mediaProperties(j).mus  = mu_s;
    mediaProperties(j).g    = g_gel;
end


%% ===== Media properties: 1650 nm =====
function mediaProperties = mediaPropertiesFunc_1650(parameters)
    % Optical properties at 1650 nm
    % Medium 1: air (outside cuboid)
    % Medium 2: gel phantom with lower μa

    mediaProperties = MCmatlab.mediumProperties;

    % "Air" (unused inside cuboid since geometry puts gel everywhere z>0,
    % but kept for consistency)
    j = 1;
    mediaProperties(j).name = 'air';
    mediaProperties(j).mua  = 1e-8;
    mediaProperties(j).mus  = 1e-8;
    mediaProperties(j).g    = 1;

    % Gel phantom @ 1650 nm
    j = 2;
    g_gel   = 0.9;
    mu_sp   = 13.0;     % [cm^-1] reduced scattering μs' (example)
    mu_s    = mu_sp / (1 - g_gel);
    mua_mix = 3.8;      % [cm^-1] lower absorption at 1650 nm

    mediaProperties(j).name = 'gel phantom @1650nm';
    mediaProperties(j).mua  = mua_mix;
    mediaProperties(j).mus  = mu_s;
    mediaProperties(j).g    = g_gel;
end
