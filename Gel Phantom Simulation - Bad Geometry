function LUT = GelPhantom_1450_1650
%% GelPhantom_1450_1650
% Monte Carlo simulation of your gel phantom sweeping μs'(1450) and
% (optionally) μa(1450) scaling.
%
% Geometry:
%   λ = 1450 nm:
%     - z = 0             : top of domain (PCB plane / PD plane)
%     - 0 <= z <= 0.10 cm : AIR layer (1 mm gap)
%     - z > 0.10 cm       : GEL phantom
%
%   λ = 1650 nm:
%     - z = 0             : top of domain (PCB/PD plane)
%     - 0 <= z <= 0.02 cm : thin AIR layer (0.2 mm)  ~ nearly flush
%     - z > 0.02 cm       : GEL phantom
%
% Phantom:
%   - 80% water, 10% gelatin, 10% Intralipid (approx.)
%   - μa(λ) obtained from Hale–Querry water absorption data, scaled by 0.8
%     for water volume fraction.
%
% Wavelength-dependent scattering:
%   - Input sweep is μs'(1450) [cm^-1] in mu_sp_list.
%   - μs'(1650) = μs'(1450) * (1650 / 1450)^(-b), with b ≈ 1.1.
%
% μa(1450) sweep hook:
%   - We allow a scaling factor of the Hale–Querry-based μa(1450) via
%     mua1450_scale_list. This lets you see how sensitive the 1450 ratio
%     is to modest μa changes.
%
% Output LUT rows:
%   [ μs'(1450),  mua1450_scale,  (R_far/R_close)_1450,  (R_far/R_close)_1650 ]

    clc;
    MCmatlab.closeMCmatlabFigures();

    %% ===== User-tunable μs'(1450) sweep =====
    % Reduced scattering coefficients μs'(1450 nm) [cm^-1]
    mu_sp_list = [5 10 15 20 25];   % EDIT as needed
    nMu        = numel(mu_sp_list);

    %% ===== μa(1450) scale sweep hook =====
    % Scale factor applied to Hale–Querry-based μa_mix(1450).
    % For example, [0.8 1.0 1.2] will do -20%, nominal, +20%.
    mua1450_scale_list = [1.0];     % EDIT this list to sweep μa(1450)
    nScale             = numel(mua1450_scale_list);

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

    %% ===== Allocate LUT: [μs'(1450), muaScale1450, R1450_ratio, R1650_ratio] =====
    LUT = zeros(nMu * nScale, 4);
    row = 0;

    %% ===== Main sweep loops =====
    for iScale = 1:nScale
        muaScale1450 = mua1450_scale_list(iScale);

        fprintf('\n############################################\n');
        fprintf('  μa(1450) scale = %.3f\n', muaScale1450);
        fprintf('############################################\n');

        for iMu = 1:nMu
            mu_sp_1450 = mu_sp_list(iMu);   % μs'(1450 nm) [cm^-1]

            fprintf('\n==============================\n');
            fprintf(' Sweeping μs''(1450) = %.2f cm^-1\n', mu_sp_1450);
            fprintf('==============================\n');

            %% ----- 1450 nm model (1 mm air gap) ----- 
            model1450 = MCmatlab.model;

            % Geometry
            model1450.G.nx = nx;
            model1450.G.ny = ny;
            model1450.G.nz = nz;
            model1450.G.Lx = Lx;
            model1450.G.Ly = Ly;
            model1450.G.Lz = Lz;

            % Air (0–0.10 cm) + gel below: 1 mm air gap
            model1450.G.geomFunc            = @geometryDefinition_1450_air_gel;
            model1450.G.mediaPropertiesFunc = @mediaPropertiesFunc_1450;
            % Pass μs'(1450) and μa scale factor as parameters
            model1450.G.mediaPropParams     = {mu_sp_1450, muaScale1450};

            if doPlotGeom
                model1450 = plot(model1450,'G');
                title(sprintf('Geometry @ 1450 nm, μs''(1450) = %.2f', mu_sp_1450));
            end

            % MC settings
            model1450.MC.simulationTimeRequested = simTime_min;
            model1450.MC.matchedInterfaces       = false;   % use n from mediaProperties
            model1450.MC.boundaryType            = 1;       % all faces escaping
            model1450.MC.wavelength              = 1450;    % [nm]

            % Light source & collector
            model1450 = configureLightSource(model1450);    % λ-dependent LED pattern
            model1450 = configureLightCollector(model1450); % PD geometry

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

            %% ----- 1650 nm model (thin air, nearly flush) ----- 
            model1650 = MCmatlab.model;

            % Geometry – thin air layer (0–0.02 cm) then gel
            model1650.G.nx = nx;
            model1650.G.ny = ny;
            model1650.G.nz = nz;
            model1650.G.Lx = Lx;
            model1650.G.Ly = Ly;
            model1650.G.Lz = Lz;

            model1650.G.geomFunc            = @geometryDefinition_1650_thin_air_gel;
            model1650.G.mediaPropertiesFunc = @mediaPropertiesFunc_1650;
            % Pass μs'(1450) and μa scale factor (only 1450 uses the latter)
            model1650.G.mediaPropParams     = {mu_sp_1450, muaScale1450};

            if doPlotGeom
                model1650 = plot(model1650,'G');
                title(sprintf('Geometry @ 1650 nm, μs''(1450) = %.2f', mu_sp_1450));
            end

            % MC settings
            model1650.MC.simulationTimeRequested = simTime_min;
            model1650.MC.matchedInterfaces       = false;
            model1650.MC.boundaryType            = 1;
            model1650.MC.wavelength              = 1650;    % [nm]

            % Light source & collector
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

            % Guard against division by zero if extremely low signal
            if R1650(1) > 0
                Rratio1650 = R1650(2) / R1650(1);
            else
                Rratio1650 = NaN;
            end

            %% ----- Store in LUT -----
            row = row + 1;
            LUT(row, :) = [mu_sp_1450, muaScale1450, Rratio1450, Rratio1650];

            fprintf(['=> μs''(1450) = %.2f, scale1450 = %.3f:  ', ...
                     'R1450_far/close = %.3f,  R1650_far/close = %.3f\n'], ...
                    mu_sp_1450, muaScale1450, Rratio1450, Rratio1650);
        end
    end

    %% ===== Print final LUT nicely =====
    fprintf('\n=========== Lookup Table (μs''(1450), μa-scale sweep) ===========\n');
    fprintf('   μs''(1450) [cm^-1]   μaScale1450   (R_far/R_close)_1450   (R_far/R_close)_1650\n');
    for r = 1:row
        fprintf('   %8.3f             %8.3f          %8.3f              %8.3f\n', ...
            LUT(r,1), LUT(r,2), LUT(r,3), LUT(r,4));
    end
    fprintf('=================================================================\n');
end


%% ===== LED-like light source (λ-dependent angular spread) =====
function model = configureLightSource(model)
    % Rectangular LED-like emitter with wavelength-dependent angular spread.
    %
    %   - sourceType = 5: X/Y factorizable rectangular
    %   - 0.3 mm × 0.3 mm emitting area (top-hat in X,Y)
    %   - Angular: cosine (Lambertian-ish) with wavelength-specific half-angle:
    %       * 1450 nm: ~120° viewing angle → half-angle ≈ 60°
    %       * 1650 nm: assume slightly narrower (e.g. ~80° total → half-angle ≈ 40°)

    model.MC.lightSource.sourceType = 5;     % X/Y factorizable rectangular

    % Beam axis centered at (0,0), pointing along +z
    model.MC.lightSource.xFocus = 0;
    model.MC.lightSource.yFocus = 0;
    model.MC.lightSource.zFocus = model.G.Lz/2;  % arbitrary focus depth
    model.MC.lightSource.theta  = 0;
    model.MC.lightSource.phi    = 0;

    % Focal-plane emitting area: 0.3 mm × 0.3 mm → 0.03 cm
    emitX = 0.03;   % [cm]
    emitY = 0.03;   % [cm]
    model.MC.lightSource.focalPlaneIntensityDistribution.XDistr = 1;  % top-hat
    model.MC.lightSource.focalPlaneIntensityDistribution.YDistr = 1;
    model.MC.lightSource.focalPlaneIntensityDistribution.XWidth = emitX;
    model.MC.lightSource.focalPlaneIntensityDistribution.YWidth = emitY;

    % Angular distribution
    lambda = model.MC.wavelength;

    if lambda <= 1500
        % 1450 nm LED: ~120° viewing angle → half-angle ~60°
        halfAngle_deg = 60;
    else
        % 1650 nm LED: assume narrower (e.g. ~80° total → half-angle ~40°)
        halfAngle_deg = 40;
    end

    halfAngle_rad = halfAngle_deg * pi/180;

    model.MC.lightSource.angularIntensityDistribution.XDistr = 2;  % cosine
    model.MC.lightSource.angularIntensityDistribution.YDistr = 2;
    model.MC.lightSource.angularIntensityDistribution.XWidth = halfAngle_rad;
    model.MC.lightSource.angularIntensityDistribution.YWidth = halfAngle_rad;
end


%% ===== Light collector (PD) =====
function model = configureLightCollector(model)
    % Configures the PD centered at (x,0) just above the top surface (z ≈ 0+).
    %
    % APX-NG011SMD (approx):
    %   - Active area ≈ 0.0095 mm² → dia ≈ 110 µm → 0.011 cm
    %   - Angle of half sensitivity: ≈ ±65°
    %     → NA ≈ sin(65°) ≈ 0.91
    %
    % No explicit mechanical aperture is modeled here, only the PD's own NA.

    model.MC.useLightCollector   = true;
    model.MC.lightCollector.f    = Inf;        % PD / single-pixel collector
    model.MC.lightCollector.y    = 0.0;        % [cm], centered in y

    % PD just above the top surface (z = 0).
    model.MC.lightCollector.z    = -1e-4;      % [cm] just above z=0

    % PD diameter = 110 µm → 0.011 cm
    model.MC.lightCollector.diam = 0.011;      % [cm]

    % PD axis normal to surface, with finite NA
    model.MC.lightCollector.theta = 0;         % facing into +z
    model.MC.lightCollector.phi   = 0;

    % Half-sensitivity angle ≈ 65° → NA ≈ sin(65°) ≈ 0.91
    model.MC.lightCollector.NA    = 0.91;

    model.MC.lightCollector.res   = 1;         % single-pixel collector
end


%% ===== Geometry: 1450 nm (air 0–0.10 cm + gel below) =====
function M = geometryDefinition_1450_air_gel(X,~,Z,~)
    % 0 <= z <= 0.10 cm → medium 1 (air)
    % z > 0.10 cm       → medium 2 (gel phantom)
    zSurface1450 = 0.10;       % [cm] air gap thickness (1 mm)
    M = ones(size(X));         % start as air (1)
    M(Z > zSurface1450) = 2;   % below zSurface → gel (2)
end


%% ===== Geometry: 1650 nm (thin air 0–0.02 cm + gel below) =====
function M = geometryDefinition_1650_thin_air_gel(X,~,Z,~)
    % 0 <= z <= 0.02 cm → medium 1 (air, thin layer ~0.2 mm)
    % z > 0.02 cm       → medium 2 (gel phantom)
    zSurface1650 = 0.02;       % [cm] thin air gap
    M = ones(size(X));         % start as air (1)
    M(Z > zSurface1650) = 2;   % below thin layer → gel (2)
end


%% ===== Media properties: 1450 nm (μs'(1450) & μa scale) =====
function mediaProperties = mediaPropertiesFunc_1450(parameters)
    % Optical properties at 1450 nm.
    % parameters{1} = μs'(1450 nm) [cm^-1]
    % parameters{2} = μa(1450) scale factor

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
        mu_sp_1450   = parameters{1};      % μs'(1450) [cm^-1]
        if numel(parameters) >= 2
            muaScale1450 = parameters{2};
        else
            muaScale1450 = 1.0;
        end
    else
        mu_sp_1450   = 15.7;
        muaScale1450 = 1.0;
    end

    % Reduced scattering to μs
    mu_s_1450 = mu_sp_1450 / (1 - g_gel);  % convert μs' -> μs

    % Hale–Querry water μa at 1450 nm, scaled by 80% water fraction
    waterFrac = 0.8;
    mu_a_water_1450 = HQ_mua_water(1450);          % pure water [cm^-1]
    mua_mix_1450    = waterFrac * mu_a_water_1450; % phantom mix (approx)

    % Apply user scale factor (sweep hook)
    mua_mix_1450 = mua_mix_1450 * muaScale1450;

    mediaProperties(j).name = 'gel phantom @1450nm';
    mediaProperties(j).mua  = mua_mix_1450;
    mediaProperties(j).mus  = mu_s_1450;
    mediaProperties(j).g    = g_gel;
    mediaProperties(j).n    = 1.33;      % approximate refractive index of gel/water
end


%% ===== Media properties: 1650 nm (μs'(λ) power law) =====
function mediaProperties = mediaPropertiesFunc_1650(parameters)
    % Optical properties at 1650 nm.
    % parameters{1} = μs'(1450 nm) [cm^-1]
    % parameters{2} = μa(1450) scale factor (unused here, kept for consistency).

    mediaProperties = MCmatlab.mediumProperties;

    % ---- Air (medium 1) ----
    % (Present in the thin-air geometry; also kept for indexing consistency.)
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
        mu_sp_1450 = parameters{1};      % μs'(1450) [cm^-1]
    else
        mu_sp_1450 = 15.7;               % default
    end

    % Power-law exponent for μs'(λ)
    b = 1.1;             % can tweak 0.7–1.5 if needed

    lambda_ref   = 1450;   % [nm]
    lambda_1650  = 1650;   % [nm]

    % Compute μs'(1650) via power law
    mu_sp_1650 = mu_sp_1450 * (lambda_1650 / lambda_ref)^(-b);
    mu_s_1650  = mu_sp_1650 / (1 - g_gel);  % convert μs' -> μs

    % Hale–Querry water μa at 1650 nm, scaled by 80% water fraction
    waterFrac        = 0.8;
    mu_a_water_1650  = HQ_mua_water(1650);          % pure water [cm^-1]
    mua_mix_1650     = waterFrac * mu_a_water_1650; % phantom mix (approx)

    mediaProperties(j).name = 'gel phantom @1650nm';
    mediaProperties(j).mua  = mua_mix_1650;
    mediaProperties(j).mus  = mu_s_1650;
    mediaProperties(j).g    = g_gel;
    mediaProperties(j).n    = 1.33;      % same gel index
end


%% ===== Hale–Querry water absorption helper =====
function muaw = HQ_mua_water(lambda_nm)
    % Returns the pure water absorption coefficient [cm^-1] at the given
    % wavelength lambda_nm (nm), interpolated from Hale–Querry data.
    %
    % You provided:
    %   λ [nm]    μa_water [cm^-1]
    %   1440      28.800
    %   1460      28.400
    %   1480      21.230
    %   1500      17.590
    %   1520      14.050
    %   1540      11.830
    %   1560       9.670
    %   1580       7.950
    %   1600       6.720
    %   1620       5.820
    %   1640       4.980
    %   1660       4.540

    persistent lambda_tab muaw_tab
    if isempty(lambda_tab)
        lambda_tab = [1440 1460 1480 1500 1520 1540 1560 1580 1600 1620 1640 1660];
        muaw_tab   = [28.800 28.400 21.230 17.590 14.050 11.830 9.670 7.950 ...
                      6.720 5.820 4.980 4.540];
    end

    muaw = interp1(lambda_tab, muaw_tab, lambda_nm, 'linear', 'extrap');
end
