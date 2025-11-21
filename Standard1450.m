%% Standard tissue example with μs' sweep, 3mm & 7mm SDS

clc;

%% ---------------- OPTICAL PROPERTIES --------------------

% Effective phantom composition ~ 88% water, 10% gel, 2% lipid
f_water = 0.88;
f_gel   = 0.10;
f_lipid = 0.02;

% Water absorption at 1450 nm (cm^-1)
mu_a_water_1450 = 28.8;

% Mix absorption (gel & lipid contributions ~0 at 1450)
mu_a_mix = f_water * mu_a_water_1450 + f_gel * 0 + f_lipid * 0;  % ≈ 25.3 cm^-1
assignin('base','mu_a_mix',mu_a_mix);

% Reduced scattering guesses μs' (cm^-1)
mu_sp_list = [3 5 8 10 12 14 16 18 20 22 25];
g = 0.9;   % anisotropy

fprintf("\n===== Optical Properties (1450 nm) =====\n");
fprintf("μa_mix = %.3f cm^-1\n", mu_a_mix);
fprintf("Sweeping μs' values: "); disp(mu_sp_list);
fprintf("========================================\n");

%% ---------------- GEOMETRY SETUP ----------------------

MCmatlab.closeMCmatlabFigures();
baseModel = MCmatlab.model;

baseModel.G.nx = 101;
baseModel.G.ny = 101;
baseModel.G.nz = 150;

baseModel.G.Lx = 7.5;   % [cm]
baseModel.G.Ly = 5.6;   % [cm]
baseModel.G.Lz = 1.5;   % [cm]

baseModel.G.mediaPropertiesFunc = @mediaPropertiesFunc;
baseModel.G.geomFunc            = @geometryDefinition;

% Plot geometry once (like the original example)
geoModel = baseModel;
geoModel = plot(geoModel,'G');

%% SDS list for detectors (cm): 3 mm & 7 mm
SDS_list = [0.3, 0.7];   % [cm]

%% ---------------- LOOP OVER μs' GUESSES -------------------------

results = zeros(length(mu_sp_list), 4);   % [μs', R_close, R_far, ratio]

for k = 1:length(mu_sp_list)

    mu_sp = mu_sp_list(k);            % μs'
    mu_s  = mu_sp / (1 - g);          % convert to μs
    assignin('base','mu_s_user',mu_s);

    fprintf("\n===========================================\n");
    fprintf(" Running simulations for μs' = %.3f cm^-1 (μs = %.1f)\n", mu_sp, mu_s);
    fprintf("===========================================\n");

    R_close = NaN;
    R_far   = NaN;

    % ---- Inner loop: run for each SDS (3mm, 7mm) ----
    for j = 1:length(SDS_list)

        sds_cm = SDS_list(j);

        % Fresh model for this run
        model = baseModel;

        % Monte Carlo settings
        model.MC.simulationTimeRequested = 1.5;   % longer run time = more photons
        model.MC.matchedInterfaces       = false; % USE real n (air vs gel)
        model.MC.boundaryType            = 1;
        model.MC.wavelength              = 1450;

        % LED model (rectangular Lambertian source)
        emitX_cm = 0.05;
        emitY_cm = 0.05;

        model.MC.lightSource.sourceType = 5;
        model.MC.lightSource.xFocus = 0;
        model.MC.lightSource.yFocus = 0;
        model.MC.lightSource.zFocus = 0.03;

        model.MC.lightSource.focalPlaneIntensityDistribution.XDistr = 0;
        model.MC.lightSource.focalPlaneIntensityDistribution.YDistr = 0;
        model.MC.lightSource.focalPlaneIntensityDistribution.XWidth = emitX_cm/2;
        model.MC.lightSource.focalPlaneIntensityDistribution.YWidth = emitY_cm/2;

        model.MC.lightSource.angularIntensityDistribution.XDistr = 2;
        model.MC.lightSource.angularIntensityDistribution.YDistr = 2;

        model.MC.lightSource.theta = 0;
        model.MC.lightSource.phi   = 0;

        % ---- Single PD collector at this SDS ----
        pd_diam_cm = 0.05;   % enlarged PD (e.g. 0.5 mm) for better statistics

        model.MC.useLightCollector  = true;
        model.MC.lightCollector.f   = Inf;
        model.MC.lightCollector.res = 1;
        model.MC.lightCollector.x   = sds_cm;   % SDS in x
        model.MC.lightCollector.y   = 0.0;
        model.MC.lightCollector.z   = -0.001;   % just above surface

        model.MC.lightCollector.theta = 0;      % facing +z (down into tissue)
        model.MC.lightCollector.phi   = 0;

        model.MC.lightCollector.diam  = pd_diam_cm;
        model.MC.lightCollector.NA    = 1.0;

        fprintf("  SDS = %.1f mm: running MC...\n", sds_cm*10);

        % Run the simulation
        model = runMonteCarlo(model);

        % *** MCmatlab visualization plots (like before) ***
        % This opens the MC plot GUI for each run.
        % Comment this out if it produces too many windows.
        model = plot(model,'MC');

        % Fraction of launched photons reaching this PD
        R_PD = model.MC.nPhotonsCollected / model.MC.nPhotons;
        fprintf("    → R_PD = %.6e\n", R_PD);

        if j == 1
            R_close = R_PD;
        else
            R_far   = R_PD;
        end

    end

    ratio = R_far / R_close;
    fprintf("Result for μs' = %.2f: R_close = %.3e, R_far = %.3e, ratio = %.4f\n", ...
            mu_sp, R_close, R_far, ratio);

    results(k,:) = [mu_sp, R_close, R_far, ratio];

end

%% ---------------- SUMMARY TABLE -------------------------

fprintf("\n=========== SUMMARY (1450 nm) ===========\n");
disp(array2table(results, ...
    'VariableNames',{'mu_s_prime','R_close_3mm','R_far_7mm','FarCloseRatio'}));

%% ---------------- SUMMARY PLOTS -------------------------

mu_s_prime = results(:,1);
R_close_3mm = results(:,2);
R_far_7mm   = results(:,3);
ratio_FR    = results(:,4);

figure;
subplot(3,1,1);
plot(mu_s_prime, R_close_3mm, '-o');
xlabel('\mu_s'' [cm^{-1}]');
ylabel('R_{close} (3 mm)');
title('Close detector reflectance vs \mu_s'' (1450 nm)');
grid on;

subplot(3,1,2);
plot(mu_s_prime, R_far_7mm, '-o');
xlabel('\mu_s'' [cm^{-1}]');
ylabel('R_{far} (7 mm)');
title('Far detector reflectance vs \mu_s'' (1450 nm)');
grid on;

subplot(3,1,3);
plot(mu_s_prime, ratio_FR, '-o');
xlabel('\mu_s'' [cm^{-1}]');
ylabel('R_{far} / R_{close}');
title('Far/Close ratio vs \mu_s'' (1450 nm)');
grid on;

%% ================= GEOMETRY FUNCTION ===========================
function M = geometryDefinition(X,Y,Z,parameters)
    zSurface = 0.03;
    M = ones(size(X));
    M(Z > zSurface) = 2;
end

%% ================= MEDIA PROPERTIES ============================
function mediaProperties = mediaPropertiesFunc(parameters)
    mediaProperties = MCmatlab.mediumProperties;

    mu_a = evalin('base','mu_a_mix');
    mu_s = evalin('base','mu_s_user');

    % Air
    j=1;
    mediaProperties(j).name = 'air';
    mediaProperties(j).mua  = 1e-8;
    mediaProperties(j).mus  = 1e-8;
    mediaProperties(j).g    = 1.0;
    mediaProperties(j).n    = 1.0;

    % Phantom
    j=2;
    mediaProperties(j).name = 'Mixture Phantom 1450nm';
    mediaProperties(j).mua  = mu_a;
    mediaProperties(j).mus  = mu_s;
    mediaProperties(j).g    = 0.9;
    mediaProperties(j).n    = 1.33;
end
