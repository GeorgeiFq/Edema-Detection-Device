%% ========================================================================
% Idealized Reflectance Simulation vs mus'
% - Pencil beam at gel surface
% - Detector at 3 mm and 7 mm SDS
% - 1450 and 1650 nm (mua from water mix)
% - Larger (idealized) detector area for better stats
% - Fixed number of photons per run
% - NOW: store and plot % photons collected vs mus' for each LED
% ========================================================================

clc; clear;
MCmatlab.closeMCmatlabFigures();
model = MCmatlab.model;

%% ---------------- GEOMETRY ----------------
model.G.nx = 101;
model.G.ny = 101;
model.G.nz = 150;

% Gel size (cm): 75 x 56 x 15 mm
model.G.Lx = 7.5;
model.G.Ly = 5.6;
model.G.Lz = 1.5;

model.G.mediaPropertiesFunc = @mediaPropertiesFunc;
model.G.geomFunc            = @geometryDefinition;

% Don't call plot(model,'G') before mediaPropParams is set

%% ---------------- ABSORPTION (FROM MIXTURE) ----------------
ConcentrationWater    = 0.8;
ExtinctionWater1450   = 28.8;  % [cm^-1]
ExtinctionWater1650   = 5.7;   % [cm^-1]

mua1450 = ConcentrationWater * ExtinctionWater1450;  % 23.04 cm^-1
mua1650 = ConcentrationWater * ExtinctionWater1650;  % 4.56 cm^-1

wavelengthList = [1450, 1650];
muaList        = [mua1450, mua1650];
wlLabels       = {'1450 nm','1650 nm'};

%% ---------------- USER-DEFINED mus' VALUES ----------------
musList = [3 5 7 10 15 20 25 30 40 50];   % [cm^-1] tweak as you like
nMus    = length(musList);

%% ---------------- STORAGE FOR FRACTIONS ----------------
% Fractions (not %) collected for each case, indexed by mus
frac1450_SDS3 = zeros(1, nMus);
frac1450_SDS7 = zeros(1, nMus);
frac1650_SDS3 = zeros(1, nMus);
frac1650_SDS7 = zeros(1, nMus);

%% ---------------- BASE MONTE CARLO SETTINGS ----------------
% MCmatlab wants positive time; we'll also cap by photon count
model.MC.simulationTimeRequested = 1;      % [min] must be > 0
model.MC.nPhotonsRequested       = 2e7;    % photons per run

model.MC.matchedInterfaces       = true;
model.MC.boundaryType            = 1;

% Pencil beam at gel surface (slightly inside tissue)
model.MC.lightSource.sourceType = 0;    % Pencil beam
model.MC.lightSource.xFocus = 0;
model.MC.lightSource.yFocus = 0;
model.MC.lightSource.zFocus = 0.03 + 1e-4;  % just below air–tissue interface
model.MC.lightSource.theta  = 0;
model.MC.lightSource.phi    = 0;

model.MC.useLightCollector = true;

% --------- COLLECTOR: working style, just larger ---------
model.MC.lightCollector.f         = 0.1;   % finite focal length (like example)
model.MC.lightCollector.diam      = 0.5;   % [cm] 5 mm diameter (bigger than APX)
model.MC.lightCollector.fieldSize = 0.5;   % [cm] imaged field diameter
model.MC.lightCollector.NA        = 0.22;  % same as example
model.MC.lightCollector.res       = 50;    % grid on detector plane

% Orientation: same as working example
model.MC.lightCollector.theta     = 0;     % [rad]
model.MC.lightCollector.phi       = pi/2;  % [rad]

%% ---------------- SDS LIST ----------------
SDSlist_cm  = [0.3, 0.7];      % 3 mm and 7 mm
SDSlabels   = {'3 mm','7 mm'};

%% ========================================================================
% MAIN LOOP OVER mus'
% ========================================================================
for m = 1:nMus
    mus = musList(m);
    fprintf('\n================ mus'' = %.1f cm^-1 ================\n', mus);

    % For each SDS (3 mm and 7 mm)
    for sdsIdx = 1:numel(SDSlist_cm)
        sds_cm   = SDSlist_cm(sdsIdx);
        SDSlabel = SDSlabels{sdsIdx};

        % Set detector position for this SDS
        model.MC.lightCollector.x = 0;
        model.MC.lightCollector.y = -sds_cm;
        model.MC.lightCollector.z = 0;   % at top boundary (air side)

        % For each wavelength (1450 & 1650 nm)
        for w = 1:2
            wl    = wavelengthList(w);
            mua   = muaList(w);
            label = wlLabels{w};

            fprintf('Running %s @ SDS %s, mus'' = %.1f ...\n', ...
                    label, SDSlabel, mus);

            % Set wavelength
            model.MC.wavelength = wl;

            % Pass mus and mua to mediaPropertiesFunc
            model.G.mediaPropParams = {mus, mua};

            % Run Monte Carlo
            model = runMonteCarlo(model);

            % Fraction of launched photon packets that hit the detector
            fracCollected = model.MC.nPhotonsCollected / model.MC.nPhotons;

            fprintf('  -> Fraction collected (%s, SDS %s): %.6g (%.6f %%)\n', ...
                    label, SDSlabel, fracCollected, fracCollected * 100);

            % -------- STORE FOR PLOTTING --------
            if w == 1  % 1450 nm
                if sdsIdx == 1      % 3 mm
                    frac1450_SDS3(m) = fracCollected;
                else                % 7 mm
                    frac1450_SDS7(m) = fracCollected;
                end
            else       % 1650 nm
                if sdsIdx == 1      % 3 mm
                    frac1650_SDS3(m) = fracCollected;
                else                % 7 mm
                    frac1650_SDS7(m) = fracCollected;
                end
            end
        end
    end
end

%% ========================================================================
% PLOTTING: % photons collected vs mus' for each LED
% ========================================================================

% Convert to percent
pct1450_SDS3 = frac1450_SDS3 * 100;
pct1450_SDS7 = frac1450_SDS7 * 100;
pct1650_SDS3 = frac1650_SDS3 * 100;
pct1650_SDS7 = frac1650_SDS7 * 100;

% ----- 1450 nm -----
figure;
plot(musList, pct1450_SDS3, '-o', 'LineWidth', 2); hold on;
plot(musList, pct1450_SDS7, '-s', 'LineWidth', 2);
grid on;
xlabel('\mu_s'' [cm^{-1}]');
ylabel('Photons collected [%]');
title('1450 nm: Collected photons vs \mu_s''');
legend('SDS = 3 mm', 'SDS = 7 mm', 'Location', 'best');

% ----- 1650 nm -----
figure;
plot(musList, pct1650_SDS3, '-o', 'LineWidth', 2); hold on;
plot(musList, pct1650_SDS7, '-s', 'LineWidth', 2);
grid on;
xlabel('\mu_s'' [cm^{-1}]');
ylabel('Photons collected [%]');
title('1650 nm: Collected photons vs \mu_s''');
legend('SDS = 3 mm', 'SDS = 7 mm', 'Location', 'best');

%% ========================================================================
% GEOMETRY FUNCTION
% ========================================================================
function M = geometryDefinition(X,Y,Z,parameters)
  zSurface = 0.03;      % [cm] air–tissue boundary
  M = ones(size(X));    % 1 = air
  M(Z > zSurface) = 2;  % 2 = tissue
end

%% ========================================================================
% MEDIA PROPERTIES FUNCTION
% ========================================================================
function mediaProperties = mediaPropertiesFunc(parameters)
  mediaProperties = MCmatlab.mediumProperties;

  % Expect parameters = {mus, mua}
  mus = parameters{1};
  mua = parameters{2};

  % Medium 1: air
  j = 1;
  mediaProperties(j).name = 'air';
  mediaProperties(j).mua  = 1e-8;
  mediaProperties(j).mus  = 1e-8;
  mediaProperties(j).g    = 1;

  % Medium 2: tissue
  j = 2;
  mediaProperties(j).name = 'tissue';
  mediaProperties(j).mua  = mua;
  mediaProperties(j).mus  = mus;
  mediaProperties(j).g    = 0.9;
end
