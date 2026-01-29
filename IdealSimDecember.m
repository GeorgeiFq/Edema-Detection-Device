%% ========================================================================
% Idealized Reflectance Simulation vs mus'
% - Pencil beam at gel surface
% - Detector at 3 mm and 7 mm SDS
% - 1450 and 1650 nm (mua from water mix)
% - Larger (idealized) detector area for better stats
% - Fixed base number of photons per run
%
% FIXES INCLUDED:
% (1) Print detected photon COUNTS (Ndet) in addition to fraction
% (2) Compute binomial uncertainty sigma_p and store for error bars
% (3) Flag low-count points (Ndet < Ndet_warn)
% (4) Adaptive re-run for low-count cases to reach a target Ndet,
%     with a hard cap on photons per run to prevent runaway runtimes
% (5) Plot with error bars (percent collected ± 1σ)
%
% LUT + SAVING:
% (6) Build a LUT structure + flat table for all conditions
% (7) Create output folder at a fixed path (your requested directory)
% (8) Save figures (.png + .fig) and save LUT (.mat + .csv)
%
% MODIFICATION (per your request):
% - For 1450 nm at 7 mm SDS ONLY: run a single pass at 2e9 photons
%   (no adaptive stepping / reruns). All other cases remain unchanged.
% ========================================================================

clc; clear;
MCmatlab.closeMCmatlabFigures();
model = MCmatlab.model;

%% ---------------- OUTPUT FOLDER (FIXED PATH) ----------------
baseOutPath = 'C:\Users\georg\Documents\MATLAB\MCmatlab-Release\Paper Data';

timestamp  = datestr(now,'yyyymmdd_HHMMSS');
outFolder  = fullfile(baseOutPath, ['MC_LUT_' timestamp]);

figFolder  = fullfile(outFolder,'figures');
lutFolder  = fullfile(outFolder,'lut');

if ~exist(baseOutPath,'dir'); mkdir(baseOutPath); end
if ~exist(outFolder,'dir');   mkdir(outFolder);   end
if ~exist(figFolder,'dir');   mkdir(figFolder);   end
if ~exist(lutFolder,'dir');   mkdir(lutFolder);   end

fprintf('\n[INFO] Saving outputs to:\n  %s\n\n', outFolder);

%% ---------------- GEOMETRY ----------------
model.G.nx = 101;
model.G.ny = 101;
model.G.nz = 150;

% Gel size (cm): 35 x 35 x 12 mm
model.G.Lx = 3.5;
model.G.Ly = 3.5;
model.G.Lz = 1.2;

model.G.mediaPropertiesFunc = @mediaPropertiesFunc;
model.G.geomFunc            = @geometryDefinition;

%% ---------------- ABSORPTION (FROM MIXTURE) ----------------
ConcentrationWater    = 0.70;
ExtinctionWater1450   = 28.8;  % [cm^-1]
ExtinctionWater1650   = 5.7;   % [cm^-1]

mua1450 = ConcentrationWater * ExtinctionWater1450;  % 23.04 cm^-1
mua1650 = ConcentrationWater * ExtinctionWater1650;  % 4.56  cm^-1

wavelengthList = [1450, 1650];
muaList        = [mua1450, mua1650];
wlLabels       = {'1450 nm','1650 nm'};

%% ---------------- USER-DEFINED mus' VALUES ----------------
g_tissue = 0.9;  % anisotropy (fixed)
%musPrimeList = [0.01 0.02 0.04 0.08 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.5 0.75 1.0 1.2 1.4 1.7 1.85 2.0];
musPrimeList_A = [ ...
  0.005 0.007 0.009 0.011 0.013 ...
  0.015 0.017 0.020 0.024 0.028 ...
  0.032 0.038 0.045 0.055 0.065 ...
  0.075 0.090 0.110 0.140 0.180 ...
];

musPrimeList_B = [ ...
  0.220 0.260 0.300 0.340 0.380 ...
  0.420 0.460 0.500 0.600 0.700 ...
  0.800 0.900 1.000 1.100 1.200 ...
  1.350 1.500 1.650 1.800 2.000 ...
];

nMus         = length(musPrimeList);

%% ---------------- SDS LIST ----------------
SDSlist_cm  = [0.3, 0.7];      % 3 mm and 7 mm
SDSlabels   = {'3 mm','7 mm'};

%% ---------------- STORAGE: FRACTIONS + COUNTS + UNCERTAINTY ----------------
% dims: (nMus x nSDS x nWavelength)
nSDS = numel(SDSlist_cm);
nW   = numel(wavelengthList);

fracCollected  = zeros(nMus, nSDS, nW);
sigmaFrac      = zeros(nMus, nSDS, nW);   % 1-sigma uncertainty on fraction
NdetCollected  = zeros(nMus, nSDS, nW);   % detected photons
Nlaunched      = zeros(nMus, nSDS, nW);   % actual simulated photons (final run)
Nrequested     = zeros(nMus, nSDS, nW);   % requested photons (final run)

%% ---------------- BASE MONTE CARLO SETTINGS ----------------
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
model.MC.lightCollector.f         = 0.1;   % finite focal length
model.MC.lightCollector.diam      = 0.5;   % [cm] 5 mm diameter
model.MC.lightCollector.fieldSize = 0.5;   % [cm] imaged field diameter
model.MC.lightCollector.NA        = 0.22;
model.MC.lightCollector.res       = 50;

% Orientation
model.MC.lightCollector.theta     = 0;
model.MC.lightCollector.phi       = pi/2;

%% ---------------- PHOTON BUDGET / ADAPTIVE SETTINGS ----------------
N_base      = 2e7;      % baseline photons per run
N_cap       = 2e9;      % hard cap (used by adaptive runs)
Ndet_target = 1000;     % target detected photons (adaptive runs)
Ndet_warn   = 500;      % warn below this
maxReruns   = 5;
scaleMax    = 20;

% --- Forced one-pass setting for 1450 nm @ 7 mm ---
N_force_1450_7mm = 2e9;

%% ========================================================================
% MAIN LOOP OVER mus'
% ========================================================================
for m = 1:nMus
    musPrime = musPrimeList(m);                 % μs' [cm^-1]
    mus      = musPrime / (1 - g_tissue);       % convert to μs for MCmatlab
    fprintf('\n================ mus'' = %.3g cm^-1 ================\n', musPrime);

    for sdsIdx = 1:nSDS
        sds_cm   = SDSlist_cm(sdsIdx);
        SDSlabel = SDSlabels{sdsIdx};

        % Set detector position for this SDS
        model.MC.lightCollector.x = 0;
        model.MC.lightCollector.y = -sds_cm;
        model.MC.lightCollector.z = 0;   % at top boundary

        for w = 1:nW
            wl    = wavelengthList(w);
            mua   = muaList(w);
            label = wlLabels{w};

            fprintf('Running %s @ SDS %s, mus'' = %.3g ...\n', ...
                label, SDSlabel, musPrime);

            % Set wavelength + medium properties
            model.MC.wavelength = wl;
            model.G.mediaPropParams = {mus, mua};

            % Identify special case: 1450 nm AND 7 mm SDS
            is1450 = (wl == 1450);
            is7mm  = (abs(sds_cm - 0.7) < 1e-12);
            isForcedCase = is1450 && is7mm;

            if isForcedCase
                % ---- ONE PASS at 2e9, no adaptive reruns ----
                model.MC.nPhotonsRequested = N_force_1450_7mm;
                modelOut = runMonteCarlo(model);

                Ndet      = modelOut.MC.nPhotonsCollected;
                NreqFinal = modelOut.MC.nPhotonsRequested;

                fprintf('  [FORCED] 1450 nm @ 7 mm: single pass at %.2e photons\n', N_force_1450_7mm);
            else
                % ---- Normal adaptive behavior for all other cases ----
                model.MC.nPhotonsRequested = N_base;
                [modelOut, Ndet, NreqFinal] = runWithTargetCounts(model, Ndet_target, N_cap, maxReruns, scaleMax);
            end

            % Fraction and uncertainty (binomial)
            Nsim = modelOut.MC.nPhotons;         % actual simulated photons
            p    = Ndet / Nsim;
            sigma_p = sqrt(max(p*(1-p)/Nsim, 0));

            % Store
            fracCollected(m, sdsIdx, w) = p;
            sigmaFrac(m, sdsIdx, w)     = sigma_p;
            NdetCollected(m, sdsIdx, w) = Ndet;
            Nlaunched(m, sdsIdx, w)     = Nsim;
            Nrequested(m, sdsIdx, w)    = NreqFinal;

            % Print
            fprintf('  -> Ndet = %d out of Nsim = %d (requested %d)\n', Ndet, Nsim, NreqFinal);
            fprintf('  -> Fraction = %.6g (%.6f %%), 1σ = %.3g (%.3g %%)\n', ...
                    p, 100*p, sigma_p, 100*sigma_p);

            if Ndet < Ndet_warn
                fprintf('  [WARN] Low-count point (Ndet < %d). Treat as noisy.\n', Ndet_warn);
            end
        end
    end
end

%% ========================================================================
% BUILD LUT + SAVE TO DISK
% ========================================================================
pct      = 100 * fracCollected;
pctSigma = 100 * sigmaFrac;

% Flatten into a tidy LUT table:
% Columns: musPrime, SDS_mm, wavelength_nm, mua_cm1, frac, fracSigma, pct, pctSigma, Ndet, Nsim, Nrequested
nRows = nMus * nSDS * nW;
musPrime_col  = zeros(nRows,1);
SDSmm_col     = zeros(nRows,1);
wl_col        = zeros(nRows,1);
mua_col       = zeros(nRows,1);
frac_col      = zeros(nRows,1);
sigma_col     = zeros(nRows,1);
pct_col       = zeros(nRows,1);
pctSigma_col  = zeros(nRows,1);
Ndet_col      = zeros(nRows,1);
Nsim_col      = zeros(nRows,1);
Nreq_col      = zeros(nRows,1);

r = 0;
for im = 1:nMus
    for is = 1:nSDS
        for iw = 1:nW
            r = r + 1;
            musPrime_col(r) = musPrimeList(im);
            SDSmm_col(r)    = 10 * SDSlist_cm(is);   % cm -> mm
            wl_col(r)       = wavelengthList(iw);
            mua_col(r)      = muaList(iw);
            frac_col(r)     = fracCollected(im,is,iw);
            sigma_col(r)    = sigmaFrac(im,is,iw);
            pct_col(r)      = pct(im,is,iw);
            pctSigma_col(r) = pctSigma(im,is,iw);
            Ndet_col(r)     = NdetCollected(im,is,iw);
            Nsim_col(r)     = Nlaunched(im,is,iw);
            Nreq_col(r)     = Nrequested(im,is,iw);
        end
    end
end

LUT_table = table( ...
    musPrime_col, SDSmm_col, wl_col, mua_col, ...
    frac_col, sigma_col, pct_col, pctSigma_col, ...
    Ndet_col, Nsim_col, Nreq_col, ...
    'VariableNames', {'musPrime_cm1','SDS_mm','wavelength_nm','mua_cm1', ...
                      'fracCollected','sigmaFrac','pctCollected','pctSigma', ...
                      'Ndet','Nsimulated','Nrequested'});

% Structured LUT (MATLAB-friendly)
LUT = struct();
LUT.meta.timestamp          = timestamp;
LUT.meta.outFolder          = outFolder;
LUT.meta.geometry           = struct('Lx_cm',model.G.Lx,'Ly_cm',model.G.Ly,'Lz_cm',model.G.Lz, ...
                                     'nx',model.G.nx,'ny',model.G.ny,'nz',model.G.nz);
LUT.meta.source             = struct('type','pencil','zFocus_cm',model.MC.lightSource.zFocus);
LUT.meta.collector          = struct('diam_cm',model.MC.lightCollector.diam,'fieldSize_cm',model.MC.lightCollector.fieldSize, ...
                                     'NA',model.MC.lightCollector.NA,'res',model.MC.lightCollector.res);
LUT.meta.boundaryType       = model.MC.boundaryType;
LUT.meta.matchedInterfaces  = model.MC.matchedInterfaces;
LUT.meta.g_tissue           = g_tissue;
LUT.meta.ConcentrationWater = ConcentrationWater;
LUT.meta.N_base             = N_base;
LUT.meta.N_cap              = N_cap;
LUT.meta.Ndet_target        = Ndet_target;
LUT.meta.Ndet_warn          = Ndet_warn;
LUT.meta.maxReruns          = maxReruns;
LUT.meta.scaleMax           = scaleMax;
LUT.meta.N_force_1450_7mm   = N_force_1450_7mm;

LUT.axes.musPrime_cm1   = musPrimeList(:);
LUT.axes.SDS_cm         = SDSlist_cm(:);
LUT.axes.SDS_mm         = 10*SDSlist_cm(:);
LUT.axes.wavelength_nm  = wavelengthList(:);
LUT.axes.mua_cm1        = muaList(:);

LUT.data.fracCollected  = fracCollected;
LUT.data.sigmaFrac      = sigmaFrac;
LUT.data.pctCollected   = pct;
LUT.data.pctSigma       = pctSigma;
LUT.data.Ndet           = NdetCollected;
LUT.data.Nsimulated     = Nlaunched;
LUT.data.Nrequested     = Nrequested;

% Save LUT files
matFile = fullfile(lutFolder, 'LUT_MC_reflectance.mat');
csvFile = fullfile(lutFolder, 'LUT_MC_reflectance.csv');

save(matFile, 'LUT', 'LUT_table');
writetable(LUT_table, csvFile);

fprintf('\n[INFO] Saved LUT:\n  %s\n  %s\n', matFile, csvFile);

%% ========================================================================
% PLOTTING WITH ERROR BARS + SAVE FIGURES
% ========================================================================
% ----- 1450 nm -----
f1 = figure('Color','w');
errorbar(musPrimeList, pct(:,1,1), pctSigma(:,1,1), '-o', 'LineWidth', 2); hold on;
errorbar(musPrimeList, pct(:,2,1), pctSigma(:,2,1), '-s', 'LineWidth', 2);
grid on;
xlabel('\mu_s'' [cm^{-1}]');
ylabel('Photons collected [%]');
title('1450 nm: Collected photons vs \mu_s'' (±1\sigma)');
legend('SDS = 3 mm', 'SDS = 7 mm', 'Location', 'best');

saveas(f1, fullfile(figFolder,'Collected_vs_musPrime_1450nm.png'));
savefig(f1, fullfile(figFolder,'Collected_vs_musPrime_1450nm.fig'));

% ----- 1650 nm -----
f2 = figure('Color','w');
errorbar(musPrimeList, pct(:,1,2), pctSigma(:,1,2), '-o', 'LineWidth', 2); hold on;
errorbar(musPrimeList, pct(:,2,2), pctSigma(:,2,2), '-s', 'LineWidth', 2);
grid on;
xlabel('\mu_s'' [cm^{-1}]');
ylabel('Photons collected [%]');
title('1650 nm: Collected photons vs \mu_s'' (±1\sigma)');
legend('SDS = 3 mm', 'SDS = 7 mm', 'Location', 'best');

saveas(f2, fullfile(figFolder,'Collected_vs_musPrime_1650nm.png'));
savefig(f2, fullfile(figFolder,'Collected_vs_musPrime_1650nm.fig'));

fprintf('[INFO] Saved figures to:\n  %s\n', figFolder);

%% ========================================================================
% OPTIONAL: QUICK TABLE PRINT
% ========================================================================
fprintf('\n\n================ SUMMARY TABLE (percent ± 1σ, and Ndet) ================\n');
for iw = 1:nW
    fprintf('\n--- %s ---\n', wlLabels{iw});
    for is = 1:nSDS
        fprintf('SDS = %s:\n', SDSlabels{is});
        for im = 1:nMus
            fprintf('  mus''=%5.1f: %8.4f ± %8.4f %%   (Ndet=%d, Nsim=%d)\n', ...
                musPrimeList(im), pct(im,is,iw), pctSigma(im,is,iw), ...
                NdetCollected(im,is,iw), Nlaunched(im,is,iw));
        end
    end
end

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

%% ========================================================================
% ADAPTIVE RUNNER: re-run with more photons until target detected counts
% ========================================================================
function [modelOut, Ndet, NrequestedFinal] = runWithTargetCounts(modelIn, NdetTarget, Ncap, maxReruns, scaleMax)
    modelOut = modelIn;

    % First run
    modelOut = runMonteCarlo(modelOut);
    Ndet = modelOut.MC.nPhotonsCollected;
    NrequestedFinal = modelOut.MC.nPhotonsRequested;

    if Ndet >= NdetTarget
        return;
    end

    % Adaptive reruns
    for k = 1:maxReruns
        Nsim = modelOut.MC.nPhotons; % actually simulated
        if Nsim <= 0
            Nsim = modelOut.MC.nPhotonsRequested;
        end

        if Ndet <= 0
            scale = 10;
        else
            scale = ceil(NdetTarget / Ndet);
        end
        scale = max(2, min(scale, scaleMax));

        Nnew = min(Ncap, scale * Nsim);
        if Nnew == Nsim
            break;
        end

        modelOut.MC.nPhotonsRequested = Nnew;
        modelOut = runMonteCarlo(modelOut);

        Ndet = modelOut.MC.nPhotonsCollected;
        NrequestedFinal = modelOut.MC.nPhotonsRequested;

        if Ndet >= NdetTarget
            break;
        end
    end
end
