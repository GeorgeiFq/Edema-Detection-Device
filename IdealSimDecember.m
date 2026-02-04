%% ========================================================================
% MC_LUT_MultiGel.m
%
% Runs MCmatlab LUT sweeps for MULTIPLE gels (e.g., two water fractions)
% in one unattended run, with optional checkpoint saving so you can resume
% if the run is interrupted.
%
% Based on your "Idealized Reflectance Simulation vs mus'" script, with:
%  - Loop over gel water fractions (fw list)
%  - User-specified mus' range (min/max/step)
%  - Per-gel output subfolders (figures + LUT + progress checkpoint)
%  - Optional resume: skips points already completed (doneMask)
%  - Optional periodic progress saves (per mus' step by default)
%
% Keeps your forced one-pass rule:
%  - 1450 nm at 7 mm SDS: single pass at 2e9 photons (no adaptive reruns)
% ========================================================================

clc; clear;
MCmatlab.closeMCmatlabFigures();
model = MCmatlab.model;

%% ===================== USER SETTINGS =====================
% ----- Output -----
baseOutPath   = 'C:\Users\georg\Documents\MATLAB\MCmatlab-Release\Paper Data';

% Set a runID so you can RESUME later by reusing the same runID.
% If empty, a timestamp-based folder is created (harder to resume).
runID         = "MC_LUT_HIRES_GELPAIR";   % e.g., "MC_LUT_HIRES_GELPAIR" or ""

resumeIfExists = true;   % if true and progress files exist, load + continue
saveEveryPoint = false;  % if true, checkpoint after EVERY point (slower)
saveEveryMus   = true;   % if true, checkpoint after each mus' index (recommended)

% ----- Gels to run (you asked for two at once; vector supports any N) -----
% These are treated as "water concentration" for mua = fw * ExtinctionWater.
% Use your true fw values you want to simulate:
gel_fw_list = [0.78, 0.8168];   % example: [0.78, 0.8168]

% ----- mus' range (cm^-1) -----
% You asked for e.g. 0 to 10. NOTE: mus'=0 is not meaningful for transport;
% we auto-bump 0 -> musPrimeEps to avoid divide-by-zero.
musPrimeMin  = 0.0;
musPrimeMax  = 10.0;
musPrimeStep = 0.10;    % high resolution (change as desired)
musPrimeEps  = 1e-4;    % used if min <= 0

% ----- SDS list (cm) -----
SDSlist_cm  = [0.3, 0.7];      % 3 mm and 7 mm
SDSlabels   = {'3 mm','7 mm'};

% ----- Wavelength + absorption model -----
wavelengthList = [1450, 1650];
wlLabels       = {'1450 nm','1650 nm'};

ExtinctionWater1450 = 28.8;  % [cm^-1]
ExtinctionWater1650 = 5.7;   % [cm^-1]

% ----- Scattering anisotropy -----
g_tissue = 0.9;  % fixed

% ----- Monte Carlo settings -----
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

% Collector (your working style, larger for better stats)
model.MC.lightCollector.f         = 0.1;   % finite focal length
model.MC.lightCollector.diam      = 0.5;   % [cm] 5 mm diameter
model.MC.lightCollector.fieldSize = 0.5;   % [cm] imaged field diameter
model.MC.lightCollector.NA        = 0.22;
model.MC.lightCollector.res       = 50;

% Orientation
model.MC.lightCollector.theta     = 0;
model.MC.lightCollector.phi       = pi/2;

% ----- Photon budget / adaptive settings -----
N_base      = 2e7;      % baseline photons per run
N_cap       = 2e9;      % hard cap (adaptive runs)
Ndet_target = 1000;     % target detected photons (adaptive runs)
Ndet_warn   = 500;      % warn below this
maxReruns   = 5;
scaleMax    = 20;

% Forced one-pass setting for 1450 nm @ 7 mm
N_force_1450_7mm = 2e9;

%% ===================== BUILD mus' LIST =====================
musPrimeList = musPrimeMin:musPrimeStep:musPrimeMax;
if isempty(musPrimeList)
    error('musPrimeList is empty. Check musPrimeMin/musPrimeMax/musPrimeStep.');
end
if musPrimeList(1) <= 0
    musPrimeList(1) = musPrimeEps;
    fprintf('[WARN] musPrimeMin <= 0. Using musPrimeList(1)=%.3g cm^-1 to avoid divide-by-zero.\n', musPrimeEps);
end
nMus = numel(musPrimeList);

%% ===================== MASTER OUTPUT FOLDER =====================
if strlength(runID) == 0
    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    outFolder = fullfile(baseOutPath, ['MC_LUT_' timestamp]);
else
    outFolder = fullfile(baseOutPath, char(runID));
end

if ~exist(baseOutPath,'dir'); mkdir(baseOutPath); end
if ~exist(outFolder,'dir');   mkdir(outFolder);   end

fprintf('\n[INFO] Saving outputs to:\n  %s\n', outFolder);

% Optional: log all console output to a diary file
diaryFile = fullfile(outFolder, 'run_log.txt');
diary(diaryFile);
diary on;

%% ===================== GEOMETRY SETUP =====================
model.G.nx = 101;
model.G.ny = 101;
model.G.nz = 150;

% Gel size (cm): 35 x 35 x 12 mm
model.G.Lx = 3.5;
model.G.Ly = 3.5;
model.G.Lz = 1.2;

model.G.mediaPropertiesFunc = @mediaPropertiesFunc;
model.G.geomFunc            = @geometryDefinition;

%% ===================== RUN MULTI-GEL LOOP =====================
nSDS = numel(SDSlist_cm);
nW   = numel(wavelengthList);

for gIdx = 1:numel(gel_fw_list)
    fw = gel_fw_list(gIdx);

    % Derive mua per wavelength
    mua1450 = fw * ExtinctionWater1450;
    mua1650 = fw * ExtinctionWater1650;
    muaList = [mua1450, mua1650];

    % Per-gel folders
    gelTag     = sprintf('gel_fw_%0.4f', fw);
    gelTag     = strrep(gelTag,'.','p');   % safer folder naming
    gelFolder  = fullfile(outFolder, gelTag);

    figFolder  = fullfile(gelFolder,'figures');
    lutFolder  = fullfile(gelFolder,'lut');

    if ~exist(gelFolder,'dir'); mkdir(gelFolder); end
    if ~exist(figFolder,'dir'); mkdir(figFolder); end
    if ~exist(lutFolder,'dir'); mkdir(lutFolder); end

    fprintf('\n============================================================\n');
    fprintf('[GEL] Starting gel fw = %.4f   (mua1450=%.3f, mua1650=%.3f cm^-1)\n', fw, mua1450, mua1650);
    fprintf('      Output: %s\n', gelFolder);
    fprintf('============================================================\n');

    % Progress file (checkpoint)
    progressFile = fullfile(lutFolder, 'progress_checkpoint.mat');

    % Allocate / resume arrays
    if resumeIfExists && exist(progressFile,'file')
        S = load(progressFile);
        state = S.state;

        % Basic compatibility checks (avoid silent mismatch)
        assert(numel(state.axes.musPrime_cm1)==nMus, 'Resume mismatch: musPrimeList length differs.');
        assert(numel(state.axes.SDS_cm)==nSDS,        'Resume mismatch: SDS list differs.');
        assert(numel(state.axes.wavelength_nm)==nW,   'Resume mismatch: wavelength list differs.');

        fracCollected = state.data.fracCollected;
        sigmaFrac     = state.data.sigmaFrac;
        NdetCollected = state.data.Ndet;
        Nlaunched     = state.data.Nsimulated;
        Nrequested    = state.data.Nrequested;
        doneMask      = state.data.doneMask;

        fprintf('[INFO] Resuming from checkpoint: %s\n', progressFile);
        fprintf('[INFO] Completed points: %d / %d\n', nnz(doneMask), numel(doneMask));
    else
        fracCollected  = zeros(nMus, nSDS, nW);
        sigmaFrac      = zeros(nMus, nSDS, nW);
        NdetCollected  = zeros(nMus, nSDS, nW);
        Nlaunched      = zeros(nMus, nSDS, nW);
        Nrequested     = zeros(nMus, nSDS, nW);
        doneMask       = false(nMus, nSDS, nW);
    end

    %% ===================== MAIN LOOP OVER mus' =====================
    for m = 1:nMus
        musPrime = musPrimeList(m);                 % μs' [cm^-1]
        mus      = musPrime / (1 - g_tissue);       % convert to μs for MCmatlab

        fprintf('\n================ mus'' = %.4g cm^-1 (gel fw=%.4f) ================\n', musPrime, fw);

        for sdsIdx = 1:nSDS
            sds_cm   = SDSlist_cm(sdsIdx);
            SDSlabel = SDSlabels{sdsIdx};

            % Set detector position for this SDS
            model.MC.lightCollector.x = 0;
            model.MC.lightCollector.y = -sds_cm;
            model.MC.lightCollector.z = 0;   % at top boundary

            for w = 1:nW
                if doneMask(m, sdsIdx, w)
                    fprintf('Skipping (already done): %s @ SDS %s, mus''=%.4g\n', wlLabels{w}, SDSlabel, musPrime);
                    continue;
                end

                wl    = wavelengthList(w);
                mua   = muaList(w);
                label = wlLabels{w};

                fprintf('Running %s @ SDS %s, mus''=%.4g, gel fw=%.4f ...\n', label, SDSlabel, musPrime, fw);

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
                doneMask(m, sdsIdx, w)      = true;

                % Print
                fprintf('  -> Ndet = %d out of Nsim = %d (requested %d)\n', Ndet, Nsim, NreqFinal);
                fprintf('  -> Fraction = %.6g (%.6f %%), 1σ = %.3g (%.3g %%)\n', ...
                        p, 100*p, sigma_p, 100*sigma_p);

                if Ndet < Ndet_warn
                    fprintf('  [WARN] Low-count point (Ndet < %d). Treat as noisy.\n', Ndet_warn);
                end

                % Real-time checkpoint (optional; slower)
                if saveEveryPoint
                    state = buildStateStruct(fw, musPrimeList, SDSlist_cm, wavelengthList, muaList, ...
                        fracCollected, sigmaFrac, NdetCollected, Nlaunched, Nrequested, doneMask, ...
                        model, g_tissue, ExtinctionWater1450, ExtinctionWater1650, ...
                        N_base, N_cap, Ndet_target, Ndet_warn, maxReruns, scaleMax, N_force_1450_7mm, outFolder);

                    save(progressFile, 'state', '-v7.3');
                end
            end
        end

        % Periodic checkpoint (recommended): save after each mus'
        if saveEveryMus
            state = buildStateStruct(fw, musPrimeList, SDSlist_cm, wavelengthList, muaList, ...
                fracCollected, sigmaFrac, NdetCollected, Nlaunched, Nrequested, doneMask, ...
                model, g_tissue, ExtinctionWater1450, ExtinctionWater1650, ...
                N_base, N_cap, Ndet_target, Ndet_warn, maxReruns, scaleMax, N_force_1450_7mm, outFolder);

            save(progressFile, 'state', '-v7.3');
            fprintf('[INFO] Checkpoint saved: %s\n', progressFile);
        end
    end

    %% ===================== FINAL LUT BUILD + SAVE =====================
    pct      = 100 * fracCollected;
    pctSigma = 100 * sigmaFrac;

    % Flatten into table
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

    % Structured LUT
    LUT = struct();
    LUT.meta.fw                 = fw;
    LUT.meta.outFolder          = gelFolder;
    LUT.meta.geometry           = struct('Lx_cm',model.G.Lx,'Ly_cm',model.G.Ly,'Lz_cm',model.G.Lz, ...
                                         'nx',model.G.nx,'ny',model.G.ny,'nz',model.G.nz);
    LUT.meta.source             = struct('type','pencil','zFocus_cm',model.MC.lightSource.zFocus);
    LUT.meta.collector          = struct('diam_cm',model.MC.lightCollector.diam,'fieldSize_cm',model.MC.lightCollector.fieldSize, ...
                                         'NA',model.MC.lightCollector.NA,'res',model.MC.lightCollector.res);
    LUT.meta.boundaryType       = model.MC.boundaryType;
    LUT.meta.matchedInterfaces  = model.MC.matchedInterfaces;
    LUT.meta.g_tissue           = g_tissue;
    LUT.meta.ExtinctionWater1450= ExtinctionWater1450;
    LUT.meta.ExtinctionWater1650= ExtinctionWater1650;
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
    LUT.data.doneMask       = doneMask;

    % Save LUT files
    matFile = fullfile(lutFolder, 'LUT_MC_reflectance.mat');
    csvFile = fullfile(lutFolder, 'LUT_MC_reflectance.csv');

    save(matFile, 'LUT', 'LUT_table', '-v7.3');
    writetable(LUT_table, csvFile);

    fprintf('\n[INFO] Saved FINAL LUT for gel fw=%.4f:\n  %s\n  %s\n', fw, matFile, csvFile);

    %% ===================== PLOTTING (per gel) =====================
    % 1450 nm plot
    f1 = figure('Color','w');
    errorbar(musPrimeList, 100*fracCollected(:,1,1), 100*sigmaFrac(:,1,1), '-o', 'LineWidth', 2); hold on;
    errorbar(musPrimeList, 100*fracCollected(:,2,1), 100*sigmaFrac(:,2,1), '-s', 'LineWidth', 2);
    grid on;
    xlabel('\mu_s'' [cm^{-1}]');
    ylabel('Photons collected [%]');
    title(sprintf('1450 nm: Collected photons vs \\mu_s'' (gel f_w=%.4f)', fw));
    legend(sprintf('SDS = %s', SDSlabels{1}), sprintf('SDS = %s', SDSlabels{2}), 'Location', 'best');

    saveas(f1, fullfile(figFolder, sprintf('Collected_vs_musPrime_1450nm_fw_%s.png', gelTag)));
    savefig(f1, fullfile(figFolder, sprintf('Collected_vs_musPrime_1450nm_fw_%s.fig', gelTag)));

    % 1650 nm plot
    f2 = figure('Color','w');
    errorbar(musPrimeList, 100*fracCollected(:,1,2), 100*sigmaFrac(:,1,2), '-o', 'LineWidth', 2); hold on;
    errorbar(musPrimeList, 100*fracCollected(:,2,2), 100*sigmaFrac(:,2,2), '-s', 'LineWidth', 2);
    grid on;
    xlabel('\mu_s'' [cm^{-1}]');
    ylabel('Photons collected [%]');
    title(sprintf('1650 nm: Collected photons vs \\mu_s'' (gel f_w=%.4f)', fw));
    legend(sprintf('SDS = %s', SDSlabels{1}), sprintf('SDS = %s', SDSlabels{2}), 'Location', 'best');

    saveas(f2, fullfile(figFolder, sprintf('Collected_vs_musPrime_1650nm_fw_%s.png', gelTag)));
    savefig(f2, fullfile(figFolder, sprintf('Collected_vs_musPrime_1650nm_fw_%s.fig', gelTag)));

    fprintf('[INFO] Saved figures to:\n  %s\n', figFolder);

    %% ===================== FINAL CHECKPOINT (done) =====================
    state = buildStateStruct(fw, musPrimeList, SDSlist_cm, wavelengthList, muaList, ...
        fracCollected, sigmaFrac, NdetCollected, Nlaunched, Nrequested, doneMask, ...
        model, g_tissue, ExtinctionWater1450, ExtinctionWater1650, ...
        N_base, N_cap, Ndet_target, Ndet_warn, maxReruns, scaleMax, N_force_1450_7mm, outFolder);

    save(progressFile, 'state', '-v7.3');
    fprintf('[INFO] Final checkpoint saved: %s\n', progressFile);
end

fprintf('\n[INFO] All gel runs complete. Master output:\n  %s\n', outFolder);
diary off;

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

%% ========================================================================
% Build checkpoint state struct (kept simple, used for resume + provenance)
% ========================================================================
function state = buildStateStruct(fw, musPrimeList, SDSlist_cm, wavelengthList, muaList, ...
    fracCollected, sigmaFrac, NdetCollected, Nlaunched, Nrequested, doneMask, ...
    model, g_tissue, Ext1450, Ext1650, N_base, N_cap, Ndet_target, Ndet_warn, maxReruns, scaleMax, N_force_1450_7mm, outFolder)

    state = struct();

    state.meta.fw                  = fw;
    state.meta.outFolder           = outFolder;
    state.meta.timestamp           = datestr(now,'yyyymmdd_HHMMSS');
    state.meta.g_tissue            = g_tissue;
    state.meta.ExtinctionWater1450 = Ext1450;
    state.meta.ExtinctionWater1650 = Ext1650;

    state.meta.N_base              = N_base;
    state.meta.N_cap               = N_cap;
    state.meta.Ndet_target         = Ndet_target;
    state.meta.Ndet_warn           = Ndet_warn;
    state.meta.maxReruns           = maxReruns;
    state.meta.scaleMax            = scaleMax;
    state.meta.N_force_1450_7mm    = N_force_1450_7mm;

    state.meta.geometry = struct('Lx_cm',model.G.Lx,'Ly_cm',model.G.Ly,'Lz_cm',model.G.Lz, ...
                                 'nx',model.G.nx,'ny',model.G.ny,'nz',model.G.nz);
    state.meta.source   = struct('type','pencil','zFocus_cm',model.MC.lightSource.zFocus);
    state.meta.collector= struct('diam_cm',model.MC.lightCollector.diam, ...
                                 'fieldSize_cm',model.MC.lightCollector.fieldSize, ...
                                 'NA',model.MC.lightCollector.NA,'res',model.MC.lightCollector.res);

    state.axes.musPrime_cm1    = musPrimeList(:);
    state.axes.SDS_cm          = SDSlist_cm(:);
    state.axes.wavelength_nm   = wavelengthList(:);
    state.axes.mua_cm1         = muaList(:);

    state.data.fracCollected   = fracCollected;
    state.data.sigmaFrac       = sigmaFrac;
    state.data.Ndet            = NdetCollected;
    state.data.Nsimulated      = Nlaunched;
    state.data.Nrequested      = Nrequested;
    state.data.doneMask        = doneMask;
end
