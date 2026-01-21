%% ========================================================================
% Reflectance LUT vs mus' WITH STL APERTURE (TUNNEL-ONLY) + HIGH-STATISTICS
% + Saves Fluence Rate + Absorbed Power plots (auto-export from MCmatlab)
%
% Geometry intent:
% - Plate (square) is flush with SENSOR plane at z = 0 (top of cuboid).
% - Tunnel extends DOWN through air and reaches the gel surface.
% - Gel begins at z = zGelTop (bottom of tunnel) and extends downward.
% - Above the gel, EVERYTHING is absorbing plastic EXCEPT the tunnel (global baffle).
%   => Collector only “sees” photons that exit gel into the tunnel.
%
% Notes:
% - MCmatlab convention: z=0 at TOP, z increases DOWN.
% - Collector MUST be outside cuboid (z < 0).
% - This script uses an adaptive photon budget per (mus', SDS, wl) to hit
%   requested statistical quality (relative SE + min Ndet), up to N_cap.
% - A single “MC plots” run is performed (configurable) and we automatically
%   save the Fluence Rate and Absorbed Power figures to disk.
% ========================================================================

clc; clear;
MCmatlab.closeMCmatlabFigures();
rng(1);

%% ---------------- USER PATHS ----------------
stlPath = 'C:\Users\georg\Documents\MATLAB\MCmatlab-Release\Aperture v2.stl';   % <-- CHANGE
voxeliseFolder = 'C:\Users\georg\Documents\MATLAB\MCmatlab-Release\Mesh_voxelisation'; % <-- CHANGE
addpath(genpath(voxeliseFolder));

baseOutPath = 'C:\Users\georg\Documents\MATLAB\MCmatlab-Release\Paper Data';

%% ---------------- QUICK PATH CHECKS ----------------
if exist('VOXELISE','file') ~= 2
    error('VOXELISE.m not found. Check voxeliseFolder and addpath(genpath(...)).');
end
if ~exist(stlPath,'file')
    error('STL not found: %s', stlPath);
end

%% ---------------- APERTURE / MATERIAL SETTINGS ----------------
AP = struct();

AP.unitScale_to_cm = 0.1;       % STL is mm -> cm
AP.flipZ           = true;      % set false if STL already oriented correctly
AP.R               = eye(3);    % optional extra rotation matrix

AP.plateTopZ_cm    = 0.0;       % sensor plane at z=0 (TOP of cuboid)

AP.tunnelDiam_mm    = 4.0;      % clear tunnel diameter
AP.collectorDiam_mm = 5.0;      % PD active diameter
AP.gelThickness_mm  = 12.0;     % gel thickness below gel surface
AP.sourceDepth_mm   = 0.20;     % source placed this far inside gel (below zGelTop)

AP.rayDirection     = 'xyz';
AP.verbose          = true;

% Plastic (baffle) optical properties
AP.plastic_mua_cm1 = 1e4;
AP.plastic_mus_cm1 = 1e-8;
AP.plastic_g       = 1.0;

%% ---------------- PHANTOM / GRID SETTINGS ----------------
% Phantom lateral size: 35 x 35 mm
Lx_cm = 3.5;
Ly_cm = 3.5;

% Voxel size target
dz_mm_desired = 0.10;                 % 0.1 mm
dz_cm_desired = dz_mm_desired / 10;   % cm

nx = 101;  % odd for centerline
ny = 101;

%% ---------------- READ STL + PLACE IT + BUILD GRID ----------------
fprintf('\n[INFO] Reading STL and computing placement...\n');

[STL, Vplaced, zGelTop_cm, stlHeight_cm] = readAndPlaceSTL(stlPath, AP);

gelThickness_cm = AP.gelThickness_mm / 10;
Lz_cm = zGelTop_cm + gelThickness_cm;

nz = max(51, round(Lz_cm/dz_cm_desired) + 1);

fprintf('  [INFO] Plate top at z=%.3f mm, Gel surface will be at z=%.3f mm\n', ...
    10*AP.plateTopZ_cm, 10*zGelTop_cm);
fprintf('[INFO] STL part height = %.3f mm, so gel surface zGelTop = %.3f mm\n', ...
    10*stlHeight_cm, 10*zGelTop_cm);
fprintf('[INFO] Set cuboid Lz = %.3f mm\n', 10*Lz_cm);
fprintf('[INFO] Grid: nx=%d, ny=%d, nz=%d\n', nx, ny, nz);

gridX = linspace(-Lx_cm/2, +Lx_cm/2, nx);
gridY = linspace(-Ly_cm/2, +Ly_cm/2, ny);
gridZ = linspace(0, Lz_cm, nz);
dz_cm = gridZ(2) - gridZ(1);
fprintf('[INFO] dz = %.4f mm\n', 10*dz_cm);

%% ---------------- VOXELIZE STL -> apMask ----------------
fprintf('[INFO] Voxelizing STL aperture into MC grid...\n');
apMask = voxelizeSTLToMask(STL.faces, Vplaced, gridX, gridY, gridZ, AP);
apMask = logical(apMask);

checkTunnelOpen(apMask, gridX, gridY, gridZ, AP, zGelTop_cm, 'BEFORE enforcement');

% Enforce: clear tunnel + global baffle above gel except tunnel
apMask = enforceTunnelAndGlobalBaffle(apMask, gridX, gridY, gridZ, AP, zGelTop_cm);

checkTunnelOpen(apMask, gridX, gridY, gridZ, AP, zGelTop_cm, 'AFTER enforcement');

%% ---------------- OUTPUT FOLDER ----------------
timestamp = datestr(now,'yyyymmdd_HHMMSS');
outFolder = fullfile(baseOutPath, ['MC_LUT_' timestamp]);
figFolder = fullfile(outFolder, 'figures');
lutFolder = fullfile(outFolder, 'lut');

if ~exist(baseOutPath,'dir'); mkdir(baseOutPath); end
if ~exist(outFolder,'dir');   mkdir(outFolder);   end
if ~exist(figFolder,'dir');   mkdir(figFolder);   end
if ~exist(lutFolder,'dir');   mkdir(lutFolder);   end

fprintf('\n[INFO] Saving outputs to:\n  %s\n\n', outFolder);

%% ---------------- BUILD MCmatlab MODEL ----------------
model = MCmatlab.model;

model.G.nx = nx; model.G.ny = ny; model.G.nz = nz;
model.G.Lx = Lx_cm; model.G.Ly = Ly_cm; model.G.Lz = Lz_cm;

model.G.mediaPropertiesFunc = @mediaPropertiesFunc;
model.G.geomFunc            = @geometryDefinition;
model.G.geomFuncParams      = {apMask, zGelTop_cm};

% MC settings
model.MC.boundaryType      = 1;        % all boundaries escaping
model.MC.matchedInterfaces = true;     % treat all refractive indices as 1

% Source: pencil beam (launched downward into +z)
model.MC.lightSource.sourceType = 0;
model.MC.lightSource.theta      = 0;
model.MC.lightSource.phi        = 0;
model.MC.lightSource.xFocus     = 0;

% Place source slightly INSIDE gel
sourceDepth_cm = AP.sourceDepth_mm / 10;
model.MC.lightSource.zFocus = zGelTop_cm + sourceDepth_cm;

% Collector
model.MC.useLightCollector = true;

model.MC.lightCollector.diam      = (AP.collectorDiam_mm/10);
model.MC.lightCollector.fieldSize = (AP.collectorDiam_mm/10);
model.MC.lightCollector.NA        = 0.22;

% IMPORTANT FIX:
% Some MCmatlab versions require lightCollector.res == 1 when lightCollector.f is Inf.
% We keep res=1 because you only need total collected photons (not imaging).
model.MC.lightCollector.res = 1;

% Put collector center OUTSIDE cuboid (z < 0)
model.MC.lightCollector.x = 0;
model.MC.lightCollector.y = 0;
model.MC.lightCollector.z = -0.5*dz_cm;

% Collector looks DOWN into the cuboid (z increases down)
model.MC.lightCollector.theta = 0;
model.MC.lightCollector.phi   = 0;

assert(model.MC.lightCollector.z < 0, 'Collector must be outside cuboid: set z < 0.');

%% ---------------- ABSORPTION / WAVELENGTHS / RESPONSIVITY ----------------
ConcentrationWater    = 0.70;
ExtinctionWater1450   = 28.8;  % [cm^-1]
ExtinctionWater1650   = 5.7;   % [cm^-1]
mua1450 = ConcentrationWater * ExtinctionWater1450;
mua1650 = ConcentrationWater * ExtinctionWater1650;

wavelengthList = [1450, 1650];
muaList        = [mua1450, mua1650];
wlLabels       = {'1450 nm','1650 nm'};
respList_AperW = [0.8, 0.4];   % optional responsivity scaling

%% ---------------- mus' SWEEP + SDS ----------------
g_tissue     = 0.9;
musPrimeList = [0.5 1 2 3 5 7 10 15 20 25 30];
nMus         = numel(musPrimeList);

SDSlist_cm = [0.3, 0.7];   % 3 mm and 7 mm
SDSlabels  = {'3 mm','7 mm'};
nSDS       = numel(SDSlist_cm);

%% ---------------- HIGH-STATISTICS SETTINGS ----------------
STATS = struct();
STATS.targetRelSE = 0.03;      % target relative standard error ~3% (p small -> Poisson)
STATS.minNdet     = 5000;      % minimum detected photons per point (if physically achievable)
STATS.N_start     = 2e6;       % starting photons per point
STATS.N_cap       = 2e9;       % hard cap (adjust if needed)
STATS.maxIter     = 7;         % max adaptive reruns per point
STATS.growFactor  = 4.0;       % minimum growth factor when increasing N

%% ---------------- PLOT GEOMETRY + VOXEL SUMS ----------------
plotGeometryOverview(Vplaced, STL.faces, zGelTop_cm, gelThickness_cm, AP, model, SDSlist_cm, figFolder);
plotVoxelSums(apMask, figFolder);

%% ---------------- DEBUG: ONE MC PLOTS RUN (FLUENCE + ABSORBED POWER SAVED) ----------------
% This produces the MCmatlab “Normalized fluence rate” and “Normalized absorbed power” plots,
% then auto-exports them to figFolder.
DEBUG = struct();
DEBUG.enable     = true;
DEBUG.wl_nm      = 1650;
DEBUG.mua_cm1    = mua1650;
DEBUG.musPrime   = 10;
DEBUG.SDS_cm     = 0.3;      % 3 mm
DEBUG.N_photons  = 5e7;      % increase if your plots still look “pixel/speckle-limited”

if DEBUG.enable
    runAndSaveMCPlotsCase(model, zGelTop_cm, sourceDepth_cm, g_tissue, AP, ...
        DEBUG.SDS_cm, DEBUG.wl_nm, DEBUG.mua_cm1, DEBUG.musPrime, DEBUG.N_photons, figFolder);
end

%% ---------------- LUT STORAGE ----------------
nW = numel(wavelengthList);

fracCollected_raw  = zeros(nMus, nSDS, nW);
sigmaFrac_raw      = zeros(nMus, nSDS, nW);

fracCollected_eff  = zeros(nMus, nSDS, nW);
sigmaFrac_eff      = zeros(nMus, nSDS, nW);

NdetCollected      = zeros(nMus, nSDS, nW);
Nlaunched          = zeros(nMus, nSDS, nW);
Nrequested         = zeros(nMus, nSDS, nW);

%% ========================================================================
% MAIN LUT LOOP (adaptive photons per point to hit STATS targets)
% ========================================================================
for m = 1:nMus
    musPrime = musPrimeList(m);
    mus      = musPrime / (1 - g_tissue);
    fprintf('\n================ musPrime = %.3g cm^-1 ================\n', musPrime);

    for sdsIdx = 1:nSDS
        sds_cm   = SDSlist_cm(sdsIdx);
        SDSlabel = SDSlabels{sdsIdx};

        model.MC.lightSource.xFocus = 0;
        model.MC.lightSource.yFocus = +sds_cm;
        model.MC.lightSource.zFocus = zGelTop_cm + sourceDepth_cm;

        for w = 1:nW
            wl   = wavelengthList(w);
            mua  = muaList(w);
            resp = respList_AperW(w);

            fprintf('Running %s @ SDS %s, musPrime=%.3g, mua=%.3g ...\n', ...
                wlLabels{w}, SDSlabel, musPrime, mua);

            model.MC.wavelength     = wl;
            model.G.mediaPropParams = {mus, mua, g_tissue, AP};

            % Run with adaptive photon budget
            [mOut, Ndet, NreqFinal] = runToTargetStats(model, STATS);

            Nsim = mOut.MC.nPhotons;
            if Nsim <= 0; Nsim = NreqFinal; end

            p_raw     = Ndet / Nsim;
            % binomial approx (Poisson regime p small -> sigma ~ sqrt(p/N))
            sigma_raw = sqrt(max(p_raw*(1-p_raw)/Nsim, 0));

            p_eff     = resp * p_raw;
            sigma_eff = resp * sigma_raw;

            fracCollected_raw(m, sdsIdx, w) = p_raw;
            sigmaFrac_raw(m, sdsIdx, w)     = sigma_raw;

            fracCollected_eff(m, sdsIdx, w) = p_eff;
            sigmaFrac_eff(m, sdsIdx, w)     = sigma_eff;

            NdetCollected(m, sdsIdx, w) = Ndet;
            Nlaunched(m, sdsIdx, w)     = Nsim;
            Nrequested(m, sdsIdx, w)    = NreqFinal;

            fprintf('  -> Ndet=%d / Nsim=%d (requested %d)\n', Ndet, Nsim, NreqFinal);
            fprintf('  -> RAW=%.6g (%.6e%%), 1sigma=%.3g (%.3e%%)\n', p_raw, 100*p_raw, sigma_raw, 100*sigma_raw);
            fprintf('  -> EFF(resp=%.3g)=%.6g (%.6e%%), 1sigma=%.3g (%.3e%%)\n', resp, p_eff, 100*p_eff, sigma_eff, 100*sigma_eff);
        end
    end
end

%% ---------------- SAVE LUT ----------------
pct_raw      = 100 * fracCollected_raw;
pctSigma_raw = 100 * sigmaFrac_raw;

pct_eff      = 100 * fracCollected_eff;
pctSigma_eff = 100 * sigmaFrac_eff;

nRows = nMus * nSDS * nW;
musPrime_col  = zeros(nRows,1);
SDSmm_col     = zeros(nRows,1);
wl_col        = zeros(nRows,1);
mua_col       = zeros(nRows,1);
resp_col      = zeros(nRows,1);

frac_raw_col     = zeros(nRows,1);
sigma_raw_col    = zeros(nRows,1);
pct_raw_col      = zeros(nRows,1);
pctSigma_raw_col = zeros(nRows,1);

frac_eff_col     = zeros(nRows,1);
sigma_eff_col    = zeros(nRows,1);
pct_eff_col      = zeros(nRows,1);
pctSigma_eff_col = zeros(nRows,1);

Ndet_col      = zeros(nRows,1);
Nsim_col      = zeros(nRows,1);
Nreq_col      = zeros(nRows,1);

r = 0;
for im = 1:nMus
    for is = 1:nSDS
        for iw = 1:nW
            r = r + 1;
            musPrime_col(r) = musPrimeList(im);
            SDSmm_col(r)    = 10 * SDSlist_cm(is);
            wl_col(r)       = wavelengthList(iw);
            mua_col(r)      = muaList(iw);
            resp_col(r)     = respList_AperW(iw);

            frac_raw_col(r)     = fracCollected_raw(im,is,iw);
            sigma_raw_col(r)    = sigmaFrac_raw(im,is,iw);
            pct_raw_col(r)      = pct_raw(im,is,iw);
            pctSigma_raw_col(r) = pctSigma_raw(im,is,iw);

            frac_eff_col(r)     = fracCollected_eff(im,is,iw);
            sigma_eff_col(r)    = sigmaFrac_eff(im,is,iw);
            pct_eff_col(r)      = pct_eff(im,is,iw);
            pctSigma_eff_col(r) = pctSigma_eff(im,is,iw);

            Ndet_col(r)      = NdetCollected(im,is,iw);
            Nsim_col(r)      = Nlaunched(im,is,iw);
            Nreq_col(r)      = Nrequested(im,is,iw);
        end
    end
end

LUT_table = table( ...
    musPrime_col, SDSmm_col, wl_col, mua_col, resp_col, ...
    frac_raw_col, sigma_raw_col, pct_raw_col, pctSigma_raw_col, ...
    frac_eff_col, sigma_eff_col, pct_eff_col, pctSigma_eff_col, ...
    Ndet_col, Nsim_col, Nreq_col, ...
    'VariableNames', {'musPrime_cm1','SDS_mm','wavelength_nm','mua_cm1','resp_AperW', ...
                      'frac_raw','sigma_raw','pct_raw','pctSigma_raw', ...
                      'frac_eff','sigma_eff','pct_eff','pctSigma_eff', ...
                      'Ndet','Nsimulated','Nrequested'} );

LUT = struct();
LUT.meta.timestamp          = timestamp;
LUT.meta.outFolder          = outFolder;
LUT.meta.apertureSTL        = stlPath;
LUT.meta.zGelTop_cm         = zGelTop_cm;
LUT.meta.gelThickness_cm    = gelThickness_cm;
LUT.meta.geometry           = struct('Lx_cm',model.G.Lx,'Ly_cm',model.G.Ly,'Lz_cm',model.G.Lz, ...
                                     'nx',model.G.nx,'ny',model.G.ny,'nz',model.G.nz);
LUT.meta.boundaryType       = model.MC.boundaryType;
LUT.meta.matchedInterfaces  = model.MC.matchedInterfaces;
LUT.meta.g_tissue           = g_tissue;
LUT.meta.ConcentrationWater = ConcentrationWater;
LUT.meta.stats              = STATS;

LUT.axes.musPrime_cm1   = musPrimeList(:);
LUT.axes.SDS_cm         = SDSlist_cm(:);
LUT.axes.SDS_mm         = 10*SDSlist_cm(:);
LUT.axes.wavelength_nm  = wavelengthList(:);
LUT.axes.mua_cm1        = muaList(:);
LUT.axes.resp_AperW     = respList_AperW(:);

LUT.data.frac_raw       = fracCollected_raw;
LUT.data.sigma_raw      = sigmaFrac_raw;
LUT.data.pct_raw        = pct_raw;
LUT.data.pctSigma_raw   = pctSigma_raw;

LUT.data.frac_eff       = fracCollected_eff;
LUT.data.sigma_eff      = sigmaFrac_eff;
LUT.data.pct_eff        = pct_eff;
LUT.data.pctSigma_eff   = pctSigma_eff;

LUT.data.Ndet           = NdetCollected;
LUT.data.Nsimulated     = Nlaunched;
LUT.data.Nrequested     = Nrequested;

matFile = fullfile(lutFolder, 'LUT_MC_reflectance_withAperture_TunnelOnly_HighStats.mat');
csvFile = fullfile(lutFolder, 'LUT_MC_reflectance_withAperture_TunnelOnly_HighStats.csv');

save(matFile, 'LUT', 'LUT_table');
writetable(LUT_table, csvFile);

fprintf('\n[INFO] Saved LUT:\n  %s\n  %s\n', matFile, csvFile);

%% ---------------- PLOTS: EFFECTIVE SIGNAL VS mus' ----------------
f1 = figure('Color','w');
errorbar(musPrimeList, squeeze(pct_eff(:,1,1)), squeeze(pctSigma_eff(:,1,1)), '-o', 'LineWidth', 2); hold on;
errorbar(musPrimeList, squeeze(pct_eff(:,2,1)), squeeze(pctSigma_eff(:,2,1)), '-s', 'LineWidth', 2);
grid on; xlabel('musPrime [cm^{-1}]'); ylabel('Effective collected signal [%]');
title('1450 nm: Effective collected vs musPrime (±1sigma)');
legend('SDS = 3 mm', 'SDS = 7 mm', 'Location', 'best');
exportgraphics(f1, fullfile(figFolder,'EffCollected_vs_musPrime_1450nm.png'), 'Resolution', 300);

f2 = figure('Color','w');
errorbar(musPrimeList, squeeze(pct_eff(:,1,2)), squeeze(pctSigma_eff(:,1,2)), '-o', 'LineWidth', 2); hold on;
errorbar(musPrimeList, squeeze(pct_eff(:,2,2)), squeeze(pctSigma_eff(:,2,2)), '-s', 'LineWidth', 2);
grid on; xlabel('musPrime [cm^{-1}]'); ylabel('Effective collected signal [%]');
title('1650 nm: Effective collected vs musPrime (±1sigma)');
legend('SDS = 3 mm', 'SDS = 7 mm', 'Location', 'best');
exportgraphics(f2, fullfile(figFolder,'EffCollected_vs_musPrime_1650nm.png'), 'Resolution', 300);

fprintf('[INFO] Saved figures to:\n  %s\n', figFolder);

%% ========================================================================
% ------------------------- LOCAL FUNCTIONS -------------------------------
% ========================================================================

function M = geometryDefinition(X,Y,Z,params)
    apMask  = params{1};
    zGelTop = params{2};

    % 1 = air, 2 = gel, 3 = plastic
    M = ones(size(X));        % air in top region
    M(Z > zGelTop) = 2;       % gel below surface
    M(apMask) = 3;            % plastic overrides
end

function mediaProperties = mediaPropertiesFunc(parameters)
    mediaProperties = MCmatlab.mediumProperties;

    mus      = parameters{1};
    mua      = parameters{2};
    g_tissue = parameters{3};
    AP       = parameters{4};

    j=1;
    mediaProperties(j).name = 'air';
    mediaProperties(j).mua  = 1e-8;
    mediaProperties(j).mus  = 1e-8;
    mediaProperties(j).g    = 1.0;

    j=2;
    mediaProperties(j).name = 'gel';
    mediaProperties(j).mua  = mua;
    mediaProperties(j).mus  = mus;
    mediaProperties(j).g    = g_tissue;

    j=3;
    mediaProperties(j).name = 'aperture_plastic';
    mediaProperties(j).mua  = AP.plastic_mua_cm1;
    mediaProperties(j).mus  = AP.plastic_mus_cm1;
    mediaProperties(j).g    = AP.plastic_g;
end

function [STL, Vplaced, zGelTop_cm, stlHeight_cm] = readAndPlaceSTL(stlPath, AP)
    STL = tryReadSTL_asFV(stlPath);
    V = double(STL.vertices);

    % mm -> cm
    V = V * AP.unitScale_to_cm;

    if AP.flipZ
        V(:,3) = -V(:,3);
    end

    V = (AP.R * V')';

    % Center X/Y
    cx = 0.5*(min(V(:,1))+max(V(:,1)));
    cy = 0.5*(min(V(:,2))+max(V(:,2)));
    V(:,1) = V(:,1) - cx;
    V(:,2) = V(:,2) - cy;

    % Shift so minZ -> plateTopZ
    zMin = min(V(:,3));
    V(:,3) = V(:,3) - zMin + AP.plateTopZ_cm;

    zGelTop_cm   = max(V(:,3));               % bottom-most point of STL = gel surface
    stlHeight_cm = max(V(:,3)) - min(V(:,3)); % should equal zGelTop if min is 0

    Vplaced = V;

    if AP.verbose
        fprintf('  [AP] STL bbox after transforms (cm):\n');
        fprintf('       X: [%.4f, %.4f] (span %.3f mm)\n', min(V(:,1)), max(V(:,1)), 10*(max(V(:,1))-min(V(:,1))));
        fprintf('       Y: [%.4f, %.4f] (span %.3f mm)\n', min(V(:,2)), max(V(:,2)), 10*(max(V(:,2))-min(V(:,2))));
        fprintf('       Z: [%.4f, %.4f] (span %.3f mm)\n', min(V(:,3)), max(V(:,3)), 10*(max(V(:,3))-min(V(:,3))));
    end
end

function FV = tryReadSTL_asFV(stlPath)
    try
        TR = stlread(stlPath);
        if isa(TR,'triangulation')
            FV.faces    = TR.ConnectivityList;
            FV.vertices = TR.Points;
            return;
        elseif isstruct(TR) && isfield(TR,'faces') && isfield(TR,'vertices')
            FV = TR;
            return;
        end
    catch
    end

    if exist('read_stl','file') == 2
        FV = read_stl(stlPath);
        if isstruct(FV) && isfield(FV,'faces') && isfield(FV,'vertices')
            return;
        end
    end

    error('Could not read STL. Ensure stlread() or read_stl() exists on path.');
end

function apMask = voxelizeSTLToMask(F, Vplaced, gridX, gridY, gridZ, AP)
    FV2 = struct('faces', double(F), 'vertices', double(Vplaced));
    try
        apMask = VOXELISE(gridX, gridY, gridZ, FV2, AP.rayDirection);
    catch ME
        error(['VOXELISE call failed. Confirm VOXELISE supports (gridX,gridY,gridZ,FV,rayDir).\n' ...
               'Original error:\n%s'], ME.message);
    end
    apMask = logical(apMask);
end

function apMask = enforceTunnelAndGlobalBaffle(apMask, gridX, gridY, gridZ, AP, zGelTop_cm)
    [Xg,Yg,Zg] = ndgrid(gridX, gridY, gridZ);

    rTunnel_cm = (AP.tunnelDiam_mm/2)/10;
    inTunnel   = (Xg.^2 + Yg.^2) <= (rTunnel_cm^2);

    % Air column: plate top down to gel surface
    inAirCol = (Zg >= AP.plateTopZ_cm) & (Zg <= zGelTop_cm);

    % Force-clear tunnel through air column
    apMask(inTunnel & inAirCol) = false;

    % Force global baffle: everything in air above gel is plastic EXCEPT tunnel
    apMask(inAirCol & ~inTunnel) = true;
end

function checkTunnelOpen(apMask, gridX, gridY, gridZ, AP, zGelTop_cm, label)
    [~,ix0] = min(abs(gridX - 0));
    [~,iy0] = min(abs(gridY - 0));
    zAir = (gridZ >= AP.plateTopZ_cm) & (gridZ <= zGelTop_cm);

    blockedCount = nnz( squeeze(apMask(ix0,iy0,zAir)) );
    if blockedCount > 0
        warning('[WARN] %s: Centerline (x=0,y=0) has %d plastic voxels in air column. Tunnel may be blocked.', label, blockedCount);
    else
        fprintf('  [OK] %s: Centerline tunnel appears open (x=0,y=0 has no plastic voxels).\n', label);
    end
end

function plotGeometryOverview(Vplaced, F, zGelTop_cm, gelThickness_cm, AP, model, SDSlist_cm, figFolder)
    f = figure('Color','w');
    hold on;

    patch('Faces',F,'Vertices',Vplaced, ...
        'FaceColor',[0.2 0.2 0.2],'FaceAlpha',0.25,'EdgeColor','none');

    Lx = model.G.Lx; Ly = model.G.Ly;
    z0 = zGelTop_cm;
    z1 = zGelTop_cm + gelThickness_cm;

    Cb = makeCuboidPatch(-Lx/2, +Lx/2, -Ly/2, +Ly/2, z0, z1);
    patch('Faces',Cb.faces,'Vertices',Cb.verts,'FaceColor',[1 0.5 0.1],'FaceAlpha',0.10,'EdgeColor','none');

    % Gel surface plane
    surf([-Lx/2 Lx/2; -Lx/2 Lx/2], [-Ly/2 -Ly/2; Ly/2 Ly/2], z0*ones(2), ...
        'FaceAlpha',0.10, 'EdgeColor','none');

    % Collector disk (visualize on z=0 plane)
    rCol_cm = (AP.collectorDiam_mm/2)/10;
    th = linspace(0,2*pi,200);
    zColPlot = AP.plateTopZ_cm;
    plot3(rCol_cm*cos(th), rCol_cm*sin(th), zColPlot*ones(size(th)), 'b', 'LineWidth', 2);
    text(0,0,zColPlot, sprintf('Collector (%.1fmm)', AP.collectorDiam_mm), 'Color','b');

    % Sources
    zSrc = zGelTop_cm + (AP.sourceDepth_mm/10);
    for k=1:numel(SDSlist_cm)
        y = SDSlist_cm(k);
        plot3(0,y,zSrc,'o','MarkerSize',8,'LineWidth',2);
        text(0,y,zSrc, sprintf(' Source @ SDS=%.1f mm', 10*SDSlist_cm(k)));
    end

    xlabel('x [cm]'); ylabel('y [cm]'); zlabel('z [cm]');
    title('Aperture (mesh) + Gel volume + Collector + Sources');
    axis equal; grid on; view(35,25);

    % Make z=0 appear on top visually
    set(gca,'ZDir','reverse');

    exportgraphics(f, fullfile(figFolder,'Geometry_Overview.png'), 'Resolution', 300);
end

function plotVoxelSums(apMask, figFolder)
    f = figure('Color','w','Position',[100 100 1200 450]);

    subplot(1,3,1);
    imagesc(squeeze(sum(apMask,1)));
    axis image; colormap(gray(256));
    title('sum over X'); xlabel('Z index'); ylabel('Y index');

    subplot(1,3,2);
    imagesc(squeeze(sum(apMask,2)));
    axis image; colormap(gray(256));
    title('sum over Y'); xlabel('Z index'); ylabel('X index');

    subplot(1,3,3);
    imagesc(squeeze(sum(apMask,3)));
    axis image; colormap(gray(256));
    title('sum over Z'); xlabel('Y index'); ylabel('X index');

    exportgraphics(f, fullfile(figFolder,'ApertureMask_VoxelSums.png'), 'Resolution', 300);
end

function [C] = makeCuboidPatch(x0,x1,y0,y1,z0,z1)
    verts = [ ...
        x0 y0 z0; x1 y0 z0; x1 y1 z0; x0 y1 z0; ...
        x0 y0 z1; x1 y0 z1; x1 y1 z1; x0 y1 z1  ...
    ];
    faces = [ ...
        1 2 3 4;
        5 6 7 8;
        1 2 6 5;
        2 3 7 6;
        3 4 8 7;
        4 1 5 8
    ];
    C.verts = verts;
    C.faces = faces;
end

function [mOut, Ndet, NreqFinal] = runToTargetStats(modelIn, STATS)
    % Runs MC adaptively to achieve:
    % - Ndet >= STATS.minNdet AND relative SE <= STATS.targetRelSE
    %   (when p>0), up to N_cap.
    %
    % If p is extremely small, it will hit N_cap and return whatever was obtained.

    modelOut = modelIn;
    Nreq = STATS.N_start;

    bestOut = [];
    bestNdet = 0;
    bestNreq = Nreq;

    for it = 1:STATS.maxIter
        modelOut.MC.nPhotonsRequested = min(Nreq, STATS.N_cap);

        modelOut = runMonteCarlo(modelOut);
        Ndet = modelOut.MC.nPhotonsCollected;
        Nsim = modelOut.MC.nPhotons;
        if Nsim <= 0; Nsim = modelOut.MC.nPhotonsRequested; end

        % Track best
        if Ndet > bestNdet
            bestOut = modelOut;
            bestNdet = Ndet;
            bestNreq = modelOut.MC.nPhotonsRequested;
        end

        % If no detections, force growth (until cap)
        if Ndet <= 0
            if modelOut.MC.nPhotonsRequested >= STATS.N_cap
                break;
            end
            Nreq = min(STATS.N_cap, max(ceil(STATS.growFactor*Nsim), Nsim+1));
            continue;
        end

        p = Ndet / Nsim;

        % Relative SE estimate (Poisson/binomial small p): relSE ~ 1/sqrt(Ndet)
        relSE = 1 / sqrt(max(Ndet,1));

        if (Ndet >= STATS.minNdet) && (relSE <= STATS.targetRelSE)
            mOut = modelOut;
            NreqFinal = modelOut.MC.nPhotonsRequested;
            return;
        end

        % Estimate needed N to reach BOTH thresholds:
        % 1) Ndet >= minNdet  => N >= minNdet / p
        % 2) relSE <= targetRelSE using Ndet ~ Np => N >= 1/(p*targetRelSE^2)
        N_need1 = STATS.minNdet / p;
        N_need2 = 1 / (p * STATS.targetRelSE^2);
        N_needed = ceil(max([N_need1, N_need2, STATS.growFactor*Nsim]));

        if modelOut.MC.nPhotonsRequested >= STATS.N_cap
            break;
        end
        Nreq = min(STATS.N_cap, max(N_needed, modelOut.MC.nPhotonsRequested+1));
    end

    % If we get here, we didn't meet targets within cap/iters. Return best achieved.
    if isempty(bestOut)
        mOut = modelOut;
        NreqFinal = modelOut.MC.nPhotonsRequested;
        Ndet = modelOut.MC.nPhotonsCollected;
    else
        mOut = bestOut;
        NreqFinal = bestNreq;
        Ndet = bestNdet;
    end
end

function runAndSaveMCPlotsCase(modelIn, zGelTop_cm, sourceDepth_cm, g_tissue, AP, ...
                              sds_cm, wl_nm, mua_cm1, musPrime, Nphot, figFolder)
    fprintf('\n[DEBUG] MC plots run (Fluence + Absorbed Power): wl=%dnm, SDS=%.1fmm, musPrime=%.3g, N=%g\n', ...
        wl_nm, 10*sds_cm, musPrime, Nphot);

    modelDbg = modelIn;

    mus = musPrime / (1 - g_tissue);

    modelDbg.MC.wavelength        = wl_nm;
    modelDbg.G.mediaPropParams    = {mus, mua_cm1, g_tissue, AP};
    modelDbg.MC.nPhotonsRequested = Nphot;

    modelDbg.MC.lightSource.xFocus = 0;
    modelDbg.MC.lightSource.yFocus = sds_cm;
    modelDbg.MC.lightSource.zFocus = zGelTop_cm + sourceDepth_cm;

    % Capture current figures
    figsBefore = findall(0,'Type','figure');

    modelDbg = runMonteCarlo(modelDbg);

    % Let MCmatlab generate its standard MC visualizations (includes FR + absorbed power)
    try
        modelDbg = plot(modelDbg,'MC');
    catch
        warning('plot(model,''MC'') failed in your MCmatlab build. Skipping auto MC plot export.');
        return;
    end

    % Capture newly created figures and export the ones we care about
    figsAfter = findall(0,'Type','figure');
    newFigs = setdiff(figsAfter, figsBefore);

    tagPrefix = sprintf('MC_%dnm_SDS%.0fmm_musPrime%.3g', wl_nm, 10*sds_cm, musPrime);

    exportMCFiguresByTitle(newFigs, figFolder, tagPrefix);
end

function exportMCFiguresByTitle(figs, figFolder, tagPrefix)
    % Saves any figure whose axes title contains:
    % - "Fluence"  -> FluenceRate
    % - "Absorbed" -> AbsorbedPower
    %
    % Also saves everything else with a generic name (optional).

    for i = 1:numel(figs)
        f = figs(i);
        axs = findall(f,'Type','axes');
        titleStr = '';

        if ~isempty(axs)
            % Use the first axes title that has content
            for k=1:numel(axs)
                t = axs(k).Title;
                if ~isempty(t) && isprop(t,'String')
                    s = t.String;
                    if iscell(s); s = strjoin(s,' '); end
                    if ~isempty(s)
                        titleStr = string(s);
                        break;
                    end
                end
            end
        end

        fname = '';
        if contains(lower(titleStr), "fluence")
            fname = sprintf('%s_FluenceRate.png', tagPrefix);
        elseif contains(lower(titleStr), "absorbed")
            fname = sprintf('%s_AbsorbedPower.png', tagPrefix);
        else
            % Optional: keep a record of other MC figures
            fname = sprintf('%s_MCfig_%02d.png', tagPrefix, i);
        end

        try
            exportgraphics(f, fullfile(figFolder, fname), 'Resolution', 300);
        catch
            try
                saveas(f, fullfile(figFolder, fname));
            catch
            end
        end
    end
end
