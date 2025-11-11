%% mc_gelphantom_dual_088_hiResGPU_1hr_v53.m
% Dual-wavelength (1450 & 1650 nm) Monte Carlo LUT for gel phantom (f_H2O ≈ 0.88).
% - GPU auto (PCT+CUDA) with CPU fallback
% - Prints progress (point/seed/percent) and wavelength tag
% - Saves CSV + PNG + FIG per wavelength in a timestamped output folder
% - Runtime tuned ~1 hour total on RTX 4070 SUPER (≈30 min per wavelength)

clear; clc;

%% ---------------- RUNTIME TARGET ----------------
target_total_minutes = 60;   % aim ~1 hour total (both wavelengths together)

% Tuned for ~30 min per lambda on 4070 SUPER:
sim_minutes_per_run  = 0.50;  % minutes per seed per μs' value
seeds_per_pt         = 3;     % seeds per μs' value (averaged)
coarse_pts           = 12;    % coarse μs' points per λ
refine_pts           = 8;     % refine μs' points per λ (around max slope)
use_refine           = true;
refine_span          = 0.20;  % +/-20% band around slope-peak μs'

% Quick sanity check of plan (approx):
% per_lambda_runs = (coarse_pts + refine_pts) * seeds_per_pt
% per_lambda_time ≈ per_lambda_runs * sim_minutes_per_run
% total_time ≈ 2 * per_lambda_time

%% ---------------- GEOMETRY & DETECTORS ----------------
% Big box + fine grid (accuracy-first), still reasonable memory
Lx = 2.8; Ly = 2.8; Lz = 3.2;   % [cm]
nx = 181; ny = 181; nz = 220;   % voxel grid
zSurface = 0.01;                % air (1) over gel (2)

% Detectors (top boundary, centered at x=±SDS, y=0)
SDS_close_cm = 0.30;    % 3 mm
SDS_far_cm   = 0.70;    % 7 mm
det_radius_cm        = 0.10;    % balance sensitivity vs SNR
surface_shell_thick  = 0.05;    % integrate 3D fluence within this top shell

% MC controls
useAllCPUs   = true;
matchedRI    = false;

%% ---------------- OUTPUT FOLDER ----------------
outdir = fullfile(pwd, ['MC_LUT_', datestr(now,'yyyymmdd_HHMMSS')]);
if ~exist(outdir,'dir'), mkdir(outdir); end
fprintf('Output folder: %s\n', outdir);

%% ---------------- GPU AUTO-SELECT ----------------
[useGPU, gpuName] = pickGPU();
if useGPU
    fprintf('GPU enabled: %s\n', gpuName);
else
    fprintf('GPU not available → using CPU (Parallel Computing Toolbox required for GPU).\n');
end

%% ---------------- CONSTANTS ----------------
% Composition (80% water, 10% gelatin, 10% IL20); effective water column ~0.88
f_water   = 0.80;  f_gelatin = 0.10;  f_IL20 = 0.10;
f_H2O_eff = 0.88;

% Optical model (van Staveren-like IL μs', simple gelatin μs'; Hale–Querry water μa)
A_IL = 239.6; b_IL = 1.155; lambda0_nm = 500;
musp_gelatin_500 = 5; beta_gel = 1.5;
g_aniso = 0.90; n_gel = 1.33; n_air = 1.00;

% Two wavelengths, in this order
lambda_list = [1450, 1650];

%% --------------- MAIN: LOOP OVER WAVELENGTHS ---------------
for LAMBDA = lambda_list

    % --- wavelength-specific optical properties
    musp_IL20      = A_IL * (LAMBDA/lambda0_nm)^(-b_IL); % 20% stock
    musp_IL_for_10 = musp_IL20 * (f_IL20/0.20);          % 10% of stock
    musp_gel       = musp_gelatin_500 * (LAMBDA/500)^(-beta_gel) * f_gelatin;

    mua_water = mua_water_HaleQuerry(LAMBDA);
    mua_gel   = 0.05;
    mua_IL    = 0.01;
    mua_mix   = f_H2O_eff*mua_water + f_gelatin*mua_gel + f_IL20*mua_IL;

    fprintf('\n============================================================\n');
    fprintf('λ = %d nm | μa_mix = %.3f cm^-1, μs''_IL(10%%) = %.2f, μs''_gel = %.2f\n', ...
        LAMBDA, mua_mix, musp_IL_for_10, musp_gel);
    fprintf('============================================================\n');

    % --- μs' TOTAL sweep (absolute totals) around a center estimate
    musp_center_est = musp_IL_for_10 + musp_gel;    % typical total μs'
    musp_min = max(5, 0.25*musp_center_est);
    musp_max = max(musp_center_est*2.25, 80);
    musp_total_coarse = linspace(musp_min, musp_max, coarse_pts);

    % --- globals for media/geometry callbacks
    clear global MC_GEL_OPTS;
    global MC_GEL_OPTS;
    MC_GEL_OPTS = struct( ...
        'zSurface', zSurface, ...
        'matchedRI', matchedRI, ...
        'g_gel', g_aniso, ...
        'n_gel', n_gel, ...
        'n_air', n_air, ...
        'mua_gel', mua_mix );

    rng('shuffle','twister');

    % --- COARSE SWEEP
    [musp_c, Rc_c, Rf_c] = run_sweep(musp_total_coarse, ...
        Lx,Ly,Lz,nx,ny,nz, SDS_close_cm,SDS_far_cm, det_radius_cm,surface_shell_thick, ...
        useAllCPUs,matchedRI,sim_minutes_per_run,seeds_per_pt,LAMBDA,g_aniso,useGPU, ...
        sprintf('Coarse λ=%d', LAMBDA));

    ratio_c = Rf_c ./ Rc_c;

    if use_refine
        % pick refine band where slope magnitude is largest
        gC = gradient(ratio_c) ./ max(1e-12, gradient(musp_c));
        [~,idx] = max(abs(gC));
        center  = musp_c(idx);
        lo = max(5, center*(1 - refine_span));
        hi = center*(1 + refine_span);
        musp_total_fine = linspace(lo, hi, refine_pts);

        [musp_f, Rc_f, Rf_f] = run_sweep(musp_total_fine, ...
            Lx,Ly,Lz,nx,ny,nz, SDS_close_cm,SDS_far_cm, det_radius_cm,surface_shell_thick, ...
            useAllCPUs,matchedRI,sim_minutes_per_run,seeds_per_pt,LAMBDA,g_aniso,useGPU, ...
            sprintf('Refine λ=%d', LAMBDA));

        musp_total_vec = [musp_c; musp_f];
        R_close = [Rc_c; Rc_f];
        R_far   = [Rf_c; Rf_f];
    else
        musp_total_vec = musp_c;
        R_close = Rc_c; R_far = Rf_c;
    end

    % --- sort and build LUT
    [musp_total_vec, order] = sort(musp_total_vec);
    R_close = R_close(order); R_far = R_far(order);
    ratio   = R_far ./ R_close;

    T = table(musp_total_vec(:), R_close(:), R_far(:), ratio(:), ...
        'VariableNames', {'mu_s_prime_total_cm1','R_close','R_far','R_far_over_R_close'});

    % --- save CSV
    csvname = fullfile(outdir, sprintf('mc_lut_%dnm_fH2O088_HIRES_GPU_v53.csv', LAMBDA));
    writetable(T, csvname);
    fprintf('Saved LUT -> %s\n', csvname);

    % --- plot & save
    f = figure('Color','w'); plot(musp_total_vec, ratio, 'o-','LineWidth',1.6); grid on;
    xlabel('\mu_s''_{total} [cm^{-1}]'); ylabel('R_{far} / R_{close}');
    title(sprintf('LUT (λ=%d nm, f_{H2O}=0.88), %d pts (seeds=%d, t=%.2f min/run)', ...
          LAMBDA, numel(musp_total_vec), seeds_per_pt, sim_minutes_per_run));
    pngname = fullfile(outdir, sprintf('mc_lut_plot_%dnm_v53.png', LAMBDA));
    figname = fullfile(outdir, sprintf('mc_lut_plot_%dnm_v53.fig', LAMBDA));
    saveas(f, pngname); saveas(f, figname); close(f);
    fprintf('Saved plots -> %s , %s\n', pngname, figname);

    % --- brief preview
    disp('----- LUT preview -----'); disp(T(1:min(12,height(T)),:));
end

fprintf('\nAll done. Files saved in: %s\n', outdir);

%% ================= Helper Functions =================
function mua = mua_water_HaleQuerry(lambda_nm)
  % Hale & Querry μa of pure water [cm^-1] (sparse nodes + pchip)
  data_nm = [400 500 600 700 800 900 1000 1100 1200 1300 1400 1450 1500 1550 1600 1650 1700 1800 1900 2000];
  data_cm1= [0.00012 0.00022 0.00040 0.00090 0.0020 0.0040 0.010 0.020 0.040 0.20 2.0 25.0 14.0 8.0 6.0 5.0 4.0 15.0 90.0 120.0];
  mua = interp1(data_nm, data_cm1, lambda_nm, 'pchip', 'extrap');
end

function mediaProperties = mediaPropertiesFunc_gel(~)
  % Global media function for 2-layer air/gel model
  global MC_GEL_OPTS;
  mediaProperties = MCmatlab.mediumProperties;
  % Air
  j=1; mediaProperties(j).name='air';
  mediaProperties(j).mua = 1e-8; mediaProperties(j).mus = 1e-8;
  mediaProperties(j).g   = 1.0;  mediaProperties(j).n   = MC_GEL_OPTS.n_air;
  % Gel (mixture)
  j=2; mediaProperties(j).name='gel';
  mediaProperties(j).mua = MC_GEL_OPTS.mua_gel;
  mediaProperties(j).mus = MC_GEL_OPTS.mus_gel;
  mediaProperties(j).g   = MC_GEL_OPTS.g_gel;
  mediaProperties(j).n   = MC_GEL_OPTS.matchedRI * MC_GEL_OPTS.n_air + ...
                           (~MC_GEL_OPTS.matchedRI) * MC_GEL_OPTS.n_gel;
end

function M = geometry_air_over_gel(X,~,Z,~)
  % media index: 1 = air, 2 = gel
  global MC_GEL_OPTS;
  M = ones(size(X));
  M(Z > MC_GEL_OPTS.zSurface) = 2;
end

function [NI2D, x, y] = grab_top_boundary_from_fig(Lx, Ly)
  % Extract top-boundary irradiance image from newest MC figure (fallback).
  NI2D = []; x = []; y = [];
  figs = flip(findobj('Type','figure'));
  for f = 1:numel(figs)
      imgs = findobj(figs(f), 'Type', 'image');
      for i = 1:numel(imgs)
          I = imgs(i);
          C = get(I, 'CData');
          if ~isnumeric(C) || ndims(C)~=2 || any(size(C)<[16 16]), continue; end
          if isempty(NI2D) || numel(C) > numel(NI2D)
              NI2D = C;
              if isprop(I,'XData') && isprop(I,'YData')
                  xd = get(I,'XData'); yd = get(I,'YData');
                  if numel(xd)==2, x = linspace(xd(1), xd(2), size(C,1)); end
                  if numel(yd)==2, y = linspace(yd(1), yd(2), size(C,2)); end
              end
          end
      end
      if ~isempty(NI2D), break; end
  end
  if ~isempty(NI2D)
      if isempty(x), x = linspace(-Lx/2, Lx/2, size(NI2D,1)); end
      if isempty(y), y = linspace(-Ly/2, Ly/2, size(NI2D,2)); end
  end
end

function [bestArr, bestPath] = find_largest_numeric_array(S, ndimsWanted)
  % Recursively scan struct/cell for largest numeric array of ndimsWanted.
  bestArr = []; bestPath = '';
  visited = containers.Map('KeyType','char','ValueType','logical');
  [bestArr, bestPath] = scan_node(S, 'model', ndimsWanted, bestArr, bestPath, visited);

  function [currBest, currPath] = scan_node(node, path, ndWanted, currBest, currPath, visited)
    key = sprintf('%s|%s', path, class(node));
    if isKey(visited, key), return; end
    visited(key) = true;

    if isnumeric(node) && ndims(node)==ndWanted && all(size(node)>1)
        if isempty(currBest) || numel(node) > numel(currBest)
            currBest = node; currPath = path;
        end
        return;
    elseif isstruct(node)
        f = fieldnames(node);
        for i=1:numel(node)
            for j=1:numel(f)
                sub = node(i).(f{j});
                [currBest, currPath] = scan_node(sub, sprintf('%s.%s(%d)', path, f{j}, i), ndWanted, currBest, currPath, visited);
            end
        end
    elseif iscell(node)
        for i=1:numel(node)
            [currBest, currPath] = scan_node(node{i}, sprintf('%s{%d}', path, i), ndWanted, currBest, currPath, visited);
        end
    end
  end
end

function [useGPU, gpuName] = pickGPU()
  % Return true and name if we can use CUDA; false otherwise.
  useGPU = false; gpuName = '';
  hasPCT = license('test','Distrib_Computing_Toolbox');
  if ~hasPCT, return; end
  if ~exist('gpuDeviceCount','file'), return; end
  try
      n = gpuDeviceCount;
      if n > 0
          g = gpuDevice(1);   % 1-based in MATLAB
          gpuName = g.Name;
          useGPU = true;
      end
  catch
      useGPU = false;
  end
end

function [musp_vec, Rclose, Rfar] = run_sweep(musp_total_vec, ...
    Lx,Ly,Lz,nx,ny,nz, SDS_close_cm,SDS_far_cm, det_radius_cm,surface_shell_thick, ...
    useAllCPUs,matchedRI,sim_minutes,seeds_per_pt,lambda_nm,g_aniso,useGPU, passName)

    K  = numel(musp_total_vec);
    TT = K * seeds_per_pt;
    Rclose = zeros(K,1); Rfar = zeros(K,1);

    counter = 0;
    for k = 1:K
        musp_total = musp_total_vec(k);
        mus_total  = musp_total / (1 - g_aniso);

        close_accum = 0; far_accum = 0;

        for rep = 1:seeds_per_pt
            counter = counter + 1;
            pct = 100*counter/TT;
            fprintf('(pt %2d/%2d, seed %d/%d) %6.2f%% [%s]  λ=%d nm, μs''_total=%6.2f cm^-1 ...\n', ...
                k, K, rep, seeds_per_pt, pct, passName, lambda_nm, musp_total);

            % ---------- Build model ----------
            model = MCmatlab.model;
            model.G.nx = nx; model.G.ny = ny; model.G.nz = nz;
            model.G.Lx = Lx; model.G.Ly = Ly; model.G.Lz = Lz;
            model.G.mediaPropertiesFunc = @mediaPropertiesFunc_gel;
            model.G.geomFunc            = @geometry_air_over_gel;

            model.MC.useAllCPUs              = useAllCPUs;
            model.MC.simulationTimeRequested = sim_minutes;
            model.MC.matchedInterfaces       = matchedRI;
            model.MC.boundaryType            = 2;
            model.MC.wavelength              = lambda_nm;

            % Prefer numeric data (reduce need for plotting fallbacks)
            if isprop(model.MC,'saveVolumeData'),   model.MC.saveVolumeData = true;   end
            if isprop(model.MC,'saveBoundaryData'), model.MC.saveBoundaryData = true; end

            % GPU toggle
            if isprop(model.MC,'useGPU'),      model.MC.useGPU = useGPU; end
            if isprop(model.MC,'GPUdevice'),   model.MC.GPUdevice = 0;   end  % first GPU

            % Version-safe seeding
            seed = randi(2^31-1);
            if isprop(model.MC,'randomSeed'), model.MC.randomSeed = seed; else, rng(seed,'twister'); end

            % LED-like source
            model.MC.lightSource.sourceType = 5;
            model.MC.lightSource.xFocus = 0; model.MC.lightSource.yFocus = 0; model.MC.lightSource.zFocus = 0;
            model.MC.lightSource.focalPlaneIntensityDistribution.XDistr = 0;
            model.MC.lightSource.focalPlaneIntensityDistribution.XWidth = 0.03;
            model.MC.lightSource.focalPlaneIntensityDistribution.YDistr = 0;
            model.MC.lightSource.focalPlaneIntensityDistribution.YWidth = 0.02;
            model.MC.lightSource.angularIntensityDistribution.XDistr = 2;
            model.MC.lightSource.angularIntensityDistribution.YDistr = 2;

            % Pass μs to media function
            global MC_GEL_OPTS; MC_GEL_OPTS.mus_gel = mus_total;

            % ---------- Run MC ----------
            model = runMonteCarlo(model);

            % ---------- Prefer numeric 3-D ----------
            [FR3, ~] = find_largest_numeric_array(model, 3);

            if ~isempty(FR3)
                [sx, sy, sz] = size(FR3);
                x = linspace(-Lx/2, Lx/2, sx);
                y = linspace(-Ly/2, Ly/2, sy);
                z = linspace(0,     Lz,    sz);
                [XX,YY,ZZ] = ndgrid(x,y,z);

                top_mask      = (ZZ <= surface_shell_thick);
                detmask_close = ((XX - SDS_close_cm).^2 + YY.^2) <= det_radius_cm^2;
                detmask_far   = ((XX - SDS_far_cm ).^2 + YY.^2) <= det_radius_cm^2;

                roi_close = top_mask & detmask_close;
                roi_far   = top_mask & detmask_far;

                dx = x(2)-x(1); dy = y(2)-y(1); dz = z(2)-z(1);
                voxvol = dx*dy*dz;

                close_accum = close_accum + sum(FR3(roi_close), 'all') * voxvol;
                far_accum   = far_accum   + sum(FR3(roi_far),   'all') * voxvol;

            else
                % ---------- Fallback: use top boundary image from figure ----------
                model = plot(model,'MC');
                [NI2D, x2d, y2d] = grab_top_boundary_from_fig(Lx, Ly);
                if isempty(NI2D)
                    [NI2D, ~] = find_largest_numeric_array(model, 2);
                    if isempty(NI2D), error('Could not obtain MC numeric arrays nor figure image.'); end
                    [sx, sy] = size(NI2D);
                    x2d = linspace(-Lx/2, Lx/2, sx);
                    y2d = linspace(-Ly/2, Ly/2, sy);
                end

                [XX2,YY2] = ndgrid(x2d, y2d);
                detmask_close = ((XX2 - SDS_close_cm).^2 + YY2.^2) <= det_radius_cm^2;
                detmask_far   = ((XX2 - SDS_far_cm ).^2 + YY2.^2) <= det_radius_cm^2;

                dx = x2d(2)-x2d(1); dy = y2d(2)-y2d(1); pixA = dx*dy;
                close_accum = close_accum + sum(NI2D(detmask_close), 'all') * pixA;
                far_accum   = far_accum   + sum(NI2D(detmask_far),   'all') * pixA;
            end
        end

        Rclose(k) = close_accum / seeds_per_pt;
        Rfar(k)   = far_accum   / seeds_per_pt;

        fprintf('... R_close = %.3e, R_far = %.3e, ratio = %.4f\n', ...
            Rclose(k), Rfar(k), Rfar(k)/Rclose(k));
    end

    musp_vec = musp_total_vec(:);
end
