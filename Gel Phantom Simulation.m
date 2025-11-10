%% mc_gelphantom_inverse.m
% MCmatlab (v4.4.9-friendly) LUT builder for gel phantom.
% Two SDS: 3 mm & 7 mm (separate). Saves CSV and plots R_far/R_close vs mu_s'.

clear; clc;

%% ---------------- USER INPUTS ----------------
lambda_nm     = 1450;

% Composition (volume fractions; sum ~ 1)
f_water       = 0.80;
f_gelatin     = 0.20;    % <-- 20% gelatin
% No intralipid in this model
f_IL20        = 0.00;    % <-- explicitly zero

% Optical model
g_aniso       = 0.90;    % Henyey–Greenstein anisotropy
n_gel         = 1.33;
n_air         = 1.00;

% (REMOVED) Intralipid mu_s' model — not used (C_IL=0)
A_IL          = 300;    b_IL = 1.0; beta_IL = 1.1; lambda0_nm = 800; %#ok<NASGU>
C_IL          = 0.0;    %#ok<NASGU>

% Optional gelatin scattering (set to 0 if unknown)
musp_gelatin_500nm = 5;   % [cm^-1] at 500 nm (tunable)
beta_gel           = 1.5;

% Component absorption coefficients at lambda (PLACEHOLDERS – replace!)
mua_water_cm  = mua_water_placeholder(lambda_nm);  % [cm^-1]
mua_gel_cm    = 0.05;                              % [cm^-1]

% mu_s' sweep
musp_grid_cm  = linspace(1.5, 10.0, 12);

% ---------------- Geometry (widened to cover 7 mm SDS) ----------------
Lz = 3.0;                 % [cm] depth
Lx = 2.5; Ly = 2.5;       % [cm] half-width = 1.25 cm > 0.70 + 0.12 + margin
nx = 151; ny = 151; nz = 180;  % grid

% Air/gel interface (z=0 top; gel below zSurface)
zSurface = 0.01;                 % [cm]

% Source: LED-like rectangular top-hat + Lambertian angles
src_type  = 5;                   % XY factorizable emitter
src_wX    = 0.03; src_wY = 0.02; % [cm] half-widths
src_theta = 0; src_phi = 0; src_psi = 0;

% Detectors (circular apertures on the top boundary)
SDS_close_cm        = 0.30;   % 3 mm
SDS_far_cm          = 0.70;   % 7 mm
det_radius_cm       = 0.12;   % enlarge a bit for SNR

% 3-D fallback near-surface shell thickness (if needed)
surface_shell_thick = 0.08;   % [cm]

% Runtime per point
sim_minutes = 0.12;
useAllCPUs  = true;
matchedRI   = false;          % false -> Fresnel (air/gel mismatch)

% OPTIONAL: measured voltages for inversion overlay
V_close_meas = NaN;   % set to your 3 mm voltage
V_far_meas   = NaN;   % set to your 7 mm voltage

%% --------- Mixed optical properties (80% H2O, 20% gelatin) ----------
% Absorption: volume-fraction-weighted (no IL)
mua_mix_cm = f_water*mua_water_cm + f_gelatin*mua_gel_cm;

% Scattering: NO intralipid; only optional gelatin term
musp_gel_cm = musp_gelatin_500nm * (lambda_nm/500)^(-beta_gel) * f_gelatin;
musp_base   = musp_gel_cm;  % informative baseline (not used in sweep)

fprintf('λ=%d nm: μa_mix=%.3f cm^-1; baseline μs''≈%.2f cm^-1 (gel %.2f; no IL)\n',...
  lambda_nm, mua_mix_cm, musp_base, musp_gel_cm);

%% --------- Global for callbacks (v4.4.9-safe) ----------
global MC_GEL_OPTS;
MC_GEL_OPTS = struct();
MC_GEL_OPTS.zSurface  = zSurface;
MC_GEL_OPTS.matchedRI = matchedRI;
MC_GEL_OPTS.g_gel     = g_aniso;
MC_GEL_OPTS.n_gel     = n_gel;
MC_GEL_OPTS.n_air     = n_air;
MC_GEL_OPTS.mua_gel   = mua_mix_cm;

%% --------- Outputs ----------
R_close = nan(numel(musp_grid_cm),1);
R_far   = nan(numel(musp_grid_cm),1);

%% --------- Sweep mu_s' ----------
for k = 1:numel(musp_grid_cm)
    musp_cm = musp_grid_cm(k);
    mus_cm  = musp_cm/(1 - g_aniso);   % convert mu_s' -> mu_s
    MC_GEL_OPTS.mus_gel = mus_cm;

    % Build model
    model = MCmatlab.model;
    model.G.nx = nx; model.G.ny = ny; model.G.nz = nz;
    model.G.Lx = Lx; model.G.Ly = Ly; model.G.Lz = Lz;
    model.G.mediaPropertiesFunc = @mediaPropertiesFunc_gel;
    model.G.geomFunc            = @geometry_air_over_gel;

    model.MC.useAllCPUs              = useAllCPUs;
    model.MC.simulationTimeRequested = sim_minutes;
    model.MC.matchedInterfaces       = matchedRI;
    model.MC.boundaryType            = 2;    % top escaping
    model.MC.wavelength              = lambda_nm;

    % Source
    model.MC.lightSource.sourceType  = src_type;
    model.MC.lightSource.xFocus = 0; model.MC.lightSource.yFocus = 0; model.MC.lightSource.zFocus = 0;
    model.MC.lightSource.focalPlaneIntensityDistribution.XDistr = 0;
    model.MC.lightSource.focalPlaneIntensityDistribution.XWidth = src_wX;
    model.MC.lightSource.focalPlaneIntensityDistribution.YDistr = 0;
    model.MC.lightSource.focalPlaneIntensityDistribution.YWidth = src_wY;
    model.MC.lightSource.angularIntensityDistribution.XDistr = 2;
    model.MC.lightSource.angularIntensityDistribution.YDistr = 2;
    model.MC.lightSource.theta = src_theta; model.MC.lightSource.phi = src_phi; model.MC.lightSource.psi = src_psi;

    % Run MC and force plotting (hydrates arrays in some builds)
    model = runMonteCarlo(model);
    model = plot(model,'MC');

    % -------- A) Try to read the top boundary image from the figure --------
    [NI2D, x2d, y2d] = grab_top_boundary_from_fig(Lx, Ly);
    used2D = ~isempty(NI2D);

    if used2D
        [XX2,YY2] = ndgrid(x2d, y2d);
        detmask_close = ( (XX2 - SDS_close_cm).^2 + (YY2 - 0).^2 ) <= det_radius_cm^2;
        detmask_far   = ( (XX2 - SDS_far_cm  ).^2 + (YY2 - 0).^2 ) <= det_radius_cm^2;

        % Guard: on-grid?
        if ~any(detmask_close(:))
            warning('Close detector mask empty (SDS_close=%.3f cm > domain half-width %.3f cm). Increase Lx,Ly.', ...
                    SDS_close_cm, Lx/2);
        end
        if ~any(detmask_far(:))
            warning('Far detector mask empty (SDS_far=%.3f cm > domain half-width %.3f cm). Increase Lx,Ly.', ...
                    SDS_far_cm, Lx/2);
        end

        dx = x2d(2)-x2d(1); dy = y2d(2)-y2d(1); pixA = dx*dy;

        R_close(k) = sum(NI2D(detmask_close), 'all') * pixA;
        R_far(k)   = sum(NI2D(detmask_far),   'all') * pixA;

    else
        % -------- B) Fallback: scan model for 2D/3D numeric arrays --------
        [arr2d, ~] = find_largest_numeric_array(model, 2);
        [arr3d, ~] = find_largest_numeric_array(model, 3);

        if ~isempty(arr2d)
            NI2D = arr2d;
            [sx, sy] = size(NI2D);
            x2d = linspace(-Lx/2, Lx/2, sx);
            y2d = linspace(-Ly/2, Ly/2, sy);
            [XX2,YY2] = ndgrid(x2d, y2d);

            detmask_close = ( (XX2 - SDS_close_cm).^2 + (YY2 - 0).^2 ) <= det_radius_cm^2;
            detmask_far   = ( (XX2 - SDS_far_cm  ).^2 + (YY2 - 0).^2 ) <= det_radius_cm^2;

            dx = x2d(2)-x2d(1); dy = y2d(2)-y2d(1); pixA = dx*dy;

            R_close(k) = sum(NI2D(detmask_close), 'all') * pixA;
            R_far(k)   = sum(NI2D(detmask_far),   'all') * pixA;

        elseif ~isempty(arr3d)
            FR3 = arr3d;
            [sx, sy, sz] = size(FR3);
            x = linspace(-Lx/2, Lx/2, sx);
            y = linspace(-Ly/2, Ly/2, sy);
            z = linspace(0,     Lz,    sz);
            [XX,YY,ZZ] = ndgrid(x,y,z);

            top_mask      = (ZZ <= surface_shell_thick);
            detmask_close = ( (XX - SDS_close_cm).^2 + (YY - 0).^2 ) <= det_radius_cm^2;
            detmask_far   = ( (XX - SDS_far_cm  ).^2 + (YY - 0).^2 ) <= det_radius_cm^2;

            roi_close = top_mask & detmask_close;
            roi_far   = top_mask & detmask_far;

            dx = x(2)-x(1); dy = y(2)-y(1); dz = z(2)-z(1);
            voxvol = dx*dy*dz;

            R_close(k) = sum(FR3(roi_close), 'all') * voxvol;
            R_far(k)   = sum(FR3(roi_far),   'all') * voxvol;

        else
            error('Could not obtain 2-D boundary image nor any 2-D/3-D array from MCmatlab.');
        end
    end

    fprintf('μs''=%5.2f cm^-1: R_close=%.3e, R_far=%.3e, ratio=%.3f\n',...
        musp_cm, R_close(k), R_far(k), safe_div(R_far(k),R_close(k)));
end

ratio = R_far ./ R_close;

%% --------- Table + CSV + Plot ----------
T = table(musp_grid_cm(:), R_close(:), R_far(:), ratio(:), ...
    'VariableNames', {'mu_s_prime_cm1','R_close','R_far','R_far_over_R_close'});
disp('----- LUT preview -----'); disp(T(1:min(10,height(T)),:));

csvname = sprintf('mc_lut_%dnm_SDS_3mm_7mm.csv', round(lambda_nm));
writetable(T, csvname);
fprintf('Saved LUT to %s\n', csvname);

figure;
plot(musp_grid_cm, ratio, 'o-','LineWidth',1.5); grid on;
xlabel('\mu_s'' [cm^{-1}]'); ylabel('R_{far} / R_{close}');
title(sprintf('Ratio vs \\mu_s'' at %d nm (SDS 3 mm vs 7 mm)', round(lambda_nm)));

% Optional overlay: measured ratio & mu_s' estimate
if ~isnan(V_close_meas) && ~isnan(V_far_meas) && V_close_meas>0
    r_meas = V_far_meas / V_close_meas;
    hold on; yline(r_meas,'--');
    musp_est = interp1(ratio, musp_grid_cm, r_meas, 'linear','extrap');
    text(musp_grid_cm(1), r_meas, sprintf('  r_{meas}=%.3f  ->  \\mu_s''≈ %.2f cm^{-1}', r_meas, musp_est), ...
        'VerticalAlignment','bottom');
    fprintf('Measured ratio r=%.4f -> estimated μs'' ≈ %.2f cm^-1\n', r_meas, musp_est);
end

%% ================= Helper Functions =================
function mediaProperties = mediaPropertiesFunc_gel(~)
  global MC_GEL_OPTS;
  mediaProperties = MCmatlab.mediumProperties;

  % Air
  j=1; mediaProperties(j).name='air';
  mediaProperties(j).mua = 1e-8; mediaProperties(j).mus = 1e-8;
  mediaProperties(j).g   = 1.0;  mediaProperties(j).n   = MC_GEL_OPTS.n_air;

  % Gel (80% water + 20% gelatin absorption baked into μa_gel)
  j=2; mediaProperties(j).name='gel';
  mediaProperties(j).mua = MC_GEL_OPTS.mua_gel;
  mediaProperties(j).mus = MC_GEL_OPTS.mus_gel;
  mediaProperties(j).g   = MC_GEL_OPTS.g_gel;
  mediaProperties(j).n   = MC_GEL_OPTS.matchedRI * MC_GEL_OPTS.n_air + (~MC_GEL_OPTS.matchedRI) * MC_GEL_OPTS.n_gel;
end

function M = geometry_air_over_gel(X,Y,Z,~)
  global MC_GEL_OPTS;
  zSurface = MC_GEL_OPTS.zSurface;
  M = ones(size(X));
  M(Z > zSurface) = 2;
end

function [NI2D, x, y] = grab_top_boundary_from_fig(Lx, Ly)
% Find the top boundary irradiance image in open figures and extract CData
% and spatial axes. Returns [] if not found.
  NI2D = []; x = []; y = [];
  % look through all images in all figures, newest first
  figs = flip(findobj('Type','figure'));
  for f = 1:numel(figs)
      imgs = findobj(figs(f), 'Type', 'image');
      for i = 1:numel(imgs)
          I = imgs(i);
          C = get(I, 'CData');
          if ~isnumeric(C) || ndims(C)~=2 || any(size(C)<[16 16]) %#ok<ISMAT>
              continue;
          end
          % Prefer the largest 2-D image
          if isempty(NI2D) || numel(C) > numel(NI2D)
              NI2D = C;
              % Axes extents: try 'XData'/'YData' on the image, else use Lx/Ly
              if isprop(I,'XData') && isprop(I,'YData')
                  xd = get(I,'XData'); yd = get(I,'YData');
                  if numel(xd)==2
                      x = linspace(xd(1), xd(2), size(C,1));
                  end
                  if numel(yd)==2
                      y = linspace(yd(1), yd(2), size(C,2));
                  end
              end
          end
      end
      if ~isempty(NI2D), break; end
  end
  % If axes weren't provided by the image, center on (0,0) with domain Lx,Ly
  if ~isempty(NI2D)
      if isempty(x), x = linspace(-Lx/2, Lx/2, size(NI2D,1)); end
      if isempty(y), y = linspace(-Ly/2, Ly/2, size(NI2D,2)); end
  end
end

function [bestArr, bestPath] = find_largest_numeric_array(S, ndimsWanted)
% Recursively scan struct/cell for largest numeric array of ndimsWanted.
% Returns [] if none found. bestPath is a string with the field path.
  bestArr = [];
  bestPath = '';
  visited = containers.Map('KeyType','char','ValueType','logical');
  [bestArr, bestPath] = scan_node(S, 'model', ndimsWanted, bestArr, bestPath, visited);

  function [currBest, currPath] = scan_node(node, path, ndWanted, currBest, currPath, visited)
    key = sprintf('%s|%s', path, class(node));
    if isKey(visited, key); return; end
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

function v = safe_div(a,b)
  if b==0, v = NaN; else, v = a/b; end
end

function mua = mua_water_placeholder(lambda_nm)
% Placeholder; replace with tabulated μa_H2O(λ) for accuracy.
  if lambda_nm < 1000
      mua = 0.01;
  elseif lambda_nm < 1400
      mua = 0.10;
  elseif lambda_nm < 1500
      mua = 0.80;
  elseif lambda_nm < 1700
      mua = 0.50;
  else
      mua = 0.20;
  end
end
