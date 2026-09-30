%% endmember_box_sensitivity.m
% Sensitivity test for PEW/PSUW endmember box choices.
%
% This script tests whether the main manuscript conclusions are robust to
% modest changes in the PEW and PSUW source boxes. It assumes that you have
% downloaded thetao and so from a domain that contains the original
% 10 x 10 degree boxes plus 2 degrees on every side.
%
% Main outputs:
%   1. sensitivity_results.csv
%      Summary metrics for each endmember-box realization.
%
%   2. sensitivity_results.mat
%      Full results, including profile-level and annual fPEW series.
%
%   3. figures/
%      Simple diagnostic figures showing sensitivity of:
%        - mean fPEW
%        - fPEW trends
%        - O2-fPEW slopes and R2
%        - low-O2 event attribution
%        - residual O'-N' slopes
%
% Required files:
%   CSOS_ctd.mat
%   P4Long_ctd.mat
%   PEW expanded thetao/so NetCDF files
%   PSUW/SAW expanded thetao/so NetCDF files
%
% Required toolboxes/functions:
%   GSW toolbox on MATLAB path.
%   mann_kendall.m on MATLAB path (same trend-significance test as manuscript).
%
% Notes:
%   - thetao from GLORYS/Copernicus is treated as potential temperature.
%   - P4 and CS-OS are evaluated on their station-specific manuscript
%     PEW-core isopycnals (26.39 and 26.48 kg m^-3, respectively).
%   - The script also stores the PEW-core density implied by each endmember
%     realization, so you can check whether the core shifts materially.

clear; close all; clc;

%% ---------------- USER SETTINGS ----------------

% Project data directory. Update this path if needed.
base_dir = '/Users/samst/Library/CloudStorage/Dropbox/Hakai/gyre/supply/data';


% Directory containing the expanded endmember downloads.
endmember_dir = fullfile(base_dir,'endmember_sensitivity');

% File patterns. These should match the files produced by the download script.
pew_pattern  = '*PEW*thetao*so*.nc';
psuw_pattern = '*PSUW*thetao*so*.nc';  % change to '*SAW*thetao*so*.nc' if needed

% Output
out_dir = fullfile(base_dir,'endmember_sensitivity_results');
fig_dir = fullfile(out_dir,'figures');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

% Analysis settings
target_sigma.CSOS = 26.48;   % station-specific manuscript PEW-core density
target_sigma.P4   = 26.39;   % station-specific manuscript PEW-core density
sgrid             = (23:0.05:27.9)';
n_random          = 100;     % 100 random + 12 deterministic = 112 total cases
rng(7);                       % reproducible random box perturbations

% Trend window for P4 fPEW
p4_trend_years = [1981 2024];

% Low-O2 event percentile
event_prctile = 5;

% Original source boxes from manuscript
orig.PEW.name = "PEW";
orig.PEW.lon  = [-113 -103];
orig.PEW.lat  = [0 10];
orig.PEW.lon0 = -108;
orig.PEW.lat0 = 5;

orig.PSUW.name = "PSUW";
orig.PSUW.lon  = [-150 -140];
orig.PSUW.lat  = [45 55];
orig.PSUW.lon0 = -145;
orig.PSUW.lat0 = 50;

%% ---------------- LOAD OBSERVATIONS ----------------
% assert(isfile(csos_mat), 'CSOS file not found: %s', csos_mat);
% assert(isfile(p4_mat),   'P4 file not found: %s', p4_mat);
load('CSOS_ctd.mat','CSOS');
load('P4Long_ctd.mat','P4');

assert(isfield(CSOS,'SA') && isfield(P4,'SA'), ...
    'CSOS and P4 must contain Absolute Salinity in field SA.');
assert(exist('mann_kendall','file') == 2, ...
    'mann_kendall.m must be on the MATLAB path.');

% Force datenum time vectors for compatibility with your existing workflow.
CSOS.mtime = local_to_datenum(CSOS.mtime);
P4.mtime   = local_to_datenum(P4.mtime);

% Basic QC consistent with the main script.
if isfield(CSOS,'oxML'),  CSOS.oxML(CSOS.oxML<10)   = CSOS.oxML(CSOS.oxML<10).*43; end
if isfield(CSOS,'oxML2'), CSOS.oxML2(CSOS.oxML2<10) = CSOS.oxML2(CSOS.oxML2<10).*43; end
if isfield(CSOS,'ox'),    CSOS.ox(CSOS.ox<10)       = CSOS.ox(CSOS.ox<10).*43; end
if isfield(CSOS,'ox2'),   CSOS.ox2(CSOS.ox2<10)     = CSOS.ox2(CSOS.ox2<10).*43; end
if isfield(P4,'ox'),      P4.ox(P4.ox>700)          = NaN; end

% Station property vectors on the station-specific manuscript target densities.
CS = local_station_on_sigma(CSOS, target_sigma.CSOS, "CSOS");
PP = local_station_on_sigma(P4,   target_sigma.P4,   "P4");

% Restrict P4 trend calculations to the manuscript trend window.
PP.pretrend_mask = PP.year < p4_trend_years(1) | PP.year > p4_trend_years(2);

%% ---------------- FIND ENDMEMBER FILES ----------------

pew_files = dir(fullfile(endmember_dir, pew_pattern));
psu_files = dir(fullfile(endmember_dir, psuw_pattern));

% Fallback in case files use SAW rather than PSUW.
if isempty(psu_files)
    psuw_pattern = '*SAW*thetao*so*.nc';
    psu_files = dir(fullfile(endmember_dir, psuw_pattern));
end

assert(~isempty(pew_files), 'No PEW files found. Check endmember_dir and pew_pattern.');
assert(~isempty(psu_files), 'No PSUW/SAW files found. Check endmember_dir and psuw_pattern.');

pew_files = fullfile({pew_files.folder},{pew_files.name});
psu_files = fullfile({psu_files.folder},{psu_files.name});

fprintf('Found %d PEW files and %d PSUW/SAW files.\n',numel(pew_files),numel(psu_files));

%% ---------------- DEFINE SENSITIVITY CASES ----------------

cases = local_build_cases(orig, n_random);

fprintf('Running %d endmember-box cases...\n',numel(cases));

%% ---------------- RUN SENSITIVITY TESTS ----------------

R = struct([]);
baseline = [];

for c = 1:numel(cases)

    fprintf('\n[%03d/%03d] %s\n', c, numel(cases), cases(c).name);

    % Source profiles for this box choice, interpolated onto sigma grid.
    PEWsrc  = local_source_profile(pew_files, cases(c).PEW,  sgrid);
    PSUWsrc = local_source_profile(psu_files, cases(c).PSUW, sgrid);

    % Profile-level fPEW on all observed depths.
    fPEW_CS_full = local_calc_fpew_profiles(CSOS, PEWsrc, PSUWsrc, sgrid);
    fPEW_P4_full = local_calc_fpew_profiles(P4,   PEWsrc, PSUWsrc, sgrid);

    % fPEW on station-specific manuscript target densities.
    fPEW_CS = local_interp_to_sigma(CSOS.s_t, fPEW_CS_full, target_sigma.CSOS);
    fPEW_P4 = local_interp_to_sigma(P4.s_t,   fPEW_P4_full, target_sigma.P4);

    % Omit early P4 values for trend consistency, but keep them in arrays.
    fPEW_P4_for_trend = fPEW_P4;
    fPEW_P4_for_trend(PP.pretrend_mask) = NaN;

    % Mean fPEW profile and implied PEW core density/depth for this case.
    coreCS = local_core_properties(CSOS, fPEW_CS_full, [26.25 26.70]);
    coreP4 = local_core_properties(P4,   fPEW_P4_full, [26.25 26.70]);

    % Main regression diagnostics.
    metCS = local_station_metrics(CS, fPEW_CS, event_prctile);
    metP4 = local_station_metrics(PP, fPEW_P4, event_prctile);

    % Annual fPEW means and trend.
    annCS = local_annual_means(CS.year, fPEW_CS);
    annP4 = local_annual_means(PP.year, fPEW_P4_for_trend);

    trCS = local_linear_trend(annCS.year, annCS.mean);
    trP4 = local_linear_trend(annP4.year, annP4.mean);

    % Propagate fPEW trend through O2-fPEW slope.
    % slope_O2 is per unit fPEW; trend is fPEW yr^-1.
    o2trendCS_dec = 10 .* metCS.slope_O2 .* trCS.slope;
    o2trendP4_dec = 10 .* metP4.slope_O2 .* trP4.slope;

    % Inter-station fPEW correlation using common years.
    [r2_inter,p_inter,n_inter] = local_common_year_corr(annP4.year,annP4.mean,annCS.year,annCS.mean);

    % Sensitivity relative to baseline.
    if c == 1
        baseline.fPEW_CS = fPEW_CS;
        baseline.fPEW_P4 = fPEW_P4;
        baseline.metCS   = metCS;
        baseline.metP4   = metP4;
    end

    dCS = local_diff_stats(fPEW_CS, baseline.fPEW_CS);
    dP4 = local_diff_stats(fPEW_P4, baseline.fPEW_P4);

    % Store compact result.
    R(c).case = cases(c).name;
    R(c).case_type = cases(c).type;
    R(c).PEW_lon = cases(c).PEW.lon;
    R(c).PEW_lat = cases(c).PEW.lat;
    R(c).PSUW_lon = cases(c).PSUW.lon;
    R(c).PSUW_lat = cases(c).PSUW.lat;
    R(c).PEWsrc = PEWsrc;
    R(c).PSUWsrc = PSUWsrc;

    R(c).CS_core_sigma = coreCS.sigma;
    R(c).CS_core_depth = coreCS.depth;
    R(c).P4_core_sigma = coreP4.sigma;
    R(c).P4_core_depth = coreP4.depth;

    R(c).CS = metCS;
    R(c).P4 = metP4;
    R(c).annCS = annCS;
    R(c).annP4 = annP4;
    R(c).trendCS = trCS;
    R(c).trendP4 = trP4;
    R(c).O2trendCS_from_fPEW_decade = o2trendCS_dec;
    R(c).O2trendP4_from_fPEW_decade = o2trendP4_dec;
    R(c).interstation_R2 = r2_inter;
    R(c).interstation_p = p_inter;
    R(c).interstation_n = n_inter;
    R(c).delta_fPEW_CS = dCS;
    R(c).delta_fPEW_P4 = dP4;
    R(c).fPEW_CS = fPEW_CS;
    R(c).fPEW_P4 = fPEW_P4;

end

%% ---------------- BUILD SUMMARY TABLE ----------------
T = local_results_table(R);

csvfile = fullfile(out_dir,'sensitivity_results.csv');
matfile = fullfile(out_dir,'sensitivity_results.mat');

writetable(T,csvfile);
save(matfile,'R','T','cases','target_sigma','orig','-v7.3');

assert(isfile(csvfile),'CSV output was not written: %s',csvfile);
assert(isfile(matfile),'MAT output was not written: %s',matfile);

disp(T(1:min(12,height(T)),:));

%% ---------------- QUICK FIGURES ----------------

local_make_figures(T, fig_dir);

fprintf('\nDone.\n');
fprintf('Wrote:\n  %s\n  %s\n',csvfile,matfile);

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function cases = local_build_cases(orig, n_random)

    cases = struct([]);
    k = 0;

    function add_case(name,type,pew_lon,pew_lat,psu_lon,psu_lat)
        k = k + 1;
        cases(k).name = string(name);
        cases(k).type = string(type);

        cases(k).PEW = orig.PEW;
        cases(k).PEW.lon = pew_lon;
        cases(k).PEW.lat = pew_lat;
        cases(k).PEW.lon0 = mean(pew_lon);
        cases(k).PEW.lat0 = mean(pew_lat);

        cases(k).PSUW = orig.PSUW;
        cases(k).PSUW.lon = psu_lon;
        cases(k).PSUW.lat = psu_lat;
        cases(k).PSUW.lon0 = mean(psu_lon);
        cases(k).PSUW.lat0 = mean(psu_lat);
    end

    % Baseline must be first.
    add_case("baseline_10x10", "deterministic", ...
        orig.PEW.lon, orig.PEW.lat, orig.PSUW.lon, orig.PSUW.lat);

    % Expanded boxes.
    add_case("both_expanded_p2deg", "deterministic", ...
        orig.PEW.lon + [-2 2], orig.PEW.lat + [-2 2], ...
        orig.PSUW.lon + [-2 2], orig.PSUW.lat + [-2 2]);

    add_case("PEW_expanded_only", "deterministic", ...
        orig.PEW.lon + [-2 2], orig.PEW.lat + [-2 2], ...
        orig.PSUW.lon, orig.PSUW.lat);

    add_case("PSUW_expanded_only", "deterministic", ...
        orig.PEW.lon, orig.PEW.lat, ...
        orig.PSUW.lon + [-2 2], orig.PSUW.lat + [-2 2]);

    % Shift PEW only by 2 degrees.
    add_case("PEW_shift_north_2deg", "deterministic", ...
        orig.PEW.lon, orig.PEW.lat + 2, orig.PSUW.lon, orig.PSUW.lat);
    add_case("PEW_shift_south_2deg", "deterministic", ...
        orig.PEW.lon, orig.PEW.lat - 2, orig.PSUW.lon, orig.PSUW.lat);
    add_case("PEW_shift_east_2deg", "deterministic", ...
        orig.PEW.lon + 2, orig.PEW.lat, orig.PSUW.lon, orig.PSUW.lat);
    add_case("PEW_shift_west_2deg", "deterministic", ...
        orig.PEW.lon - 2, orig.PEW.lat, orig.PSUW.lon, orig.PSUW.lat);

    % Shift PSUW only by 2 degrees.
    add_case("PSUW_shift_north_2deg", "deterministic", ...
        orig.PEW.lon, orig.PEW.lat, orig.PSUW.lon, orig.PSUW.lat + 2);
    add_case("PSUW_shift_south_2deg", "deterministic", ...
        orig.PEW.lon, orig.PEW.lat, orig.PSUW.lon, orig.PSUW.lat - 2);
    add_case("PSUW_shift_east_2deg", "deterministic", ...
        orig.PEW.lon, orig.PEW.lat, orig.PSUW.lon + 2, orig.PSUW.lat);
    add_case("PSUW_shift_west_2deg", "deterministic", ...
        orig.PEW.lon, orig.PEW.lat, orig.PSUW.lon - 2, orig.PSUW.lat);

    % Random boundary perturbations within +/-2 degrees of each original
    % box boundary. These boxes are all contained in the expanded downloads.
    for ii = 1:n_random

        pew_lon = sort(orig.PEW.lon + (rand(1,2)*4 - 2));
        pew_lat = sort(orig.PEW.lat + (rand(1,2)*4 - 2));
        psu_lon = sort(orig.PSUW.lon + (rand(1,2)*4 - 2));
        psu_lat = sort(orig.PSUW.lat + (rand(1,2)*4 - 2));

        % Avoid extremely narrow boxes from unlucky random draws.
        if diff(pew_lon) < 4, pew_lon = orig.PEW.lon; end
        if diff(pew_lat) < 4, pew_lat = orig.PEW.lat; end
        if diff(psu_lon) < 4, psu_lon = orig.PSUW.lon; end
        if diff(psu_lat) < 4, psu_lat = orig.PSUW.lat; end

        add_case(sprintf("random_%03d",ii), "random", pew_lon, pew_lat, psu_lon, psu_lat);
    end
end

function src = local_source_profile(files, box, sgrid)

    sumT = [];
    sumS = [];
    cntT = [];
    cntS = [];
    depth = [];

    for f = 1:numel(files)

        fn = files{f};

        lon = local_read_coord(fn, ["longitude","lon"]);
        lat = local_read_coord(fn, ["latitude","lat"]);
        dep = local_read_coord(fn, ["depth","deptht","lev"]);

        if isempty(depth), depth = double(dep(:)); end

        ilon = local_lon_indices(lon, box.lon);
        ilat = find(lat >= min(box.lat) & lat <= max(box.lat));

        if isempty(ilon) || isempty(ilat)
            error('No grid cells found for %s in file %s.', box.name, fn);
        end

        theta = local_read_subset(fn, 'thetao', ilon, ilat);
        salt  = local_read_subset(fn, 'so',     ilon, ilat);

        [profT_sum, profT_cnt] = local_sum_profile(theta, fn, 'thetao');
        [profS_sum, profS_cnt] = local_sum_profile(salt,  fn, 'so');

        if isempty(sumT)
            sumT = zeros(size(profT_sum));
            sumS = zeros(size(profS_sum));
            cntT = zeros(size(profT_cnt));
            cntS = zeros(size(profS_cnt));
        end

        sumT = sumT + profT_sum;
        sumS = sumS + profS_sum;
        cntT = cntT + profT_cnt;
        cntS = cntS + profS_cnt;
    end

    thetao = sumT ./ cntT;
    so     = sumS ./ cntS;

    thetao(cntT==0) = NaN;
    so(cntS==0)     = NaN;

    % Convert Practical Salinity to Absolute Salinity using pressure, then
    % calculate sigma0. GLORYS thetao is treated as potential temperature.
    depth = depth(:);
    thetao = thetao(:);
    so = so(:);

    p = gsw_p_from_z(-depth, box.lat0);
    SA = gsw_SA_from_SP(so, p, ...
        box.lon0.*ones(size(depth)), box.lat0.*ones(size(depth)));
    CT = gsw_CT_from_pt(SA, thetao);
    sig = gsw_sigma0(SA, CT);

    % Sort by density before interpolation so small density inversions do not
    % violate interp1's requirement for monotonic sample points.
    msk = isfinite(sig) & isfinite(thetao) & isfinite(SA);
    sig_tmp = sig(msk);
    t_tmp = thetao(msk);
    SA_tmp = SA(msk);
    [sig_tmp, ii] = sort(sig_tmp);
    t_tmp = t_tmp(ii);
    SA_tmp = SA_tmp(ii);
    [sigU, ia] = unique(sig_tmp);
    t_tmp = t_tmp(ia);
    SA_tmp = SA_tmp(ia);

    src.name = box.name;
    src.box = box;
    src.depth = depth;
    src.thetao_z = thetao;
    src.SP_z = so;
    src.SA_z = SA;
    src.sigma_z = sig;
    src.thetao = interp1(sigU, t_tmp,  sgrid, 'linear', NaN);
    src.SA     = interp1(sigU, SA_tmp, sgrid, 'linear', NaN);
    src.sgrid  = sgrid;
end

function X = local_read_subset(fn, varname, ilon, ilat)
    % Assumes standard Copernicus variable dims include longitude and latitude.
    info = ncinfo(fn, varname);
    dimnames = string({info.Dimensions.Name});

    start = ones(1,numel(dimnames));
    count = [info.Dimensions.Length];

    lonDim = find(ismember(lower(dimnames), ["longitude","lon"]),1);
    latDim = find(ismember(lower(dimnames), ["latitude","lat"]),1);

    start(lonDim) = min(ilon);
    count(lonDim) = numel(ilon);
    start(latDim) = min(ilat);
    count(latDim) = numel(ilat);

    % This assumes ilon/ilat are contiguous, which they are for these boxes.
    X = double(ncread(fn, varname, start, count));

    % If longitude wraps and indices are non-contiguous, fall back to full read.
    if any(diff(ilon) ~= 1) || any(diff(ilat) ~= 1)
        Xall = double(ncread(fn, varname));
        idx = repmat({':'},1,ndims(Xall));
        idx{lonDim} = ilon;
        idx{latDim} = ilat;
        X = Xall(idx{:});
    end
end

function [prof_sum, prof_cnt] = local_sum_profile(X, fn, varname)

    info = ncinfo(fn,varname);
    dimnames = string({info.Dimensions.Name});
    depDim = find(ismember(lower(dimnames), ["depth","deptht","lev"]),1);

    % Move depth to first dimension, collapse all other dimensions.
    order = [depDim setdiff(1:ndims(X), depDim)];
    Xp = permute(X, order);
    Xp = reshape(Xp, size(Xp,1), []);

    prof_sum = sum(Xp, 2, 'omitnan');
    prof_cnt = sum(isfinite(Xp), 2);
end

function idx = local_lon_indices(lon, boxlon)
    lon = double(lon(:));
    b = double(boxlon(:))';

    if max(lon) > 180
        b = mod(b,360);
    end

    if b(1) <= b(2)
        idx = find(lon >= b(1) & lon <= b(2));
    else
        idx = find(lon >= b(1) | lon <= b(2));
    end
end

function coord = local_read_coord(fn, names)
    coord = [];
    info = ncinfo(fn);
    vnames = string({info.Variables.Name});
    for n = names
        ii = find(lower(vnames) == lower(n),1);
        if ~isempty(ii)
            coord = double(ncread(fn, vnames(ii)));
            return
        end
    end
    error('Could not find coordinate variable: %s', strjoin(names,", "));
end

function fPEW = local_calc_fpew_profiles(STA, PEWsrc, PSUWsrc, sgrid)

    [nz,np] = size(STA.temp);
    fPEW = NaN(nz,np);

    for i = 1:nz
        for j = 1:np

            st = STA.s_t(i,j);
            if ~isfinite(st), continue; end

            psu_t = interp1(sgrid, PSUWsrc.thetao, st, 'linear', NaN);
            psu_s = interp1(sgrid, PSUWsrc.SA,     st, 'linear', NaN);
            pew_t = interp1(sgrid, PEWsrc.thetao,  st, 'linear', NaN);
            pew_s = interp1(sgrid, PEWsrc.SA,      st, 'linear', NaN);

            if ~all(isfinite([psu_t psu_s pew_t pew_s])), continue; end

            T = STA.temp(i,j);
            S = STA.SA(i,j);
            if ~all(isfinite([T S])), continue; end

            dT   = psu_t - pew_t;
            dS   = psu_s - pew_s;
            dT_s = T - pew_t;
            dS_s = S - pew_s;

            denom = dT.^2 + dS.^2;
            if ~isfinite(denom) || denom <= eps
                continue
            end

            alphaPSUW = (dT_s*dT + dS_s*dS) ./ denom;
            alphaPSUW = max(min(alphaPSUW,1),0);

            fPEW(i,j) = 1 - alphaPSUW;
        end
    end
end

function y = local_interp_to_sigma(sig, field, target_sigma)

    [~,np] = size(sig);
    y = NaN(1,np);

    for j = 1:np
        x = sig(:,j);
        v = field(:,j);
        m = isfinite(x) & isfinite(v);

        if sum(m) >= 2
            x_tmp = x(m);
            v_tmp = v(m);
            [x_tmp, ii] = sort(x_tmp);
            v_tmp = v_tmp(ii);
            [xu, ia] = unique(x_tmp);
            vu = v_tmp(ia);
            if numel(xu) >= 2
                y(j) = interp1(xu, vu, target_sigma, 'linear', NaN);
            end
        end
    end
end

function S = local_station_on_sigma(STA, target_sigma, stationName)

    S.station = string(stationName);
    S.time = local_to_datenum(STA.mtime);
    S.year = year(datetime(S.time,'ConvertFrom','datenum'));
    S.month = month(datetime(S.time,'ConvertFrom','datenum'));

    S.O2 = local_interp_oxygen(STA, target_sigma);
    S.Nox = local_interp_field(STA.s_t, STA.no3, target_sigma);

    if isfield(STA,'po4')
        S.PO4 = local_interp_field(STA.s_t, STA.po4, target_sigma);
    else
        S.PO4 = NaN(size(S.O2));
    end

    if isfield(STA,'spice')
        S.spice = local_interp_field(STA.s_t, STA.spice, target_sigma);
    else
        S.spice = NaN(size(S.O2));
    end

    S.NO = S.O2 + 8.625 .* S.Nox;
end

function O2 = local_interp_oxygen(STA, target_sigma)

    [~,np] = size(STA.s_t);
    O2field = NaN(size(STA.s_t));

    if isfield(STA,'oxML') && isfield(STA,'oxML2') && isfield(STA,'ox') && isfield(STA,'ox2')
        O2field = mean(cat(3,STA.oxML,STA.oxML2,STA.ox,STA.ox2),3,'omitnan');
    elseif isfield(STA,'ox')
        O2field = STA.ox;
    else
        error('No recognized oxygen field found.');
    end

    O2 = NaN(1,np);
    for j = 1:np
        O2(j) = local_interp_vector(STA.s_t(:,j), O2field(:,j), target_sigma);
    end
end

function y = local_interp_field(sig, field, target_sigma)
    [~,np] = size(sig);
    y = NaN(1,np);
    for j = 1:np
        y(j) = local_interp_vector(sig(:,j), field(:,j), target_sigma);
    end
end

function yi = local_interp_vector(x, y, xi)
    m = isfinite(x) & isfinite(y);
    yi = NaN;
    if sum(m) >= 2
        x_tmp = x(m);
        y_tmp = y(m);
        [x_tmp, ii] = sort(x_tmp);
        y_tmp = y_tmp(ii);
        [xu, ia] = unique(x_tmp);
        yu = y_tmp(ia);
        if numel(xu) >= 2
            yi = interp1(xu, yu, xi, 'linear', NaN);
        end
    end
end

function core = local_core_properties(STA, fPEW_full, core_range)

    % Reproduce the manuscript PEW-core definition as closely as possible:
    % profile-specific maximum fPEW within sigma_theta = 26.25--26.70,
    % sufficient profile depth, CTD profiles only where profile_type exists,
    % followed by a 3-scaled-MAD filter on core density. The median of the
    % retained profile-specific core densities is reported.

    nzmax = min(500,size(fPEW_full,1));
    np = size(fPEW_full,2);

    core_sigma = NaN(1,np);
    core_depth = NaN(1,np);
    core_fPEW  = NaN(1,np);

    for j = 1:np

        last_valid = find(isfinite(fPEW_full(1:nzmax,j)),1,'last');
        if isempty(last_valid) || last_valid <= 175
            continue
        end

        if isfield(STA,'profile_type') && STA.profile_type(j) ~= 1
            continue
        end

        st = STA.s_t(1:nzmax,j);
        ff = fPEW_full(1:nzmax,j);

        m = isfinite(st) & isfinite(ff) & ...
            st >= core_range(1) & st <= core_range(2);

        if ~any(m)
            continue
        end

        inds = find(m);
        [~,jj] = max(ff(m));
        idx = inds(jj);

        core_sigma(j) = st(idx);
        core_fPEW(j) = ff(idx);

        % Use pressure as the vertical coordinate where available; otherwise
        % fall back to row index (the manuscript station grids are ~1 m/dbar).
        if isfield(STA,'pr') && size(STA.pr,1) >= idx && size(STA.pr,2) >= j
            core_depth(j) = STA.pr(idx,j);
        else
            core_depth(j) = idx;
        end
    end

    med_sig = median(core_sigma,'omitnan');
    scaled_mad = 1.4826 .* median(abs(core_sigma-med_sig),'omitnan');

    if isfinite(scaled_mad) && scaled_mad > 0
        bad = abs(core_sigma-med_sig) > 3.*scaled_mad;
        core_sigma(bad) = NaN;
        core_depth(bad) = NaN;
        core_fPEW(bad) = NaN;
    end

    core.sigma = median(core_sigma,'omitnan');
    core.depth = median(core_depth,'omitnan');
    core.fPEW  = median(core_fPEW,'omitnan');
    core.n     = sum(isfinite(core_sigma));
    core.sigma_all = core_sigma;
    core.depth_all = core_depth;
    core.fPEW_all  = core_fPEW;
end

function M = local_station_metrics(S, fPEW, event_prctile)

    M = struct();

    M.mean_fPEW = mean(fPEW,'omitnan');
    M.std_fPEW  = std(fPEW,'omitnan');
    M.n_fPEW    = sum(isfinite(fPEW));

    [M.slope_O2, M.int_O2, M.R2_O2, M.p_O2, M.n_O2] = local_regress(fPEW, S.O2);
    [M.slope_Nox,M.int_Nox,M.R2_Nox,M.p_Nox,M.n_Nox] = local_regress(fPEW, S.Nox);
    [M.slope_NO, M.int_NO, M.R2_NO, M.p_NO, M.n_NO] = local_regress(fPEW, S.NO);

    M.slope_O2_per_0p1  = 0.1 .* M.slope_O2;
    M.slope_Nox_per_0p1 = 0.1 .* M.slope_Nox;
    M.slope_NO_per_0p1  = 0.1 .* M.slope_NO;

    % Residual O'-N' slope.
    Ohat = polyval([M.slope_O2 M.int_O2], fPEW);
    Nhat = polyval([M.slope_Nox M.int_Nox], fPEW);
    Ores = S.O2  - Ohat;
    Nres = S.Nox - Nhat;

    [M.slope_ON_resid, M.int_ON_resid, M.R2_ON_resid, M.p_ON_resid, M.n_ON_resid] = ...
        local_regress(Nres, Ores);

    % Low-O2 event attribution. The O2 threshold is defined from all valid
    % station O2 values, while the decomposition uses a common event sample
    % with both O2 and fPEW available.
    Ocrit = prctile(S.O2(isfinite(S.O2)), event_prctile);
    evt = S.O2 < Ocrit & isfinite(S.O2) & isfinite(fPEW);

    Omean = mean(S.O2,'omitnan');
    O2_anom_total = S.O2 - Omean;
    O2_anom_PEW   = Ohat - Omean;
    O2_anom_resid = Ores;

    M.lowO_threshold = Ocrit;
    M.lowO_n = sum(evt);
    M.lowO_event_anom = mean(O2_anom_total(evt),'omitnan');
    M.lowO_PEW_comp   = mean(O2_anom_PEW(evt),'omitnan');
    M.lowO_resid_comp = mean(O2_anom_resid(evt),'omitnan');
    M.lowO_frac_PEW   = M.lowO_PEW_comp ./ M.lowO_event_anom;
    M.lowO_frac_resid = M.lowO_resid_comp ./ M.lowO_event_anom;
end

function [slope,intcpt,R2,p,n] = local_regress(x,y)

    x = x(:); y = y(:);
    m = isfinite(x) & isfinite(y);
    n = sum(m);

    slope = NaN; intcpt = NaN; R2 = NaN; p = NaN;

    if n >= 3 && numel(unique(x(m))) > 1
        pp = polyfit(x(m), y(m), 1);
        slope = pp(1);
        intcpt = pp(2);

        [R,P] = corrcoef(x(m),y(m));
        R2 = R(1,2).^2;
        p = P(1,2);
    end
end

function A = local_annual_means(years, x)

    years = years(:);
    x = x(:);

    valid_years = years(isfinite(years));
    if isempty(valid_years)
        A.year = []; A.mean = []; A.n = [];
        return
    end

    uy = (min(valid_years):max(valid_years))';
    mn = NaN(size(uy));
    nn = NaN(size(uy));

    for k = 1:numel(uy)
        idx = years == uy(k);
        mn(k) = mean(x(idx),'omitnan');
        nn(k) = sum(isfinite(x(idx)));
    end

    A.year = uy;
    A.mean = mn;
    A.n = nn;
end

function tr = local_linear_trend(years, x)

    m = isfinite(years) & isfinite(x);
    tr.slope = NaN;
    tr.intercept = NaN;
    tr.R2 = NaN;
    tr.p = NaN;
    tr.n = sum(m);

    if sum(m) >= 3
        pp = polyfit(years(m), x(m), 1);
        tr.slope = pp(1);
        tr.intercept = pp(2);
        [R,~] = corrcoef(years(m), x(m));
        tr.R2 = R(1,2).^2;

        % Match the manuscript: OLS for trend magnitude, Mann--Kendall for
        % significance of a monotonic trend.
        [~,tr.p] = mann_kendall(x(m),0.1);
    end
end

function [R2,p,n] = local_common_year_corr(y1,x1,y2,x2)

    [yc,ia,ib] = intersect(y1,y2);
    x = x1(ia);
    y = x2(ib);

    m = isfinite(x) & isfinite(y);
    n = sum(m);
    R2 = NaN; p = NaN;

    if n >= 3
        [R,P] = corrcoef(x(m),y(m));
        R2 = R(1,2).^2;
        p = P(1,2);
    end
end

function D = local_diff_stats(x, xb)

    m = isfinite(x) & isfinite(xb);
    if ~any(m)
        D.bias = NaN; D.mad = NaN; D.rmse = NaN; D.maxabs = NaN;
        return
    end

    dx = x(m)-xb(m);
    D.bias = mean(dx,'omitnan');
    D.mad  = mean(abs(dx),'omitnan');
    D.rmse = sqrt(mean(dx.^2,'omitnan'));
    D.maxabs = max(abs(dx),[],'omitnan');
end

function T = local_results_table(R)

    n = numel(R);

    T = table;
    for i = 1:n
        T.case(i,1) = string(R(i).case);
        T.case_type(i,1) = string(R(i).case_type);

        T.PEW_lon_min(i,1) = R(i).PEW_lon(1);
        T.PEW_lon_max(i,1) = R(i).PEW_lon(2);
        T.PEW_lat_min(i,1) = R(i).PEW_lat(1);
        T.PEW_lat_max(i,1) = R(i).PEW_lat(2);

        T.PSUW_lon_min(i,1) = R(i).PSUW_lon(1);
        T.PSUW_lon_max(i,1) = R(i).PSUW_lon(2);
        T.PSUW_lat_min(i,1) = R(i).PSUW_lat(1);
        T.PSUW_lat_max(i,1) = R(i).PSUW_lat(2);

        T.CS_core_sigma(i,1) = R(i).CS_core_sigma;
        T.P4_core_sigma(i,1) = R(i).P4_core_sigma;

        % fPEW propagation.
        T.CS_mean_fPEW(i,1) = R(i).CS.mean_fPEW;
        T.P4_mean_fPEW(i,1) = R(i).P4.mean_fPEW;
        T.CS_delta_fPEW_bias(i,1) = R(i).delta_fPEW_CS.bias;
        T.P4_delta_fPEW_bias(i,1) = R(i).delta_fPEW_P4.bias;
        T.CS_delta_fPEW_rmse(i,1) = R(i).delta_fPEW_CS.rmse;
        T.P4_delta_fPEW_rmse(i,1) = R(i).delta_fPEW_P4.rmse;

        % O2-fPEW relationship.
        T.CS_O2_slope_per_0p1(i,1) = R(i).CS.slope_O2_per_0p1;
        T.P4_O2_slope_per_0p1(i,1) = R(i).P4.slope_O2_per_0p1;
        T.CS_O2_R2(i,1) = R(i).CS.R2_O2;
        T.P4_O2_R2(i,1) = R(i).P4.R2_O2;
        T.CS_O2_p(i,1) = R(i).CS.p_O2;
        T.P4_O2_p(i,1) = R(i).P4.p_O2;

        % Nutrient and NO relationships.
        T.CS_Nox_slope_per_0p1(i,1) = R(i).CS.slope_Nox_per_0p1;
        T.P4_Nox_slope_per_0p1(i,1) = R(i).P4.slope_Nox_per_0p1;
        T.CS_Nox_R2(i,1) = R(i).CS.R2_Nox;
        T.P4_Nox_R2(i,1) = R(i).P4.R2_Nox;

        T.CS_NO_slope_per_0p1(i,1) = R(i).CS.slope_NO_per_0p1;
        T.P4_NO_slope_per_0p1(i,1) = R(i).P4.slope_NO_per_0p1;
        T.CS_NO_R2(i,1) = R(i).CS.R2_NO;
        T.P4_NO_R2(i,1) = R(i).P4.R2_NO;

        % Trend and trend propagation.
        T.CS_fPEW_trend_decade(i,1) = 10 .* R(i).trendCS.slope;
        T.P4_fPEW_trend_decade(i,1) = 10 .* R(i).trendP4.slope;
        T.CS_fPEW_trend_p(i,1) = R(i).trendCS.p;
        T.P4_fPEW_trend_p(i,1) = R(i).trendP4.p;
        T.CS_O2trend_from_fPEW_decade(i,1) = R(i).O2trendCS_from_fPEW_decade;
        T.P4_O2trend_from_fPEW_decade(i,1) = R(i).O2trendP4_from_fPEW_decade;

        % Regional coherence.
        T.P4_CS_fPEW_annual_R2(i,1) = R(i).interstation_R2;
        T.P4_CS_fPEW_annual_p(i,1) = R(i).interstation_p;
        T.P4_CS_fPEW_annual_n(i,1) = R(i).interstation_n;

        % Low-O2 event attribution.
        T.CS_lowO_frac_PEW(i,1) = R(i).CS.lowO_frac_PEW;
        T.P4_lowO_frac_PEW(i,1) = R(i).P4.lowO_frac_PEW;
        T.CS_lowO_PEW_comp(i,1) = R(i).CS.lowO_PEW_comp;
        T.P4_lowO_PEW_comp(i,1) = R(i).P4.lowO_PEW_comp;
        T.CS_lowO_resid_comp(i,1) = R(i).CS.lowO_resid_comp;
        T.P4_lowO_resid_comp(i,1) = R(i).P4.lowO_resid_comp;

        % Residual O'-N' relationship.
        T.CS_ON_resid_slope(i,1) = R(i).CS.slope_ON_resid;
        T.P4_ON_resid_slope(i,1) = R(i).P4.slope_ON_resid;
        T.CS_ON_resid_R2(i,1) = R(i).CS.R2_ON_resid;
        T.P4_ON_resid_R2(i,1) = R(i).P4.R2_ON_resid;
    end
end

function local_make_figures(T, fig_dir)

    isBase = T.case == "baseline_10x10";
    isRand = T.case_type == "random";
    isDet  = T.case_type == "deterministic";

    % 1. O2-fPEW slopes.
    f = figure('Color','w','Units','centimeters','Position',[2 2 16 8]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on;
    histogram(T.P4_O2_slope_per_0p1(isRand), 20);
    xline(T.P4_O2_slope_per_0p1(isBase),'k','LineWidth',1.5);
    scatter(T.P4_O2_slope_per_0p1(isDet), zeros(sum(isDet),1), 30, 'filled');
    xlabel('P4 O_2 slope per 0.1 f_{PEW} (\mumol kg^{-1})');
    ylabel('Count'); grid on; box on;

    nexttile; hold on;
    histogram(T.CS_O2_slope_per_0p1(isRand), 20);
    xline(T.CS_O2_slope_per_0p1(isBase),'k','LineWidth',1.5);
    scatter(T.CS_O2_slope_per_0p1(isDet), zeros(sum(isDet),1), 30, 'filled');
    xlabel('CS-OS O_2 slope per 0.1 f_{PEW} (\mumol kg^{-1})');
    ylabel('Count'); grid on; box on;

    exportgraphics(f, fullfile(fig_dir,'sensitivity_O2_fPEW_slopes.png'), 'Resolution',300);

    % 2. R2 values.
    f = figure('Color','w','Units','centimeters','Position',[2 2 16 8]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on;
    histogram(T.P4_O2_R2(isRand), 20);
    xline(T.P4_O2_R2(isBase),'k','LineWidth',1.5);
    xlabel('P4 O_2-f_{PEW} R^2');
    ylabel('Count'); grid on; box on;

    nexttile; hold on;
    histogram(T.CS_O2_R2(isRand), 20);
    xline(T.CS_O2_R2(isBase),'k','LineWidth',1.5);
    xlabel('CS-OS O_2-f_{PEW} R^2');
    ylabel('Count'); grid on; box on;

    exportgraphics(f, fullfile(fig_dir,'sensitivity_O2_fPEW_R2.png'), 'Resolution',300);

    % 3. P4 trend propagation.
    f = figure('Color','w','Units','centimeters','Position',[2 2 16 8]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on;
    histogram(T.P4_fPEW_trend_decade(isRand), 20);
    xline(T.P4_fPEW_trend_decade(isBase),'k','LineWidth',1.5);
    xlabel('P4 f_{PEW} trend (decade^{-1})');
    ylabel('Count'); grid on; box on;

    nexttile; hold on;
    histogram(T.P4_O2trend_from_fPEW_decade(isRand), 20);
    xline(T.P4_O2trend_from_fPEW_decade(isBase),'k','LineWidth',1.5);
    xlabel('P4 O_2 trend implied by f_{PEW} (\mumol kg^{-1} decade^{-1})');
    ylabel('Count'); grid on; box on;

    exportgraphics(f, fullfile(fig_dir,'sensitivity_P4_trend_propagation.png'), 'Resolution',300);

    % 4. Low-O2 event attribution.
    f = figure('Color','w','Units','centimeters','Position',[2 2 16 8]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on;
    histogram(100*T.P4_lowO_frac_PEW(isRand), 20);
    xline(100*T.P4_lowO_frac_PEW(isBase),'k','LineWidth',1.5);
    xlabel('P4 low-O_2 anomaly attributed to PEW (%)');
    ylabel('Count'); grid on; box on;

    nexttile; hold on;
    histogram(100*T.CS_lowO_frac_PEW(isRand), 20);
    xline(100*T.CS_lowO_frac_PEW(isBase),'k','LineWidth',1.5);
    xlabel('CS-OS low-O_2 anomaly attributed to PEW (%)');
    ylabel('Count'); grid on; box on;

    exportgraphics(f, fullfile(fig_dir,'sensitivity_lowO_attribution.png'), 'Resolution',300);

    % 5. Residual O'-N' slopes.
    f = figure('Color','w','Units','centimeters','Position',[2 2 16 8]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on;
    histogram(T.P4_ON_resid_slope(isRand), 20);
    xline(T.P4_ON_resid_slope(isBase),'k','LineWidth',1.5);
    xline(-8.625,'--k');
    xlabel('P4 residual O''-N'' slope');
    ylabel('Count'); grid on; box on;

    nexttile; hold on;
    histogram(T.CS_ON_resid_slope(isRand), 20);
    xline(T.CS_ON_resid_slope(isBase),'k','LineWidth',1.5);
    xline(-8.625,'--k');
    xlabel('CS-OS residual O''-N'' slope');
    ylabel('Count'); grid on; box on;

    exportgraphics(f, fullfile(fig_dir,'sensitivity_residual_ON_slopes.png'), 'Resolution',300);
end

function dn = local_to_datenum(t)
    if isdatetime(t)
        dn = datenum(t);
    else
        dn = double(t);
    end
end
