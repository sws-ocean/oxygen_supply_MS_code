%% Recreate TnK using 5-year sliding-window endmembers

addpath(genpath('/Users/samst/Dropbox/Hakai/'));
addpath(genpath('/Users/samst/Dropbox/UBC/GSW/'));
addpath(genpath('/Users/samst/Dropbox/UBC/m_map/'));
addpath(genpath('/Users/samst/Dropbox/UBC/Misc/'));

%%
clear


%% Load ORAS5 source-water data
% Monthly spatially averaged T/S profiles
%
% PEW:
% Longitude: 113 W to 103 W
% Latitude:  0 to 10 N
%
% PSUW:
% Longitude: 150 W to 140 W
% Latitude:  45 to 55 N

oras_dir = '/Users/samst/Library/CloudStorage/Dropbox/Hakai/gyre/supply/data/ORAS5_OPENDAP';

% Approximate centres of source-water boxes
PEW_lon = -108;
PEW_lat = 5;

SAW_lon = -145;
SAW_lat = 50;

% ============================================================
% PEW
% ============================================================

temp_files = dir(fullfile(oras_dir,'PEW','votemper','PEW_votemper_*.nc'));
salt_files = dir(fullfile(oras_dir,'PEW','vosaline','PEW_vosaline_*.nc'));

% Extract YYYYMM from filenames
temp_dates = regexp({temp_files.name},'(\d{6})','tokens','once');
temp_dates = string(cellfun(@(x) x{1},temp_dates,'UniformOutput',false));

salt_dates = regexp({salt_files.name},'(\d{6})','tokens','once');
salt_dates = string(cellfun(@(x) x{1},salt_dates,'UniformOutput',false));

% Keep only months where both T and S exist
common_dates = intersect(temp_dates,salt_dates);
common_dates = sort(common_dates);

nmonths = length(common_dates);

% Depth grid
PEWz = double(ncread( ...
    fullfile(temp_files(1).folder,temp_files(1).name), ...
    'deptht'));

nz = length(PEWz);

PEWt_monthly = nan(nz,nmonths);
PEWs_monthly = nan(nz,nmonths);
PEW.time = NaT(nmonths,1);

for i = 1:nmonths

    datecode = common_dates(i);

    tfile = temp_files(temp_dates == datecode);
    sfile = salt_files(salt_dates == datecode);

    yr = str2double(extractBetween(datecode,1,4));
    mo = str2double(extractBetween(datecode,5,6));

    PEW.time(i) = datetime(yr,mo,15);

    T = double(ncread( ...
        fullfile(tfile.folder,tfile.name), ...
        'votemper'));

    T(abs(T) > 1e20) = NaN;

    PEWt_monthly(:,i) = squeeze( ...
        mean(T,[1 2],'omitnan'));

    S = double(ncread( ...
        fullfile(sfile.folder,sfile.name), ...
        'vosaline'));

    S(abs(S) > 1e20) = NaN;

    PEWs_monthly(:,i) = squeeze( ...
        mean(S,[1 2],'omitnan'));

end

% ============================================================
% PSUW
% ============================================================

temp_files = dir(fullfile(oras_dir,'PSUW','votemper','PSUW_votemper_*.nc'));
salt_files = dir(fullfile(oras_dir,'PSUW','vosaline','PSUW_vosaline_*.nc'));

temp_dates = regexp({temp_files.name},'(\d{6})','tokens','once');
temp_dates = string(cellfun(@(x) x{1},temp_dates,'UniformOutput',false));

salt_dates = regexp({salt_files.name},'(\d{6})','tokens','once');
salt_dates = string(cellfun(@(x) x{1},salt_dates,'UniformOutput',false));

common_dates = intersect(temp_dates,salt_dates);
common_dates = sort(common_dates);

nmonths = length(common_dates);

SAWz = double(ncread( ...
    fullfile(temp_files(1).folder,temp_files(1).name), ...
    'deptht'));

nz = length(SAWz);

SAWt_monthly = nan(nz,nmonths);
SAWs_monthly = nan(nz,nmonths);
SAW.time = NaT(nmonths,1);

for i = 1:nmonths

    datecode = common_dates(i);

    tfile = temp_files(temp_dates == datecode);
    sfile = salt_files(salt_dates == datecode);

    yr = str2double(extractBetween(datecode,1,4));
    mo = str2double(extractBetween(datecode,5,6));

    SAW.time(i) = datetime(yr,mo,15);

    T = double(ncread( ...
        fullfile(tfile.folder,tfile.name), ...
        'votemper'));

    T(abs(T) > 1e20) = NaN;

    SAWt_monthly(:,i) = squeeze( ...
        mean(T,[1 2],'omitnan'));

    S = double(ncread( ...
        fullfile(sfile.folder,sfile.name), ...
        'vosaline'));

    S(abs(S) > 1e20) = NaN;

    SAWs_monthly(:,i) = squeeze( ...
        mean(S,[1 2],'omitnan'));

end

%% Raw ORAS5 grid-point fPEW at the PEW-core density
%
% Purpose:
%   Calculate fPEW independently at every ORAS5 horizontal grid point
%   in the PEW and PSUW source boxes for every month.
%
% Important:
%   These values are used to show within-box spatial variability only.
%   The primary source-water time series remains based on the monthly
%   spatially averaged source-water profiles calculated above.
%
% Output:
%   PEW_grid_fPEW   [grid point x month]
%   SAW_grid_fPEW   [grid point x month]
%
%   PEW_grid_time   datetime vector corresponding to columns
%   SAW_grid_time   datetime vector corresponding to columns

% fprintf('\nCalculating raw ORAS5 grid-point fPEW...\n');
% 
% core_st = 26.50;
% 
% % ============================================================
% % Fixed PEW--PSUW mixing line evaluated at core density
% % ============================================================
% 
% SAW_t_core = interp1(sgrid,SAWt,core_st);
% SAW_s_core = interp1(sgrid,SAWs,core_st);
% 
% PEW_t_core = interp1(sgrid,PEWt,core_st);
% PEW_s_core = interp1(sgrid,PEWs,core_st);
% 
% dT_core = SAW_t_core - PEW_t_core;
% dS_core = SAW_s_core - PEW_s_core;
% 
% denom_core = dT_core.^2 + dS_core.^2;
% 
% if ~isfinite(denom_core) || denom_core == 0
%     error('Invalid fixed PEW--PSUW mixing vector at sigma=%.2f.',core_st);
% end
% 
% 
% % ============================================================
% % PEW source region
% % ============================================================
% 
% temp_files_PEW = dir(fullfile(oras_dir,'PEW','votemper','PEW_votemper_*.nc'));
% salt_files_PEW = dir(fullfile(oras_dir,'PEW','vosaline','PEW_vosaline_*.nc'));
% 
% temp_dates_PEW = regexp({temp_files_PEW.name},'(\d{6})','tokens','once');
% temp_dates_PEW = string(cellfun(@(x) x{1},temp_dates_PEW,'UniformOutput',false));
% 
% salt_dates_PEW = regexp({salt_files_PEW.name},'(\d{6})','tokens','once');
% salt_dates_PEW = string(cellfun(@(x) x{1},salt_dates_PEW,'UniformOutput',false));
% 
% common_dates_PEW = sort(intersect(temp_dates_PEW,salt_dates_PEW));
% 
% nmonths_PEW = numel(common_dates_PEW);
% 
% % Read first file to determine horizontal grid
% firstfile = fullfile(temp_files_PEW(1).folder,temp_files_PEW(1).name);
% 
% [PEW_grid_lon,PEW_grid_lat] = get_oras_lonlat(firstfile);
% 
% PEW_grid_lon = double(PEW_grid_lon);
% PEW_grid_lat = double(PEW_grid_lat);
% 
% % Ensure lon/lat are 2-D arrays matching horizontal model grid
% if isvector(PEW_grid_lon) && isvector(PEW_grid_lat)
%     [PEW_grid_lon,PEW_grid_lat] = ndgrid(PEW_grid_lon,PEW_grid_lat);
% end
% 
% nxy_PEW = numel(PEW_grid_lon);
% 
% PEW_grid_fPEW = nan(nxy_PEW,nmonths_PEW);
% PEW_grid_time = NaT(nmonths_PEW,1);
% 
% PEWz_raw = double(ncread(firstfile,'deptht'));
% 
% for im = 1:nmonths_PEW
% 
%     datecode = common_dates_PEW(im);
% 
%     tfile = temp_files_PEW(temp_dates_PEW == datecode);
%     sfile = salt_files_PEW(salt_dates_PEW == datecode);
% 
%     yr = str2double(extractBetween(datecode,1,4));
%     mo = str2double(extractBetween(datecode,5,6));
% 
%     PEW_grid_time(im) = datetime(yr,mo,15);
% 
%     T = double(ncread(fullfile(tfile.folder,tfile.name),'votemper'));
%     S = double(ncread(fullfile(sfile.folder,sfile.name),'vosaline'));
% 
%     T(abs(T) > 1e20) = NaN;
%     S(abs(S) > 1e20) = NaN;
% 
%     % Expected dimensions are x x y x depth.
%     % Flatten horizontal dimensions while retaining depth.
%     T2 = reshape(T,[],size(T,ndims(T)));
%     S2 = reshape(S,[],size(S,ndims(S)));
% 
%     % If depth is not the final dimension, move it there
%     if size(T2,2) ~= numel(PEWz_raw)
% 
%         depth_dim = find(size(T) == numel(PEWz_raw),1,'last');
% 
%         ord = 1:ndims(T);
%         ord([depth_dim end]) = ord([end depth_dim]);
% 
%         T = permute(T,ord);
%         S = permute(S,ord);
% 
%         T2 = reshape(T,[],numel(PEWz_raw));
%         S2 = reshape(S,[],numel(PEWz_raw));
% 
%     end
% 
%     lonvec = PEW_grid_lon(:);
%     latvec = PEW_grid_lat(:);
% 
%     for ig = 1:nxy_PEW
% 
%         tprof = T2(ig,:)';
%         sprof = S2(ig,:)';
% 
%         if sum(isfinite(tprof) & isfinite(sprof)) < 2
%             continue
%         end
% 
%         lon0 = lonvec(ig);
%         lat0 = latvec(ig);
% 
%         if ~isfinite(lon0) || ~isfinite(lat0)
%             continue
%         end
% 
%         p = gsw_p_from_z(-PEWz_raw,lat0);
% 
%         SA = gsw_SA_from_SP( ...
%             sprof, ...
%             p, ...
%             lon0 .* ones(size(PEWz_raw)), ...
%             lat0 .* ones(size(PEWz_raw)));
% 
%         CT = gsw_CT_from_pt(SA,tprof);
% 
%         st = gsw_sigma0(SA,CT);
% 
%         msk = ...
%             isfinite(st) & ...
%             isfinite(tprof) & ...
%             isfinite(SA);
% 
%         if sum(msk) < 2
%             continue
%         end
% 
%         st_tmp = st(msk);
%         t_tmp  = tprof(msk);
%         s_tmp  = SA(msk);
% 
%         [st_tmp,ii] = sort(st_tmp);
% 
%         t_tmp = t_tmp(ii);
%         s_tmp = s_tmp(ii);
% 
%         [st_tmp,ia] = unique(st_tmp);
% 
%         t_tmp = t_tmp(ia);
%         s_tmp = s_tmp(ia);
% 
%         if numel(st_tmp) < 2 || ...
%            core_st < min(st_tmp) || ...
%            core_st > max(st_tmp)
%             continue
%         end
% 
%         t_core = interp1(st_tmp,t_tmp,core_st);
%         s_core = interp1(st_tmp,s_tmp,core_st);
% 
%         if ~isfinite(t_core) || ~isfinite(s_core)
%             continue
%         end
% 
%         % Project onto fixed PEW -> PSUW mixing vector
%         dT_s = t_core - PEW_t_core;
%         dS_s = s_core - PEW_s_core;
% 
%         alpha_PSUW = ...
%             (dT_s.*dT_core + dS_s.*dS_core) ./ ...
%             denom_core;
% 
%         % Raw, unclipped PEW-axis coordinate
%         PEW_grid_fPEW(ig,im) = 1 - alpha_PSUW;
% 
%     end
% 
%     fprintf('PEW %s: %d valid grid points\n', ...
%         datestr(PEW_grid_time(im),'yyyy-mm'), ...
%         sum(isfinite(PEW_grid_fPEW(:,im))));
% 
% end
% 
% 
% % ============================================================
% % PSUW source region
% % ============================================================
% 
% temp_files_SAW = dir(fullfile(oras_dir,'PSUW','votemper','PSUW_votemper_*.nc'));
% salt_files_SAW = dir(fullfile(oras_dir,'PSUW','vosaline','PSUW_vosaline_*.nc'));
% 
% temp_dates_SAW = regexp({temp_files_SAW.name},'(\d{6})','tokens','once');
% temp_dates_SAW = string(cellfun(@(x) x{1},temp_dates_SAW,'UniformOutput',false));
% 
% salt_dates_SAW = regexp({salt_files_SAW.name},'(\d{6})','tokens','once');
% salt_dates_SAW = string(cellfun(@(x) x{1},salt_dates_SAW,'UniformOutput',false));
% 
% common_dates_SAW = sort(intersect(temp_dates_SAW,salt_dates_SAW));
% 
% nmonths_SAW = numel(common_dates_SAW);
% 
% firstfile = fullfile(temp_files_SAW(1).folder,temp_files_SAW(1).name);
% 
% [SAW_grid_lon,SAW_grid_lat] = get_oras_lonlat(firstfile);
% 
% SAW_grid_lon = double(SAW_grid_lon);
% SAW_grid_lat = double(SAW_grid_lat);
% 
% if isvector(SAW_grid_lon) && isvector(SAW_grid_lat)
%     [SAW_grid_lon,SAW_grid_lat] = ndgrid(SAW_grid_lon,SAW_grid_lat);
% end
% 
% nxy_SAW = numel(SAW_grid_lon);
% 
% SAW_grid_fPEW = nan(nxy_SAW,nmonths_SAW);
% SAW_grid_time = NaT(nmonths_SAW,1);
% 
% SAWz_raw = double(ncread(firstfile,'deptht'));
% 
% for im = 1:nmonths_SAW
% 
%     datecode = common_dates_SAW(im);
% 
%     tfile = temp_files_SAW(temp_dates_SAW == datecode);
%     sfile = salt_files_SAW(salt_dates_SAW == datecode);
% 
%     yr = str2double(extractBetween(datecode,1,4));
%     mo = str2double(extractBetween(datecode,5,6));
% 
%     SAW_grid_time(im) = datetime(yr,mo,15);
% 
%     T = double(ncread(fullfile(tfile.folder,tfile.name),'votemper'));
%     S = double(ncread(fullfile(sfile.folder,sfile.name),'vosaline'));
% 
%     T(abs(T) > 1e20) = NaN;
%     S(abs(S) > 1e20) = NaN;
% 
%     T2 = reshape(T,[],size(T,ndims(T)));
%     S2 = reshape(S,[],size(S,ndims(S)));
% 
%     if size(T2,2) ~= numel(SAWz_raw)
% 
%         depth_dim = find(size(T) == numel(SAWz_raw),1,'last');
% 
%         ord = 1:ndims(T);
%         ord([depth_dim end]) = ord([end depth_dim]);
% 
%         T = permute(T,ord);
%         S = permute(S,ord);
% 
%         T2 = reshape(T,[],numel(SAWz_raw));
%         S2 = reshape(S,[],numel(SAWz_raw));
% 
%     end
% 
%     lonvec = SAW_grid_lon(:);
%     latvec = SAW_grid_lat(:);
% 
%     for ig = 1:nxy_SAW
% 
%         tprof = T2(ig,:)';
%         sprof = S2(ig,:)';
% 
%         if sum(isfinite(tprof) & isfinite(sprof)) < 2
%             continue
%         end
% 
%         lon0 = lonvec(ig);
%         lat0 = latvec(ig);
% 
%         if ~isfinite(lon0) || ~isfinite(lat0)
%             continue
%         end
% 
%         p = gsw_p_from_z(-SAWz_raw,lat0);
% 
%         SA = gsw_SA_from_SP( ...
%             sprof, ...
%             p, ...
%             lon0 .* ones(size(SAWz_raw)), ...
%             lat0 .* ones(size(SAWz_raw)));
% 
%         CT = gsw_CT_from_pt(SA,tprof);
% 
%         st = gsw_sigma0(SA,CT);
% 
%         msk = ...
%             isfinite(st) & ...
%             isfinite(tprof) & ...
%             isfinite(SA);
% 
%         if sum(msk) < 2
%             continue
%         end
% 
%         st_tmp = st(msk);
%         t_tmp  = tprof(msk);
%         s_tmp  = SA(msk);
% 
%         [st_tmp,ii] = sort(st_tmp);
% 
%         t_tmp = t_tmp(ii);
%         s_tmp = s_tmp(ii);
% 
%         [st_tmp,ia] = unique(st_tmp);
% 
%         t_tmp = t_tmp(ia);
%         s_tmp = s_tmp(ia);
% 
%         if numel(st_tmp) < 2 || ...
%            core_st < min(st_tmp) || ...
%            core_st > max(st_tmp)
%             continue
%         end
% 
%         t_core = interp1(st_tmp,t_tmp,core_st);
%         s_core = interp1(st_tmp,s_tmp,core_st);
% 
%         if ~isfinite(t_core) || ~isfinite(s_core)
%             continue
%         end
% 
%         dT_s = t_core - PEW_t_core;
%         dS_s = s_core - PEW_s_core;
% 
%         alpha_PSUW = ...
%             (dT_s.*dT_core + dS_s.*dS_core) ./ ...
%             denom_core;
% 
%         SAW_grid_fPEW(ig,im) = 1 - alpha_PSUW;
% 
%     end
% 
%     fprintf('PSUW %s: %d valid grid points\n', ...
%         datestr(SAW_grid_time(im),'yyyy-mm'), ...
%         sum(isfinite(SAW_grid_fPEW(:,im))));
% 
% end
% 
% 
% % ============================================================
% % Quick diagnostics
% % ============================================================
% 
% fprintf('\nRaw grid-point fPEW ranges at sigma=%.2f:\n',core_st);
% 
% fprintf('PEW source:  %.3f to %.3f\n', ...
%     min(PEW_grid_fPEW,[],'all','omitnan'), ...
%     max(PEW_grid_fPEW,[],'all','omitnan'));
% 
% fprintf('PSUW source: %.3f to %.3f\n', ...
%     min(SAW_grid_fPEW,[],'all','omitnan'), ...
%     max(SAW_grid_fPEW,[],'all','omitnan'));

% save('/Users/samst/Library/CloudStorage/Dropbox/Hakai/gyre/supply/data/ORAS5_gridpoint_fPEW.mat', ...
%     'PEW_grid_fPEW','SAW_grid_fPEW', ...
%     'PEW_grid_time','SAW_grid_time', ...
%     'PEW_grid_lon','PEW_grid_lat', ...
%     'SAW_grid_lon','SAW_grid_lat', ...
%     'core_st','-v7.3');

load('/Users/samst/Library/CloudStorage/Dropbox/Hakai/gyre/supply/data/ORAS5_gridpoint_fPEW.mat');

%% Create monthly PEW and PSUW candidate profiles
%
% Each monthly source-water profile provides one candidate endmember.
%
% Stationary endmembers are selected using spiciness at
% sigma-theta =26.47 kg m^-3:
%
% PSUW = least spicy monthly profile
% PEW  = most spicy monthly profile
%
% The entire monthly profile corresponding to the extreme
% spiciness value at sigma-theta =26.47 is retained.

core_st =26.5;

% Interpolate monthly candidate profiles onto regular density grid

sgrid = [23:0.05:27.9]';

% ============================================================
% PSUW monthly candidates
% ============================================================

nSAW = length(SAW.time);

SAWt_monthly_grid = NaN(length(sgrid),nSAW);
SAWs_monthly_grid = NaN(length(sgrid),nSAW);
SAWz_monthly_grid = NaN(length(sgrid),nSAW);

SAW_spice_core = NaN(1,nSAW);

p = gsw_p_from_z(-SAWz,SAW_lat);

for i = 1:nSAW

    if all(isnan(SAWt_monthly(:,i))) || ...
            all(isnan(SAWs_monthly(:,i)))
        continue
    end

    % Monthly ORAS5 salinity is Practical Salinity
    SA = gsw_SA_from_SP( ...
        SAWs_monthly(:,i), ...
        p, ...
        SAW_lon .* ones(size(SAWz)), ...
        SAW_lat .* ones(size(SAWz)));

    % Monthly ORAS5 temperature is potential temperature
    CT = gsw_CT_from_pt( ...
        SA, ...
        SAWt_monthly(:,i));

    st = gsw_sigma0(SA,CT);

    msk = ...
        isfinite(st) & ...
        isfinite(SAWt_monthly(:,i)) & ...
        isfinite(SA);

    if sum(msk) < 2
        continue
    end

    [st_tmp,idx] = sort(st(msk));

    t_tmp = SAWt_monthly(msk,i);
    s_tmp = SA(msk);
    z_tmp = SAWz(msk);

    t_tmp = t_tmp(idx);
    s_tmp = s_tmp(idx);
    z_tmp = z_tmp(idx);

    [st_tmp,ia] = unique(st_tmp);

    t_tmp = t_tmp(ia);
    s_tmp = s_tmp(ia);
    z_tmp = z_tmp(ia);

    if length(st_tmp) < 2
        continue
    end

    % Interpolate complete monthly profile onto common density grid
    SAWt_monthly_grid(:,i) = ...
        interp1(st_tmp,t_tmp,sgrid);

    SAWs_monthly_grid(:,i) = ...
        interp1(st_tmp,s_tmp,sgrid);

    SAWz_monthly_grid(:,i) = ...
        interp1(st_tmp,z_tmp,sgrid);

    % Evaluate spiciness at PEW-core density
    t_core = interp1( ...
        st_tmp,t_tmp,core_st);

    s_core = interp1( ...
        st_tmp,s_tmp,core_st);

    if isfinite(t_core) && isfinite(s_core)

        CT_core = gsw_CT_from_pt( ...
            s_core,t_core);

        SAW_spice_core(i) = ...
            gsw_spiciness0(s_core,CT_core);

    end

end

% ============================================================
% PEW monthly candidates
% ============================================================

nPEW = length(PEW.time);

PEWt_monthly_grid = NaN(length(sgrid),nPEW);
PEWs_monthly_grid = NaN(length(sgrid),nPEW);
PEWz_monthly_grid = NaN(length(sgrid),nPEW);

PEW_spice_core = NaN(1,nPEW);

p = gsw_p_from_z(-PEWz,PEW_lat);

for i = 1:nPEW

    if all(isnan(PEWt_monthly(:,i))) || ...
            all(isnan(PEWs_monthly(:,i)))
        continue
    end

    % Monthly ORAS5 salinity is Practical Salinity
    SA = gsw_SA_from_SP( ...
        PEWs_monthly(:,i), ...
        p, ...
        PEW_lon .* ones(size(PEWz)), ...
        PEW_lat .* ones(size(PEWz)));

    % Monthly ORAS5 temperature is potential temperature
    CT = gsw_CT_from_pt( ...
        SA, ...
        PEWt_monthly(:,i));

    st = gsw_sigma0(SA,CT);

    msk = ...
        isfinite(st) & ...
        isfinite(PEWt_monthly(:,i)) & ...
        isfinite(SA);

    if sum(msk) < 2
        continue
    end

    [st_tmp,idx] = sort(st(msk));

    t_tmp = PEWt_monthly(msk,i);
    s_tmp = SA(msk);
    z_tmp = PEWz(msk);

    t_tmp = t_tmp(idx);
    s_tmp = s_tmp(idx);
    z_tmp = z_tmp(idx);

    [st_tmp,ia] = unique(st_tmp);

    t_tmp = t_tmp(ia);
    s_tmp = s_tmp(ia);
    z_tmp = z_tmp(ia);

    if length(st_tmp) < 2
        continue
    end

    % Interpolate complete monthly profile onto common density grid
    PEWt_monthly_grid(:,i) = ...
        interp1(st_tmp,t_tmp,sgrid);

    PEWs_monthly_grid(:,i) = ...
        interp1(st_tmp,s_tmp,sgrid);

    PEWz_monthly_grid(:,i) = ...
        interp1(st_tmp,z_tmp,sgrid);

    % Evaluate spiciness at PEW-core density
    t_core = interp1( ...
        st_tmp,t_tmp,core_st);

    s_core = interp1( ...
        st_tmp,s_tmp,core_st);

    if isfinite(t_core) && isfinite(s_core)

        CT_core = gsw_CT_from_pt( ...
            s_core,t_core);

        PEW_spice_core(i) = ...
            gsw_spiciness0(s_core,CT_core);

    end

end

% Convert monthly source-water salinity to Absolute Salinity

% ============================================================
% PSUW
% ============================================================

SAWs_monthly_SA = nan(size(SAWs_monthly));

p_SAW = gsw_p_from_z(-SAWz,SAW_lat);

for i = 1:length(SAW.time)

    SAWs_monthly_SA(:,i) = gsw_SA_from_SP( ...
        SAWs_monthly(:,i), ...
        p_SAW, ...
        SAW_lon .* ones(size(SAWz)), ...
        SAW_lat .* ones(size(SAWz)));

end


% ============================================================
% PEW
% ============================================================

PEWs_monthly_SA = nan(size(PEWs_monthly));

p_PEW = gsw_p_from_z(-PEWz,PEW_lat);

for i = 1:length(PEW.time)

    PEWs_monthly_SA(:,i) = gsw_SA_from_SP( ...
        PEWs_monthly(:,i), ...
        p_PEW, ...
        PEW_lon .* ones(size(PEWz)), ...
        PEW_lat .* ones(size(PEWz)));

end

%% Select stationary monthly endmember profiles
%
% PSUW = least-spicy monthly profile at sigma-theta =26.47
% PEW  = most-spicy monthly profile at sigma-theta =26.47
%
% The spiciness criterion selects an entire monthly profile.

[~,iSAW] = min(SAW_spice_core,[],'omitnan');
[~,iPEW] = max(PEW_spice_core,[],'omitnan');

if isempty(iSAW) || ~isfinite(SAW_spice_core(iSAW))
    error('Could not identify a valid PSUW stationary endmember.')
end

if isempty(iPEW) || ~isfinite(PEW_spice_core(iPEW))
    error('Could not identify a valid PEW stationary endmember.')
end

% Fixed PSUW profile
SAWt = SAWt_monthly_grid(:,iSAW);
SAWs = SAWs_monthly_grid(:,iSAW);
SAWz_grid = SAWz_monthly_grid(:,iSAW);

% Fixed PEW profile
PEWt = PEWt_monthly_grid(:,iPEW);
PEWs = PEWs_monthly_grid(:,iPEW);
PEWz_grid = PEWz_monthly_grid(:,iPEW);

fprintf('\nSelected stationary endmembers:\n');

fprintf(['PSUW: %s monthly profile, ' ...
    'spiciness at sigma=%.2f = %.4f kg m^-3\n'], ...
    datestr(SAW.time(iSAW),'mmm yyyy'), ...
    core_st, ...
    SAW_spice_core(iSAW));

fprintf(['PEW: %s monthly profile, ' ...
    'spiciness at sigma=%.2f = %.4f kg m^-3\n\n'], ...
    datestr(PEW.time(iPEW),'mmm yyyy'), ...
    core_st, ...
    PEW_spice_core(iPEW));

%% Create stationary PEW-PSUW mixtures
%
% Dimensions:
% density x mixture fraction

mixGridT = NaN(length(sgrid),9);
mixGridS = NaN(length(sgrid),9);

for i = 1:9

    % i/10 = PSUW fraction
    % 1-i/10 = PEW fraction

    mixGridT(:,i) = ...
        SAWt.*(i/10) + ...
        PEWt.*(1-(i/10));

    mixGridS(:,i) = ...
        SAWs.*(i/10) + ...
        PEWs.*(1-(i/10));

end

%% Find fractions using stationary linear mixing model along isopycnals

load CSOS_ctd.mat
load P4Long_ctd.mat

% ============================================================
% CS-OS
% ============================================================

[nz_cs,nprof_cs] = size(CSOS.temp);

alpha_SAW_CSOS = nan(nz_cs,nprof_cs);
alpha_PEW_CSOS = nan(nz_cs,nprof_cs);

for j = 1:nprof_cs

    for i = 1:nz_cs

        st_val = CSOS.s_t(i,j);

        if ~isfinite(st_val)
            continue
        end

        % Fixed endmembers evaluated at density of observation
        SAW_t = interp1( ...
            sgrid,SAWt,st_val);

        SAW_s = interp1( ...
            sgrid,SAWs,st_val);

        PEW_t = interp1( ...
            sgrid,PEWt,st_val);

        PEW_s = interp1( ...
            sgrid,PEWs,st_val);

        % Observed temperature and Absolute Salinity
        T_val = CSOS.temp(i,j);
        S_val = CSOS.SA(i,j);

        if ~isfinite(T_val) || ...
           ~isfinite(S_val) || ...
           ~isfinite(SAW_t) || ...
           ~isfinite(SAW_s) || ...
           ~isfinite(PEW_t) || ...
           ~isfinite(PEW_s)

            continue
        end

        % PEW -> PSUW mixing vector
        dT = SAW_t - PEW_t;
        dS = SAW_s - PEW_s;

        % PEW -> observation vector
        dT_s = T_val - PEW_t;
        dS_s = S_val - PEW_s;

        denom = dT^2 + dS^2;

        if denom == 0
            continue
        end

        alpha = ...
            (dT_s*dT + dS_s*dS) / denom;

        % Restrict fraction to 0-1
        alpha = max(min(alpha,1),0);

        alpha_SAW_CSOS(i,j) = alpha;
        alpha_PEW_CSOS(i,j) = 1 - alpha;

    end

end

% ============================================================
% P4
% ============================================================

[nz_p4,nprof_p4] = size(P4.temp);

alpha_SAW_P4 = nan(nz_p4,nprof_p4);
alpha_PEW_P4 = nan(nz_p4,nprof_p4);

for j = 1:nprof_p4

    for i = 1:nz_p4

        st_val = P4.s_t(i,j);

        if ~isfinite(st_val)
            continue
        end

        % Fixed endmembers evaluated at density of observation
        SAW_t = interp1( ...
            sgrid,SAWt,st_val);

        SAW_s = interp1( ...
            sgrid,SAWs,st_val);

        PEW_t = interp1( ...
            sgrid,PEWt,st_val);

        PEW_s = interp1( ...
            sgrid,PEWs,st_val);

        % Observed temperature and Absolute Salinity
        T_val = P4.temp(i,j);
        S_val = P4.SA(i,j);

        if ~isfinite(T_val) || ...
           ~isfinite(S_val) || ...
           ~isfinite(SAW_t) || ...
           ~isfinite(SAW_s) || ...
           ~isfinite(PEW_t) || ...
           ~isfinite(PEW_s)

            continue
        end

        % PEW -> PSUW mixing vector
        dT = SAW_t - PEW_t;
        dS = SAW_s - PEW_s;

        % PEW -> observation vector
        dT_s = T_val - PEW_t;
        dS_s = S_val - PEW_s;

        denom = dT^2 + dS^2;

        if denom == 0
            continue
        end

        alpha = ...
            (dT_s*dT + dS_s*dS) / denom;

        % Restrict fraction to 0-1
        alpha = max(min(alpha,1),0);

        alpha_SAW_P4(i,j) = alpha;
        alpha_PEW_P4(i,j) = 1 - alpha;

    end

end


% Monthly PSUW source-water fPEW
%
% Project each monthly PSUW source-water profile onto the SAME
% stationary PEW--PSUW mixing line used for CS-OS and P4.
%
% Raw fractions are retained because individual monthly source-water
% states may lie outside the stationary endmember range.

alpha_PEW_SAW_source_raw = ...
    nan(length(sgrid),length(SAW.time));

alpha_PEW_SAW_source = ...
    nan(length(sgrid),length(SAW.time));

for j = 1:length(SAW.time)

    % Monthly PSUW profile is already available on sgrid
    T_month = SAWt_monthly_grid(:,j);
    S_month = SAWs_monthly_grid(:,j);

    for i = 1:length(sgrid)

        if ~isfinite(T_month(i)) || ...
           ~isfinite(S_month(i)) || ...
           ~isfinite(SAWt(i)) || ...
           ~isfinite(SAWs(i)) || ...
           ~isfinite(PEWt(i)) || ...
           ~isfinite(PEWs(i))

            continue
        end

        % PEW -> PSUW mixing vector
        dT = SAWt(i) - PEWt(i);
        dS = SAWs(i) - PEWs(i);

        % PEW -> monthly PSUW-source observation vector
        dT_s = T_month(i) - PEWt(i);
        dS_s = S_month(i) - PEWs(i);

        denom = dT^2 + dS^2;

        if denom == 0
            continue
        end

        % alpha = PSUW fraction
        alpha = ...
            (dT_s*dT + dS_s*dS) / denom;

        % PEW fraction before clipping
        alpha_PEW_SAW_source_raw(i,j) = 1 - alpha;

        % PEW fraction constrained to physical 0-1 range
        alpha = max(min(alpha,1),0);

        alpha_PEW_SAW_source(i,j) = 1 - alpha;

    end

end

% Monthly PEW source-water fPEW

alpha_PEW_PEW_source_raw = ...
    nan(length(sgrid),length(PEW.time));

alpha_PEW_PEW_source = ...
    nan(length(sgrid),length(PEW.time));

for j = 1:length(PEW.time)

    % Monthly PEW profile is already available on sgrid
    T_month = PEWt_monthly_grid(:,j);
    S_month = PEWs_monthly_grid(:,j);

    for i = 1:length(sgrid)

        if ~isfinite(T_month(i)) || ...
           ~isfinite(S_month(i)) || ...
           ~isfinite(SAWt(i)) || ...
           ~isfinite(SAWs(i)) || ...
           ~isfinite(PEWt(i)) || ...
           ~isfinite(PEWs(i))

            continue
        end

        % PEW -> PSUW mixing vector
        dT = SAWt(i) - PEWt(i);
        dS = SAWs(i) - PEWs(i);

        % PEW -> monthly PEW-source observation vector
        dT_s = T_month(i) - PEWt(i);
        dS_s = S_month(i) - PEWs(i);

        denom = dT^2 + dS^2;

        if denom == 0
            continue
        end

        % alpha = PSUW fraction
        alpha = ...
            (dT_s*dT + dS_s*dS) / denom;

        % PEW fraction before clipping
        alpha_PEW_PEW_source_raw(i,j) = 1 - alpha;

        % PEW fraction constrained to physical 0-1 range
        alpha = max(min(alpha,1),0);

        alpha_PEW_PEW_source(i,j) = 1 - alpha;

    end

end

% Extract source-water fPEW at PEW-core density

[~,i_core] = min(abs(sgrid-core_st));

fPEW_SAW_source_core = ...
    alpha_PEW_SAW_source(i_core,:);

fPEW_PEW_source_core = ...
    alpha_PEW_PEW_source(i_core,:);

fPEW_SAW_source_core_raw = ...
    alpha_PEW_SAW_source_raw(i_core,:);

fPEW_PEW_source_core_raw = ...
    alpha_PEW_PEW_source_raw(i_core,:);

fprintf('\nMonthly source-water fPEW at sigma=%.2f:\n', ...
    sgrid(i_core));

fprintf(['PSUW source: mean = %.3f, ' ...
    'range = %.3f to %.3f\n'], ...
    mean(fPEW_SAW_source_core,'omitnan'), ...
    min(fPEW_SAW_source_core,[],'omitnan'), ...
    max(fPEW_SAW_source_core,[],'omitnan'));

fprintf(['PEW source:  mean = %.3f, ' ...
    'range = %.3f to %.3f\n'], ...
    mean(fPEW_PEW_source_core,'omitnan'), ...
    min(fPEW_PEW_source_core,[],'omitnan'), ...
    max(fPEW_PEW_source_core,[],'omitnan'));



%% PEW core depth and sigma-theta variability diagnostics

csst = mean(CSOS.s_t,2,'omitnan');
p4st = mean(P4.s_t,2,'omitnan');

% Characteristic PEW-core density range
PEW_st_min = 26.25;
PEW_st_max = 26.7;

% Mean-profile PEW maxima
% CS-OS
mean_alpha_CS = mean(alpha_PEW_CSOS(1:500,:),2,'omitnan');

msk = ...
    csst(1:500) >= PEW_st_min & ...
    csst(1:500) <= PEW_st_max & ...
    isfinite(mean_alpha_CS);

inds = find(msk);

[CS_PEW_max,jj] = max(mean_alpha_CS(msk));
CS_z_max = inds(jj);
CS_st_max = csst(CS_z_max);


% P4
mean_alpha_P4 = mean(alpha_PEW_P4(1:500,:),2,'omitnan');

msk = ...
    p4st(1:500) >= PEW_st_min & ...
    p4st(1:500) <= PEW_st_max & ...
    isfinite(mean_alpha_P4);

inds = find(msk);

[P4_PEW_max,jj] = max(mean_alpha_P4(msk));
P4_z_max = inds(jj);
P4_st_max = p4st(P4_z_max);

% Profile-by-profile PEW max depth/index and sigma-theta variability
% Maximum fPEW is only searched for within the characteristic
% PEW-core density range.

% -----------------
% CS-OS
% -----------------

CS_z_all  = NaN(1,size(alpha_PEW_CSOS,2));
CS_st_all = NaN(1,size(alpha_PEW_CSOS,2));

for i = 1:size(alpha_PEW_CSOS,2)

    idx = find(isfinite(alpha_PEW_CSOS(1:500,i)),1,'last');

    if ~isempty(idx) && idx > 175 && CSOS.profile_type(i) == 1

        st = CSOS.s_t(1:500,i);

        msk = ...
            isfinite(alpha_PEW_CSOS(1:500,i)) & ...
            isfinite(st) & ...
            st >= PEW_st_min & ...
            st <= PEW_st_max;

        if any(msk)

            inds = find(msk);

            [~,jj] = max(alpha_PEW_CSOS(msk,i));

            CS_z_all(i) = inds(jj);
            CS_st_all(i) = CSOS.s_t(CS_z_all(i),i);

        end
    end
end

% % Remove deep maxima unrelated to CUC
% CS_st_all(CS_z_all > 300) = NaN;
% CS_z_all(CS_z_all > 300) = NaN;

nout = 3;

CS_st_med = median(CS_st_all,'omitnan');
CS_st_sig = 1.4826 .* median(abs(CS_st_all - CS_st_med),'omitnan');

CS_out = abs(CS_st_all - CS_st_med) > nout .* CS_st_sig;

CS_st_all(CS_out) = NaN;
CS_z_all(CS_out)  = NaN;


CS_z_max_mean = mean(CS_z_all,'omitnan');
CS_z_max_std  = std(CS_z_all,'omitnan');

CS_st_mean   = mean(CS_st_all,'omitnan');
CS_st_median = median(CS_st_all,'omitnan');
CS_st_std    = std(CS_st_all,'omitnan');


% -----------------
% P4
% -----------------

P4_z_all  = NaN(1,size(alpha_PEW_P4,2));
P4_st_all = NaN(1,size(alpha_PEW_P4,2));

for i = 1:size(alpha_PEW_P4,2)

    idx = find(isfinite(alpha_PEW_P4(1:500,i)),1,'last');

    if ~isempty(idx) && idx > 175 && P4.profile_type(i) == 1

        st = P4.s_t(1:500,i);

        msk = ...
            isfinite(alpha_PEW_P4(1:500,i)) & ...
            isfinite(st) & ...
            st >= PEW_st_min & ...
            st <= PEW_st_max;

        if any(msk)

            inds = find(msk);

            [~,jj] = max(alpha_PEW_P4(msk,i));

            P4_z_all(i) = inds(jj);
            P4_st_all(i) = P4.s_t(P4_z_all(i),i);

        end
    end
end


% P4_st_all(P4_z_all > 300) = NaN;
% P4_z_all(P4_z_all > 300) = NaN;

P4_st_med = median(P4_st_all,'omitnan');
P4_st_sig = 1.4826 .* median(abs(P4_st_all - P4_st_med),'omitnan');

P4_out = abs(P4_st_all - P4_st_med) > nout .* P4_st_sig;

P4_st_all(P4_out) = NaN;
P4_z_all(P4_out)  = NaN;

P4_z_max_mean = mean(P4_z_all,'omitnan');
P4_z_max_std  = std(P4_z_all,'omitnan');

P4_st_mean   = mean(P4_st_all,'omitnan');
P4_st_median = median(P4_st_all,'omitnan');
P4_st_std    = std(P4_st_all,'omitnan');


% fPEW at profile-specific core

CS_PEW_all = NaN(1,size(alpha_PEW_CSOS,2));

for i = 1:size(alpha_PEW_CSOS,2)

    if isfinite(CS_z_all(i))

        CS_PEW_all(i) = ...
            alpha_PEW_CSOS(CS_z_all(i),i);

    end
end


P4_PEW_all = NaN(1,size(alpha_PEW_P4,2));

for i = 1:size(alpha_PEW_P4,2)

    if isfinite(P4_z_all(i))

        P4_PEW_all(i) = ...
            alpha_PEW_P4(P4_z_all(i),i);

    end
end


%  PSUW source-water PEW-core diagnostics

% Mean-profile PEW maximum
mean_alpha_SAW = mean(alpha_PEW_SAW_source,2,'omitnan');

msk = ...
    sgrid >= PEW_st_min & ...
    sgrid <= PEW_st_max & ...
    isfinite(mean_alpha_SAW);

inds = find(msk);

[SAW_PEW_max,jj] = max(mean_alpha_SAW(msk));

SAW_i_max  = inds(jj);
SAW_st_max = sgrid(SAW_i_max);

% Mean depth of this density surface
SAW_z_max = mean( ...
    SAWz_monthly_grid(SAW_i_max,:), ...
    'omitnan');


% -----------------
% Profile-by-profile PEW maximum
% -----------------

SAW_i_all   = NaN(1,size(alpha_PEW_SAW_source,2));
SAW_z_all   = NaN(1,size(alpha_PEW_SAW_source,2));
SAW_st_all  = NaN(1,size(alpha_PEW_SAW_source,2));
SAW_PEW_all = NaN(1,size(alpha_PEW_SAW_source,2));

for i = 1:size(alpha_PEW_SAW_source,2)

    msk = ...
        isfinite(alpha_PEW_SAW_source(:,i)) & ...
        sgrid >= PEW_st_min & ...
        sgrid <= PEW_st_max;

    if any(msk)

        inds = find(msk);

        [~,jj] = max(alpha_PEW_SAW_source(msk,i));

        SAW_i_all(i)   = inds(jj);
        SAW_st_all(i)  = sgrid(SAW_i_all(i));
        SAW_z_all(i)   = SAWz_monthly_grid(SAW_i_all(i),i);
        SAW_PEW_all(i) = alpha_PEW_SAW_source(SAW_i_all(i),i);

    end
end


% -----------------
% Remove density outliers
% -----------------

SAW_st_med = median(SAW_st_all,'omitnan');

SAW_st_sig = ...
    1.4826 .* ...
    median(abs(SAW_st_all-SAW_st_med),'omitnan');

SAW_out = ...
    abs(SAW_st_all-SAW_st_med) > ...
    nout .* SAW_st_sig;

SAW_i_all(SAW_out)   = NaN;
SAW_z_all(SAW_out)   = NaN;
SAW_st_all(SAW_out)  = NaN;
SAW_PEW_all(SAW_out) = NaN;


% -----------------
% CS-OS summary statistics
% -----------------

% Keep existing mean/SD diagnostics
CS_z_max_mean = mean(CS_z_all,'omitnan');
CS_z_max_std  = std(CS_z_all,'omitnan');

CS_st_mean   = mean(CS_st_all,'omitnan');
CS_st_median = median(CS_st_all,'omitnan');
CS_st_std    = std(CS_st_all,'omitnan');

% Median and IQR
tmp = CS_z_all(isfinite(CS_z_all));
CS_z_median = median(tmp);
CS_z_q25 = prctile(tmp,25);
CS_z_q75 = prctile(tmp,75);
CS_z_iqr = CS_z_q75 - CS_z_q25;

tmp = CS_st_all(isfinite(CS_st_all));
CS_st_median = median(tmp);
CS_st_q25 = prctile(tmp,25);
CS_st_q75 = prctile(tmp,75);
CS_st_iqr = CS_st_q75 - CS_st_q25;


% -----------------
% P4 summary statistics
% -----------------

% Keep existing mean/SD diagnostics
P4_z_max_mean = mean(P4_z_all,'omitnan');
P4_z_max_std  = std(P4_z_all,'omitnan');

P4_st_mean   = mean(P4_st_all,'omitnan');
P4_st_median = median(P4_st_all,'omitnan');
P4_st_std    = std(P4_st_all,'omitnan');

% Median and IQR
tmp = P4_z_all(isfinite(P4_z_all));
P4_z_median = median(tmp);
P4_z_q25 = prctile(tmp,25);
P4_z_q75 = prctile(tmp,75);
P4_z_iqr = P4_z_q75 - P4_z_q25;

tmp = P4_st_all(isfinite(P4_st_all));
P4_st_median = median(tmp);
P4_st_q25 = prctile(tmp,25);
P4_st_q75 = prctile(tmp,75);
P4_st_iqr = P4_st_q75 - P4_st_q25;


% -----------------
% PSUW source summary statistics
% -----------------

% Keep mean/SD diagnostics
SAW_z_mean = mean(SAW_z_all,'omitnan');
SAW_z_std  = std(SAW_z_all,'omitnan');

SAW_st_mean   = mean(SAW_st_all,'omitnan');
SAW_st_median = median(SAW_st_all,'omitnan');
SAW_st_std    = std(SAW_st_all,'omitnan');

SAW_PEW_mean = mean(SAW_PEW_all,'omitnan');
SAW_PEW_std  = std(SAW_PEW_all,'omitnan');

% Median and IQR
tmp = SAW_z_all(isfinite(SAW_z_all));
SAW_z_median = median(tmp);
SAW_z_q25 = prctile(tmp,25);
SAW_z_q75 = prctile(tmp,75);
SAW_z_iqr = SAW_z_q75 - SAW_z_q25;

tmp = SAW_st_all(isfinite(SAW_st_all));
SAW_st_median = median(tmp);
SAW_st_q25 = prctile(tmp,25);
SAW_st_q75 = prctile(tmp,75);
SAW_st_iqr = SAW_st_q75 - SAW_st_q25;

tmp = SAW_PEW_all(isfinite(SAW_PEW_all));
SAW_PEW_median = median(tmp);
SAW_PEW_q25 = prctile(tmp,25);
SAW_PEW_q75 = prctile(tmp,75);
SAW_PEW_iqr = SAW_PEW_q75 - SAW_PEW_q25;


% -----------------
% fPEW at profile-specific core
% -----------------

CS_PEW_all = NaN(1,size(alpha_PEW_CSOS,2));
for i = 1:size(alpha_PEW_CSOS,2)
    if isfinite(CS_z_all(i))
        CS_PEW_all(i) = alpha_PEW_CSOS(CS_z_all(i),i);
    end
end

P4_PEW_all = NaN(1,size(alpha_PEW_P4,2));
for i = 1:size(alpha_PEW_P4,2)
    if isfinite(P4_z_all(i))
        P4_PEW_all(i) = alpha_PEW_P4(P4_z_all(i),i);
    end
end

% Keep existing mean/SD
CS_PEW_mean = mean(CS_PEW_all,'omitnan');
CS_PEW_std  = std(CS_PEW_all,'omitnan');

P4_PEW_mean = mean(P4_PEW_all,'omitnan');
P4_PEW_std  = std(P4_PEW_all,'omitnan');

% Median and IQR
tmp = CS_PEW_all(isfinite(CS_PEW_all));
CS_PEW_median = median(tmp);
CS_PEW_q25 = prctile(tmp,25);
CS_PEW_q75 = prctile(tmp,75);
CS_PEW_iqr = CS_PEW_q75 - CS_PEW_q25;

tmp = P4_PEW_all(isfinite(P4_PEW_all));
P4_PEW_median = median(tmp);
P4_PEW_q25 = prctile(tmp,25);
P4_PEW_q75 = prctile(tmp,75);
P4_PEW_iqr = P4_PEW_q75 - P4_PEW_q25;


% -----------------
% Print summaries
% -----------------

fprintf(['CS-OS PEW core depth: %.1f m [IQR %.1f--%.1f]; ' ...
    'sigma-theta: %.3f [IQR %.3f--%.3f]; ' ...
    'fPEW: %.3f [IQR %.3f--%.3f]\n'], ...
    CS_z_median, CS_z_q25, CS_z_q75, ...
    CS_st_median, CS_st_q25, CS_st_q75, ...
    CS_PEW_median, CS_PEW_q25, CS_PEW_q75);

fprintf(['P4 PEW core depth: %.1f m [IQR %.1f--%.1f]; ' ...
    'sigma-theta: %.3f [IQR %.3f--%.3f]; ' ...
    'fPEW: %.3f [IQR %.3f--%.3f]\n'], ...
    P4_z_median, P4_z_q25, P4_z_q75, ...
    P4_st_median, P4_st_q25, P4_st_q75, ...
    P4_PEW_median, P4_PEW_q25, P4_PEW_q75);

fprintf(['PSUW source PEW core depth: %.1f m [IQR %.1f--%.1f]; ' ...
    'sigma-theta: %.3f [IQR %.3f--%.3f]; ' ...
    'fPEW: %.3f [IQR %.3f--%.3f]\n'], ...
    SAW_z_median, SAW_z_q25, SAW_z_q75, ...
    SAW_st_median, SAW_st_q25, SAW_st_q75, ...
    SAW_PEW_median, SAW_PEW_q25, SAW_PEW_q75);

%% Calculate time variability along isopycnal

% Target isopycnal
ptargetCS = CS_st_median;
ptargetP4 = P4_st_median;

% CSOS: interpolate alpha 
alpha_SAW_CSOS_26p5 = nan(1, nprof_cs);
alpha_PEW_CSOS_26p5 = nan(1, nprof_cs);

for j = 1:nprof_cs
  st_prof  = CSOS.s_t(:,j);
  saw_prof = alpha_SAW_CSOS(:,j);
  pew_prof = alpha_PEW_CSOS(:,j);

  valid = ~isnan(st_prof) & ~isnan(saw_prof);
  if sum(valid) >= 2
    % linear interpolation at ptarget, out-of-range → NaN
    alpha_SAW_CSOS_26p5(j) = interp1(st_prof(valid), saw_prof(valid), ptargetCS, 'linear', NaN);
    alpha_PEW_CSOS_26p5(j) = interp1(st_prof(valid), pew_prof(valid), ptargetCS, 'linear', NaN);
  end
end

% P4: interpolate alpha
alpha_SAW_P4_26p5 = nan(1, nprof_p4);
alpha_PEW_P4_26p5 = nan(1, nprof_p4);

for j = 1:nprof_p4
  st_prof  = P4.s_t(:,j);
  saw_prof = alpha_SAW_P4(:,j);
  pew_prof = alpha_PEW_P4(:,j);

  valid = ~isnan(st_prof) & ~isnan(saw_prof);
  if sum(valid) >= 2
    alpha_SAW_P4_26p5(j) = interp1(st_prof(valid), saw_prof(valid), ptargetP4, 'linear', NaN);
    alpha_PEW_P4_26p5(j) = interp1(st_prof(valid), pew_prof(valid), ptargetP4, 'linear', NaN);
  end
end

%% TS figure
%
% Annual-mean source-water envelopes are shown as shaded regions.
% Dashed blue/red profiles show the selected stationary endmembers:
%
%   PSUW = least-spicy annual mean at sigma-theta =26.47
%   PEW  = most-spicy annual mean at sigma-theta =26.47
%
% Grey dashed lines show the stationary PEW--PSUW mixing lines.

Sgrid = 32:0.01:35.5;
Tgrid = 2:0.1:20;

[SS,TT] = meshgrid(Sgrid,Tgrid);

SA_grid = gsw_SA_from_SP(SS,0,-142.5,50);
CT_grid = gsw_CT_from_t(SA_grid,TT,0);

sigma0 = gsw_sigma0(SA_grid,CT_grid);
spice  = gsw_spiciness0(SA_grid,CT_grid);

f1 = figure('units','centimeters','outerposition',...
    [0.01 0.01 14 18],'color','w');

set(gcf,'Renderer','opengl');

% ax1 = axes('Position',[0.075 0.075 0.6 0.9]);
ax1=axes;
hold on;

[h1,c1] = contour(SS,TT,spice,'k:');
[h2,c2] = contour(SS,TT,sigma0,'k-');

% ============================================================
% Source waters
% ============================================================

% Monthly PSUW source-water points
scatter( ...
    SAWs_monthly_SA, ...
    SAWt_monthly, ...
    1, ...
    'filled', ...
    'markerfacecolor',rgb_x('greyish purple'), ...
    'markerfacealpha',0.05, ...
    'HandleVisibility','off');

% Monthly PEW source-water points
scatter( ...
    PEWs_monthly_SA, ...
    PEWt_monthly, ...
    1, ...
    'filled', ...
    'markerfacecolor',rgb_x('burnt orange'), ...
    'markerfacealpha',0.05, ...
    'HandleVisibility','off');

% -----------------
% Selected stationary PSUW endmember
% -----------------

% White underlay
plot( ...
    SAWs, ...
    SAWt, ...
    'w', ...
    'linewidth',3, ...
    'HandleVisibility','off');

ll1 = plot( ...
    SAWs, ...
    SAWt,'-', ...
    'linewidth',2, ...
    'color',rgb_x('greyish purple'));

% -----------------
% Selected stationary PEW endmember
% -----------------

% White underlay
plot( ...
    PEWs, ...
    PEWt, ...
    'w', ...
    'linewidth',3, ...
    'HandleVisibility','off');

ll2 = plot( ...
    PEWs, ...
    PEWt, '-', ...
    'linewidth',2, ...
    'color',rgb_x('burnt orange'));

% ============================================================
% Observations
% ============================================================

scatter( ...
    CSOS.SA, ...
    CSOS.temp, ...
    1, ...
    'filled', ...
    'markerfacecolor',rgb_x('rose'), ...
    'markerfacealpha',0.05, ...
    'HandleVisibility','off');

scatter( ...
    P4.SA, ...
    P4.temp, ...
    1, ...
    'filled', ...
    'markerfacecolor',rgb_x('steel blue'), ...
    'markerfacealpha',0.05, ...
    'HandleVisibility','off');

% -----------------
% Mean P4 profile
% -----------------

plot( ...
    mean(P4.SA,2,'omitnan'), ...
    mean(P4.temp,2,'omitnan'), ...
    'w', ...
    'linewidth',3, ...
    'HandleVisibility','off');

l3 = plot( ...
    mean(P4.SA,2,'omitnan'), ...
    mean(P4.temp,2,'omitnan'), ...
    'linewidth',2, ...
    'color',rgb_x('steel blue'));

% -----------------
% Mean CS-OS profile
% -----------------

plot( ...
    mean(CSOS.SA,2,'omitnan'), ...
    mean(CSOS.temp,2,'omitnan'), ...
    'w', ...
    'linewidth',3, ...
    'HandleVisibility','off');

l4 = plot( ...
    mean(CSOS.SA,2,'omitnan'), ...
    mean(CSOS.temp,2,'omitnan'), ...
    'linewidth',2, ...
    'color',rgb_x('rose'));

% ============================================================
% Stationary mixing lines
% ============================================================

for i = 1:size(mixGridT,2)

    ml=plot( ...
        mixGridS(:,i), ...
        mixGridT(:,i), ...
        '--', ...
        'Color',[0.75 0.75 0.75], ...
        'LineWidth',1.5, ...
        'HandleVisibility','off');

end

% ============================================================
% Main-axis formatting
% ============================================================

xlabel('Absolute Salinity (g kg^{-1})');
ylabel('Potential Temperature, \theta (°C)');

box on;

set(gca,'FontSize',10);

xlim([32.75 35.25]);
ylim([2 14]);

clabel( ...
    h1, ...
    'FontSize',8, ...
    'Color',[0.4 0.4 0.4], ...
    'fontname','Latin Modern Roman');

clabel( ...
    h2, ...
    'FontSize',8, ...
    'Color',[0 0 0], ...
    'fontname','Latin Modern Roman');

% Legend labels identifying selected annual endmembers
legend( ...
    [c1,c2,ll1,ll2,l3,l4,ml], ...
    {'Spiciness', ...
     '\sigma_\theta', ...
     'PSUW endmember (Feb. 2013)', ...
     'PEW endmember (Jun. 1997)','P4', ...
     'CS-OS', 'Endmember mixtures (10% intervals)' ...
     },...
    'Location','northwest', ...
    'FontSize',10, ...
    'NumColumns',1);

set(findall(gcf,'-property','fontname'),...
    'fontname','Latin Modern Roman');
set(findall(gcf,'-property','fontsize'),'fontsize',8);
%%
% export_fig /Users/samst/Library/CloudStorage/Dropbox/Hakai/gyre/supply/figs/TS.jpg -m5

%% Isopycnal ox/spice

% -----------------
% Unit cleanup
% -----------------
CSOS.oxML(CSOS.oxML<10)   = CSOS.oxML(CSOS.oxML<10).*43;
CSOS.oxML2(CSOS.oxML2<10) = CSOS.oxML2(CSOS.oxML2<10).*43;
CSOS.ox(CSOS.ox<10)       = CSOS.ox(CSOS.ox<10).*43;
CSOS.ox2(CSOS.ox2<10)     = CSOS.ox2(CSOS.ox2<10).*43;
P4.ox(P4.ox>700)   = NaN;

% -----------------
% CSOS isopycnals
% -----------------
CSOS_sig_thresh = [26.2 CS_st_median 26.2:0.05:24.05];

niso = length(CSOS_sig_thresh);
nprof = size(CSOS.ox,2);

CSOS.isoDO    = NaN(niso,nprof);
CSOS.isoAOU   = NaN(niso,nprof);
CSOS.isoT     = NaN(niso,nprof);
CSOS.isoNO3   = NaN(niso,nprof);
CSOS.isoPO4   = NaN(niso,nprof);
CSOS.isoSR    = NaN(niso,nprof);
CSOS.isoSpice = NaN(niso,nprof);
CSOS.isoN2    = NaN(niso,nprof);
CSOS.isoZ     = NaN(niso,nprof);

for i = 1:niso
    for ii = 1:nprof

        st = CSOS.s_t(:,ii);

        ox  = mean([CSOS.oxML(:,ii), CSOS.oxML2(:,ii), ...
                    CSOS.ox(:,ii),   CSOS.ox2(:,ii)], 2, 'omitnan');
        no3 = CSOS.no3(:,ii);
        po4 = CSOS.po4(:,ii);

        % Oxygen-linked fields
        if sum(~isnan(ox)) > 1 && sum(~isnan(st)) > 1

            PT      = gsw_pt0_from_t(CSOS.SA(:,ii),CSOS.temp(:,ii),CSOS.pr(:,ii));
            tmp_aou = aou(CSOS.sal(:,ii),PT,ox);

            msk = ~isnan(st) & ~isnan(ox);

            CSOS.isoDO(i,ii)    = interp1(st(msk), ox(msk),              CSOS_sig_thresh(i));
            CSOS.isoAOU(i,ii)   = interp1(st(msk), tmp_aou(msk),         CSOS_sig_thresh(i));
            CSOS.isoT(i,ii)     = interp1(st(msk), CSOS.temp(msk,ii),    CSOS_sig_thresh(i));
            CSOS.isoSR(i,ii)    = interp1(st(msk), CSOS.SR(msk,ii),      CSOS_sig_thresh(i));
            CSOS.isoSpice(i,ii) = interp1(st(msk), CSOS.spice(msk,ii),   CSOS_sig_thresh(i));
            CSOS.isoZ(i,ii)     = interp1(st(msk), CSOS.pr(msk,ii),      CSOS_sig_thresh(i));

            if isfield(CSOS,'N2')
                CSOS.isoN2(i,ii) = interp1(st(msk), CSOS.N2(msk,ii), CSOS_sig_thresh(i));
            end
        end

        % Nitrate
        if sum(~isnan(no3)) > 1 && sum(~isnan(st)) > 1
            msk = ~isnan(st) & ~isnan(no3);
            CSOS.isoNO3(i,ii) = interp1(st(msk), no3(msk), CSOS_sig_thresh(i));
        end

        % Phosphate
        if sum(~isnan(po4)) > 1 && sum(~isnan(st)) > 1
            msk = ~isnan(st) & ~isnan(po4);
            CSOS.isoPO4(i,ii) = interp1(st(msk), po4(msk), CSOS_sig_thresh(i));
        end
    end
end

CSOS.isoDO(CSOS.isoDO > 1000) = NaN;

% -----------------
% P4 isopycnals
% -----------------
sig_thresh = [26.2 P4_st_median 26.2:0.05:24.05];

P4.isoDO    = NaN(length(sig_thresh),size(P4.ox,2));
P4.isoAOU   = NaN(length(sig_thresh),size(P4.ox,2));
P4.isoT     = NaN(length(sig_thresh),size(P4.ox,2));
P4.isoNO3   = NaN(length(sig_thresh),size(P4.ox,2));
P4.isoPO4   = NaN(length(sig_thresh),size(P4.ox,2));
P4.isoSR    = NaN(length(sig_thresh),size(P4.ox,2));
P4.isoSpice = NaN(length(sig_thresh),size(P4.ox,2));
P4.isoN2    = NaN(length(sig_thresh),size(P4.ox,2));
P4.isoZ     = NaN(length(sig_thresh),size(P4.ox,2));

for i = 1:length(sig_thresh)
    for ii = 1:size(P4.ox,2)
        ox  = P4.ox(:,ii);
        no3 = P4.no3(:,ii);
        po4 = P4.po4(:,ii);

        if sum(~isnan(ox))>1
            tmp_aou = aou(P4.sal(:,ii),P4.temp(:,ii),ox);
            t       = P4.temp(:,ii);
            SR      = P4.SR(:,ii);
            spice   = P4.spice(:,ii);
            msk     = ~isnan(P4.s_t(:,ii)) & ~isnan(ox);

            P4.isoDO(i,ii)    = interp1(P4.s_t(msk,ii),ox(msk),sig_thresh(i));
            P4.isoAOU(i,ii)   = interp1(P4.s_t(msk,ii),tmp_aou(msk),sig_thresh(i));
            P4.isoT(i,ii)     = interp1(P4.s_t(msk,ii),t(msk),sig_thresh(i));
            P4.isoSR(i,ii)    = interp1(P4.s_t(msk,ii),SR(msk),sig_thresh(i));
            P4.isoSpice(i,ii) = interp1(P4.s_t(msk,ii),spice(msk),sig_thresh(i));
            P4.isoZ(i,ii)     = interp1(P4.s_t(msk,ii),P4.pr(msk,ii),sig_thresh(i));
        end

        if sum(~isnan(no3))>1
            msk = ~isnan(P4.s_t(:,ii)) & ~isnan(no3);
            P4.isoNO3(i,ii) = interp1(P4.s_t(msk,ii),no3(msk),sig_thresh(i));
        end
        
        if sum(~isnan(po4))>1
            msk = ~isnan(P4.s_t(:,ii)) & ~isnan(po4);
            P4.isoPO4(i,ii) = interp1(P4.s_t(msk,ii),po4(msk),sig_thresh(i));
        end
    end
end
P4.isoPO4(2,(P4.isoPO4(2,:)<1.5)) = NaN;

% -----------------
% Masks for extreme hypoxia years
% -----------------
P4msk    = year(P4.mtime)    == 2021 & month(P4.mtime)    >= 5 & month(P4.mtime)    <= 10;
CSmsk    = year(CSOS.mtime)  == 2023 & month(CSOS.mtime)  >= 5 & month(CSOS.mtime)  <= 10;


%% fPEW vertical structure + interannual PEW-core variability

f1 = figure('units','centimeters','outerposition',[0.01 0.01 11 22],'color','w');

% PSUW mean profiles
mean_alpha_SAW = mean(alpha_PEW_SAW_source,2,'omitnan');
mean_z_SAW = mean(SAWz_monthly_grid,2,'omitnan');

% ============================================================
% Annual PEW-core diagnostics for lower panels
% ============================================================

Y_cs = year(CSOS.mtime);
Y_p4 = year(P4.mtime);

% Full year grids so gaps remain visible
years_cs = min(Y_cs(isfinite(Y_cs))):max(Y_cs(isfinite(Y_cs)));
years_p4 = 1980:max(Y_p4(isfinite(Y_p4)));

mean_z_cs  = NaN(size(years_cs));
mean_st_cs = NaN(size(years_cs));
std_z_cs   = NaN(size(years_cs));
std_st_cs  = NaN(size(years_cs));

mean_z_p4  = NaN(size(years_p4));
mean_st_p4 = NaN(size(years_p4));
std_z_p4   = NaN(size(years_p4));
std_st_p4  = NaN(size(years_p4));

% -----------------
% CS-OS
% -----------------

for k = 1:length(years_cs)

    idx = Y_cs == years_cs(k);

    if any(idx)

        if sum(isfinite(CS_z_all(idx))) > 1
            mean_z_cs(k) = mean(CS_z_all(idx),'omitnan');
            std_z_cs(k)  = std(CS_z_all(idx),'omitnan');
        end

        if sum(isfinite(CS_st_all(idx))) > 1
            mean_st_cs(k) = mean(CS_st_all(idx),'omitnan');
            std_st_cs(k)  = std(CS_st_all(idx),'omitnan');
        end

    end
end

% -----------------
% P4
% -----------------

for k = 1:length(years_p4)

    idx = Y_p4 == years_p4(k);

    if any(idx)

        if sum(isfinite(P4_z_all(idx))) > 1
            mean_z_p4(k) = mean(P4_z_all(idx),'omitnan');
            std_z_p4(k)  = std(P4_z_all(idx),'omitnan');
        end

        if sum(isfinite(P4_st_all(idx))) > 1
            mean_st_p4(k) = mean(P4_st_all(idx),'omitnan');
            std_st_p4(k)  = std(P4_st_all(idx),'omitnan');
        end

    end
end


% ============================================================
% Manual panel positions
%
% Approximate equivalent of:
% subplot(6,1,1:2)
% subplot(6,1,3:4)
% subplot(6,1,5)
% subplot(6,1,6)
% ============================================================

left = 0.1;
width = 0.79;

ax_sub1 = axes('Position',[left 0.695 width 0.285]);  % depth profile
ax_sub2 = axes('Position',[left 0.395 width 0.285]);  % density profile
ax_sub3 = axes('Position',[left 0.205 width 0.145]);  % depth time series
ax_sub4 = axes('Position',[left 0.045 width 0.145]);  % density time series


% ============================================================
% Panel 1: depth vs fPEW
% ============================================================

axes(ax_sub1);
hold on;

plot(alpha_PEW_P4,1:1000,'color',[rgb_x('steel blue') 0.02],'linewidth',0.75);
plot(alpha_PEW_CSOS,1:1000,'color',[rgb_x('rose') 0.05],'linewidth',0.75);
plot(alpha_PEW_SAW_source,SAWz_monthly_grid,'color',[rgb_x('greyish purple') 0.03],'linewidth',0.75);

% -----------------
% Mean P4
% -----------------

plot(mean(alpha_PEW_P4,2,'omitnan'),1:1000,'w','linewidth',4);
lp4 = plot(mean(alpha_PEW_P4,2,'omitnan'),1:1000,'color',rgb_x('steel blue'),'linewidth',2);

% -----------------
% Mean CS-OS
% -----------------

plot(mean(alpha_PEW_CSOS,2,'omitnan'),1:1000,'w','linewidth',4);
lcs = plot(mean(alpha_PEW_CSOS,2,'omitnan'),1:1000,'color',rgb_x('rose'),'linewidth',2);

% -----------------
% Mean PSUW source
% -----------------

plot(mean_alpha_SAW,mean_z_SAW,'w','linewidth',4);
lsaw = plot(mean_alpha_SAW,mean_z_SAW,'color',rgb_x('greyish purple'),'linewidth',2);

axis ij;
ylabel('Depth (m)');
grid on;


% ============================================================
% 0.2 sigma-theta spacing
% ============================================================

for i = 25:0.2:27

    % -----------------
    % CS-OS
    % -----------------

    y = interp1(csst,1:1000,i);
    x = interp1(1:1000,mean(alpha_PEW_CSOS,2,'omitnan'),y);

    scatter(x,y,10,'filled','markerfacecolor','w','markeredgecolor',rgb_x('rose'));

    % -----------------
    % P4
    % -----------------

    y = interp1(p4st,1:1000,i);
    x = interp1(1:1000,mean(alpha_PEW_P4,2,'omitnan'),y);

    s1 = scatter(x,y,10,'filled','markerfacecolor','w','markeredgecolor',rgb_x('steel blue'));

    % -----------------
    % PSUW source
    % -----------------

    msk_SAW = isfinite(mean_z_SAW) & isfinite(sgrid) & isfinite(mean_alpha_SAW);

    z_tmp = mean_z_SAW(msk_SAW);
    st_tmp = sgrid(msk_SAW);
    a_tmp = mean_alpha_SAW(msk_SAW);

    [z_tmp,ii] = sort(z_tmp);

    st_tmp = st_tmp(ii);
    a_tmp = a_tmp(ii);

    [z_tmp,ia] = unique(z_tmp);

    st_tmp = st_tmp(ia);
    a_tmp = a_tmp(ia);

    if length(z_tmp) >= 2

        y = interp1(st_tmp,z_tmp,i);
        x = interp1(st_tmp,a_tmp,i);

        if isfinite(x) && isfinite(y)
            scatter(x,y,10,'filled','markerfacecolor','w','markeredgecolor',rgb_x('greyish purple'));
        end
    end

end


% ============================================================
% PEW cores
% ============================================================

% -----------------
% CS-OS
% -----------------

tmp = mean(alpha_PEW_CSOS,2,'omitnan');

errorbar(tmp(CS_z_max),CS_z_max,CS_z_max_std,CS_z_max_std,...
    'LineStyle','none','Color',rgb_x('rose'),'LineWidth',1.3,...
    'CapSize',4,'HandleVisibility','off');

scatter(tmp(CS_z_max),CS_z_max,30,'filled',...
    'markerfacecolor',rgb_x('light grey'),...
    'markeredgecolor',rgb_x('rose'));

% -----------------
% P4
% -----------------

tmp = mean(alpha_PEW_P4,2,'omitnan');

errorbar(tmp(P4_z_max),P4_z_max,P4_z_max_std,P4_z_max_std,...
    'LineStyle','none','Color',rgb_x('steel blue'),'LineWidth',1.3,...
    'CapSize',4,'HandleVisibility','off');

scatter(tmp(P4_z_max),P4_z_max,30,'filled',...
    'markerfacecolor',rgb_x('light grey'),...
    'markeredgecolor',rgb_x('steel blue'));

% -----------------
% PSUW source
% -----------------

errorbar(mean_alpha_SAW(SAW_i_max),SAW_z_mean,SAW_z_std,SAW_z_std,...
    'LineStyle','none','Color',rgb_x('greyish purple'),'LineWidth',1.3,...
    'CapSize',4,'HandleVisibility','off');

scatter(mean_alpha_SAW(SAW_i_max),SAW_z_mean,30,'filled',...
    'markerfacecolor',rgb_x('light grey'),...
    'markeredgecolor',rgb_x('greyish purple'));

s2 = scatter(NaN,NaN,30,'filled',...
    'markerfacecolor',rgb_x('light grey'),...
    'markeredgecolor',rgb_x('black'));

axis tight;
ylim([75 500]);

box on;

% ============================================================
% Labels above mean profiles
% ============================================================

x_p4_lab  = mean(alpha_PEW_P4,2,'omitnan');
x_p4_lab  = x_p4_lab(75);

x_cs_lab  = mean(alpha_PEW_CSOS,2,'omitnan');
x_cs_lab  = x_cs_lab(75);

x_saw_lab = mean(mean_alpha_SAW,2,'omitnan');
x_saw_lab = x_saw_lab(75);

xl = xlim;

x_p4_norm  = (x_p4_lab-xl(1))  ./ diff(xl);
x_cs_norm  = (x_cs_lab-xl(1))  ./ diff(xl);
x_saw_norm = (x_saw_lab-xl(1)) ./ diff(xl);

text(x_p4_norm,1.015,'P4', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'Color',rgb_x('steel blue'), ...
    'FontSize',8,'fontweight','bold');

text(x_cs_norm,1.015,'CS-OS', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'Color',rgb_x('rose'), ...
    'FontSize',8,'fontweight','bold');

text(x_saw_norm,1.015,'PSUW', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'Color',rgb_x('greyish purple'), ...
    'FontSize',8,'fontweight','bold');

legend([s1,s2], ...
    {'0.2 \sigma_\theta spacing','PEW core'}, ...
    'Location','southeast', ...
    'FontSize',8);

xticklabels([]);


% ============================================================
% Panel 2: density vs fPEW
% ============================================================

axes(ax_sub2);
hold on;

plot(alpha_PEW_P4,P4.s_t,'color',[rgb_x('steel blue') 0.02],'linewidth',0.75);
plot(alpha_PEW_CSOS,CSOS.s_t,'color',[rgb_x('rose') 0.05],'linewidth',0.75);
plot(alpha_PEW_SAW_source,sgrid,'color',[rgb_x('greyish purple') 0.03],'linewidth',0.75);

% -----------------
% Mean P4
% -----------------

plot(mean(alpha_PEW_P4,2,'omitnan'),p4st,'w','linewidth',4);
lp4_st = plot(mean(alpha_PEW_P4,2,'omitnan'),p4st,...
    'color',rgb_x('steel blue'),'linewidth',2);

% -----------------
% Mean CS-OS
% -----------------

plot(mean(alpha_PEW_CSOS,2,'omitnan'),csst,'w','linewidth',4);
lcs_st = plot(mean(alpha_PEW_CSOS,2,'omitnan'),csst,...
    'color',rgb_x('rose'),'linewidth',2);

% -----------------
% Mean PSUW source
% -----------------

plot(mean_alpha_SAW,sgrid,'w','linewidth',4);
lsaw_st = plot(mean_alpha_SAW,sgrid,...
    'color',rgb_x('greyish purple'),'linewidth',2);

xlabel('\it{f_{PEW}}');
ylabel('\sigma_\theta (kg m^{-3})');

grid on;
box on;


% ============================================================
% 50 m spacing markers
% ============================================================

for i = 50:50:500

    % -----------------
    % CS-OS
    % -----------------

    scatter(mean(alpha_PEW_CSOS(i,:),2,'omitnan'),csst(i),10,'filled',...
        'markerfacecolor','w','markeredgecolor',rgb_x('rose'));

    % -----------------
    % P4
    % -----------------

    s1b = scatter(mean(alpha_PEW_P4(i,:),2,'omitnan'),p4st(i),10,'filled',...
        'markerfacecolor','w','markeredgecolor',rgb_x('steel blue'));

    % -----------------
    % PSUW source
    % -----------------

    msk_SAW = isfinite(mean_z_SAW) & isfinite(sgrid) & isfinite(mean_alpha_SAW);

    z_tmp = mean_z_SAW(msk_SAW);
    st_tmp = sgrid(msk_SAW);
    a_tmp = mean_alpha_SAW(msk_SAW);

    [z_tmp,ii] = sort(z_tmp);

    st_tmp = st_tmp(ii);
    a_tmp = a_tmp(ii);

    [z_tmp,ia] = unique(z_tmp);

    st_tmp = st_tmp(ia);
    a_tmp = a_tmp(ia);

    if length(z_tmp) >= 2

        st_i = interp1(z_tmp,st_tmp,i);

        if isfinite(st_i)

            x = interp1(st_tmp,a_tmp,st_i);

            if isfinite(x)
                scatter(x,st_i,10,'filled',...
                    'markerfacecolor','w',...
                    'markeredgecolor',rgb_x('greyish purple'));
            end
        end
    end

end


% ============================================================
% PEW maxima
% ============================================================

% -----------------
% CS-OS
% -----------------

tmp = mean(alpha_PEW_CSOS,2,'omitnan');

errorbar(tmp(CS_z_max),CS_st_median,CS_st_std,CS_st_std,...
    'LineStyle','none','Color',rgb_x('rose'),'LineWidth',1.3,...
    'CapSize',4,'HandleVisibility','off');

scatter(tmp(CS_z_max),CS_st_median,30,'filled',...
    'markerfacecolor',rgb_x('light grey'),...
    'markeredgecolor',rgb_x('rose'));

% -----------------
% P4
% -----------------

tmp = mean(alpha_PEW_P4,2,'omitnan');

errorbar(tmp(P4_z_max),P4_st_median,P4_st_std,P4_st_std,...
    'LineStyle','none','Color',rgb_x('steel blue'),'LineWidth',1.3,...
    'CapSize',4,'HandleVisibility','off');

scatter(tmp(P4_z_max),P4_st_median,30,'filled',...
    'markerfacecolor',rgb_x('light grey'),...
    'markeredgecolor',rgb_x('steel blue'));

% -----------------
% PSUW source
% -----------------

errorbar(mean_alpha_SAW(SAW_i_max),SAW_st_median,SAW_st_std,SAW_st_std,...
    'LineStyle','none','Color',rgb_x('greyish purple'),'LineWidth',1.3,...
    'CapSize',4,'HandleVisibility','off');

scatter(mean_alpha_SAW(SAW_i_max),SAW_st_median,30,'filled',...
    'markerfacecolor',rgb_x('light grey'),...
    'markeredgecolor',rgb_x('greyish purple'));

s2b = scatter(NaN,NaN,30,'filled',...
    'markerfacecolor',rgb_x('light grey'),...
    'markeredgecolor',rgb_x('black'));

set(gca,'YDir','reverse');

axis tight;
ylim([26 27]);

legend([s1b,s2b],...
    {'50 m spacing','PEW core'},...
    'Location','southeast',...
    'FontSize',8);


% ============================================================
% Panel 3: PEW-core depth time series
% ============================================================

axes(ax_sub3);
hold on;

% Individual profiles
scatter(CSOS.mtime,CS_z_all,8,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.10,...
    'MarkerEdgeColor','none');

p4mask = P4.mtime > datenum(1980,1,1);

scatter(P4.mtime(p4mask),P4_z_all(p4mask),8,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.075,...
    'MarkerEdgeColor','none');

% Annual ±1 SD
errorbar(datenum(years_cs,1,1),mean_z_cs,std_z_cs,...
    'LineStyle','none',...
    'Color',rgb_x('rose'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');

errorbar(datenum(years_p4,1,1),mean_z_p4,std_z_p4,...
    'LineStyle','none',...
    'Color',rgb_x('steel blue'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');

% Annual means
lcs_z = plot(datenum(years_cs,1,1),mean_z_cs,'-',...
    'color',rgb_x('rose'),'LineWidth',1.25);

scatter(datenum(years_cs,1,1),mean_z_cs,10,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.8,...
    'MarkerEdgeColor',rgb_x('rose'),...
    'HandleVisibility','off');

lp4_z = plot(datenum(years_p4,1,1),mean_z_p4,'-',...
    'color',rgb_x('steel blue'),'LineWidth',1.25);

scatter(datenum(years_p4,1,1),mean_z_p4,10,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.8,...
    'MarkerEdgeColor',rgb_x('steel blue'),...
    'HandleVisibility','off');

ylabel('PEW core depth (m)');

set(gca,'YDir','reverse');

xlim([datenum(1980,1,1) datenum(2025,1,1)]);

datetick('x','yyyy','keeplimits');
set(gca,'XTickLabel',[]);

grid on;
box on;

% legend([lp4_z,lcs_z],{'P4','CS-OS'},'Location','best','FontSize',7);


% ============================================================
% Panel 4: PEW-core sigma-theta time series
% ============================================================

axes(ax_sub4);
hold on;

scatter(CSOS.mtime,CS_st_all,8,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.10,...
    'MarkerEdgeColor','none');

scatter(P4.mtime(p4mask),P4_st_all(p4mask),8,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.075,...
    'MarkerEdgeColor','none');

% Annual ±1 SD
errorbar(datenum(years_cs,1,1),mean_st_cs,std_st_cs,...
    'LineStyle','none',...
    'Color',rgb_x('rose'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');

errorbar(datenum(years_p4,1,1),mean_st_p4,std_st_p4,...
    'LineStyle','none',...
    'Color',rgb_x('steel blue'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');

% Annual means
op2=plot(datenum(years_cs,1,1),mean_st_cs,'-',...
    'color',rgb_x('rose'),'LineWidth',1.25);

scatter(datenum(years_cs,1,1),mean_st_cs,10,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.8,...
    'MarkerEdgeColor',rgb_x('rose'),...
    'HandleVisibility','off');

op1=plot(datenum(years_p4,1,1),mean_st_p4,'-',...
    'color',rgb_x('steel blue'),'LineWidth',1.25);

scatter(datenum(years_p4,1,1),mean_st_p4,10,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.8,...
    'MarkerEdgeColor',rgb_x('steel blue'),...
    'HandleVisibility','off');

ylabel('PEW core \sigma_\theta (kg m^{-3})');
xlabel('Year');

xlim([datenum(1980,1,1) datenum(2025,1,1)]);

datetick('x','yyyy','keeplimits');

grid on;
box on;

legend([op1 op2], 'P4', 'CS-OS', 'location', 'northeast');


% ============================================================
% Link time-series x axes
% ============================================================

linkaxes([ax_sub3 ax_sub4],'x');


% ============================================================
% Formatting
% ============================================================

set(findall(gcf,'-property','fontsize'),'fontsize',8);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');


% Inter-station correlations: PEW-core depth and sigma-theta

[~,ia,ib] = intersect(years_p4,years_cs);

% PEW-core depth
x = mean_z_p4(ia);
y = mean_z_cs(ib);
msk = isfinite(x) & isfinite(y);

[R,P] = corrcoef(x(msk),y(msk));

fprintf('P4 vs CS-OS PEW-core depth: r = %.3f, R2 = %.3f, p = %.3g, N = %d\n', ...
    R(1,2),R(1,2)^2,P(1,2),sum(msk));

% PEW-core sigma-theta
x = mean_st_p4(ia);
y = mean_st_cs(ib);
msk = isfinite(x) & isfinite(y);

[R,P] = corrcoef(x(msk),y(msk));

fprintf('P4 vs CS-OS PEW-core sigma-theta: r = %.3f, R2 = %.3f, p = %.3g, N = %d\n', ...
    R(1,2),R(1,2)^2,P(1,2),sum(msk));
% keyboard

%% fPEW vertical structure + interannual PEW-core variability
f1 = figure('units','centimeters','outerposition',[0.01 0.01 11 22],'color','w');
% PSUW mean profiles
mean_alpha_SAW = mean(alpha_PEW_SAW_source,2,'omitnan');
mean_z_SAW = mean(SAWz_monthly_grid,2,'omitnan');
% ============================================================
% Annual PEW-core diagnostics for lower panels
% ============================================================
Y_cs = year(CSOS.mtime);
Y_p4 = year(P4.mtime);
% Full year grids so gaps remain visible
years_cs = min(Y_cs(isfinite(Y_cs))):max(Y_cs(isfinite(Y_cs)));
years_p4 = 1980:max(Y_p4(isfinite(Y_p4)));
mean_z_cs  = NaN(size(years_cs));
mean_st_cs = NaN(size(years_cs));
std_z_cs   = NaN(size(years_cs));
std_st_cs  = NaN(size(years_cs));
mean_z_p4  = NaN(size(years_p4));
mean_st_p4 = NaN(size(years_p4));
std_z_p4   = NaN(size(years_p4));
std_st_p4  = NaN(size(years_p4));
% -----------------
% CS-OS
% -----------------
for k = 1:length(years_cs)
    idx = Y_cs == years_cs(k);
    if any(idx)
        if sum(isfinite(CS_z_all(idx))) > 1
            mean_z_cs(k) = mean(CS_z_all(idx),'omitnan');
            std_z_cs(k)  = std(CS_z_all(idx),'omitnan');
        end
        if sum(isfinite(CS_st_all(idx))) > 1
            mean_st_cs(k) = mean(CS_st_all(idx),'omitnan');
            std_st_cs(k)  = std(CS_st_all(idx),'omitnan');
        end
    end
end
% -----------------
% P4
% -----------------
for k = 1:length(years_p4)
    idx = Y_p4 == years_p4(k);
    if any(idx)
        if sum(isfinite(P4_z_all(idx))) > 1
            mean_z_p4(k) = mean(P4_z_all(idx),'omitnan');
            std_z_p4(k)  = std(P4_z_all(idx),'omitnan');
        end
        if sum(isfinite(P4_st_all(idx))) > 1
            mean_st_p4(k) = mean(P4_st_all(idx),'omitnan');
            std_st_p4(k)  = std(P4_st_all(idx),'omitnan');
        end
    end
end
% ============================================================
% Manual panel positions
%
% Approximate equivalent of:
% subplot(6,1,1:2)
% subplot(6,1,3:4)
% subplot(6,1,5)
% subplot(6,1,6)
% ============================================================
left = 0.1;
width = 0.79;
ax_sub1 = axes('Position',[left 0.695 width 0.285]);  % depth profile
ax_sub2 = axes('Position',[left 0.395 width 0.285]);  % density profile
ax_sub3 = axes('Position',[left 0.205 width 0.145]);  % depth time series
ax_sub4 = axes('Position',[left 0.045 width 0.145]);  % density time series
% ============================================================
% Panel 1: depth vs fPEW
% ============================================================
axes(ax_sub1);
hold on;
plot(alpha_PEW_P4,1:1000,'color',[rgb_x('steel blue') 0.02],'linewidth',0.75);
plot(alpha_PEW_CSOS,1:1000,'color',[rgb_x('rose') 0.05],'linewidth',0.75);
plot(alpha_PEW_SAW_source,SAWz_monthly_grid,'color',[rgb_x('greyish purple') 0.03],'linewidth',0.75);
% -----------------
% Mean P4
% -----------------
plot(mean(alpha_PEW_P4,2,'omitnan'),1:1000,'w','linewidth',4);
lp4 = plot(mean(alpha_PEW_P4,2,'omitnan'),1:1000,'color',rgb_x('steel blue'),'linewidth',2);
% -----------------
% Mean CS-OS
% -----------------
plot(mean(alpha_PEW_CSOS,2,'omitnan'),1:1000,'w','linewidth',4);
lcs = plot(mean(alpha_PEW_CSOS,2,'omitnan'),1:1000,'color',rgb_x('rose'),'linewidth',2);
% -----------------
% Mean PSUW source
% -----------------
plot(mean_alpha_SAW,mean_z_SAW,'w','linewidth',4);
lsaw = plot(mean_alpha_SAW,mean_z_SAW,'color',rgb_x('greyish purple'),'linewidth',2);
axis ij;
ylabel('Depth (m)');
grid on;
% ============================================================
% 0.2 sigma-theta spacing
% ============================================================
for i = 25:0.2:27
    % -----------------
    % CS-OS
    % -----------------
    y = interp1(csst,1:1000,i);
    x = interp1(1:1000,mean(alpha_PEW_CSOS,2,'omitnan'),y);
    scatter(x,y,10,'filled','markerfacecolor','w','markeredgecolor',rgb_x('rose'));
    % -----------------
    % P4
    % -----------------
    y = interp1(p4st,1:1000,i);
    x = interp1(1:1000,mean(alpha_PEW_P4,2,'omitnan'),y);
    s1 = scatter(x,y,10,'filled','markerfacecolor','w','markeredgecolor',rgb_x('steel blue'));
    % -----------------
    % PSUW source
    % -----------------
    msk_SAW = isfinite(mean_z_SAW) & isfinite(sgrid) & isfinite(mean_alpha_SAW);
    z_tmp = mean_z_SAW(msk_SAW);
    st_tmp = sgrid(msk_SAW);
    a_tmp = mean_alpha_SAW(msk_SAW);
    [z_tmp,ii] = sort(z_tmp);
    st_tmp = st_tmp(ii);
    a_tmp = a_tmp(ii);
    [z_tmp,ia] = unique(z_tmp);
    st_tmp = st_tmp(ia);
    a_tmp = a_tmp(ia);
    if length(z_tmp) >= 2
        y = interp1(st_tmp,z_tmp,i);
        x = interp1(st_tmp,a_tmp,i);
        if isfinite(x) && isfinite(y)
            scatter(x,y,10,'filled','markerfacecolor','w','markeredgecolor',rgb_x('greyish purple'));
        end
    end
end
% ============================================================
% PEW cores
% ============================================================
% -----------------
% CS-OS
% -----------------
tmp = mean(alpha_PEW_CSOS,2,'omitnan');
errorbar(tmp(CS_z_max),CS_z_max,CS_z_max_std,CS_z_max_std,...
    'LineStyle','none','Color',rgb_x('rose'),'LineWidth',1.3,...
    'CapSize',4,'HandleVisibility','off');
scatter(tmp(CS_z_max),CS_z_max,30,'filled',...
    'markerfacecolor',rgb_x('light grey'),...
    'markeredgecolor',rgb_x('rose'));
% -----------------
% P4
% -----------------
tmp = mean(alpha_PEW_P4,2,'omitnan');
errorbar(tmp(P4_z_max),P4_z_max,P4_z_max_std,P4_z_max_std,...
    'LineStyle','none','Color',rgb_x('steel blue'),'LineWidth',1.3,...
    'CapSize',4,'HandleVisibility','off');
scatter(tmp(P4_z_max),P4_z_max,30,'filled',...
    'markerfacecolor',rgb_x('light grey'),...
    'markeredgecolor',rgb_x('steel blue'));
% -----------------
% PSUW source
% -----------------
errorbar(mean_alpha_SAW(SAW_i_max),SAW_z_mean,SAW_z_std,SAW_z_std,...
    'LineStyle','none','Color',rgb_x('greyish purple'),'LineWidth',1.3,...
    'CapSize',4,'HandleVisibility','off');
scatter(mean_alpha_SAW(SAW_i_max),SAW_z_mean,30,'filled',...
    'markerfacecolor',rgb_x('light grey'),...
    'markeredgecolor',rgb_x('greyish purple'));
s2 = scatter(NaN,NaN,30,'filled',...
    'markerfacecolor',rgb_x('light grey'),...
    'markeredgecolor',rgb_x('black'));
axis tight;
ylim([75 500]);
box on;
% ============================================================
% Labels above mean profiles
% ============================================================
x_p4_lab  = mean(alpha_PEW_P4,2,'omitnan');
x_p4_lab  = x_p4_lab(75);
x_cs_lab  = mean(alpha_PEW_CSOS,2,'omitnan');
x_cs_lab  = x_cs_lab(75);
x_saw_lab = mean(mean_alpha_SAW,2,'omitnan');
x_saw_lab = x_saw_lab(75);
xl = xlim;
x_p4_norm  = (x_p4_lab-xl(1))  ./ diff(xl);
x_cs_norm  = (x_cs_lab-xl(1))  ./ diff(xl);
x_saw_norm = (x_saw_lab-xl(1)) ./ diff(xl);
text(x_p4_norm,1.015,'P4', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'Color',rgb_x('steel blue'), ...
    'FontSize',8,'fontweight','bold');
text(x_cs_norm,1.015,'CS-OS', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'Color',rgb_x('rose'), ...
    'FontSize',8,'fontweight','bold');
text(x_saw_norm,1.015,'PSUW', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'Color',rgb_x('greyish purple'), ...
    'FontSize',8,'fontweight','bold');
legend([s1,s2], ...
    {'0.2 \sigma_\theta spacing','PEW core'}, ...
    'Location','southeast', ...
    'FontSize',8);
xticklabels([]);
% ============================================================
% Panel 2: density vs fPEW
% ============================================================
axes(ax_sub2);
hold on;
plot(alpha_PEW_P4,P4.s_t,'color',[rgb_x('steel blue') 0.02],'linewidth',0.75);
plot(alpha_PEW_CSOS,CSOS.s_t,'color',[rgb_x('rose') 0.05],'linewidth',0.75);
plot(alpha_PEW_SAW_source,sgrid,'color',[rgb_x('greyish purple') 0.03],'linewidth',0.75);
% -----------------
% Mean P4
% -----------------
plot(mean(alpha_PEW_P4,2,'omitnan'),p4st,'w','linewidth',4);
lp4_st = plot(mean(alpha_PEW_P4,2,'omitnan'),p4st,...
    'color',rgb_x('steel blue'),'linewidth',2);
% -----------------
% Mean CS-OS
% -----------------
plot(mean(alpha_PEW_CSOS,2,'omitnan'),csst,'w','linewidth',4);
lcs_st = plot(mean(alpha_PEW_CSOS,2,'omitnan'),csst,...
    'color',rgb_x('rose'),'linewidth',2);
% -----------------
% Mean PSUW source
% -----------------
plot(mean_alpha_SAW,sgrid,'w','linewidth',4);
lsaw_st = plot(mean_alpha_SAW,sgrid,...
    'color',rgb_x('greyish purple'),'linewidth',2);
xlabel('\it{f_{PEW}}');
ylabel('\sigma_\theta (kg m^{-3})');
grid on;
box on;
% ============================================================
% 50 m spacing markers
% ============================================================
for i = 50:50:500
    % -----------------
    % CS-OS
    % -----------------
    scatter(mean(alpha_PEW_CSOS(i,:),2,'omitnan'),csst(i),10,'filled',...
        'markerfacecolor','w','markeredgecolor',rgb_x('rose'));
    % -----------------
    % P4
    % -----------------
    s1b = scatter(mean(alpha_PEW_P4(i,:),2,'omitnan'),p4st(i),10,'filled',...
        'markerfacecolor','w','markeredgecolor',rgb_x('steel blue'));
    % -----------------
    % PSUW source
    % -----------------
    msk_SAW = isfinite(mean_z_SAW) & isfinite(sgrid) & isfinite(mean_alpha_SAW);
    z_tmp = mean_z_SAW(msk_SAW);
    st_tmp = sgrid(msk_SAW);
    a_tmp = mean_alpha_SAW(msk_SAW);
    [z_tmp,ii] = sort(z_tmp);
    st_tmp = st_tmp(ii);
    a_tmp = a_tmp(ii);
    [z_tmp,ia] = unique(z_tmp);
    st_tmp = st_tmp(ia);
    a_tmp = a_tmp(ia);
    if length(z_tmp) >= 2
        st_i = interp1(z_tmp,st_tmp,i);
        if isfinite(st_i)
            x = interp1(st_tmp,a_tmp,st_i);
            if isfinite(x)
                scatter(x,st_i,10,'filled',...
                    'markerfacecolor','w',...
                    'markeredgecolor',rgb_x('greyish purple'));
            end
        end
    end
end
% ============================================================
% PEW maxima
% ============================================================
% -----------------
% CS-OS
% -----------------
tmp = mean(alpha_PEW_CSOS,2,'omitnan');
errorbar(tmp(CS_z_max),CS_st_median,CS_st_std,CS_st_std,...
    'LineStyle','none','Color',rgb_x('rose'),'LineWidth',1.3,...
    'CapSize',4,'HandleVisibility','off');
scatter(tmp(CS_z_max),CS_st_median,30,'filled',...
    'markerfacecolor',rgb_x('light grey'),...
    'markeredgecolor',rgb_x('rose'));
% -----------------
% P4
% -----------------
tmp = mean(alpha_PEW_P4,2,'omitnan');
errorbar(tmp(P4_z_max),P4_st_median,P4_st_std,P4_st_std,...
    'LineStyle','none','Color',rgb_x('steel blue'),'LineWidth',1.3,...
    'CapSize',4,'HandleVisibility','off');
scatter(tmp(P4_z_max),P4_st_median,30,'filled',...
    'markerfacecolor',rgb_x('light grey'),...
    'markeredgecolor',rgb_x('steel blue'));
% -----------------
% PSUW source
% -----------------
errorbar(mean_alpha_SAW(SAW_i_max),SAW_st_median,SAW_st_std,SAW_st_std,...
    'LineStyle','none','Color',rgb_x('greyish purple'),'LineWidth',1.3,...
    'CapSize',4,'HandleVisibility','off');
scatter(mean_alpha_SAW(SAW_i_max),SAW_st_median,30,'filled',...
    'markerfacecolor',rgb_x('light grey'),...
    'markeredgecolor',rgb_x('greyish purple'));
s2b = scatter(NaN,NaN,30,'filled',...
    'markerfacecolor',rgb_x('light grey'),...
    'markeredgecolor',rgb_x('black'));
set(gca,'YDir','reverse');
axis tight;
ylim([26 27]);
legend([s1b,s2b],...
    {'50 m spacing','PEW core'},...
    'Location','southeast',...
    'FontSize',8);
% ============================================================
% Panel 3: PEW-core depth time series
% ============================================================
axes(ax_sub3);
hold on;
% Individual profiles
scatter(CSOS.mtime,CS_z_all,8,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.10,...
    'MarkerEdgeColor','none');
p4mask = P4.mtime > datenum(1980,1,1);
scatter(P4.mtime(p4mask),P4_z_all(p4mask),8,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.075,...
    'MarkerEdgeColor','none');
% Annual ±1 SD
errorbar(datenum(years_cs,1,1),mean_z_cs,std_z_cs,...
    'LineStyle','none',...
    'Color',rgb_x('rose'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');
errorbar(datenum(years_p4,1,1),mean_z_p4,std_z_p4,...
    'LineStyle','none',...
    'Color',rgb_x('steel blue'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');
% Annual means
lcs_z = plot(datenum(years_cs,1,1),mean_z_cs,'-',...
    'color',rgb_x('rose'),'LineWidth',1.25);
scatter(datenum(years_cs,1,1),mean_z_cs,10,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.8,...
    'MarkerEdgeColor',rgb_x('rose'),...
    'HandleVisibility','off');
lp4_z = plot(datenum(years_p4,1,1),mean_z_p4,'-',...
    'color',rgb_x('steel blue'),'LineWidth',1.25);
scatter(datenum(years_p4,1,1),mean_z_p4,10,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.8,...
    'MarkerEdgeColor',rgb_x('steel blue'),...
    'HandleVisibility','off');
ylabel('PEW core depth (m)');
set(gca,'YDir','reverse');
xlim([datenum(1980,1,1) datenum(2025,1,1)]);
datetick('x','yyyy','keeplimits');
set(gca,'XTickLabel',[]);
grid on;
box on;
% legend([lp4_z,lcs_z],{'P4','CS-OS'},'Location','best','FontSize',7);
% ============================================================
% Panel 4: PEW-core sigma-theta time series
% ============================================================
axes(ax_sub4);
hold on;
scatter(CSOS.mtime,CS_st_all,8,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.10,...
    'MarkerEdgeColor','none');
scatter(P4.mtime(p4mask),P4_st_all(p4mask),8,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.075,...
    'MarkerEdgeColor','none');
% Annual ±1 SD
errorbar(datenum(years_cs,1,1),mean_st_cs,std_st_cs,...
    'LineStyle','none',...
    'Color',rgb_x('rose'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');
errorbar(datenum(years_p4,1,1),mean_st_p4,std_st_p4,...
    'LineStyle','none',...
    'Color',rgb_x('steel blue'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');
% Annual means
op2=plot(datenum(years_cs,1,1),mean_st_cs,'-',...
    'color',rgb_x('rose'),'LineWidth',1.25);
scatter(datenum(years_cs,1,1),mean_st_cs,10,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.8,...
    'MarkerEdgeColor',rgb_x('rose'),...
    'HandleVisibility','off');
op1=plot(datenum(years_p4,1,1),mean_st_p4,'-',...
    'color',rgb_x('steel blue'),'LineWidth',1.25);
scatter(datenum(years_p4,1,1),mean_st_p4,10,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.8,...
    'MarkerEdgeColor',rgb_x('steel blue'),...
    'HandleVisibility','off');
ylabel('PEW core \sigma_\theta (kg m^{-3})');
xlabel('Year');
xlim([datenum(1980,1,1) datenum(2025,1,1)]);
datetick('x','yyyy','keeplimits');
grid on;
box on;
legend([op1 op2], 'P4', 'CS-OS', 'location', 'northeast');
% ============================================================
% Link time-series x axes
% ============================================================
linkaxes([ax_sub3 ax_sub4],'x');
% ============================================================
% Formatting
% ============================================================
set(findall(gcf,'-property','fontsize'),'fontsize',8);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
% Inter-station correlations: PEW-core depth and sigma-theta
[~,ia,ib] = intersect(years_p4,years_cs);
% PEW-core depth
x = mean_z_p4(ia);
y = mean_z_cs(ib);
msk = isfinite(x) & isfinite(y);
[R,P] = corrcoef(x(msk),y(msk));
fprintf('P4 vs CS-OS PEW-core depth: r = %.3f, R2 = %.3f, p = %.3g, N = %d\n', ...
    R(1,2),R(1,2)^2,P(1,2),sum(msk));
% PEW-core sigma-theta
x = mean_st_p4(ia);
y = mean_st_cs(ib);
msk = isfinite(x) & isfinite(y);
[R,P] = corrcoef(x(msk),y(msk));
fprintf('P4 vs CS-OS PEW-core sigma-theta: r = %.3f, R2 = %.3f, p = %.3g, N = %d\n', ...
    R(1,2),R(1,2)^2,P(1,2),sum(msk));

%% PEW-core depth and density trends
trend_years = [];   % [] uses all valid years in the annual series above
nboot = 2000;
% Depth uses the existing core row indices, as in the plotted time series.
% Interpret depth slopes as m/yr only if row i corresponds to i metres.

% P4 core depth trend
msk = isfinite(years_p4) & isfinite(mean_z_p4);

if ~isempty(trend_years)
    msk = msk & years_p4 >= trend_years(1) & years_p4 <= trend_years(2);
end

if sum(msk) >= 3
    pfit = polyfit(years_p4(msk),mean_z_p4(msk),1);
    [R,~] = corrcoef(years_p4(msk),mean_z_p4(msk));
    [~,p_mk] = mann_kendall(mean_z_p4(msk),0.1);
    berr = bootstrap_trend_error(years_p4(msk),mean_z_p4(msk),nboot);

    fprintf('P4 core depth trend %.0f-%.0f: %.4g +/- %.4g grid rows yr^-1, r=%.2f, MK p=%.3f\n', ...
        min(years_p4(msk)),max(years_p4(msk)),pfit(1),berr,R(1,2),p_mk);

    if p_mk < 0.1
        plot(ax_sub3,datenum(years_p4(msk),6,1),...
            polyval(pfit,years_p4(msk)),...
            '--','color',rgb_x('steel blue'),...
            'linewidth',1.5,...
            'HandleVisibility','off');
    end
end

% CS-OS core depth trend
msk = isfinite(years_cs) & isfinite(mean_z_cs);

if ~isempty(trend_years)
    msk = msk & years_cs >= trend_years(1) & years_cs <= trend_years(2);
end

if sum(msk) >= 3
    pfit = polyfit(years_cs(msk),mean_z_cs(msk),1);
    [R,~] = corrcoef(years_cs(msk),mean_z_cs(msk));
    [~,p_mk] = mann_kendall(mean_z_cs(msk),0.1);
    berr = bootstrap_trend_error(years_cs(msk),mean_z_cs(msk),nboot);

    fprintf('CS-OS core depth trend %.0f-%.0f: %.4g +/- %.4g grid rows yr^-1, r=%.2f, MK p=%.3f\n', ...
        min(years_cs(msk)),max(years_cs(msk)),pfit(1),berr,R(1,2),p_mk);

    if p_mk < 0.1
        plot(ax_sub3,datenum(years_cs(msk),6,1),...
            polyval(pfit,years_cs(msk)),...
            '--','color',rgb_x('rose'),...
            'linewidth',1.5,...
            'HandleVisibility','off');
    end
end

% P4 core density trend
msk = isfinite(years_p4) & isfinite(mean_st_p4);

if ~isempty(trend_years)
    msk = msk & years_p4 >= trend_years(1) & years_p4 <= trend_years(2);
end

if sum(msk) >= 3
    pfit = polyfit(years_p4(msk),mean_st_p4(msk),1);
    [R,~] = corrcoef(years_p4(msk),mean_st_p4(msk));
    [~,p_mk] = mann_kendall(mean_st_p4(msk),0.1);
    berr = bootstrap_trend_error(years_p4(msk),mean_st_p4(msk),nboot);

    fprintf('P4 core density trend %.0f-%.0f: %.4g +/- %.4g kg m^-3 yr^-1, r=%.2f, MK p=%.3f\n', ...
        min(years_p4(msk)),max(years_p4(msk)),pfit(1),berr,R(1,2),p_mk);

    if p_mk < 0.1
        plot(ax_sub4,datenum(years_p4(msk),6,1),...
            polyval(pfit,years_p4(msk)),...
            '--','color',rgb_x('steel blue'),...
            'linewidth',1.5,...
            'HandleVisibility','off');
    end
end

% CS-OS core density trend
msk = isfinite(years_cs) & isfinite(mean_st_cs);

if ~isempty(trend_years)
    msk = msk & years_cs >= trend_years(1) & years_cs <= trend_years(2);
end

if sum(msk) >= 3
    pfit = polyfit(years_cs(msk),mean_st_cs(msk),1);
    [R,~] = corrcoef(years_cs(msk),mean_st_cs(msk));
    [~,p_mk] = mann_kendall(mean_st_cs(msk),0.1);
    berr = bootstrap_trend_error(years_cs(msk),mean_st_cs(msk),nboot);

    fprintf('CS-OS core density trend %.0f-%.0f: %.4g +/- %.4g kg m^-3 yr^-1, r=%.2f, MK p=%.3f\n', ...
        min(years_cs(msk)),max(years_cs(msk)),pfit(1),berr,R(1,2),p_mk);

    if p_mk < 0.1
        plot(ax_sub4,datenum(years_cs(msk),6,1),...
            polyval(pfit,years_cs(msk)),...
            '--','color',rgb_x('rose'),...
            'linewidth',1.5,...
            'HandleVisibility','off');
    end
end

%%
% export_fig /Users/samst/Library/CloudStorage/Dropbox/Hakai/gyre/supply/figs/PEWvsSAUWvertical.jpg -m5

%% fPEW timeseries
%
% Source-water background:
%   Individual ORAS5 horizontal grid-point fPEW values evaluated
%   independently at sigma-theta = 26.50 kg m^-3.
%
% Source-water primary time series:
%   Monthly spatially averaged source-water profiles, projected onto
%   the same fixed PEW--PSUW mixing axis.
%
% P4 and CS-OS:
%   Existing profile-based fPEW evaluated at the local PEW core.
%
% Annual means are shown with +/-1 SD of the values contributing
% to each annual mean.
%
% Trends are calculated from annual-mean fPEW.
% Dashed lines indicate significant trends (Mann-Kendall p < 0.1).

core_st = 26.50;
nboot = 2000;

% ============================================================
% Primary source-water fractions and station-core fractions
% ============================================================

fPEW_PEW = fPEW_PEW_source_core_raw(:);
fPEW_SAW = fPEW_SAW_source_core_raw(:);

fPEW_P4_core   = alpha_PEW_P4_26p5';
fPEW_CSOS_core = alpha_PEW_CSOS_26p5';

P4_time_core = datetime(P4.mtime(:),'ConvertFrom','datenum');
CSOS_time_core = datetime(CSOS.mtime(:),'ConvertFrom','datenum');


% ============================================================
% Annual-mean fPEW and annual spread
% ============================================================

% -----------------
% PEW source
% -----------------

PEW_year = year(PEW.time(:));
PEW_years = unique(PEW_year);

PEW_fPEW_annual = nan(size(PEW_years));
PEW_fPEW_annual_std = nan(size(PEW_years));

for i = 1:length(PEW_years)

    msk = ...
        PEW_year == PEW_years(i) & ...
        isfinite(fPEW_PEW);

    if any(msk)

        PEW_fPEW_annual(i) = ...
            mean(fPEW_PEW(msk),'omitnan');

        PEW_fPEW_annual_std(i) = ...
            std(fPEW_PEW(msk),'omitnan');

    end

end

PEW_fPEW_anom = ...
    PEW_fPEW_annual - ...
    mean(PEW_fPEW_annual,'omitnan');


% -----------------
% PSUW source
% -----------------

SAW_year = year(SAW.time(:));
SAW_years = unique(SAW_year);

SAW_fPEW_annual = nan(size(SAW_years));
SAW_fPEW_annual_std = nan(size(SAW_years));

for i = 1:length(SAW_years)

    msk = ...
        SAW_year == SAW_years(i) & ...
        isfinite(fPEW_SAW);

    if any(msk)

        SAW_fPEW_annual(i) = ...
            mean(fPEW_SAW(msk),'omitnan');

        SAW_fPEW_annual_std(i) = ...
            std(fPEW_SAW(msk),'omitnan');

    end

end

SAW_fPEW_anom = ...
    SAW_fPEW_annual - ...
    mean(SAW_fPEW_annual,'omitnan');


% -----------------
% P4
% -----------------

P4_year = year(P4_time_core);
P4_years = unique(P4_year);

P4_fPEW_annual = nan(size(P4_years));
P4_fPEW_annual_std = nan(size(P4_years));

for i = 1:length(P4_years)

    msk = ...
        P4_year == P4_years(i) & ...
        isfinite(fPEW_P4_core);

    if sum(msk) > 1 && P4_years(i) > 1980

        P4_fPEW_annual(i) = ...
            mean(fPEW_P4_core(msk),'omitnan');

        P4_fPEW_annual_std(i) = ...
            std(fPEW_P4_core(msk),'omitnan');

    end

end

P4_fPEW_anom = ...
    P4_fPEW_annual - ...
    mean(P4_fPEW_annual,'omitnan');


% -----------------
% CS-OS
% -----------------

CSOS_year = year(CSOS_time_core);
CSOS_years = unique(CSOS_year);

CSOS_fPEW_annual = nan(size(CSOS_years));
CSOS_fPEW_annual_std = nan(size(CSOS_years));

for i = 1:length(CSOS_years)

    msk = ...
        CSOS_year == CSOS_years(i) & ...
        isfinite(fPEW_CSOS_core);

    if sum(msk) > 1

        CSOS_fPEW_annual(i) = ...
            mean(fPEW_CSOS_core(msk),'omitnan');

        CSOS_fPEW_annual_std(i) = ...
            std(fPEW_CSOS_core(msk),'omitnan');

    end

end

CSOS_fPEW_anom = ...
    CSOS_fPEW_annual - ...
    mean(CSOS_fPEW_annual,'omitnan');


% ============================================================
% Annual-mean trend analysis
% ============================================================

% -----------------
% PEW source
% -----------------

msk = isfinite(PEW_years) & isfinite(PEW_fPEW_annual);

x = PEW_years(msk);
y = PEW_fPEW_annual(msk);

if sum(msk) > 2

    PEW_fPEW_pfit = polyfit(x,y,1);

    [R,~] = corrcoef(x,y);
    PEW_fPEW_r = R(1,2);

    [~,PEW_fPEW_p] = mann_kendall(y,0.1);

    btrend = nan(nboot,1);

    for b = 1:nboot

        ii = randi(length(x),length(x),1);

        if numel(unique(x(ii))) > 1

            pp = polyfit(x(ii),y(ii),1);
            btrend(b) = pp(1)*10;

        end

    end

    PEW_fPEW_booterr = std(btrend,'omitnan');

else

    PEW_fPEW_pfit = [NaN NaN];
    PEW_fPEW_r = NaN;
    PEW_fPEW_p = NaN;
    PEW_fPEW_booterr = NaN;

end

fprintf(['PEW source annual-mean fPEW trend %.0f-%.0f: ' ...
    '%.4g +/- %.4g decade^-1, r=%.2f, MK p=%.3g\n'], ...
    min(x),max(x),PEW_fPEW_pfit(1)*10,...
    PEW_fPEW_booterr,PEW_fPEW_r,PEW_fPEW_p);


% -----------------
% PSUW source
% -----------------

msk = isfinite(SAW_years) & isfinite(SAW_fPEW_annual);

x = SAW_years(msk);
y = SAW_fPEW_annual(msk);

if sum(msk) > 2

    SAW_fPEW_pfit = polyfit(x,y,1);

    [R,~] = corrcoef(x,y);
    SAW_fPEW_r = R(1,2);

    [~,SAW_fPEW_p] = mann_kendall(y,0.1);

    btrend = nan(nboot,1);

    for b = 1:nboot

        ii = randi(length(x),length(x),1);

        if numel(unique(x(ii))) > 1

            pp = polyfit(x(ii),y(ii),1);
            btrend(b) = pp(1)*10;

        end

    end

    SAW_fPEW_booterr = std(btrend,'omitnan');

else

    SAW_fPEW_pfit = [NaN NaN];
    SAW_fPEW_r = NaN;
    SAW_fPEW_p = NaN;
    SAW_fPEW_booterr = NaN;

end

fprintf(['PSUW source annual-mean fPEW trend %.0f-%.0f: ' ...
    '%.4g +/- %.4g decade^-1, r=%.2f, MK p=%.3g\n'], ...
    min(x),max(x),SAW_fPEW_pfit(1)*10,...
    SAW_fPEW_booterr,SAW_fPEW_r,SAW_fPEW_p);


% -----------------
% P4
% -----------------

msk = isfinite(P4_years) & isfinite(P4_fPEW_annual);

x = P4_years(msk);
y = P4_fPEW_annual(msk);

if sum(msk) > 2

    P4_fPEW_pfit = polyfit(x,y,1);

    [R,~] = corrcoef(x,y);
    P4_fPEW_r = R(1,2);

    [~,P4_fPEW_p] = mann_kendall(y,0.1);

    btrend = nan(nboot,1);

    for b = 1:nboot

        ii = randi(length(x),length(x),1);

        if numel(unique(x(ii))) > 1

            pp = polyfit(x(ii),y(ii),1);
            btrend(b) = pp(1)*10;

        end

    end

    P4_fPEW_booterr = std(btrend,'omitnan');

else

    P4_fPEW_pfit = [NaN NaN];
    P4_fPEW_r = NaN;
    P4_fPEW_p = NaN;
    P4_fPEW_booterr = NaN;

end

fprintf(['P4 annual-mean fPEW trend %.0f-%.0f: ' ...
    '%.4g +/- %.4g decade^-1, r=%.2f, MK p=%.3g\n'], ...
    min(x),max(x),P4_fPEW_pfit(1)*10,...
    P4_fPEW_booterr,P4_fPEW_r,P4_fPEW_p);


% -----------------
% CS-OS
% -----------------

msk = isfinite(CSOS_years) & isfinite(CSOS_fPEW_annual);

x = CSOS_years(msk);
y = CSOS_fPEW_annual(msk);

if sum(msk) > 2

    CSOS_fPEW_pfit = polyfit(x,y,1);

    [R,~] = corrcoef(x,y);
    CSOS_fPEW_r = R(1,2);

    [~,CSOS_fPEW_p] = mann_kendall(y,0.1);

    btrend = nan(nboot,1);

    for b = 1:nboot

        ii = randi(length(x),length(x),1);

        if numel(unique(x(ii))) > 1

            pp = polyfit(x(ii),y(ii),1);
            btrend(b) = pp(1)*10;

        end

    end

    CSOS_fPEW_booterr = std(btrend,'omitnan');

else

    CSOS_fPEW_pfit = [NaN NaN];
    CSOS_fPEW_r = NaN;
    CSOS_fPEW_p = NaN;
    CSOS_fPEW_booterr = NaN;

end

fprintf(['CS-OS annual-mean fPEW trend %.0f-%.0f: ' ...
    '%.4g +/- %.4g decade^-1, r=%.2f, MK p=%.3g\n'], ...
    min(x),max(x),CSOS_fPEW_pfit(1)*10,...
    CSOS_fPEW_booterr,CSOS_fPEW_r,CSOS_fPEW_p);


% ============================================================
% Figure
% ============================================================

f2 = figure( ...
    'units','centimeters', ...
    'outerposition',[0.01 0.01 13 15], ...
    'color','w');


% ============================================================
% Main panel
% ============================================================

ax1 = axes('Position',[0.125 0.3 0.86 0.59]);
hold on;


% ============================================================
% Source-water spatial envelopes: monthly 10th-90th percentiles
% ============================================================

PEW_q10 = nan(size(PEW_grid_time));
PEW_q90 = nan(size(PEW_grid_time));

for i = 1:numel(PEW_grid_time)
    tmp = PEW_grid_fPEW(:,i);
    tmp = tmp(isfinite(tmp));

    if ~isempty(tmp)
        q = prctile(tmp,[10 90]);
        PEW_q10(i) = q(1);
        PEW_q90(i) = q(2);
    end
end

SAW_q10 = nan(size(SAW_grid_time));
SAW_q90 = nan(size(SAW_grid_time));

for i = 1:numel(SAW_grid_time)
    tmp = SAW_grid_fPEW(:,i);
    tmp = tmp(isfinite(tmp));

    if ~isempty(tmp)
        q = prctile(tmp,[10 90]);
        SAW_q10(i) = q(1);
        SAW_q90(i) = q(2);
    end
end


% -----------------
% PEW envelope
% -----------------

msk = isfinite(PEW_q10) & isfinite(PEW_q90);

fill( ...
    [PEW_grid_time(msk); flipud(PEW_grid_time(msk))], ...
    [PEW_q10(msk); flipud(PEW_q90(msk))], ...
    rgb_x('burnt orange'), ...
    'FaceAlpha',0.1, ...
    'EdgeColor','none', ...
    'HandleVisibility','off');


% -----------------
% PSUW envelope
% -----------------

msk = isfinite(SAW_q10) & isfinite(SAW_q90);

fill( ...
    [SAW_grid_time(msk); flipud(SAW_grid_time(msk))], ...
    [SAW_q10(msk); flipud(SAW_q90(msk))], ...
    rgb_x('greyish purple'), ...
    'FaceAlpha',0.1, ...
    'EdgeColor','none', ...
    'HandleVisibility','off');


% ============================================================
% Optional monthly spatially averaged source-water values
% ============================================================

% These are retained very subtly to show the actual monthly values
% from which the source-water annual means are calculated.
% 
% scatter( ...
%     PEW.time(:), ...
%     fPEW_PEW, ...
%     5, ...
%     rgb_x('burnt orange'), ...
%     'filled', ...
%     'MarkerFaceAlpha',0.10, ...
%     'MarkerEdgeAlpha',0.10, ...
%     'HandleVisibility','off');
% 
% scatter( ...
%     SAW.time(:), ...
%     fPEW_SAW, ...
%     5, ...
%     rgb_x('greyish purple'), ...
%     'filled', ...
%     'MarkerFaceAlpha',0.10, ...
%     'MarkerEdgeAlpha',0.10, ...
%     'HandleVisibility','off');


% ============================================================
% Source-water annual means +/-1 SD
% ============================================================

errorbar( ...
    datetime(PEW_years,7,1), ...
    PEW_fPEW_annual, ...
    PEW_fPEW_annual_std, ...
    'LineStyle','none', ...
    'Color',rgb_x('burnt orange'), ...
    'LineWidth',0.25, ...
    'CapSize',2, ...
    'HandleVisibility','off');

errorbar( ...
    datetime(SAW_years,7,1), ...
    SAW_fPEW_annual, ...
    SAW_fPEW_annual_std, ...
    'LineStyle','none', ...
    'Color',rgb_x('greyish purple'), ...
    'LineWidth',0.25, ...
    'CapSize',2, ...
    'HandleVisibility','off');


% Annual source-water means
% Annual source-water means
la1 = plot(datetime(PEW_years,7,1),PEW_fPEW_annual, ...
    '-s','Color',rgb_x('burnt orange'), ...
    'LineWidth',1.2, ...
    'MarkerSize',2, ...
    'MarkerFaceColor',rgb_x('burnt orange'), ...
    'MarkerEdgeColor',rgb_x('burnt orange'));

la2 = plot(datetime(SAW_years,7,1),SAW_fPEW_annual, ...
    '-^','Color',rgb_x('greyish purple'), ...
    'LineWidth',1.2, ...
    'MarkerSize',2, ...
    'MarkerFaceColor',rgb_x('greyish purple'), ...
    'MarkerEdgeColor',rgb_x('greyish purple'));

% la1 = plot(datetime(PEW_years,7,1),PEW_fPEW_annual, ...
%     '-','Color',rgb_x('burnt orange'), ...
%     'LineWidth',1.2);
% 
% scatter(datetime(PEW_years,7,1),PEW_fPEW_annual, ...
%     14,rgb_x('burnt orange'),'filled','marker','s', ...
%     'MarkerFaceAlpha',0.85, ...
%     'MarkerEdgeColor',rgb_x('burnt orange'), ...
%     'HandleVisibility','off');
% 
% la2 = plot(datetime(SAW_years,7,1),SAW_fPEW_annual, ...
%     '-','Color',rgb_x('greyish purple'), ...
%     'LineWidth',1.2);
% 
% scatter(datetime(SAW_years,7,1),SAW_fPEW_annual, ...
%     14,rgb_x('greyish purple'),'filled','marker','s', ...
%     'MarkerFaceAlpha',0.85, ...
%     'MarkerEdgeColor',rgb_x('greyish purple'), ...
%     'HandleVisibility','off');


% ============================================================
% Significant source-water trends
% ============================================================

if PEW_fPEW_p < 0.1

    msk = ...
        isfinite(PEW_years) & ...
        isfinite(PEW_fPEW_annual);

    xtrend = linspace( ...
        min(PEW_years(msk)), ...
        max(PEW_years(msk)), ...
        200);

    ttrend = ...
        datetime(floor(xtrend),1,1) + ...
        calmonths(round(12*(xtrend-floor(xtrend))));

    ytrend = polyval(PEW_fPEW_pfit,xtrend);

    plot( ...
        ttrend,ytrend,'--', ...
        'color',rgb_x('burnt orange'), ...
        'linewidth',1.5, ...
        'HandleVisibility','off');

end


if SAW_fPEW_p < 0.1

    msk = ...
        isfinite(SAW_years) & ...
        isfinite(SAW_fPEW_annual);

    xtrend = linspace( ...
        min(SAW_years(msk)), ...
        max(SAW_years(msk)), ...
        200);

    ttrend = ...
        datetime(floor(xtrend),1,1) + ...
        calmonths(round(12*(xtrend-floor(xtrend))));

    ytrend = polyval(SAW_fPEW_pfit,xtrend);

    plot( ...
        ttrend,ytrend,'-', ...
        'color',rgb_x('white'), ...
        'linewidth',2.5, ...
        'HandleVisibility','off');

    plot( ...
        ttrend,ytrend,'--', ...
        'color',rgb_x('greyish purple'), ...
        'linewidth',1.5, ...
        'HandleVisibility','off');

end


% ============================================================
% P4 and CS-OS observations
% ============================================================

scatter( ...
    P4_time_core, ...
    fPEW_P4_core, ...
    4, ...
    rgb_x('steel blue'), ...
    'filled', ...
    'MarkerFaceAlpha',0.15, ...
    'MarkerEdgeAlpha',0.25, ...
    'HandleVisibility','off');

scatter( ...
    CSOS_time_core, ...
    fPEW_CSOS_core, ...
    4, ...
    rgb_x('rose'), ...
    'filled', ...
    'MarkerFaceAlpha',0.15, ...
    'MarkerEdgeAlpha',0.25, ...
    'HandleVisibility','off');


% ============================================================
% P4 / CS-OS annual +/-1 SD
% ============================================================

errorbar( ...
    datetime(P4_years,7,1), ...
    P4_fPEW_annual, ...
    P4_fPEW_annual_std, ...
    'LineStyle','none', ...
    'Color',rgb_x('steel blue'), ...
    'LineWidth',0.5, ...
    'CapSize',2, ...
    'HandleVisibility','off');

errorbar( ...
    datetime(CSOS_years,7,1), ...
    CSOS_fPEW_annual, ...
    CSOS_fPEW_annual_std, ...
    'LineStyle','none', ...
    'Color',rgb_x('rose'), ...
    'LineWidth',0.5, ...
    'CapSize',2, ...
    'HandleVisibility','off');


% ============================================================
% P4 / CS-OS annual means
% ============================================================

plot(datetime(P4_years,7,1),P4_fPEW_annual, ...
    'w','LineWidth',2.5,'HandleVisibility','off');

plot(datetime(CSOS_years,7,1),CSOS_fPEW_annual, ...
    'w','LineWidth',2.5,'HandleVisibility','off');

% la3 = plot(datetime(P4_years,7,1),P4_fPEW_annual, ...
%     '-','Color',rgb_x('steel blue'), ...
%     'LineWidth',1.2);
% 
% la4 = plot(datetime(CSOS_years,7,1),CSOS_fPEW_annual, ...
%     '-','Color',rgb_x('rose'), ...
%     'LineWidth',1.2);

% scatter(datetime(P4_years,7,1),P4_fPEW_annual, ...
%     14,rgb_x('steel blue'),'filled', ...
%     'MarkerFaceAlpha',0.85, ...
%     'MarkerEdgeColor',rgb_x('steel blue'), ...
%     'HandleVisibility','off');
% 
% scatter(datetime(CSOS_years,7,1),CSOS_fPEW_annual, ...
%     14,rgb_x('rose'),'filled', ...
%     'MarkerFaceAlpha',0.85, ...
%     'MarkerEdgeColor',rgb_x('rose'), ...
%     'HandleVisibility','off');
la3 = plot(datetime(P4_years,7,1),P4_fPEW_annual, ...
    '-o','Color',rgb_x('steel blue'), ...
    'LineWidth',1.2, ...
    'MarkerSize',2, ...
    'MarkerFaceColor',rgb_x('steel blue'), ...
    'MarkerEdgeColor',rgb_x('steel blue'));

la4 = plot(datetime(CSOS_years,7,1),CSOS_fPEW_annual, ...
    '-o','Color',rgb_x('rose'), ...
    'LineWidth',1.2, ...
    'MarkerSize',2, ...
    'MarkerFaceColor',rgb_x('rose'), ...
    'MarkerEdgeColor',rgb_x('rose'));


% ============================================================
% Significant station trends
% ============================================================

if P4_fPEW_p < 0.1

    msk = ...
        isfinite(P4_years) & ...
        isfinite(P4_fPEW_annual);

    xtrend = linspace( ...
        min(P4_years(msk)), ...
        max(P4_years(msk)), ...
        200);

    ttrend = ...
        datetime(floor(xtrend),1,1) + ...
        calmonths(round(12*(xtrend-floor(xtrend))));

    ytrend = polyval(P4_fPEW_pfit,xtrend);

    plot( ...
        ttrend,ytrend,'--', ...
        'color',rgb_x('steel blue'), ...
        'linewidth',1.5, ...
        'HandleVisibility','off');

end


if CSOS_fPEW_p < 0.1

    msk = ...
        isfinite(CSOS_years) & ...
        isfinite(CSOS_fPEW_annual);

    xtrend = linspace( ...
        min(CSOS_years(msk)), ...
        max(CSOS_years(msk)), ...
        200);

    ttrend = ...
        datetime(floor(xtrend),1,1) + ...
        calmonths(round(12*(xtrend-floor(xtrend))));

    ytrend = polyval(CSOS_fPEW_pfit,xtrend);

    plot( ...
        ttrend,ytrend,'-', ...
        'color',rgb_x('white'), ...
        'linewidth',2.5, ...
        'HandleVisibility','off');

    plot( ...
        ttrend,ytrend,'--', ...
        'color',rgb_x('rose'), ...
        'linewidth',1.5, ...
        'HandleVisibility','off');

end


% ============================================================
% Limits and labels
% ============================================================

all_times = [ ...
    PEW_grid_time(:); ...
    SAW_grid_time(:); ...
    P4_time_core(:); ...
    CSOS_time_core(:)];

xlims = [min(all_times) max(all_times)];

xlim(xlims);


all_fPEW = [ ...
    PEW_grid_fPEW(isfinite(PEW_grid_fPEW)); ...
    SAW_grid_fPEW(isfinite(SAW_grid_fPEW)); ...
    fPEW_P4_core(isfinite(fPEW_P4_core)); ...
    fPEW_CSOS_core(isfinite(fPEW_CSOS_core)); ...
    0; 1];

yrange = range(all_fPEW);

axis tight

ylabel('$f_{\mathrm{PEW}}$','Interpreter','latex');

grid on;
box on;

ax1.XTickLabel = [];

legend( ...
    [la1,la2,la3,la4], ...
    {'PEW source','PSUW source','P4','CS-OS'}, ...
    'Location','best', ...
    'FontSize',8, ...
    'NumColumns',1);

% ============================================================
% Lower panel: annual-mean fPEW anomaly
% ============================================================

ax2 = axes('Position',[0.125 0.12 0.86 0.16]);
hold on;

% -----------------
% Annual +/-1 SD
% -----------------

% errorbar( ...
%     datetime(PEW_years,7,1), ...
%     PEW_fPEW_anom, ...
%     PEW_fPEW_annual_std, ...
%     'LineStyle','none', ...
%     'Color',rgb_x('burnt orange'), ...
%     'LineWidth',0.25, ...
%     'CapSize',2, ...
%     'HandleVisibility','off');
% 
% errorbar( ...
%     datetime(SAW_years,7,1), ...
%     SAW_fPEW_anom, ...
%     SAW_fPEW_annual_std, ...
%     'LineStyle','none', ...
%     'Color',rgb_x('greyish purple'), ...
%     'LineWidth',0.25, ...
%     'CapSize',2, ...
%     'HandleVisibility','off');
% 
% errorbar( ...
%     datetime(P4_years,7,1), ...
%     P4_fPEW_anom, ...
%     P4_fPEW_annual_std, ...
%     'LineStyle','none', ...
%     'Color',rgb_x('steel blue'), ...
%     'LineWidth',0.25, ...
%     'CapSize',2, ...
%     'HandleVisibility','off');
% 
% errorbar( ...
%     datetime(CSOS_years,7,1), ...
%     CSOS_fPEW_anom, ...
%     CSOS_fPEW_annual_std, ...
%     'LineStyle','none', ...
%     'Color',rgb_x('rose'), ...
%     'LineWidth',0.25, ...
%     'CapSize',2, ...
%     'HandleVisibility','off');


% -----------------
% PEW source
% -----------------

pPEW = plot( ...
    datetime(PEW_years,7,1), ...
    PEW_fPEW_anom, ...
    '-', ...
    'color',rgb_x('burnt orange'), ...
    'LineWidth',1);

scatter( ...
    datetime(PEW_years,7,1), ...
    PEW_fPEW_anom, ...
    10, ...
    rgb_x('burnt orange'), ...
    'filled','marker','s',...
    'MarkerFaceAlpha',0.75, ...
    'MarkerEdgeColor',rgb_x('burnt orange'), ...
    'HandleVisibility','off');


% -----------------
% PSUW source
% -----------------

pSAW = plot( ...
    datetime(SAW_years,7,1), ...
    SAW_fPEW_anom, ...
    '-', ...
    'color',rgb_x('greyish purple'), ...
    'LineWidth',1);

scatter( ...
    datetime(SAW_years,7,1), ...
    SAW_fPEW_anom, ...
    10, ...
    rgb_x('greyish purple'), ...
    'filled', ...
    'MarkerFaceAlpha',0.75,'marker','^', ...
    'MarkerEdgeColor',rgb_x('greyish purple'), ...
    'HandleVisibility','off');


% -----------------
% P4
% -----------------

pP4 = plot( ...
    datetime(P4_years,7,1), ...
    P4_fPEW_anom, ...
    '-', ...
    'color',rgb_x('steel blue'), ...
    'LineWidth',1);

scatter( ...
    datetime(P4_years,7,1), ...
    P4_fPEW_anom, ...
    10, ...
    rgb_x('steel blue'), ...
    'filled', ...
    'MarkerFaceAlpha',0.75, ...
    'MarkerEdgeColor',rgb_x('steel blue'), ...
    'HandleVisibility','off');


% -----------------
% CS-OS
% -----------------

pCS = plot( ...
    datetime(CSOS_years,7,1), ...
    CSOS_fPEW_anom, ...
    '-', ...
    'color',rgb_x('rose'), ...
    'LineWidth',1);

scatter( ...
    datetime(CSOS_years,7,1), ...
    CSOS_fPEW_anom, ...
    10, ...
    rgb_x('rose'), ...
    'filled', ...
    'MarkerFaceAlpha',0.75, ...
    'MarkerEdgeColor',rgb_x('rose'), ...
    'HandleVisibility','off');


xlabel('Year');

ylabel( ...
    '$f_{\mathrm{PEW}} - \overline{f_{\mathrm{PEW}}}$', ...
    'Interpreter','latex');

xlim(xlims);

grid on;
box on;
yticks([-0.1:0.05:0.1]);


% ============================================================
% Formatting
% ============================================================

linkaxes([ax1 ax2],'x');

set(findall(f2,'-property','fontsize'), ...
    'fontsize',10);

set(findall(f2,'-property','fontname'), ...
    'fontname','Latin Modern Roman');


% ============================================================
% Interannual variability of annual-mean fPEW
% ============================================================

PEW_fPEW_ann_var = ...
    var(PEW_fPEW_annual,'omitnan');

PEW_fPEW_ann_std = ...
    std(PEW_fPEW_annual,'omitnan');

SAW_fPEW_ann_var = ...
    var(SAW_fPEW_annual,'omitnan');

SAW_fPEW_ann_std = ...
    std(SAW_fPEW_annual,'omitnan');

P4_fPEW_ann_var = ...
    var(P4_fPEW_annual,'omitnan');

P4_fPEW_ann_std = ...
    std(P4_fPEW_annual,'omitnan');

CSOS_fPEW_ann_var = ...
    var(CSOS_fPEW_annual,'omitnan');

CSOS_fPEW_ann_std = ...
    std(CSOS_fPEW_annual,'omitnan');


fprintf('\nInterannual variability in annual-mean fPEW:\n');

fprintf( ...
    'PEW source: variance = %.4f, SD = %.3f\n', ...
    PEW_fPEW_ann_var, ...
    PEW_fPEW_ann_std);

fprintf( ...
    'PSUW source: variance = %.4f, SD = %.3f\n', ...
    SAW_fPEW_ann_var, ...
    SAW_fPEW_ann_std);

fprintf( ...
    'P4: variance = %.4f, SD = %.3f\n', ...
    P4_fPEW_ann_var, ...
    P4_fPEW_ann_std);

fprintf( ...
    'CS-OS: variance = %.4f, SD = %.3f\n', ...
    CSOS_fPEW_ann_var, ...
    CSOS_fPEW_ann_std);

%%
% export_fig /Users/samst/Library/CloudStorage/Dropbox/Hakai/gyre/supply/figs/PEWvsSAUWinterannual.jpg -m5

%% Sensitivity of fPEW variability to core density +/- 0.15 sigma-theta

core0 = [26.50 26.55 26.39 26.48];   % PEW, PSUW, P4, CS-OS
dc = -0.15:0.01:0.15;

names = ["PEW source","PSUW source","P4","CS-OS"];

times = { ...
    PEW.time(:), ...
    SAW.time(:), ...
    P4_time_core(:), ...
    CSOS_time_core(:)};

S = nan(4,length(dc));
trend = nan(4,length(dc));
pval = nan(4,length(dc));

for k = 1:length(dc)

    st = core0 + dc(k);

    % Source waters
    f{1} = interp1( ...
        sgrid, ...
        alpha_PEW_PEW_source_raw, ...
        st(1))';

    f{2} = interp1( ...
        sgrid, ...
        alpha_PEW_SAW_source_raw, ...
        st(2))';

    % P4
    f{3} = nan(size(P4_time_core));

    for j = 1:numel(f{3})

        x = P4.s_t(:,j);
        y = alpha_PEW_P4(:,j);

        m = isfinite(x) & isfinite(y);

        if sum(m) > 1

            [x,ii] = unique(x(m));
            yy = y(m);

            f{3}(j) = interp1( ...
                x, ...
                yy(ii), ...
                st(3));

        end
    end

    % CS-OS
    f{4} = nan(size(CSOS_time_core));

    for j = 1:numel(f{4})

        x = CSOS.s_t(:,j);
        y = alpha_PEW_CSOS(:,j);

        m = isfinite(x) & isfinite(y);

        if sum(m) > 1

            [x,ii] = unique(x(m));
            yy = y(m);

            f{4}(j) = interp1( ...
                x, ...
                yy(ii), ...
                st(4));

        end
    end

    % Annual variability and trends
    for n = 1:4

        yr = year(times{n});
        yrs = unique(yr);

        ann = nan(size(yrs));

        for i = 1:numel(yrs)

            m = ...
                yr == yrs(i) & ...
                isfinite(f{n});

            if any(m)

                ann(i) = ...
                    mean(f{n}(m),'omitnan');

            end
        end

        m = isfinite(ann);

        % Interannual variability
        S(n,k) = std(ann(m));

        % Decadal trend
        pp = polyfit( ...
            yrs(m), ...
            ann(m), ...
            1);

        trend(n,k) = pp(1)*10;

        [~,pval(n,k)] = mann_kendall(ann(m),0.1);

    end
end

% ============================================================
% Print full sensitivity table
% ============================================================

fprintf('\nCore-density sensitivity of annual fPEW\n');
fprintf('--------------------------------------------------------------------------\n');
fprintf('%-12s %7s %9s %10s %10s\n', ...
    'Record','dSigma','SD','Trend/dec','MK p');
fprintf('--------------------------------------------------------------------------\n');

for n = 1:4

    for k = 1:length(dc)

        fprintf('%-12s %+7.2f %9.4f %10.4f %10.3g\n', ...
            names(n), ...
            dc(k), ...
            S(n,k), ...
            trend(n,k), ...
            pval(n,k));

    end
end

fprintf('--------------------------------------------------------------------------\n');

% ============================================================
% Sensitivity relative to nominal core density
% ============================================================

[~,i0] = min(abs(dc));

SD_change = ...
    100 .* ...
    (S-S(:,i0)) ./ ...
    S(:,i0);

trend_change = ...
    100 .* ...
    (trend-trend(:,i0)) ./ ...
    abs(trend(:,i0));

fprintf('\nMaximum sensitivity relative to nominal core density:\n');

for n = 1:4

    fprintf('%-12s SD: %.1f%%   trend: %.1f%%\n', ...
        names(n), ...
        max(abs(SD_change(n,:)),[],'omitnan'), ...
        max(abs(trend_change(n,:)),[],'omitnan'));

end

% ============================================================
% Range across all tested core densities
% ============================================================

fprintf('\nRange across core-density sensitivity test:\n');

for n = 1:4

    fprintf(['%-12s SD = %.4f--%.4f; ' ...
        'trend = %.4f--%.4f decade^-1\n'], ...
        names(n), ...
        min(S(n,:),[],'omitnan'), ...
        max(S(n,:),[],'omitnan'), ...
        min(trend(n,:),[],'omitnan'), ...
        max(trend(n,:),[],'omitnan'));

end

%% Interannual PEW core biogeochemical signal
clc

f1 = figure('units','centimeters','outerposition',...
    [0.01 0.01, 20, 15],'color','w');

xlims = [datenum(1992,1,1) datenum(2025,1,1)];

trend_years = [];   % set to [] to use all available years
nboot = 2000;

CSOS.isoNO = CSOS.isoDO + 9.2 .* CSOS.isoNO3;
P4.isoNO   = P4.isoDO   + 9.2 .* P4.isoNO3;

Y_cs = year(CSOS.mtime);
Y_p4 = year(P4.mtime);

% Full year grids so missing years stay NaN and plotted lines are not joined
years_cs = min(Y_cs(isfinite(Y_cs))):max(Y_cs(isfinite(Y_cs)));
years_p4 = 1980:max(Y_p4(isfinite(Y_p4)));


% ============================================================
% Annual means and standard deviations: CS-OS
% ============================================================

mean_DO_cs  = NaN(size(years_cs));
mean_N_cs   = NaN(size(years_cs));
mean_NO_cs  = NaN(size(years_cs));
mean_z_cs   = NaN(size(years_cs));
mean_st_cs  = NaN(size(years_cs));
mean_PEW_cs = NaN(size(years_cs));
mean_AOU_cs = NaN(size(years_cs));

std_DO_cs   = NaN(size(years_cs));
std_N_cs    = NaN(size(years_cs));
std_NO_cs   = NaN(size(years_cs));
std_AOU_cs  = NaN(size(years_cs));

for k = 1:length(years_cs)

    idx = Y_cs == years_cs(k);

    if any(idx)

        if sum(~isnan(CSOS.isoDO(2,idx))) > 2
            mean_DO_cs(k) = mean(CSOS.isoDO(2,idx),'omitnan');
            std_DO_cs(k)  = std(CSOS.isoDO(2,idx),'omitnan');
        end

        if sum(~isnan(CSOS.isoNO3(2,idx))) > 1
            mean_N_cs(k) = mean(CSOS.isoNO3(2,idx),'omitnan');
            std_N_cs(k)  = std(CSOS.isoNO3(2,idx),'omitnan');
        end

        if sum(~isnan(CSOS.isoNO(2,idx))) > 1
            mean_NO_cs(k) = mean(CSOS.isoNO(2,idx),'omitnan');
            std_NO_cs(k)  = std(CSOS.isoNO(2,idx),'omitnan');
        end

        if sum(~isnan(CS_z_all(idx))) > 1
            mean_z_cs(k) = mean(CS_z_all(idx),'omitnan');
        end

        if sum(~isnan(CS_st_all(idx))) > 1
            mean_st_cs(k) = mean(CS_st_all(idx),'omitnan');
        end

        if sum(~isnan(alpha_PEW_CSOS_26p5(idx))) > 1
            mean_PEW_cs(k) = mean(alpha_PEW_CSOS_26p5(idx),'omitnan');
        end

        if sum(~isnan(CSOS.isoAOU(2,idx))) > 1
            mean_AOU_cs(k) = mean(CSOS.isoAOU(2,idx),'omitnan');
            std_AOU_cs(k)  = std(CSOS.isoAOU(2,idx),'omitnan');
        end

    end
end


% ============================================================
% Annual means and standard deviations: P4
% ============================================================

mean_DO_p4  = NaN(size(years_p4));
mean_N_p4   = NaN(size(years_p4));
mean_NO_p4  = NaN(size(years_p4));
mean_z_p4   = NaN(size(years_p4));
mean_st_p4  = NaN(size(years_p4));
mean_PEW_p4 = NaN(size(years_p4));
mean_AOU_p4 = NaN(size(years_p4));

std_DO_p4   = NaN(size(years_p4));
std_N_p4    = NaN(size(years_p4));
std_NO_p4   = NaN(size(years_p4));
std_AOU_p4  = NaN(size(years_p4));

for k = 1:length(years_p4)

    idx = Y_p4 == years_p4(k);

    if any(idx)

        if sum(~isnan(P4.isoDO(2,idx))) > 1
            mean_DO_p4(k) = mean(P4.isoDO(2,idx),'omitnan');
            std_DO_p4(k)  = std(P4.isoDO(2,idx),'omitnan');
        end

        if sum(~isnan(P4.isoNO3(2,idx))) > 1
            mean_N_p4(k) = mean(P4.isoNO3(2,idx),'omitnan');
            std_N_p4(k)  = std(P4.isoNO3(2,idx),'omitnan');
        end

        if sum(~isnan(P4.isoNO(2,idx))) > 1
            mean_NO_p4(k) = mean(P4.isoNO(2,idx),'omitnan');
            std_NO_p4(k)  = std(P4.isoNO(2,idx),'omitnan');
        end

        if sum(~isnan(P4_z_all(idx))) > 1
            mean_z_p4(k) = mean(P4_z_all(idx),'omitnan');
        end

        if sum(~isnan(P4_st_all(idx))) > 1
            mean_st_p4(k) = mean(P4_st_all(idx),'omitnan');
        end

        if sum(~isnan(alpha_PEW_P4_26p5(idx))) > 1
            mean_PEW_p4(k) = mean(alpha_PEW_P4_26p5(idx),'omitnan');
        end

        if sum(~isnan(P4.isoAOU(2,idx))) > 1
            mean_AOU_p4(k) = mean(P4.isoAOU(2,idx),'omitnan');
            std_AOU_p4(k)  = std(P4.isoAOU(2,idx),'omitnan');
        end

    end
end


% ============================================================
% Manual axes positions
% ============================================================

f1.Position = [0.01 0.01 9 13];

left = 0.13;
width = 0.82;
h = 0.175;

ax1 = axes('Position',[left 0.80-0.05 width h]);   % O2
ax2 = axes('Position',[left 0.59-0.05 width h]);   % Nox
ax3 = axes('Position',[left 0.38-0.05 width h]);   % NO
ax4 = axes('Position',[left 0.17-0.05 width h]);   % AOU


% ============================================================
% 1. O2 variability
% ============================================================

axes(ax1);
hold on;

% Individual observations
scatter(CSOS.mtime,CSOS.isoDO(2,:),10,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.1);

p4msk = P4.mtime > datenum(1980,1,1);

scatter(P4.mtime(p4msk),P4.isoDO(2,p4msk),10,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.075);

% Annual +/-1 SD
errorbar(datenum(years_cs,1,1),mean_DO_cs,std_DO_cs,...
    'LineStyle','none',...
    'Color',rgb_x('rose'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');

errorbar(datenum(years_p4,1,1),mean_DO_p4,std_DO_p4,...
    'LineStyle','none',...
    'Color',rgb_x('steel blue'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');

% Annual means
lcs_o = plot(datenum(years_cs,1,1),mean_DO_cs,...
    'color',rgb_x('rose'),'LineWidth',1.5);

scatter(datenum(years_cs,1,1),mean_DO_cs,10,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.7,...
    'MarkerEdgeColor',rgb_x('rose'),...
    'HandleVisibility','off');

lp4_o = plot(datenum(years_p4,1,1),mean_DO_p4,...
    'color',rgb_x('steel blue'),'LineWidth',1.5);

scatter(datenum(years_p4,1,1),mean_DO_p4,10,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.7,...
    'MarkerEdgeColor',rgb_x('steel blue'),...
    'HandleVisibility','off');

% CS-OS O2 trend
msk = isfinite(years_cs) & isfinite(mean_DO_cs);

if ~isempty(trend_years)
    msk = msk & years_cs >= trend_years(1) & years_cs <= trend_years(2);
end

pfit = polyfit(years_cs(msk),mean_DO_cs(msk),1);
[R,~] = corrcoef(years_cs(msk),mean_DO_cs(msk));
[~,p_mk] = mann_kendall(mean_DO_cs(msk),0.1);
berr = bootstrap_trend_error(years_cs(msk),mean_DO_cs(msk),nboot);

fprintf('CS-OS O2 trend %.0f-%.0f: %.4g +/- %.4g umol kg^-1 yr^-1, r=%.2f, MK p=%.3f\n', ...
    min(years_cs(msk)),max(years_cs(msk)),pfit(1),berr,R(1,2),p_mk);

if p_mk < 0.1
    plot(datenum(years_cs(msk),6,1),...
        polyval(pfit,years_cs(msk)),...
        '--','color',rgb_x('rose'),...
        'linewidth',1.5,...
        'HandleVisibility','off');
end

% P4 O2 trend
msk = isfinite(years_p4) & isfinite(mean_DO_p4);

if ~isempty(trend_years)
    msk = msk & years_p4 >= trend_years(1) & years_p4 <= trend_years(2);
end

pfit = polyfit(years_p4(msk),mean_DO_p4(msk),1);
[R,~] = corrcoef(years_p4(msk),mean_DO_p4(msk));
[~,p_mk] = mann_kendall(mean_DO_p4(msk),0.1);
berr = bootstrap_trend_error(years_p4(msk),mean_DO_p4(msk),nboot);

fprintf('P4 O2 trend %.0f-%.0f: %.4g +/- %.4g umol kg^-1 yr^-1, r=%.2f, MK p=%.3f\n', ...
    min(years_p4(msk)),max(years_p4(msk)),pfit(1),berr,R(1,2),p_mk);

if p_mk < 0.1
    plot(datenum(years_p4(msk),6,1),...
        polyval(pfit,years_p4(msk)),...
        '--','color',rgb_x('steel blue'),...
        'linewidth',1.5,...
        'HandleVisibility','off');
end

% Inter-station O2 correlation
[~,ia,ib] = intersect(years_p4,years_cs);

x = mean_DO_p4(ia);
y = mean_DO_cs(ib);

msk = isfinite(x) & isfinite(y);

[R,P] = corrcoef(x(msk),y(msk));

fprintf('P4 vs CS-OS O2: r^2=%.2f, p=%.3f, N=%d\n',...
    R(1,2)^2,P(1,2),sum(msk));

legend([lp4_o lcs_o],{'P4','CS-OS'},...
    'location','northwest','fontsize',8);

ylabel('O_2 (\mumol kg^{-1})');

axis tight;
xlim(xlims);

datetick('x','yyyy','keeplimits');
set(gca,'XTickLabel',[]);

box on;
grid on;


% ============================================================
% 2. Nox variability
% ============================================================

axes(ax2);
hold on;

% Individual observations
scatter(CSOS.mtime,CSOS.isoNO3(2,:),10,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.1);

scatter(P4.mtime(p4msk),P4.isoNO3(2,p4msk),10,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.075);

% Annual +/-1 SD
errorbar(datenum(years_cs,1,1),mean_N_cs,std_N_cs,...
    'LineStyle','none',...
    'Color',rgb_x('rose'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');

errorbar(datenum(years_p4,1,1),mean_N_p4,std_N_p4,...
    'LineStyle','none',...
    'Color',rgb_x('steel blue'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');

% Annual means
plot(datenum(years_cs,1,1),mean_N_cs,...
    'color',rgb_x('rose'),'LineWidth',1.5);

scatter(datenum(years_cs,1,1),mean_N_cs,10,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.7,...
    'MarkerEdgeColor',rgb_x('rose'),...
    'HandleVisibility','off');

plot(datenum(years_p4,1,1),mean_N_p4,...
    'color',rgb_x('steel blue'),'LineWidth',1.5);

scatter(datenum(years_p4,1,1),mean_N_p4,10,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.7,...
    'MarkerEdgeColor',rgb_x('steel blue'),...
    'HandleVisibility','off');

% CS-OS Nox trend
msk = isfinite(years_cs) & isfinite(mean_N_cs);

if ~isempty(trend_years)
    msk = msk & years_cs >= trend_years(1) & years_cs <= trend_years(2);
end

pfit = polyfit(years_cs(msk),mean_N_cs(msk),1);
[R,~] = corrcoef(years_cs(msk),mean_N_cs(msk));
[~,p_mk] = mann_kendall(mean_N_cs(msk),0.1);
berr = bootstrap_trend_error(years_cs(msk),mean_N_cs(msk),nboot);

fprintf('CS-OS Nox trend %.0f-%.0f: %.4g +/- %.4g umol kg^-1 yr^-1, r=%.2f, MK p=%.3f\n', ...
    min(years_cs(msk)),max(years_cs(msk)),pfit(1),berr,R(1,2),p_mk);

if p_mk < 0.1
    plot(datenum(years_cs(msk),6,1),...
        polyval(pfit,years_cs(msk)),...
        '--','color',rgb_x('rose'),...
        'linewidth',1.5,...
        'HandleVisibility','off');
end

% P4 Nox trend
msk = isfinite(years_p4) & isfinite(mean_N_p4);

if ~isempty(trend_years)
    msk = msk & years_p4 >= trend_years(1) & years_p4 <= trend_years(2);
end

pfit = polyfit(years_p4(msk),mean_N_p4(msk),1);
[R,~] = corrcoef(years_p4(msk),mean_N_p4(msk));
[~,p_mk] = mann_kendall(mean_N_p4(msk),0.1);
berr = bootstrap_trend_error(years_p4(msk),mean_N_p4(msk),nboot);

fprintf('P4 Nox trend %.0f-%.0f: %.4g +/- %.4g umol kg^-1 yr^-1, r=%.2f, MK p=%.3f\n', ...
    min(years_p4(msk)),max(years_p4(msk)),pfit(1),berr,R(1,2),p_mk);

if p_mk < 0.1
    plot(datenum(years_p4(msk),6,1),...
        polyval(pfit,years_p4(msk)),...
        '--','color',rgb_x('steel blue'),...
        'linewidth',1.5,...
        'HandleVisibility','off');
end

% Inter-station N correlation
[~,ia,ib] = intersect(years_p4,years_cs);

x = mean_N_p4(ia);
y = mean_N_cs(ib);

msk = isfinite(x) & isfinite(y);

[R,P] = corrcoef(x(msk),y(msk));

fprintf('P4 vs CS-OS Nox: r^2=%.2f, p=%.3f, N=%d\n',...
    R(1,2)^2,P(1,2),sum(msk));

ylabel('$N_{\mathrm{ox}}$ ($\mu$mol kg$^{-1}$)',...
    'interpreter','latex');

axis tight;
xlim(xlims);

datetick('x','yyyy','keeplimits');
set(gca,'XTickLabel',[]);

box on;
grid on;


% ============================================================
% 3. NO variability
% ============================================================

axes(ax3);
hold on;

% Individual observations
scatter(CSOS.mtime,CSOS.isoNO(2,:),10,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.1);

scatter(P4.mtime(p4msk),P4.isoNO(2,p4msk),10,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.075);

% Annual +/-1 SD
errorbar(datenum(years_cs,1,1),mean_NO_cs,std_NO_cs,...
    'LineStyle','none',...
    'Color',rgb_x('rose'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');

errorbar(datenum(years_p4,1,1),mean_NO_p4,std_NO_p4,...
    'LineStyle','none',...
    'Color',rgb_x('steel blue'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');

% Annual means
plot(datenum(years_cs,1,1),mean_NO_cs,...
    'color',rgb_x('rose'),'LineWidth',1.5);

scatter(datenum(years_cs,1,1),mean_NO_cs,10,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.7,...
    'MarkerEdgeColor',rgb_x('rose'),...
    'HandleVisibility','off');

plot(datenum(years_p4,1,1),mean_NO_p4,...
    'color',rgb_x('steel blue'),'LineWidth',1.5);

scatter(datenum(years_p4,1,1),mean_NO_p4,10,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.7,...
    'MarkerEdgeColor',rgb_x('steel blue'),...
    'HandleVisibility','off');

% CS-OS NO trend
msk = isfinite(years_cs) & isfinite(mean_NO_cs);

if ~isempty(trend_years)
    msk = msk & years_cs >= trend_years(1) & years_cs <= trend_years(2);
end

pfit = polyfit(years_cs(msk),mean_NO_cs(msk),1);
[R,~] = corrcoef(years_cs(msk),mean_NO_cs(msk));
[~,p_mk] = mann_kendall(mean_NO_cs(msk),0.1);
berr = bootstrap_trend_error(years_cs(msk),mean_NO_cs(msk),nboot);

fprintf('CS-OS NO trend %.0f-%.0f: %.4g +/- %.4g umol kg^-1 yr^-1, r=%.2f, MK p=%.3f\n', ...
    min(years_cs(msk)),max(years_cs(msk)),pfit(1),berr,R(1,2),p_mk);

if p_mk < 0.103
    plot(datenum(years_cs(msk),6,1),...
        polyval(pfit,years_cs(msk)),...
        '--','color',rgb_x('rose'),...
        'linewidth',1.5,...
        'HandleVisibility','off');
end

% P4 NO trend
msk = isfinite(years_p4) & isfinite(mean_NO_p4);

if ~isempty(trend_years)
    msk = msk & years_p4 >= trend_years(1) & years_p4 <= trend_years(2);
end

pfit = polyfit(years_p4(msk),mean_NO_p4(msk),1);
[R,~] = corrcoef(years_p4(msk),mean_NO_p4(msk));
[~,p_mk] = mann_kendall(mean_NO_p4(msk),0.1);
berr = bootstrap_trend_error(years_p4(msk),mean_NO_p4(msk),nboot);

fprintf('P4 NO trend %.0f-%.0f: %.4g +/- %.4g umol kg^-1 yr^-1, r=%.2f, MK p=%.3f\n', ...
    min(years_p4(msk)),max(years_p4(msk)),pfit(1),berr,R(1,2),p_mk);

if p_mk < 0.103
    plot(datenum(years_p4(msk),6,1),...
        polyval(pfit,years_p4(msk)),...
        '--','color',rgb_x('steel blue'),...
        'linewidth',1.5,...
        'HandleVisibility','off');
end

% Inter-station NO correlation
[~,ia,ib] = intersect(years_p4,years_cs);

x = mean_NO_p4(ia);
y = mean_NO_cs(ib);

msk = isfinite(x) & isfinite(y);

[R,P] = corrcoef(x(msk),y(msk));

fprintf('P4 vs CS-OS NO: r^2=%.2f, p=%.3f, N=%d\n',...
    R(1,2)^2,P(1,2),sum(msk));

ylabel('\it{NO} \rm{(\mumol kg^{-1})}');

axis tight;
xlim(xlims);

datetick('x','yyyy','keeplimits');
set(gca,'XTickLabel',[]);

box on;
grid on;


% ============================================================
% 4. AOU variability
% ============================================================

axes(ax4);
hold on;

% Individual observations
scatter(CSOS.mtime,CSOS.isoAOU(2,:),10,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.1);

scatter(P4.mtime(p4msk),P4.isoAOU(2,p4msk),10,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.075);

% Annual +/-1 SD
errorbar(datenum(years_cs,1,1),mean_AOU_cs,std_AOU_cs,...
    'LineStyle','none',...
    'Color',rgb_x('rose'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');

errorbar(datenum(years_p4,1,1),mean_AOU_p4,std_AOU_p4,...
    'LineStyle','none',...
    'Color',rgb_x('steel blue'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');

% Annual means
plot(datenum(years_cs,1,1),mean_AOU_cs,...
    'color',rgb_x('rose'),'LineWidth',1.5);

scatter(datenum(years_cs,1,1),mean_AOU_cs,10,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.7,...
    'MarkerEdgeColor',rgb_x('rose'),...
    'HandleVisibility','off');

plot(datenum(years_p4,1,1),mean_AOU_p4,...
    'color',rgb_x('steel blue'),'LineWidth',1.5);

scatter(datenum(years_p4,1,1),mean_AOU_p4,10,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.7,...
    'MarkerEdgeColor',rgb_x('steel blue'),...
    'HandleVisibility','off');

% CS-OS AOU trend
msk = isfinite(years_cs) & isfinite(mean_AOU_cs);

if ~isempty(trend_years)
    msk = msk & years_cs >= trend_years(1) & years_cs <= trend_years(2);
end

pfit = polyfit(years_cs(msk),mean_AOU_cs(msk),1);
[R,~] = corrcoef(years_cs(msk),mean_AOU_cs(msk));
[~,p_mk] = mann_kendall(mean_AOU_cs(msk),0.1);
berr = bootstrap_trend_error(years_cs(msk),mean_AOU_cs(msk),nboot);

fprintf('CS-OS AOU trend %.0f-%.0f: %.4g +/- %.4g umol kg^-1 yr^-1, r=%.2f, MK p=%.3f\n', ...
    min(years_cs(msk)),max(years_cs(msk)),pfit(1),berr,R(1,2),p_mk);

if p_mk < 0.1
    plot(datenum(years_cs(msk),6,1),...
        polyval(pfit,years_cs(msk)),...
        '--','color',rgb_x('rose'),...
        'linewidth',1.5,...
        'HandleVisibility','off');
end

% P4 AOU trend
msk = isfinite(years_p4) & isfinite(mean_AOU_p4);

if ~isempty(trend_years)
    msk = msk & years_p4 >= trend_years(1) & years_p4 <= trend_years(2);
end

pfit = polyfit(years_p4(msk),mean_AOU_p4(msk),1);
[R,~] = corrcoef(years_p4(msk),mean_AOU_p4(msk));
[~,p_mk] = mann_kendall(mean_AOU_p4(msk),0.1);
berr = bootstrap_trend_error(years_p4(msk),mean_AOU_p4(msk),nboot);

fprintf('P4 AOU trend %.0f-%.0f: %.4g +/- %.4g umol kg^-1 yr^-1, r=%.2f, MK p=%.3f\n', ...
    min(years_p4(msk)),max(years_p4(msk)),pfit(1),berr,R(1,2),p_mk);

if p_mk < 0.1
    plot(datenum(years_p4(msk),6,1),...
        polyval(pfit,years_p4(msk)),...
        '--','color',rgb_x('steel blue'),...
        'linewidth',1.5,...
        'HandleVisibility','off');
end

% Inter-station AOU correlation
[~,ia,ib] = intersect(years_p4,years_cs);

x = mean_AOU_p4(ia);
y = mean_AOU_cs(ib);

msk = isfinite(x) & isfinite(y);

[R,P] = corrcoef(x(msk),y(msk));

fprintf('P4 vs CS-OS AOU: r^2=%.2f, p=%.3f, N=%d\n',...
    R(1,2)^2,P(1,2),sum(msk));

ylabel('AOU (\mumol kg^{-1})');
xlabel('Year');

axis tight;
xlim(xlims);

datetick('x','yyyy','keeplimits');

box on;
grid on;


% ============================================================
% Formatting
% ============================================================

set(findall(gcf,'-property','fontsize'),'fontsize',8);

set(findall(gcf,'-property','fontname'),...
    'fontname','Latin Modern Roman');


%%
% export_fig /Users/samst/Library/CloudStorage/Dropbox/Hakai/gyre/supply/figs/interannual.jpg -m5

%% Annual diagnostics

vars = { ...
    'fPEW', ...
    'depth', ...
    'sigma-theta', ...
    'O2', ...
    'Nox', ...
    'NO', ...
    'AOU'};

P4_data = { ...
    mean_PEW_p4, ...
    mean_z_p4, ...
    mean_st_p4, ...
    mean_DO_p4, ...
    mean_N_p4, ...
    mean_NO_p4, ...
    mean_AOU_p4};

CS_data = { ...
    mean_PEW_cs, ...
    mean_z_cs, ...
    mean_st_cs, ...
    mean_DO_cs, ...
    mean_N_cs, ...
    mean_NO_cs, ...
    mean_AOU_cs};

for v = 1:length(vars)

    [P4_min, i1] = min(P4_data{v}, [], 'omitnan');
    [P4_max, i2] = max(P4_data{v}, [], 'omitnan');

    [CS_min, j1] = min(CS_data{v}, [], 'omitnan');
    [CS_max, j2] = max(CS_data{v}, [], 'omitnan');

    fprintf('\n%s\n', vars{v});

    fprintf('P4    min %.3f in %d; max %.3f in %d; range %.3f\n', ...
        P4_min, years_p4(i1), ...
        P4_max, years_p4(i2), ...
        P4_max-P4_min);

    fprintf('CS-OS min %.3f in %d; max %.3f in %d; range %.3f\n', ...
        CS_min, years_cs(j1), ...
        CS_max, years_cs(j2), ...
        CS_max-CS_min);

end


%% Plot O2, Nox, NO, AOU, and spiciness vs fPEW
% Prints correlation/regression statistics to Command Window

f1 = figure( ...
    'units','centimeters', ...
    'outerposition',[0.01 0.01 14 20], ...
    'color','w');

% ============================================================
% Manual axes positions
% ============================================================

w  = 0.39;
h  = 0.29;

x1 = 0.09;
x2 = 0.58;

y1 = 0.705;
y2 = 0.407;
y3 = 0.109;

ax1 = axes('Position',[x1 y1 w h]); % O2
ax2 = axes('Position',[x2 y1 w h]); % Nox
ax3 = axes('Position',[x1 y2 w h]); % NO
ax4 = axes('Position',[x2 y2 w h]); % AOU
ax5 = axes('Position',[x1 y3 w h]); % Spice

fprintf('\n============================================================\n');
fprintf('fPEW CORRELATION / REGRESSION STATISTICS\n');
fprintf('============================================================\n');


% ===================== O2 vs PEW =====================

axes(ax1);
hold on

fprintf('\n---------------- O2 vs fPEW ----------------\n');

% ----- CSOS -----

ll(1) = scatter(alpha_PEW_CSOS_26p5,CSOS.isoDO(2,:),10,'filled', ...
    'markerfacecolor',rgb_x('rose'),'markerfacealpha',0.4);

grid on
axis tight

msk = ~isnan(alpha_PEW_CSOS_26p5) & ~isnan(CSOS.isoDO(2,:));

[C,p] = corrcoef(alpha_PEW_CSOS_26p5(msk),CSOS.isoDO(2,msk));

r_CSOS_DO  = C(1,2);
r2_CSOS_DO = r_CSOS_DO^2;
p_CSOS_DO  = p(1,2);

coeffs = polyfit(alpha_PEW_CSOS_26p5(msk),CSOS.isoDO(2,msk),1);

slope_CSOS_DO = coeffs(1);
delta01_CSOS_DO = 0.1*slope_CSOS_DO;

fprintf(['CS-OS O2:\n' ...
         '   n = %d\n' ...
         '   r = %.3f\n' ...
         '   R2 = %.3f\n' ...
         '   p = %.3g\n' ...
         '   slope = %.3f\n' ...
         '   0.1 change in fPEW = %+.3f umol/kg change in O2\n'], ...
         sum(msk),r_CSOS_DO,r2_CSOS_DO,p_CSOS_DO, ...
         slope_CSOS_DO,delta01_CSOS_DO);

x = linspace(min(alpha_PEW_CSOS_26p5(msk)), ...
             max(alpha_PEW_CSOS_26p5(msk)));

fittedy = polyval(coeffs,x);

plot(x,fittedy,'linewidth',3,'color',rgb_x('white'));
plot(x,fittedy,'-.','linewidth',2,'color',rgb_x('rose'));


% ----- P4 -----

ll(2) = scatter(alpha_PEW_P4_26p5,P4.isoDO(2,:),10,'filled', ...
    'markerfacecolor',rgb_x('steel blue'),'markerfacealpha',0.4);

msk = ~isnan(alpha_PEW_P4_26p5) & ~isnan(P4.isoDO(2,:));

[C,p] = corrcoef(alpha_PEW_P4_26p5(msk),P4.isoDO(2,msk));

r_P4_DO  = C(1,2);
r2_P4_DO = r_P4_DO^2;
p_P4_DO  = p(1,2);

coeffs = polyfit(alpha_PEW_P4_26p5(msk),P4.isoDO(2,msk),1);

slope_P4_DO = coeffs(1);
delta01_P4_DO = 0.1*slope_P4_DO;

fprintf(['P4 O2:\n' ...
         '   n = %d\n' ...
         '   r = %.3f\n' ...
         '   R2 = %.3f\n' ...
         '   p = %.3g\n' ...
         '   slope = %.3f\n' ...
         '   0.1 change in fPEW = %+.3f umol/kg change in O2\n'], ...
         sum(msk),r_P4_DO,r2_P4_DO,p_P4_DO, ...
         slope_P4_DO,delta01_P4_DO);

x = linspace(min(alpha_PEW_P4_26p5(msk)), ...
             max(alpha_PEW_P4_26p5(msk)));

fittedy = polyval(coeffs,x);

plot(x,fittedy,'linewidth',3,'color',rgb_x('white'));
plot(x,fittedy,'-.','linewidth',2,'color',rgb_x('steel blue'));

xticks(0.1:0.1:1);
set(gca,'XTickLabel',[]);

line([min(alpha_PEW_CSOS_26p5) max(alpha_PEW_P4_26p5)], ...
    [109 109],'color','w','linewidth',1);

line([min(alpha_PEW_CSOS_26p5) max(alpha_PEW_P4_26p5)], ...
    [109 109],'color',rgb_x('rose'),'linewidth',0.5);

line([min(alpha_PEW_CSOS_26p5) max(alpha_PEW_P4_26p5)], ...
    [94 94],'color','w','linewidth',1);

line([min(alpha_PEW_CSOS_26p5) max(alpha_PEW_P4_26p5)], ...
    [94 94],'color',rgb_x('steel blue'),'linewidth',0.5);

ylabel('$\mathrm{O_2}\;(\mu\mathrm{mol}\,\mathrm{kg}^{-1})$', ...
       'Interpreter','latex');

box on

text(0.98,0.96, sprintf('CS-OS R^2 = %.2f', r2_CSOS_DO), ...
    'units','normalized','fontsize',8,'horizontalalignment','right', ...
    'verticalalignment','top','color',rgb_x('rose'));

text(0.98,0.90, sprintf('P4 R^2 = %.2f', r2_P4_DO), ...
    'units','normalized','fontsize',8,'horizontalalignment','right', ...
    'verticalalignment','top','color',rgb_x('steel blue'));


% ===================== NO3 vs PEW =====================

axes(ax2);
hold on

fprintf('\n---------------- Nox vs fPEW ----------------\n');

scatter(alpha_PEW_P4_26p5,P4.isoNO3(2,:),10,'filled', ...
    'markerfacecolor',rgb_x('steel blue'),'markerfacealpha',0.5);

scatter(alpha_PEW_CSOS_26p5,CSOS.isoNO3(2,:),10,'filled', ...
    'markerfacecolor',rgb_x('rose'),'markerfacealpha',0.5);

grid on
axis tight
box on


% ----- P4 -----

msk = ~isnan(alpha_PEW_P4_26p5) & ~isnan(P4.isoNO3(2,:));

[C,p] = corrcoef(alpha_PEW_P4_26p5(msk),P4.isoNO3(2,msk));

r_P4_NO3  = C(1,2);
r2_P4_NO3 = r_P4_NO3^2;
p_P4_NO3  = p(1,2);

coeffs = polyfit(alpha_PEW_P4_26p5(msk),P4.isoNO3(2,msk),1);

slope_P4_NO3 = coeffs(1);
delta01_P4_NO3 = 0.1*slope_P4_NO3;

fprintf(['P4 Nox:\n' ...
         '   n = %d\n' ...
         '   r = %.3f\n' ...
         '   R2 = %.3f\n' ...
         '   p = %.3g\n' ...
         '   slope = %.3f\n' ...
         '   0.1 change in fPEW = %+.3f umol/kg change in Nox\n'], ...
         sum(msk),r_P4_NO3,r2_P4_NO3,p_P4_NO3, ...
         slope_P4_NO3,delta01_P4_NO3);

x = linspace(min(alpha_PEW_P4_26p5(msk)), ...
             max(alpha_PEW_P4_26p5(msk)));

fittedy = polyval(coeffs,x);

plot(x,fittedy,'linewidth',3,'color','w');
plot(x,fittedy,'-.','linewidth',2,'color',rgb_x('steel blue'));


% ----- CSOS -----

msk = ~isnan(alpha_PEW_CSOS_26p5) & ~isnan(CSOS.isoNO3(2,:));

[C,p] = corrcoef(alpha_PEW_CSOS_26p5(msk),CSOS.isoNO3(2,msk));

r_CSOS_NO3  = C(1,2);
r2_CSOS_NO3 = r_CSOS_NO3^2;
p_CSOS_NO3  = p(1,2);

coeffs = polyfit(alpha_PEW_CSOS_26p5(msk),CSOS.isoNO3(2,msk),1);

slope_CSOS_NO3 = coeffs(1);
delta01_CSOS_NO3 = 0.1*slope_CSOS_NO3;

fprintf(['CS-OS Nox:\n' ...
         '   n = %d\n' ...
         '   r = %.3f\n' ...
         '   R2 = %.3f\n' ...
         '   p = %.3g\n' ...
         '   slope = %.3f\n' ...
         '   0.1 change in fPEW = %+.3f umol/kg change in Nox\n'], ...
         sum(msk),r_CSOS_NO3,r2_CSOS_NO3,p_CSOS_NO3, ...
         slope_CSOS_NO3,delta01_CSOS_NO3);

x = linspace(min(alpha_PEW_CSOS_26p5(msk)), ...
             max(alpha_PEW_CSOS_26p5(msk)));

fittedy = polyval(coeffs,x);

plot(x,fittedy,'linewidth',3,'color','w');
plot(x,fittedy,'-.','linewidth',2,'color',rgb_x('rose'));

ylabel('$N_{\mathrm{ox}}$ ($\mu$mol kg$^{-1}$)','interpreter','latex');

xticks(0.1:0.1:1);
set(gca,'XTickLabel',[]);

text(0.5,0.96, sprintf('CS-OS R^2 = %.2f', r2_CSOS_NO3), ...
    'units','normalized','fontsize',8,'horizontalalignment','right', ...
    'verticalalignment','top','color',rgb_x('rose'));

text(0.5,0.90, sprintf('P4 R^2 = %.2f', r2_P4_NO3), ...
    'units','normalized','fontsize',8,'horizontalalignment','right', ...
    'verticalalignment','top','color',rgb_x('steel blue'));


% ===================== NO vs PEW =====================

axes(ax3);
hold on

fprintf('\n---------------- NO vs fPEW ----------------\n');

% CSOS.isoNO = CSOS.isoDO + 8.625 .* CSOS.isoNO3;
% CSOS.isoPO = CSOS.isoDO + 138   .* CSOS.isoPO4;
% 
% CSOS.isoNO(CSOS.isoNO < 360) = NaN;
% CSOS.isoPO(CSOS.isoPO < 360) = NaN;
% 
% P4.isoNO = P4.isoDO + 8.625 .* P4.isoNO3;
% P4.isoPO = P4.isoDO + 138   .* P4.isoPO4;


% ----- CSOS -----

ll(1) = scatter(alpha_PEW_CSOS_26p5,CSOS.isoNO(2,:),10,'filled', ...
    'markerfacecolor',rgb_x('rose'),'markerfacealpha',0.5);

grid on
axis tight
box on

msk = ~isnan(alpha_PEW_CSOS_26p5) & ~isnan(CSOS.isoNO(2,:));

[C,p] = corrcoef(alpha_PEW_CSOS_26p5(msk),CSOS.isoNO(2,msk));

r_CSOS_NO  = C(1,2);
r2_CSOS_NO = r_CSOS_NO^2;
p_CSOS_NO  = p(1,2);

coeffs = polyfit(alpha_PEW_CSOS_26p5(msk),CSOS.isoNO(2,msk),1);

slope_CSOS_NO = coeffs(1);
delta01_CSOS_NO = 0.1*slope_CSOS_NO;

fprintf(['CS-OS NO:\n' ...
         '   n = %d\n' ...
         '   r = %.3f\n' ...
         '   R2 = %.3f\n' ...
         '   p = %.3g\n' ...
         '   slope = %.3f\n' ...
         '   0.1 change in fPEW = %+.3f umol/kg change in NO\n'], ...
         sum(msk),r_CSOS_NO,r2_CSOS_NO,p_CSOS_NO, ...
         slope_CSOS_NO,delta01_CSOS_NO);

x = linspace(min(alpha_PEW_CSOS_26p5(msk)), ...
             max(alpha_PEW_CSOS_26p5(msk)));

fittedy = polyval(coeffs,x);

plot(x,fittedy,'linewidth',3,'color',rgb_x('white'));
plot(x,fittedy,'-.','linewidth',2,'color',rgb_x('rose'));


% ----- P4 -----

ll(2) = scatter(alpha_PEW_P4_26p5,P4.isoNO(2,:),10,'filled', ...
    'markerfacecolor',rgb_x('steel blue'),'markerfacealpha',0.5);

msk = ~isnan(alpha_PEW_P4_26p5) & ~isnan(P4.isoNO(2,:));

[C,p] = corrcoef(alpha_PEW_P4_26p5(msk),P4.isoNO(2,msk));

r_P4_NO  = C(1,2);
r2_P4_NO = r_P4_NO^2;
p_P4_NO  = p(1,2);

coeffs = polyfit(alpha_PEW_P4_26p5(msk),P4.isoNO(2,msk),1);

slope_P4_NO = coeffs(1);
delta01_P4_NO = 0.1*slope_P4_NO;

fprintf(['P4 NO:\n' ...
         '   n = %d\n' ...
         '   r = %.3f\n' ...
         '   R2 = %.3f\n' ...
         '   p = %.3g\n' ...
         '   slope = %.3f\n' ...
         '   0.1 change in fPEW = %+.3f umol/kg change in NO\n'], ...
         sum(msk),r_P4_NO,r2_P4_NO,p_P4_NO, ...
         slope_P4_NO,delta01_P4_NO);

x = linspace(min(alpha_PEW_P4_26p5(msk)), ...
             max(alpha_PEW_P4_26p5(msk)));

fittedy = polyval(coeffs,x);

plot(x,fittedy,'linewidth',3,'color',rgb_x('white'));
plot(x,fittedy,'-.','linewidth',2,'color',rgb_x('steel blue'));

ylabel('$NO$ ($\mu$mol kg$^{-1}$)','interpreter','latex');

xticks(0.1:0.1:1);
set(gca,'XTickLabel',[]);

text(0.98,0.96, sprintf('CS-OS R^2 = %.2f', r2_CSOS_NO), ...
    'units','normalized','fontsize',8,'horizontalalignment','right', ...
    'verticalalignment','top','color',rgb_x('rose'));

text(0.98,0.90, sprintf('P4 R^2 = %.2f', r2_P4_NO), ...
    'units','normalized','fontsize',8,'horizontalalignment','right', ...
    'verticalalignment','top','color',rgb_x('steel blue'));


% ===================== AOU vs PEW =====================

axes(ax4);
hold on

fprintf('\n---------------- AOU vs fPEW ----------------\n');


% ----- CSOS -----

msk = ~isnan(alpha_PEW_CSOS_26p5) & ~isnan(CSOS.isoAOU(2,:));

[C,p] = corrcoef(alpha_PEW_CSOS_26p5(msk),CSOS.isoAOU(2,msk));

r_CSOS_AOU  = C(1,2);
r2_CSOS_AOU = r_CSOS_AOU^2;
p_CSOS_AOU  = p(1,2);

coeffs = polyfit(alpha_PEW_CSOS_26p5(msk),CSOS.isoAOU(2,msk),1);

slope_CSOS_AOU = coeffs(1);
delta01_CSOS_AOU = 0.1*slope_CSOS_AOU;

fprintf(['CS-OS AOU:\n' ...
         '   n = %d\n' ...
         '   r = %.3f\n' ...
         '   R2 = %.3f\n' ...
         '   p = %.3g\n' ...
         '   slope = %.3f\n' ...
         '   0.1 change in fPEW = %+.3f umol/kg change in AOU\n'], ...
         sum(msk),r_CSOS_AOU,r2_CSOS_AOU,p_CSOS_AOU, ...
         slope_CSOS_AOU,delta01_CSOS_AOU);

x = linspace(min(alpha_PEW_CSOS_26p5(msk)), ...
             max(alpha_PEW_CSOS_26p5(msk)));

fittedy = polyval(coeffs,x);

plot(x,fittedy,'linewidth',3,'color',rgb_x('white'));
plot(x,fittedy,'-.','linewidth',2,'color',rgb_x('rose'));

scatter(alpha_PEW_CSOS_26p5,CSOS.isoAOU(2,:),10,'filled', ...
    'markerfacecolor',rgb_x('rose'),'markerfacealpha',0.4);


% ----- P4 -----

msk = ~isnan(alpha_PEW_P4_26p5) & ~isnan(P4.isoAOU(2,:));

[C,p] = corrcoef(alpha_PEW_P4_26p5(msk),P4.isoAOU(2,msk));

r_P4_AOU  = C(1,2);
r2_P4_AOU = r_P4_AOU^2;
p_P4_AOU  = p(1,2);

coeffs = polyfit(alpha_PEW_P4_26p5(msk),P4.isoAOU(2,msk),1);

slope_P4_AOU = coeffs(1);
delta01_P4_AOU = 0.1*slope_P4_AOU;

fprintf(['P4 AOU:\n' ...
         '   n = %d\n' ...
         '   r = %.3f\n' ...
         '   R2 = %.3f\n' ...
         '   p = %.3g\n' ...
         '   slope = %.3f\n' ...
         '   0.1 change in fPEW = %+.3f umol/kg change in AOU\n'], ...
         sum(msk),r_P4_AOU,r2_P4_AOU,p_P4_AOU, ...
         slope_P4_AOU,delta01_P4_AOU);

x = linspace(min(alpha_PEW_P4_26p5(msk)), ...
             max(alpha_PEW_P4_26p5(msk)));

fittedy = polyval(coeffs,x);

plot(x,fittedy,'linewidth',3,'color',rgb_x('white'));
plot(x,fittedy,'-.','linewidth',2,'color',rgb_x('steel blue'));

scatter(alpha_PEW_P4_26p5,P4.isoAOU(2,:),10,'filled', ...
    'markerfacecolor',rgb_x('steel blue'),'markerfacealpha',0.4);

xlabel('\it{f_{PEW}}');
ylabel('AOU (\mumol kg^{-1})');

xticks(0.1:0.1:1);

grid on
box on
axis tight

text(0.5,0.96, sprintf('CS-OS R^2 = %.2f', r2_CSOS_AOU), ...
    'units','normalized','fontsize',8,'horizontalalignment','right', ...
    'verticalalignment','top','color',rgb_x('rose'));

text(0.5,0.90, sprintf('P4 R^2 = %.2f', r2_P4_AOU), ...
    'units','normalized','fontsize',8,'horizontalalignment','right', ...
    'verticalalignment','top','color',rgb_x('steel blue'));


% ===================== Spice vs PEW =====================

axes(ax5);
hold on

fprintf('\n---------------- Spiciness vs fPEW ----------------\n');


% ----- CSOS -----

msk = ~isnan(alpha_PEW_CSOS_26p5) & ~isnan(CSOS.isoSpice(2,:));

[C,p] = corrcoef(alpha_PEW_CSOS_26p5(msk),CSOS.isoSpice(2,msk));

r_CSOS_spice  = C(1,2);
r2_CSOS_spice = r_CSOS_spice^2;
p_CSOS_spice  = p(1,2);

coeffs = polyfit(alpha_PEW_CSOS_26p5(msk),CSOS.isoSpice(2,msk),1);

slope_CSOS_spice = coeffs(1);
delta01_CSOS_spice = 0.1*slope_CSOS_spice;

fprintf(['CS-OS Spiciness:\n' ...
         '   n = %d\n' ...
         '   r = %.3f\n' ...
         '   R2 = %.3f\n' ...
         '   p = %.3g\n' ...
         '   slope = %.3f\n' ...
         '   0.1 change in fPEW = %+.3f kg/m^3 change in spiciness\n'], ...
         sum(msk),r_CSOS_spice,r2_CSOS_spice,p_CSOS_spice, ...
         slope_CSOS_spice,delta01_CSOS_spice);

x = linspace(min(alpha_PEW_CSOS_26p5(msk)), ...
             max(alpha_PEW_CSOS_26p5(msk)));

fittedy = polyval(coeffs,x);

plot(x,fittedy,'linewidth',1.5,'color',rgb_x('white'));
plot(x,fittedy,'-.','linewidth',1,'color',rgb_x('rose'));

scatter(alpha_PEW_CSOS_26p5,CSOS.isoSpice(2,:),10,'filled', ...
    'markerfacecolor',rgb_x('rose'),'markerfacealpha',0.4);


% ----- P4 -----

msk = ~isnan(alpha_PEW_P4_26p5) & ~isnan(P4.isoSpice(2,:));

[C,p] = corrcoef(alpha_PEW_P4_26p5(msk),P4.isoSpice(2,msk));

r_P4_spice  = C(1,2);
r2_P4_spice = r_P4_spice^2;
p_P4_spice  = p(1,2);

coeffs = polyfit(alpha_PEW_P4_26p5(msk),P4.isoSpice(2,msk),1);

slope_P4_spice = coeffs(1);
delta01_P4_spice = 0.1*slope_P4_spice;

fprintf(['P4 Spiciness:\n' ...
         '   n = %d\n' ...
         '   r = %.3f\n' ...
         '   R2 = %.3f\n' ...
         '   p = %.3g\n' ...
         '   slope = %.3f\n' ...
         '   0.1 change in fPEW = %+.3f kg/m^3 change in spiciness\n'], ...
         sum(msk),r_P4_spice,r2_P4_spice,p_P4_spice, ...
         slope_P4_spice,delta01_P4_spice);

x = linspace(min(alpha_PEW_P4_26p5(msk)), ...
             max(alpha_PEW_P4_26p5(msk)));

fittedy = polyval(coeffs,x);

plot(x,fittedy,'linewidth',1.5,'color',rgb_x('white'));
plot(x,fittedy,'-.','linewidth',1,'color',rgb_x('steel blue'));

scatter(alpha_PEW_P4_26p5,P4.isoSpice(2,:),10,'filled', ...
    'markerfacecolor',rgb_x('steel blue'),'markerfacealpha',0.4);

xlabel('\it{f_{PEW}}');
xticks(0.1:0.1:1);

ylabel('Spiciness (kg m^{-3})');

grid on
box on
axis tight

text(0.55,0.96, sprintf('CS-OS R^2 = %.2f', r2_CSOS_spice), ...
    'units','normalized','fontsize',8,'horizontalalignment','right', ...
    'verticalalignment','top','color',rgb_x('rose'));

text(0.55,0.90, sprintf('P4 R^2 = %.2f', r2_P4_spice), ...
    'units','normalized','fontsize',8,'horizontalalignment','right', ...
    'verticalalignment','top','color',rgb_x('steel blue'));


% ============================================================
% Common formatting
% ============================================================

set(findall(gcf,'-property','fontname'), ...
    'fontname','Latin Modern Roman');

fprintf('\n============================================================\n');
fprintf('Done.\n');
fprintf('============================================================\n\n');

%%
export_fig /Users/samst/Library/CloudStorage/Dropbox/Hakai/gyre/supply/figs/DO_N_NO.jpg -m5

%% Calculate O*, N', and P'

cal_col = [230 159 0]/255; % muted mustard

% -----------------
% CSOS residuals
% -----------------
msk = ~isnan(alpha_PEW_CSOS_26p5) & ~isnan(CSOS.isoDO(2,:));
coeffs = polyfit(alpha_PEW_CSOS_26p5(msk), CSOS.isoDO(2,msk), 1);
O_fit = polyval(coeffs, alpha_PEW_CSOS_26p5);
CSOS.Ostar = CSOS.isoDO(2,:) - O_fit;    

msk = ~isnan(alpha_PEW_CSOS_26p5) & ~isnan(CSOS.isoNO3(2,:));
coeffs = polyfit(alpha_PEW_CSOS_26p5(msk), CSOS.isoNO3(2,msk), 1);
N_fit = polyval(coeffs, alpha_PEW_CSOS_26p5);
CSOS.Nres = CSOS.isoNO3(2,:) - N_fit;    
CSOS.Nres(CSOS.Nres < -4) = NaN;

% -----------------
% P4 residuals
% -----------------
msk = ~isnan(alpha_PEW_P4_26p5) & ~isnan(P4.isoDO(2,:));
coeffs = polyfit(alpha_PEW_P4_26p5(msk), P4.isoDO(2,msk), 1);
O_fit = polyval(coeffs, alpha_PEW_P4_26p5);
P4.Ostar = P4.isoDO(2,:) - O_fit;

msk = ~isnan(alpha_PEW_P4_26p5) & ~isnan(P4.isoNO3(2,:));
coeffs = polyfit(alpha_PEW_P4_26p5(msk), P4.isoNO3(2,msk), 1);
N_fit = polyval(coeffs, alpha_PEW_P4_26p5);
P4.Nres = P4.isoNO3(2,:) - N_fit;

%% O' and N' residual time series + combined O'/N' relationship

f1 = figure('units','centimeters','outerposition',[0.01 0.01 18 10],'color','w');

% ============================================================
% Manual axes positions
% [left bottom width height]
% ============================================================

left_x  = 0.085;
left_w  = 0.39;

ts_h    = 0.4;
top_y   = 0.55;
bot_y   = 0.11;

right_x = 0.54;
right_y = 0.11;
right_w = 0.42;
right_h = 0.84;

ax1 = axes('Position',[left_x  top_y   left_w  ts_h]);   % O'
ax2 = axes('Position',[left_x  bot_y   left_w  ts_h]);   % N'
ax3 = axes('Position',[right_x right_y right_w right_h]); % O'/N'

% ============================================================
% Annual means + SD for O'
% ============================================================

[Y_cs,~,~] = datevec(CSOS.mtime);
years_cs = unique(Y_cs);

mean_Ostar_CSOS = nan(size(years_cs));
std_Ostar_CSOS  = nan(size(years_cs));

for k = 1:numel(years_cs)

    idx = Y_cs == years_cs(k) & isfinite(CSOS.Ostar);

    if sum(idx) > 1
        mean_Ostar_CSOS(k) = mean(CSOS.Ostar(idx),'omitnan');
        std_Ostar_CSOS(k)  = std(CSOS.Ostar(idx),'omitnan');
    end
end


[Y_p4,~,~] = datevec(P4.mtime);
years_p4 = unique(Y_p4);

mean_Ostar_P4 = nan(size(years_p4));
std_Ostar_P4  = nan(size(years_p4));

for k = 1:numel(years_p4)

    idx = Y_p4 == years_p4(k) & isfinite(P4.Ostar);

    if sum(idx) > 1
        mean_Ostar_P4(k) = mean(P4.Ostar(idx),'omitnan');
        std_Ostar_P4(k)  = std(P4.Ostar(idx),'omitnan');
    end
end


% ============================================================
% Annual means + SD for N'
% ============================================================

[YcsN,~,~] = datevec(CSOS.mtime);
years_csN = unique(YcsN);

mean_Nres_CSOS = nan(size(years_csN));
std_Nres_CSOS  = nan(size(years_csN));

for k = 1:numel(years_csN)

    idx = YcsN == years_csN(k) & isfinite(CSOS.Nres);

    if sum(idx) > 1
        mean_Nres_CSOS(k) = mean(CSOS.Nres(idx),'omitnan');
        std_Nres_CSOS(k)  = std(CSOS.Nres(idx),'omitnan');
    end
end


[Yp4N,~,~] = datevec(P4.mtime);
years_p4N = unique(Yp4N);

mean_Nres_P4 = nan(size(years_p4N));
std_Nres_P4  = nan(size(years_p4N));

for k = 1:numel(years_p4N)

    idx = Yp4N == years_p4N(k) & isfinite(P4.Nres);

    if sum(idx) > 1
        mean_Nres_P4(k) = mean(P4.Nres(idx),'omitnan');
        std_Nres_P4(k)  = std(P4.Nres(idx),'omitnan');
    end
end
% ============================================================
% O' time series
% ============================================================

axes(ax1); hold on;

scatter(CSOS.mtime,CSOS.Ostar,10,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.10,...
    'markeredgecolor','none');

scatter(P4.mtime,P4.Ostar,10,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.10,...
    'markeredgecolor','none');

% Annual +/-1 SD
errorbar(datenum(years_cs,1,1),mean_Ostar_CSOS,std_Ostar_CSOS,...
    'LineStyle','none',...
    'Color',rgb_x('rose'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');

errorbar(datenum(years_p4,1,1),mean_Ostar_P4,std_Ostar_P4,...
    'LineStyle','none',...
    'Color',rgb_x('steel blue'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');

% Annual means
l_cs = plot(datenum(years_cs,1,1),mean_Ostar_CSOS,...
    'color',rgb_x('rose'),'LineWidth',0.75);

l_p4 = plot(datenum(years_p4,1,1),mean_Ostar_P4,...
    'color',rgb_x('steel blue'),'LineWidth',0.75);

scatter(datenum(years_cs,1,1),mean_Ostar_CSOS,12,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.8,...
    'MarkerEdgeColor',rgb_x('rose'),...
    'HandleVisibility','off');

scatter(datenum(years_p4,1,1),mean_Ostar_P4,12,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.8,...
    'MarkerEdgeColor',rgb_x('steel blue'),...
    'HandleVisibility','off');

ylabel('O'' (\mumol kg^{-1})');

box on;
grid on;
axis tight;

xlim([datenum(1992,1,1) datenum(2024,1,1)]);
xticks(datenum(1995:5:2020,1,1));
datetick('x','yyyy','keepticks','keeplimits');
yticks(-50:25:50);

ax1.XTickLabel = [];

legend([l_p4 l_cs],{'P4','CS-OS'},...
    'location','southwest','fontsize',8);
% ============================================================
% N' time series
% ============================================================

axes(ax2); hold on;

scatter(CSOS.mtime,CSOS.Nres,10,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.10,...
    'markeredgecolor','none');

scatter(P4.mtime,P4.Nres,10,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.10,...
    'markeredgecolor','none');

% Annual +/-1 SD
errorbar(datenum(years_csN,1,1),mean_Nres_CSOS,std_Nres_CSOS,...
    'LineStyle','none',...
    'Color',rgb_x('rose'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');

errorbar(datenum(years_p4N,1,1),mean_Nres_P4,std_Nres_P4,...
    'LineStyle','none',...
    'Color',rgb_x('steel blue'),...
    'LineWidth',0.25,...
    'CapSize',2,...
    'HandleVisibility','off');

% Annual means
ln_cs = plot(datenum(years_csN,1,1),mean_Nres_CSOS,...
    'color',rgb_x('rose'),'LineWidth',0.75);

ln_p4 = plot(datenum(years_p4N,1,1),mean_Nres_P4,...
    'color',rgb_x('steel blue'),'LineWidth',0.75);

scatter(datenum(years_csN,1,1),mean_Nres_CSOS,12,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'MarkerFaceAlpha',0.8,...
    'MarkerEdgeColor',rgb_x('rose'),...
    'HandleVisibility','off');

scatter(datenum(years_p4N,1,1),mean_Nres_P4,12,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'MarkerFaceAlpha',0.8,...
    'MarkerEdgeColor',rgb_x('steel blue'),...
    'HandleVisibility','off');

xlabel('Year');
ylabel('N'' (\mumol kg^{-1})');

box on;
grid on;
axis tight;

xlim([datenum(1992,1,1) datenum(2024,1,1)]);
xticks(datenum(1995:5:2020,1,1));
datetick('x','yyyy','keepticks','keeplimits');

% ============================================================
% O'/N' and O/N anomaly relationship
% ============================================================

mN_ref = -8.625;
clearvars l

% Residual O' vs N'
okN_p4 = isfinite(P4.Ostar) & isfinite(P4.Nres);
okN_cs = isfinite(CSOS.Ostar) & isfinite(CSOS.Nres);

O_N_p4  = P4.Ostar(okN_p4);
Nres_p4 = P4.Nres(okN_p4);

O_N_cs  = CSOS.Ostar(okN_cs);
Nres_cs = CSOS.Nres(okN_cs);

% Raw O vs N anomalies
okON_p4 = isfinite(P4.isoDO(2,:)) & isfinite(P4.isoNO3(2,:));
okON_cs = isfinite(CSOS.isoDO(2,:)) & isfinite(CSOS.isoNO3(2,:));

Nraw_p4  = P4.isoNO3(2,okON_p4);
Oraw_p4  = P4.isoDO(2,okON_p4);
Nrawa_p4 = Nraw_p4 - mean(Nraw_p4,'omitnan');
Orawa_p4 = Oraw_p4 - mean(Oraw_p4,'omitnan');

Nraw_cs  = CSOS.isoNO3(2,okON_cs);
Oraw_cs  = CSOS.isoDO(2,okON_cs);
Nrawa_cs = Nraw_cs - mean(Nraw_cs,'omitnan');
Orawa_cs = Oraw_cs - mean(Oraw_cs,'omitnan');

% Fits
pN_p4  = polyfit(Nres_p4,O_N_p4,1);
pN_cs  = polyfit(Nres_cs,O_N_cs,1);

pON_p4 = polyfit(Nrawa_p4,Orawa_p4,1);
pON_cs = polyfit(Nrawa_cs,Orawa_cs,1);

Ohat_p4    = polyval(pN_p4,Nres_p4);
Ohat_cs    = polyval(pN_cs,Nres_cs);
Ohat_ON_p4 = polyval(pON_p4,Nrawa_p4);
Ohat_ON_cs = polyval(pON_cs,Nrawa_cs);

R2_p4 = 1 - sum((O_N_p4 - Ohat_p4).^2) ./ sum((O_N_p4 - mean(O_N_p4)).^2);
R2_cs = 1 - sum((O_N_cs - Ohat_cs).^2) ./ sum((O_N_cs - mean(O_N_cs)).^2);

R2_ON_p4 = 1 - sum((Orawa_p4 - Ohat_ON_p4).^2) ./ sum((Orawa_p4 - mean(Orawa_p4)).^2);
R2_ON_cs = 1 - sum((Orawa_cs - Ohat_ON_cs).^2) ./ sum((Orawa_cs - mean(Orawa_cs)).^2);

fprintf('OLS O''~N'' (P4):   slope = %.3f, intercept = %.3f, R^2 = %.3f (Redfield %.3f), n=%d\n',...
    pN_p4(1),pN_p4(2),R2_p4,mN_ref,numel(O_N_p4));

fprintf('OLS O''~N'' (CSOS): slope = %.3f, intercept = %.3f, R^2 = %.3f (Redfield %.3f), n=%d\n',...
    pN_cs(1),pN_cs(2),R2_cs,mN_ref,numel(O_N_cs));

fprintf('OLS O~N anomalies (P4):   slope = %.3f, intercept = %.3f, R^2 = %.3f (Redfield %.3f), n=%d\n',...
    pON_p4(1),pON_p4(2),R2_ON_p4,mN_ref,sum(okON_p4));

fprintf('OLS O~N anomalies (CSOS): slope = %.3f, intercept = %.3f, R^2 = %.3f (Redfield %.3f), n=%d\n',...
    pON_cs(1),pON_cs(2),R2_ON_cs,mN_ref,sum(okON_cs));

axes(ax3); hold on;

% Raw O vs N anomalies
scatter(Nrawa_cs,Orawa_cs,16,'x',...
    'markeredgecolor',rgb_x('rose'),...
    'markeredgealpha',0.12,...
    'linewidth',0.75);

scatter(Nrawa_p4,Orawa_p4,16,'x',...
    'markeredgecolor',rgb_x('steel blue'),...
    'markeredgealpha',0.12,...
    'linewidth',0.75);

% Residual O' vs N'
scatter(Nres_cs,O_N_cs,12,'filled',...
    'markerfacecolor',rgb_x('rose'),...
    'markerfacealpha',0.18,...
    'markeredgecolor','none');

scatter(Nres_p4,O_N_p4,12,'filled',...
    'markerfacecolor',rgb_x('steel blue'),...
    'markerfacealpha',0.18,...
    'markeredgecolor','none');

xx = linspace(-6,6,200);

% Raw anomaly fits
l(1) = plot(xx,polyval(pON_cs,xx),':',...
    'color',rgb_x('rose'),'linewidth',1.5);

l(2) = plot(xx,polyval(pON_p4,xx),':',...
    'color',rgb_x('steel blue'),'linewidth',1.5);

% Residual fits
l(3) = plot(xx,polyval(pN_cs,xx),'-',...
    'color',rgb_x('rose'),'linewidth',2);

l(4) = plot(xx,polyval(pN_p4,xx),'-',...
    'color',rgb_x('steel blue'),'linewidth',2);

% Redfield line
l(5) = plot(xx,mN_ref*xx,'--k','linewidth',1);
xlabel('$N_{\mathrm{ox}}$ ($\mu$mol kg$^{-1}$)','interpreter','latex');
ylabel('O$_2$ ($\mu$mol kg$^{-1}$)','interpreter','latex');

xticks(-6:3:6);
yticks(-75:25:75);
xlim([-6 6]);
ylim([-75 75]);

grid on;
box on;

legend(l,{...
    sprintf('CS-OS O_2/N_{ox} = %.2f',pON_cs(1)),...
    sprintf('P4 O_2/N_{ox} = %.2f',pON_p4(1)),...
    sprintf('CS-OS O''/N'' = %.2f',pN_cs(1)),...
    sprintf('P4 O''/N'' = %.2f',pN_p4(1)),...
    sprintf('Redfield = %.2f',mN_ref)},...
    'location','northeast',...
    'fontsize',8);
legend boxoff

set(findall(gcf,'-property','fontsize'),'fontsize',8);
set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');



%%
 % export_fig /Users/samst/Library/CloudStorage/Dropbox/Hakai/gyre/supply/figs/Ostar_stoich.jpg -m5


%% Annual fPEW correlations with ENSO and NPGO

load('aENSO.mat')       % aENSO(:,1)=year, aENSO(:,2)=index
load('aNPGO2024.mat')   % NPGOa(:,1)=year, NPGOa(:,2)=index

idx  = {'ENSO','NPGO'};
Iyr  = {aENSO(:,1), NPGOa(:,1)};
Ival = {aENSO(:,2), NPGOa(:,2)};

stations = {'P4','CS-OS'};
Syr0  = {years_p4(:), years_cs(:)};
Sf0   = {mean_PEW_p4(:), mean_PEW_cs(:)};

for ii = 1:2
    for ss = 1:2

        n0 = min(numel(Syr0{ss}), numel(Sf0{ss}));
        Syr = Syr0{ss}(1:n0);
        Sf  = Sf0{ss}(1:n0);

        [~,ia,ib] = intersect(Syr, Iyr{ii});
        x = Ival{ii}(ib);
        y = Sf(ia);

        m = isfinite(x) & isfinite(y);

        [R,P] = corrcoef(x(m), y(m));
        b = polyfit(x(m), y(m), 1);

        fprintf('\n%s: fPEW vs %s\n', stations{ss}, idx{ii});
        fprintf('n = %d, r = %.3f, r2 = %.3f, p = %.4f, slope = %.4f fPEW/index unit\n', ...
            sum(m), R(1,2), R(1,2)^2, P(1,2), b(1));
    end
end

%% Lagged correlations among annual fPEW records

clear fPEWLag

lags = 4:-1:-1;

names = ["PEW source","PSUW source","P4","CS-OS"];

yrs = { ...
    PEW_years(:), ...
    SAW_years(:), ...
    P4_years(:), ...
    CSOS_years(:)};

vals = { ...
    PEW_fPEW_annual(:), ...
    SAW_fPEW_annual(:), ...
    P4_fPEW_annual(:), ...
    CSOS_fPEW_annual(:)};

pairs = [ ...
    1 2; ...
    1 3; ...
    1 4; ...
    2 3; ...
    2 4; ...
    3 4];

row = 0;

for ip = 1:size(pairs,1)

    i1 = pairs(ip,1);
    i2 = pairs(ip,2);

    for lag = lags

        x = nan(size(yrs{i2}));
        y = vals{i2};

        for k = 1:numel(yrs{i2})

            jj = find( ...
                yrs{i1} == yrs{i2}(k)-lag, ...
                1);

            if ~isempty(jj)
                x(k) = vals{i1}(jj);
            end

        end

        msk = isfinite(x) & isfinite(y);

        row = row + 1;

        fPEWLag.pair(row,1) = ...
            names(i1) + " - " + names(i2);

        fPEWLag.lag_yr(row,1) = lag;
        fPEWLag.n(row,1) = sum(msk);

        if sum(msk) >= 3

            [R,P] = corrcoef(x(msk),y(msk));

            fPEWLag.r(row,1)  = R(1,2);
            fPEWLag.R2(row,1) = R(1,2)^2;
            fPEWLag.p(row,1)  = P(1,2);

        else

            fPEWLag.r(row,1)  = NaN;
            fPEWLag.R2(row,1) = NaN;
            fPEWLag.p(row,1)  = NaN;

        end

    end

end

% Plot lag correlations

figure( ...
    'units','centimeters', ...
    'outerposition',[0.01 0.01 12 7], ...
    'color','w');

hold on

pair_names = unique(fPEWLag.pair,'stable');

cols = { ...
    rgb_x('dark grey'), ...
    rgb_x('burnt orange'), ...
    rgb_x('dark orange'), ...
    rgb_x('greyish purple'), ...
    rgb_x('purple'), ...
    rgb_x('steel blue')};

ll = gobjects(numel(pair_names),1);

for i = 1:numel(pair_names)

    idx = fPEWLag.pair == pair_names(i);

    ll(i) = plot( ...
        fPEWLag.lag_yr(idx), ...
        fPEWLag.R2(idx), ...
        '-', ...
        'color',cols{i}, ...
        'linewidth',1.5);

    sig = idx & fPEWLag.p < 0.05;
    nsig = idx & fPEWLag.p >= 0.05;

    scatter( ...
        fPEWLag.lag_yr(sig), ...
        fPEWLag.R2(sig), ...
        28, ...
        'filled', ...
        'markerfacecolor',cols{i}, ...
        'markeredgecolor','k', ...
        'linewidth',0.5);

    scatter( ...
        fPEWLag.lag_yr(nsig), ...
        fPEWLag.R2(nsig), ...
        28, ...
        'markerfacecolor','w', ...
        'markeredgecolor',cols{i}, ...
        'linewidth',1.1);

end

xline(0,'--','color',[0.5 0.5 0.5]);

xlabel('Lag (years; positive = first-listed record leads)');
ylabel('R^2');

ylim([0 1]);

legend( ...
    ll, ...
    pair_names, ...
    'location','best', ...
    'fontsize',7);

grid on
box on

set(findall(gcf,'-property','fontname'), ...
    'fontname','Latin Modern Roman');

set(findall(gcf,'-property','fontsize'), ...
    'fontsize',8);


% ============================================================
% Print significant lag correlations only
% ============================================================

sig = fPEWLag.p < 0.05;

fprintf('\nSignificant lagged annual fPEW correlations (p < 0.05)\n');
fprintf('---------------------------------------------------------------\n');
fprintf('%-24s %5s %5s %8s %8s %10s\n', ...
    'Pair','Lag','n','r','R^2','p');
fprintf('---------------------------------------------------------------\n');

for i = find(sig)'

    fprintf('%-24s %5.0f %5.0f %8.2f %8.2f %10.3g\n', ...
        fPEWLag.pair(i), ...
        fPEWLag.lag_yr(i), ...
        fPEWLag.n(i), ...
        fPEWLag.r(i), ...
        fPEWLag.R2(i), ...
        fPEWLag.p(i));

end

fprintf('---------------------------------------------------------------\n');
%% Low-oxygen event attribution: PEW vs remineralization history
% Assumes these already exist:
% CSOS.isoDO, CSOS.isoNO3, CSOS.isoNO, CSOS.isoAOU, CSOS.isoZ, CSOS.isoSpice
% P4.isoDO,   P4.isoNO3,   P4.isoNO,   P4.isoAOU,   P4.isoZ,   P4.isoSpice
% alpha_PEW_CSOS_26p5, alpha_PEW_P4_26p5
% CSOS.Ostar, CSOS.Nres, P4.Ostar, P4.Nres

clear lowO

event_prctile = 5;   % lowest 10% oxygen events
target_row    = 2;    % target isopycnal row

% ===================== CSOS event table =====================

CSOS_tbl = table;
CSOS_tbl.station = repmat("CSOS", numel(CSOS.mtime), 1);
CSOS_tbl.mtime   = CSOS.mtime(:);
CSOS_tbl.year    = year(datetime(CSOS.mtime(:),'ConvertFrom','datenum'));
CSOS_tbl.month   = month(datetime(CSOS.mtime(:),'ConvertFrom','datenum'));

CSOS_tbl.PEW     = alpha_PEW_CSOS_26p5(:);
CSOS_tbl.O2      = CSOS.isoDO(target_row,:)';
CSOS_tbl.NO3     = CSOS.isoNO3(target_row,:)';
CSOS_tbl.NO      = CSOS.isoNO(target_row,:)';
CSOS_tbl.AOU     = CSOS.isoAOU(target_row,:)';
CSOS_tbl.Z       = CSOS.isoZ(target_row,:)';
CSOS_tbl.Spice   = CSOS.isoSpice(target_row,:)';
CSOS_tbl.Ostar   = CSOS.Ostar(:);
CSOS_tbl.Nres    = CSOS.Nres(:);

CSOS_tbl.season = strings(height(CSOS_tbl),1);
CSOS_tbl.season(ismember(CSOS_tbl.month,[12 1 2])) = "winter";
CSOS_tbl.season(ismember(CSOS_tbl.month,[3 4 5]))  = "spring";
CSOS_tbl.season(ismember(CSOS_tbl.month,[6 7 8]))  = "summer";
CSOS_tbl.season(ismember(CSOS_tbl.month,[9 10 11]))= "fall";

Ocrit_CSOS = prctile(CSOS_tbl.O2(isfinite(CSOS_tbl.O2)), event_prctile);
CSOS_tbl.lowO_event = CSOS_tbl.O2 <= Ocrit_CSOS;


% ===================== P4 event table =====================

P4_tbl = table;
P4_tbl.station = repmat("P4", numel(P4.mtime), 1);
P4_tbl.mtime   = P4.mtime(:);
P4_tbl.year    = year(datetime(P4.mtime(:),'ConvertFrom','datenum'));
P4_tbl.month   = month(datetime(P4.mtime(:),'ConvertFrom','datenum'));

P4_tbl.PEW     = alpha_PEW_P4_26p5(:);
P4_tbl.O2      = P4.isoDO(target_row,:)';
P4_tbl.NO3     = P4.isoNO3(target_row,:)';
P4_tbl.NO      = P4.isoNO(target_row,:)';
P4_tbl.AOU     = P4.isoAOU(target_row,:)';
P4_tbl.Z       = P4.isoZ(target_row,:)';
P4_tbl.Spice   = P4.isoSpice(target_row,:)';
P4_tbl.Ostar   = P4.Ostar(:);
P4_tbl.Nres    = P4.Nres(:);

P4_tbl.season = strings(height(P4_tbl),1);
P4_tbl.season(ismember(P4_tbl.month,[12 1 2])) = "winter";
P4_tbl.season(ismember(P4_tbl.month,[3 4 5]))  = "spring";
P4_tbl.season(ismember(P4_tbl.month,[6 7 8]))  = "summer";
P4_tbl.season(ismember(P4_tbl.month,[9 10 11]))= "fall";

Ocrit_P4 = prctile(P4_tbl.O2(isfinite(P4_tbl.O2)), event_prctile);
P4_tbl.lowO_event = P4_tbl.O2 <= Ocrit_P4;


% ===================== Combine table =====================

T = [CSOS_tbl; P4_tbl];

% Remove rows missing core variables
T = T(isfinite(T.O2) & isfinite(T.PEW),:);


% ===================== Fit O2 vs PEW source-water relationships =====================

% CSOS fit
idx = T.station=="CSOS" & isfinite(T.PEW) & isfinite(T.O2);
pO_CSOS = polyfit(T.PEW(idx), T.O2(idx), 1);

% P4 fit
idx = T.station=="P4" & isfinite(T.PEW) & isfinite(T.O2);
pO_P4 = polyfit(T.PEW(idx), T.O2(idx), 1);

% Add fitted source-water oxygen and residual anomaly
T.O2_fit_PEW = NaN(height(T),1);

idx = T.station=="CSOS";
T.O2_fit_PEW(idx) = polyval(pO_CSOS, T.PEW(idx));

idx = T.station=="P4";
T.O2_fit_PEW(idx) = polyval(pO_P4, T.PEW(idx));

T.O2_resid_from_PEW = T.O2 - T.O2_fit_PEW;


% ===================== Event attribution =====================
% Total low-O anomaly is relative to each station mean.
% PEW component is fitted O2 from PEW minus station mean O2.
% Residual component is observed minus fitted O2.

T.O2_station_mean = NaN(height(T),1);
T.O2_anom_total   = NaN(height(T),1);
T.O2_anom_PEW     = NaN(height(T),1);
T.O2_anom_resid   = NaN(height(T),1);

idx = T.station=="CSOS";
Omean = mean(T.O2(idx),'omitnan');
T.O2_station_mean(idx) = Omean;
T.O2_anom_total(idx)   = T.O2(idx) - Omean;
T.O2_anom_PEW(idx)     = T.O2_fit_PEW(idx) - Omean;
T.O2_anom_resid(idx)   = T.O2(idx) - T.O2_fit_PEW(idx);

idx = T.station=="P4";
Omean = mean(T.O2(idx),'omitnan');
T.O2_station_mean(idx) = Omean;
T.O2_anom_total(idx)   = T.O2(idx) - Omean;
T.O2_anom_PEW(idx)     = T.O2_fit_PEW(idx) - Omean;
T.O2_anom_resid(idx)   = T.O2(idx) - T.O2_fit_PEW(idx);

% Fractional attribution for low-O events only
T.frac_PEW   = NaN(height(T),1);
T.frac_resid = NaN(height(T),1);

idx = T.lowO_event & T.O2_anom_total < 0;
T.frac_PEW(idx)   = T.O2_anom_PEW(idx)   ./ T.O2_anom_total(idx);
T.frac_resid(idx) = T.O2_anom_resid(idx) ./ T.O2_anom_total(idx);


% ===================== Remineralization diagnostics =====================

% O:N residual ratio. Negative values near -8.625 imply Redfield-like remineralization.
T.ON_resid_ratio = T.Ostar ./ T.Nres;

% Keep only sensible finite ratios for summary; do not overwrite raw values.
T.ON_resid_ratio_QC = T.ON_resid_ratio;
T.ON_resid_ratio_QC(abs(T.ON_resid_ratio_QC) > 50) = NaN;


% ===================== Global summaries =====================

fprintf('\n=== Low-O event thresholds ===\n');
fprintf('CS-OS lowest %.0f%% O2 threshold: %.2f umol/kg\n', event_prctile, Ocrit_CSOS);
fprintf('P4   lowest %.0f%% O2 threshold: %.2f umol/kg\n', event_prctile, Ocrit_P4);

fprintf('\n=== Event attribution by station ===\n');

for st = ["CSOS","P4"]
    idx_all = T.station==st;
    idx_evt = T.station==st & T.lowO_event;

    fprintf('\n%s\n', st);
    fprintf('n all = %d, n events = %d\n', sum(idx_all), sum(idx_evt));
    fprintf('mean all O2     = %.2f\n', mean(T.O2(idx_all),'omitnan'));
    fprintf('mean event O2   = %.2f\n', mean(T.O2(idx_evt),'omitnan'));
    fprintf('mean all PEW    = %.2f\n', mean(T.PEW(idx_all),'omitnan'));
    fprintf('mean event PEW  = %.2f\n', mean(T.PEW(idx_evt),'omitnan'));
    fprintf('mean event O2 anomaly total = %.2f\n', mean(T.O2_anom_total(idx_evt),'omitnan'));
    fprintf('mean event PEW component    = %.2f\n', mean(T.O2_anom_PEW(idx_evt),'omitnan'));
    fprintf('mean event residual comp.   = %.2f\n', mean(T.O2_anom_resid(idx_evt),'omitnan'));
    fprintf('mean event frac PEW         = %.2f\n', mean(T.frac_PEW(idx_evt),'omitnan'));
    fprintf('mean event frac residual    = %.2f\n', mean(T.frac_resid(idx_evt),'omitnan'));
    fprintf('mean event Ostar            = %.2f\n', mean(T.Ostar(idx_evt),'omitnan'));
    fprintf('mean event Nres             = %.2f\n', mean(T.Nres(idx_evt),'omitnan'));
    fprintf('mean event Ostar/Nres       = %.2f\n', mean(T.ON_resid_ratio_QC(idx_evt),'omitnan'));
    fprintf('mean event NO               = %.2f\n', mean(T.NO(idx_evt),'omitnan'));
    fprintf('mean event AOU              = %.2f\n', mean(T.AOU(idx_evt),'omitnan'));
end



% ===================== Figure 1: O2 vs PEW with low-O events =====================

figure('units','centimeters','outerposition',[0.01 0.01 14 7],'color','w');
tiledlayout(1,2,'padding','compact','tilespacing','compact');

for kk = 1:2
    if kk == 1
        st = "CSOS";
        col = rgb_x('rose');
        pO = pO_CSOS;
    else
        st = "P4";
        col = rgb_x('steel blue');
        pO = pO_P4;
    end

    nexttile; hold on;

    idx = T.station==st;
    idx_evt = idx & T.lowO_event;

    scatter(T.PEW(idx), T.O2(idx), 14, 'filled', ...
        'markerfacecolor', col, 'markerfacealpha', 0.25);

    scatter(T.PEW(idx_evt), T.O2(idx_evt), 36, 'filled', ...
        'markerfacecolor', col, 'markeredgecolor', 'k', ...
        'linewidth', 0.5);

    x = linspace(min(T.PEW(idx),[],'omitnan'), max(T.PEW(idx),[],'omitnan'), 200);
    y = polyval(pO,x);
    plot(x,y,'k-','linewidth',1.5);

    xlabel('\it{f_{PEW}}');
    ylabel('O_2 (\mumol kg^{-1})');
    title(st);
    grid on; box on;
end


% ===================== Figure 2: Ostar vs Nres with low-O events =====================

figure('units','centimeters','outerposition',[0.01 0.01 14 7],'color','w');
tiledlayout(1,2,'padding','compact','tilespacing','compact');

m_ref = -8.625;

for kk = 1:2
    if kk == 1
        st = "CSOS";
        col = rgb_x('rose');
    else
        st = "P4";
        col = rgb_x('steel blue');
    end

    nexttile; hold on;

    idx = T.station==st & isfinite(T.Ostar) & isfinite(T.Nres);
    idx_evt = idx & T.lowO_event;

    scatter(T.Nres(idx), T.Ostar(idx), 14, 'filled', ...
        'markerfacecolor', col, 'markerfacealpha', 0.25);

    scatter(T.Nres(idx_evt), T.Ostar(idx_evt), 36, 'filled', ...
        'markerfacecolor', col, 'markeredgecolor', 'k', ...
        'linewidth', 0.5);

    p = polyfit(T.Nres(idx), T.Ostar(idx), 1);
    xx = linspace(min(T.Nres(idx),[],'omitnan'), max(T.Nres(idx),[],'omitnan'), 200);

    plot(xx, polyval(p,xx), '-', 'color', col, 'linewidth', 1.8);

    x0 = mean(T.Nres(idx),'omitnan');
    y0 = mean(T.Ostar(idx),'omitnan');
    plot(xx, m_ref*(xx-x0)+y0, '--k', 'linewidth', 1.2);

    xlabel('N'' (\mumol kg^{-1})');
    ylabel('O'' (\mumol kg^{-1})');
    title(sprintf('%s: O'' vs N''', st));
    grid on; box on;
end


% ===================== Figure 3: Event attribution bars =====================

figure('units','centimeters','outerposition',[0.01 0.01 10 7],'color','w');
hold on;

stn = ["CSOS","P4"];
bar_data = NaN(2,2);

for k = 1:2
    idx_evt = T.station==stn(k) & T.lowO_event;
    bar_data(k,1) = mean(T.O2_anom_PEW(idx_evt),'omitnan');
    bar_data(k,2) = mean(T.O2_anom_resid(idx_evt),'omitnan');
end

bar(bar_data,'stacked');
xticks(1:2);
xticklabels(stn);
ylabel('Mean low-O_2 anomaly contribution (\mumol kg^{-1})');
legend({'PEW/source-water component','Residual O'' component'},'location','best');
grid on; box on;

% 
% % ===================== Figure 4: seasonal low-O event frequency =====================
% figure('units','centimeters','outerposition',[0.01 0.01 12 7],'color','w');
% tiledlayout(1,2,'padding','compact','tilespacing','compact');
% 
% for kk = 1:2
%     st = stn(kk);
%     nexttile; hold on;
% 
%     event_frac = NaN(numel(seasons),1);
% 
%     for s = 1:numel(seasons)
%         idx = T.station==st & T.season==seasons(s);
%         event_frac(s) = sum(T.lowO_event(idx)) ./ sum(idx);
%     end
% 
%     bar(event_frac);
%     xticks(1:numel(seasons));
%     xticklabels(seasons);
%     ylabel('Low-O event fraction');
%     title(st);
%     ylim([0 max(event_frac)*1.2]);
%     grid on; box on;
% end



% ===================== Save event table =====================

save('/Users/samst/Library/CloudStorage/Dropbox/Hakai/gyre/supply/data/lowO_event_attribution.mat', ...
    'T','pO_CSOS','pO_P4','Ocrit_CSOS','Ocrit_P4');

writetable(T, ...
    '/Users/samst/Library/CloudStorage/Dropbox/Hakai/gyre/supply/data/lowO_veent_attribution.csv');




%% Local function

function err = bootstrap_trend_error(x,y,nboot)

x = x(:);
y = y(:);

m = isfinite(x) & isfinite(y);
x = x(m);
y = y(m);

btrend = nan(nboot,1);

for b = 1:nboot

    ii = randi(numel(x),numel(x),1);

    % Extremely rare resample containing only one unique year
    if numel(unique(x(ii))) > 1
        pp = polyfit(x(ii),y(ii),1);
        btrend(b) = pp(1);
    end

end

err = std(btrend,'omitnan');

end

function [lon,lat] = get_oras_lonlat(ncfile)

info = ncinfo(ncfile);
vnames = string({info.Variables.Name});

lon_candidates = [ ...
    "nav_lon", ...
    "longitude", ...
    "lon", ...
    "glamt", ...
    "x"];

lat_candidates = [ ...
    "nav_lat", ...
    "latitude", ...
    "lat", ...
    "gphit", ...
    "y"];

lon_name = "";
lat_name = "";

for k = 1:numel(lon_candidates)
    if any(strcmpi(vnames,lon_candidates(k)))
        lon_name = vnames(find(strcmpi(vnames,lon_candidates(k)),1));
        break
    end
end

for k = 1:numel(lat_candidates)
    if any(strcmpi(vnames,lat_candidates(k)))
        lat_name = vnames(find(strcmpi(vnames,lat_candidates(k)),1));
        break
    end
end

if strlength(lon_name)==0 || strlength(lat_name)==0

    fprintf('\nVariables in %s:\n',ncfile);
    disp(vnames');

    error(['Could not automatically identify longitude/latitude ' ...
        'variables in ORAS5 file.']);

end

lon = squeeze(ncread(ncfile,lon_name));
lat = squeeze(ncread(ncfile,lat_name));

end


%% NO vs temperature, salinity, and fPEW at PEW-core isopycnals

stations = {CSOS,P4};
names = {'CS-OS','P4'};
cols = {rgb_x('rose'),rgb_x('steel blue')};
fractions = {alpha_PEW_CSOS_26p5,alpha_PEW_P4_26p5};
vars = {'isoT','isoSR','fPEW'};
xlabels = {'Potential Temperature (°C)', ...
           'Reference Salinity (g kg^{-1})','f_{PEW}'};

figure('units','centimeters','outerposition',[0.01 0.01 25 9],'color','w');

for j = 1:3
    subplot(1,3,j); hold on;

    for k = 1:2
        y = stations{k}.isoNO(2,:)';

        if j < 3
            x = stations{k}.(vars{j})(2,:)';
        else
            x = fractions{k}(:);
        end

        msk = isfinite(x) & isfinite(y);
        if sum(msk) < 3, continue; end

        [R,P] = corrcoef(x(msk),y(msk));
        coeffs = polyfit(x(msk),y(msk),1);
        xx = linspace(min(x(msk)),max(x(msk)),100);

        scatter(x(msk),y(msk),12,'filled', ...
            'MarkerFaceColor',cols{k},'MarkerFaceAlpha',0.3, ...
            'DisplayName',sprintf('%s (R^2 = %.2f)',names{k},R(1,2)^2));

        plot(xx,polyval(coeffs,xx),'--','Color',cols{k}, ...
            'LineWidth',1.5,'HandleVisibility','off');

        fprintf('%s NO vs %s: N=%d, R2=%.3f, p=%.3g, slope=%.3f\n', ...
            names{k},vars{j},sum(msk),R(1,2)^2,P(1,2),coeffs(1));
    end

    xlabel(xlabels{j});
    ylabel('NO (\mumol kg^{-1})');
    legend('Location','best');
    grid on; box on;
end

set(findall(gcf,'-property','fontname'),'fontname','Latin Modern Roman');
set(findall(gcf,'-property','fontsize'),'fontsize',9);