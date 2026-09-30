clear
load CSOS_ctd.mat 
load P4_ctd.mat
gO=load('globalOxI.mat');

%% maps
col=[255 214 140]/255; % YELLOW!
lat_lim=[50 53]; lon_lim=[-132 -127];
fname='/Users/samst/Dropbox/UBC/Misc/british_columbia_3_msl_2013.nc';
lat=ncread(fname,'lat');    
lon=ncread(fname,'lon');
ilon=lon>=lon_lim(1) & lon<=lon_lim(2);
ilat=lat>=lat_lim(1) & lat<=lat_lim(2);
Z=ncread(fname,...
    'Band1',[ find((ilon),1,'first') find((ilat),1,'first')],...
    [ sum(ilon) sum(ilat)],[1 1]);

figure('units','centimeters','outerposition',[0 0 18 21],'color','w');

% ax1=axes('position',[0.525 0.5-0.05 0.425 0.425]);
% t = tiledlayout(2,1,'TileSpacing', 'Compact', 'Padding', 'Compact');

%%%%%%%
axGlob=axes('position',[0.1 0.55 0.75 0.4], 'clipping' , 'on');
% axGlob=nexttile;
% Add small inset map highlighting Queen Charlotte Sound
m_proj('equidistant', 'lon', [120 260], 'lat', [0 60]); % North Pacific bounds
% m_proj('ortho', 'lon', -179, 'lat', 20); 
hold on;

[gO.lat_grid,gO.lon_grid]=meshgrid(gO.lat,gO.lon);
lon_grid_corrected= mod(gO.lon_grid + 360, 360);
cmap = cmocean('tempo'); cmap(1:30,:)=[];
cmap = [[0.5 0.5 0.5]; cmap];  % Add grey for NaNs at the start
colormap(gca, cmap);

% Set NaN values to a specific value outside the contour range
data = gO.oxI(:,:,1); 
msk=gO.lon_grid<-100 | gO.lon_grid>120;
data(~msk)=NaN;
data(isnan(data)) = -1;  % Replace NaNs with a value outside the plotted range

% Plot the data with adjusted NaN handling
[CS, CH] = m_contourf(lon_grid_corrected, gO.lat_grid(:,:,1), data, [0:30:360 -1], 'linestyle', 'none');

% Set the color axis to include the NaN color
caxis([-2 360]);  % Ensure the colormap includes the NaN value (-1)
% 
m_contour(lon_grid_corrected,gO.lat_grid(:,:,1),gO.geoI,-20:1.25:0,'linecolor',[0.3 0.3 0.3]);

m_coast('patch', col, 'edgecolor', 'none'); % Light gray continents

lbls={'BC','CA','OR','WA'};

count=0;
Mp = m_shaperead('ne_50m_admin_1_states_provinces'); 
for k = [38 54 87 97]
    lonC = mod(Mp.ncst{k}(:,1) + 360, 360);  % Convert longitudes to 0-360
    latC = Mp.ncst{k}(:,2);

    % Find NaN indices to break the segments
    nanIdx = isnan(lonC) | isnan(latC);
    segmentStart = [1; find(nanIdx) + 1];  % Start of each segment
    segmentEnd = [find(nanIdx) - 1; length(lonC)];  % End of each segment

    % Remove invalid indices (where start > end)
    validSeg = segmentStart <= segmentEnd;
    segmentStart = segmentStart(validSeg);
    segmentEnd = segmentEnd(validSeg);

    % Loop through segments and plot them separately
    for j = 1:length(segmentStart)
        m_plot(lonC(segmentStart(j):segmentEnd(j)), ...
             latC(segmentStart(j):segmentEnd(j)), 'color',[0.6 0 0], 'LineWidth', 0.1);
    end
    count=count+1;
    m_text(mean(lonC,'omitnan'),mean(latC,'omitnan'), lbls{count});
end

Mp = m_shaperead('ne_50m_admin_0_countries'); 
for k = [17 56 76 90 203]
    lonC = mod(Mp.ncst{k}(:,1) + 360, 360);  % Convert longitudes to 0-360
    latC = Mp.ncst{k}(:,2);

    % Find NaN indices to break the segments
    nanIdx = isnan(lonC) | isnan(latC);
    segmentStart = [1; find(nanIdx) + 1];  % Start of each segment
    segmentEnd = [find(nanIdx) - 1; length(lonC)];  % End of each segment

    % Remove invalid indices (where start > end)
    validSeg = segmentStart <= segmentEnd;
    segmentStart = segmentStart(validSeg);
    segmentEnd = segmentEnd(validSeg);

    % Loop through segments and plot them separately
    for j = 1:length(segmentStart)
        m_plot(lonC(segmentStart(j):segmentEnd(j)), ...
             latC(segmentStart(j):segmentEnd(j)), 'color',[0.6 0.6 0.6], 'LineWidth', 0.2);
    end
end


% Extents
BClon_lim = [-132 -124];
BClat_lim = [48   52];

% Highlight Queen Charlotte Sound region
m_patch(mod(BClon_lim([1 2 2 1 1]) + 360, 360), BClat_lim([1 1 2 2 1]), 'r', ...
    'facealpha', 0, 'edgecolor', 'k','linewidth',2); % Semi-transparent red box

m_text(247.97,53.24,'Canada','fontweight','bold')
m_text(250.51,42.28,'U.S.A.','fontweight','bold')

% =========================================================
% Source-water / ventilation boxes
% =========================================================

% ---- PSUW source region (10° x 10° box centered at 50°N, 145°W)
PSUW_lon_lim = [-150 -140];
PSUW_lat_lim = [45 55];

% ---- PEW source region (10° x 10° box centered at 5°N, 108°W)
PEW_lon_lim  = [-113 -103];
PEW_lat_lim  = [1 10];

% ---- Subarctic ventilation zone
% Fill these in with the exact bounds you want
VENT_lon_lim = [-140 140];   % e.g. [___ ___]
VENT_lat_lim = [30 59];   % e.g. [___ ___]

% Convert all longitudes to 0–360
PSUW_lon_lim_0360 = mod(PSUW_lon_lim + 360, 360);
PEW_lon_lim_0360  = mod(PEW_lon_lim  + 360, 360);

% ---- Plot PSUW box
m_patch(PSUW_lon_lim_0360([1 2 2 1 1]), ...
        PSUW_lat_lim([1 1 2 2 1]), ...
        'b', 'facealpha', 0, 'edgecolor', rgb_x('greyish purple'), 'linewidth', 2);

% ---- Plot PEW box
m_patch(PEW_lon_lim_0360([1 2 2 1 1]), ...
        PEW_lat_lim([1 1 2 2 1]), ...
        'r', 'facealpha', 0, 'edgecolor',rgb_x('burnt orange'), 'linewidth', 2);
% 
% % ---- Plot ventilation box only if bounds are filled
% VENT_lon_lim_0360 = mod(VENT_lon_lim + 360, 360);
% 
% m_patch(VENT_lon_lim_0360([1 2 2 1 1]), ...
%     VENT_lat_lim([1 1 2 2 1]), ...
%     [0 0.5 0], 'facealpha', 0, ...
%     'edgecolor', [0.2 0.2 0.2], 'linewidth', 2, 'linestyle', '--');

% Add grid lines for the orthographic map
m_grid('xtick',4 , 'ytick', 5, 'linestyle', 'none'); % Hide grid lines
text(0.05, 0.95,'a)','fontsize',9,'units','normalized');
%
ax=m_contfbar(1.05,[.2 .8],CS,CH);
ylabel(ax,'O$_2$ ($\mu$mol kg$^{-1}$)','Interpreter','latex');


% Station lon/lat (re-use existing arrays; avoid new vars)
P4_lon = mode(P4.lon);  P4_lat = mode(P4.lat);
msk = strcmp('IOS_CS00',CSOS.station); CS00_lon = mode(CSOS.lon(msk)); CS00_lat = mode(CSOS.lat(msk));
msk = strcmp('IOS_CS01',CSOS.station); CS01_lon = mode(CSOS.lon(msk)); CS01_lat = mode(CSOS.lat(msk));
msk = strcmp('IOS_CS02',CSOS.station); CS02_lon = mode(CSOS.lon(msk)); CS02_lat = mode(CSOS.lat(msk));

% Land color to match your style
col_land = [255 214 140]/255;

% ===== Main map =====
ax1 = axes('position',[0.25 0.15 0.5 0.40]); %#ok<LAXES>
m_proj('equidistant','lon',BClon_lim,'lat',BClat_lim);

% Bathy + coast (same workflow you used elsewhere)
[CS,CH]=m_etopo2('contourf',-4000:250:0,'edgecolor','none'); 
caxis([-4000 0]); colormap(ax1,m_colmap('blues'));
m_gshhs_i('patch',col_land,'edgecolor',[0.3 0.3 0.3],'linewidth',0.3);

% Grid (call BEFORE m_ruler per docs)
m_grid('tickdir','out','linestyle','none','box','on');

% Stations
m_scatter([P4_lon],[P4_lat], ...
          50,'filled','markerfacecolor',rgb_x('steel blue'), 'markeredgecolor','k','linew',0.8);
m_scatter([CS00_lon CS01_lon CS02_lon], ...
          [CS00_lat CS01_lat CS02_lat], ...
          50,'filled','markerfacecolor',rgb_x('rose'),'markeredgecolor','k','linew',0.8);
m_text(P4_lon+0.02,  P4_lat , '  P4'    ,'color','k');
m_text(CS00_lon+0.02,CS00_lat-0.01,'  CS00'  ,'color','k');
m_text(CS01_lon+0.02,CS01_lat,'  CS01'  ,'color','k');
m_text(CS02_lon+0.02,CS02_lat+0.01,'  CS02'  ,'color','k');

m_ruler([0.05 0.3], 0.9,3);

ax=m_contfbar(1.05,[.2 .8],CS,CH,'endpiece','no');
ylabel(ax,{'Depth (m)'}); % Move up by inserting a blank line

set(findall(gcf, '-property', 'fontname'), 'fontname', 'Latin Modern Roman');
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 10);

%%
export_fig /Users/samst/Library/CloudStorage/Dropbox/Hakai/gyre/supply/figs/map_v2.jpg -m5


