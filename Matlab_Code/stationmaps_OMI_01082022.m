%% Plot global stations, colored based on what MJO phase exhibits highest wet bulb

addpath /home/ivanov/matlab/Code/m_map

OMI_prep
%% Indices
% Global lat/lon
lat = [-90 90];
lon = [-180 180];

% Regional lat/lon
rlat = [0 36];
rlon = [30 130];

latlon = readmatrix('finalstnlatlons.txt','Delimiter',',');

% Identify regional stations that pass test
regional = latlon(latlon(:,1) > rlat(1) & latlon(:,1) < rlat(2) & latlon(:,2) > rlon(1) & latlon(:,2) < rlon(2),:);

% Find station indices  
reg_lats = ismember(stnlats, regional(:,1));
reg_lons = ismember(stnlons, regional(:,2));
reg_ind = find(reg_lats + reg_lons == 2);

%% Plot seasonal maps

% Plot MJ 

season = 1;

% Initialize array, MJ
OMI_wbt_MJJA_no0 = squeeze(OMI_wbt_MJJA(2:9,season,:));
OMI_highest_phase = 0*stnelevs;

% Find highest WBT
for ii = 1:length(stnelevs)
    for kk = 1:8
        stn_data = OMI_wbt_MJJA_no0(:,ii);
      
        if OMI_wbt_MJJA_no0(kk,ii) == max(stn_data)
            OMI_highest_phase(ii) = kk; 
        end
    end
end

% Plot stations by highest phase

CME_lon = [45 60 60 45 45];
CME_lat = [20 20 36 36 20];

NW_lon = [68 78 78 68 68];
NW_lat = [22 22 32 32 22];

SE_lon = [78 90 90 78 78];
SE_lat = [8 8 22 22 8];

addpath /home/ivanov/matlab/Code/cbrewer/cbrewer/

sz = 55;      

c = OMI_highest_phase(reg_ind);
fig = figure('Position', [10 10 1200 700]);
set(gca,'DefaultTextFontSize',10)
subplot(2,1,1)

set(0,'DefaultAxesColorOrder',cbrewer('qual','Set2',8))
set(0,'DefaultLineLineWidth',1.2)
m_proj('Equidistant Cylindrical','lon',rlon,'lat',rlat);
m_coast('patch',[.85 .85 .85],'edgecolor','black');
m_grid('linestyle','none','tickdir','out');
hold on
m_scatter(stnlons(reg_ind),stnlats(reg_ind),sz, c, 'filled')
hold on

% Add subregion boxes
m_line(CME_lon,CME_lat,'linewi',2,'color',[0.6 0.6 0.6]);
m_line(NW_lon,NW_lat,'linewi',2,'color',[0.6 0.6 0.6]);
m_line(SE_lon,SE_lat,'linewi',2,'color',[0.6 0.6 0.6]);


% Plot JA
season = 2;

% Initialize array, MJ
OMI_wbt_MJJA_no0 = squeeze(OMI_wbt_MJJA(2:9,season,:));
OMI_highest_phase = 0*stnelevs;

% Find highest WBT
for ii = 1:length(stnelevs)
    for kk = 1:8
        stn_data = OMI_wbt_MJJA_no0(:,ii);
      
        if OMI_wbt_MJJA_no0(kk,ii) == max(stn_data)
            OMI_highest_phase(ii) = kk; 
        end
    end
end

% Plot stations by highest phase

CME_lon = [45 60 60 45 45];
CME_lat = [20 20 36 36 20];

NW_lon = [68 78 78 68 68];
NW_lat = [22 22 32 32 22];

SE_lon = [78 90 90 78 78];
SE_lat = [8 8 22 22 8];

addpath /home/ivanov/matlab/Code/cbrewer/cbrewer/

sz = 55;      

c = OMI_highest_phase(reg_ind);

s2 = subplot(2,1,2);

set(0,'DefaultAxesColorOrder',cbrewer('qual','Set2',8))
set(0,'DefaultLineLineWidth',1.2)  
m_proj('Equidistant Cylindrical','lon',rlon,'lat',rlat);
m_coast('patch',[.85 .85 .85],'edgecolor','black');
m_grid('linestyle','none','tickdir','out');
hold on
m_scatter(stnlons(reg_ind),stnlats(reg_ind),sz, c, 'filled')

pos = get(s2, 'Position');
posnew = pos;
posnew(2) = posnew(2) + 0.06;
set(s2, 'Position', posnew);

hold on

% Add subregion boxes
m_line(CME_lon,CME_lat,'linewi',2,'color',[0.6 0.6 0.6]);
m_line(NW_lon,NW_lat,'linewi',2,'color',[0.6 0.6 0.6]);
m_line(SE_lon,SE_lat,'linewi',2,'color',[0.6 0.6 0.6]);

h = axes(fig,'visible','off'); 
c = colorbar(h,'southoutside','position',[0.37 0.1 0.3 0.02]);
colormap(cbrewer('div','Spectral',8))
caxis(h,[1,9]);
set(c, 'Ticks', [1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5], 'FontSize',10);
set(c, 'TickLabels', {'1', '2', '3','4','5','6','7','8',' '});
ylabel(c, 'Phase', 'Fontsize',10)