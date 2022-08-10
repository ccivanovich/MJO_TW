%% Plot global stations, colored based on what MJO phase exhibits highest wet bulb

addpath /home/ivanov/matlab/Code/m_map

addpath /home/ivanov/matlab/Code/MJO_TW/Manuscript_prep/MJJA

OMI_prep

% Global lat/lon
lat = [-90 90];
lon = [-180 180];

% Regional lat/lon
rlat = [0 36];
rlon = [30 130];

%% MJJA

% Initialize array
OMI_wbt_MJJA_no0 = OMI_wbt_MJJA(2:9,:);
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

%% Plot stations by highest phase

addpath /home/ivanov/matlab/Code/cbrewer/cbrewer/

% Select only stations that pass threshold criterion
goodstns = readmatrix('stnreqts_50pct_mjjas_Fig1box.csv');

goodlats = stnlats(goodstns);
goodlons = stnlons(goodstns);
c = OMI_highest_phase(goodstns);

sz = 55;      

figure('Position', [10 10 1200 700])
set(0,'DefaultAxesColorOrder',cbrewer('qual','Set2',8))
set(0,'DefaultLineLineWidth',1.2)
set(0,'DefaultFigureColormap',cbrewer('seq','YlGnBu',8));    
m_proj('Equidistant Cylindrical','lon',rlon,'lat',rlat);
m_coast('patch',[.85 .85 .85],'edgecolor','black');
m_grid('linestyle','none','tickdir','out');
hold on
m_scatter(goodlons,goodlats,sz, c, 'filled')
colormap(cbrewer('seq','YlGnBu',8));
hCbar = colorbar('eastoutside');
set(hCbar, 'TickLabels', {'1', '2', '3','4','5','6','7','8'});
ylabel(hCbar, 'Phase')

%% MJ

% Initialize array
OMI_wbt_MJJA_no0 = OMI_wbt_MJJA(2:9,1,:);
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

%% Plot stations by highest phase

addpath /home/ivanov/matlab/Code/cbrewer/cbrewer/

sz = 55;      

c = OMI_highest_phase;
figure('Position', [10 10 1200 700])
set(0,'DefaultAxesColorOrder',cbrewer('qual','Set2',8))
set(0,'DefaultLineLineWidth',1.2)
set(0,'DefaultFigureColormap',cbrewer('seq','YlGnBu',8));    
m_proj('Equidistant Cylindrical','lon',rlon,'lat',rlat);
m_coast('patch',[.85 .85 .85],'edgecolor','black');
m_grid('linestyle','none','tickdir','out');
hold on
m_scatter(stnlons,stnlats,sz, c, 'filled')
colormap(cbrewer('seq','YlGnBu',8));
hCbar = colorbar('eastoutside');
set(hCbar, 'TickLabels', {'1', '2', '3','4','5','6','7','8'});
ylabel(hCbar, 'Phase')

%% JA

% Initialize array
OMI_wbt_MJJA_no0 = OMI_wbt_MJJA(2:9,1,:);
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

%% Plot stations by highest phase

addpath /home/ivanov/matlab/Code/cbrewer/cbrewer/

sz = 55;      

c = OMI_highest_phase;
figure('Position', [10 10 1200 700])
set(0,'DefaultAxesColorOrder',cbrewer('qual','Set2',8))
set(0,'DefaultLineLineWidth',1.2)
set(0,'DefaultFigureColormap',cbrewer('seq','YlGnBu',8));    
m_proj('Equidistant Cylindrical','lon',rlon,'lat',rlat);
m_coast('patch',[.85 .85 .85],'edgecolor','black');
m_grid('linestyle','none','tickdir','out');
hold on
m_scatter(stnlons,stnlats,sz, c, 'filled')
colormap(cbrewer('seq','YlGnBu',8));
hCbar = colorbar('eastoutside');
set(hCbar, 'TickLabels', {'1', '2', '3','4','5','6','7','8'});
ylabel(hCbar, 'Phase')