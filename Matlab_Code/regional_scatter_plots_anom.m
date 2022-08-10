%% Reliable Coastal ME
% Load data
addpath /home/ivanov/matlab/Code/MJO_TW/

OMI_TWanom_01112022
OMI_stats_anom_01112022

%%
fid = fopen('persiangulfstns.txt','rt');
P = textscan(fid, '%f%f');
fclose(fid);
PG_locations = cell2mat(P);
phases = 0:1:8;
season = ["DJF", "MAM", "JJA", "SON"];

% Identify coastal middle east stations
CME_lon = [45 60];
CME_lat = [20 36];
coastal = PG_locations(PG_locations(:,1) > CME_lat(1) & PG_locations(:,1) < CME_lat(2) & PG_locations(:,2) > CME_lon(1) & PG_locations(:,2) < CME_lon(2)-6,:);

% Find station indices  
coast_lats = ismember(stnlats, coastal(:,1));
coast_lons = ismember(stnlons, coastal(:,2));
coastal_ind = find(coast_lats + coast_lons == 2);

% Convert text file to matrix of lat/lon
fid = fopen('southasiastns.txt','rt');
S = textscan(fid, '%f%f');
fclose(fid);
SA_locations = cell2mat(S);

season = 2;

% Select CME anomaly data
CME_anom = mean(squeeze(OMI_wbt(:,season,coastal_ind)),2);

% Calculate scrambled 95th and 5th percentiles
scram_CME = squeeze(OMI_scrambled(:,season,1:length(coastal_ind),:));
scram_95_CME = mean(prctile(squeeze(mean(scram_CME,2)),95,2));
scram_5_CME = mean(prctile(squeeze(mean(scram_CME,2)),5,2));

%Create background shading
subplot(3,1,1)
v = [0 scram_5_CME; 8 scram_5_CME; 8 scram_95_CME; 0 scram_95_CME];
f = [1 2 3 4];
grey = [0.95 0.95 0.95];
patch('Faces', f, 'Vertices', v, 'FaceColor',grey, 'EdgeColor', 'none')
hold on

%Indicate 95th and 5th percentiles
line([0,8],[scram_5_CME,scram_5_CME],'Color','blue')
hold on

line([0,8],[scram_95_CME,scram_95_CME],'Color','red')
hold on

%Plot scatter
a = scatter(phases,CME_anom,1000,'k','.');
xlabel('Phase')
ylabel({'Wet Bulb Anomaly (C)'})
xticks(phases)

uistack(a,'top')
hold off

season = 2;

% Nothern South Asia
NW_lon = [68 78];
NW_lat = [22 32];

% Find station indices  
NW_lats = ismember(stnlats, SA_locations(:,1));
NW_lons = ismember(stnlons, SA_locations(:,2));
NW_ind = find(NW_lats + NW_lons == 2 & stnlats > NW_lat(1) & stnlats < NW_lat(2) & stnlons > NW_lon(1) & stnlons < NW_lon(2));

% Select NW anomaly data
NW_anom = mean(squeeze(OMI_wbt(:,season,NW_ind)),2);

% Calculate scrambled 95th and 5th percentiles
scram_NW = squeeze(OMI_scrambled(:,season,length(coastal_ind)+1:length(coastal_ind)+ length(NW_ind),:));
scram_95_NW = mean(prctile(squeeze(mean(scram_NW,2)),95,2));
scram_5_NW = mean(prctile(squeeze(mean(scram_NW,2)),5,2));

%Create background shading
subplot(3,1,2)
v = [0 scram_5_NW; 8 scram_5_NW; 8 scram_95_NW; 0 scram_95_NW];
f = [1 2 3 4];
grey = [0.95 0.95 0.95];
patch('Faces', f, 'Vertices', v, 'FaceColor',grey, 'EdgeColor', 'none')
hold on

%Indicate 95th and 5th percentiles
line([0,8],[scram_5_NW,scram_5_NW],'Color','blue')
hold on

line([0,8],[scram_95_NW,scram_95_NW],'Color','red')
hold on

%Plot scatter
a = scatter(phases,NW_anom,1000,'k','.');
xlabel('Phase')
ylabel({'Wet Bulb Anomaly (C)'})
xticks(phases)

uistack(a,'top')
hold off

season = 1;

% SE South Asia
SE_lon = [78 90];
SE_lat = [8 22];

% Find station indices  
SE_lats = ismember(stnlats, SA_locations(:,1));
SE_lons = ismember(stnlons, SA_locations(:,2));
SE_ind = find(SE_lats + SE_lons == 2 & stnlats > SE_lat(1) & stnlats < SE_lat(2) & stnlons > SE_lon(1) & stnlons < SE_lon(2));

% Select NW anomaly data
SE_anom = nanmean(squeeze(OMI_wbt(:,season,SE_ind)),2);

% Calculate scrambled 95th and 5th percentiles
scram_SE = squeeze(OMI_scrambled(:,season,length(coastal_ind)+ length(NW_ind)+1:end,:));
scram_95_SE = nanmean(prctile(squeeze(nanmean(scram_SE,2)),95,2));
scram_5_SE = nanmean(prctile(squeeze(nanmean(scram_SE,2)),5,2));

%Create background shading
subplot(3,1,3)
v = [0 scram_5_SE; 8 scram_5_SE; 8 scram_95_SE; 0 scram_95_SE];
f = [1 2 3 4];
grey = [0.95 0.95 0.95];
patch('Faces', f, 'Vertices', v, 'FaceColor',grey, 'EdgeColor', 'none')
hold on

%Indicate 95th and 5th percentiles
line([0,8],[scram_5_SE,scram_5_SE],'Color','blue')
hold on

line([0,8],[scram_95_SE,scram_95_SE],'Color','red')
hold on

%Plot scatter
a = scatter(phases,SE_anom,1000,'k','.');
xlabel('Phase')
ylabel({'Wet Bulb Anomaly (C)'})
xticks(phases)

uistack(a,'top')
hold off