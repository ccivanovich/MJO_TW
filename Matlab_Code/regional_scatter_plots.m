%% Reliable Coastal ME
% Load data
addpath /home/ivanov/matlab/Code/MJO_TW

OMI_prep_01062022
OMI_stats_absthresh_01062022
%% Regional Indices
fid = fopen('persiangulfstns.txt','rt');
P = textscan(fid, '%f%f');
fclose(fid);
PG_locations = cell2mat(P);
phases = 0:1:8;
j = 3; % set TW threshold

% Identify coastal middle east stations
CME_lon = [45 60];
CME_lat = [20 36];
coastal = PG_locations(PG_locations(:,1) > CME_lat(1) & PG_locations(:,1) < CME_lat(2) & PG_locations(:,2) > CME_lon(1) & PG_locations(:,2) < CME_lon(2)-6,:);

% Find station indices  
coast_lats = ismember(stnlats, coastal(:,1));
coast_lons = ismember(stnlons, coastal(:,2));
coastal_ind = find(coast_lats + coast_lons == 2);

% NEW: Only use three top stations
seasonal_exceed = squeeze(sum(OMI_thresh_1std(:,j,2,:),1));
coastal_exceed = squeeze(sum(OMI_thresh_1std(:,j,2,coastal_ind),1));
[sorted, sortedInd] = sort(coastal_exceed,'descend');
top3 = sorted(1);
coastal_3 = find(seasonal_exceed >= top3 & stnlats > CME_lat(1) & stnlats < CME_lat(2) & stnlons > CME_lon(1) & stnlons < CME_lon(2) & coast_lats + coast_lons == 2);

% Convert text file to matrix of lat/lon
fid = fopen('southasiastns.txt','rt');
S = textscan(fid, '%f%f');
fclose(fid);
SA_locations = cell2mat(S);

% Northwestern South Asia
NW_lon = [68 78];
NW_lat = [22 32];

% Find station indices  
NW_lats = ismember(stnlats, SA_locations(:,1));
NW_lons = ismember(stnlons, SA_locations(:,2));
NW_ind = find(NW_lats + NW_lons == 2 & stnlats > NW_lat(1) & stnlats < NW_lat(2) & stnlons > NW_lon(1) & stnlons < NW_lon(2));

% NEW: Only use three top stations
seasonal_exceed = squeeze(sum(OMI_thresh_1std(:,j,2,:),1));
NW_exceed = squeeze(sum(OMI_thresh_1std(:,j,2,NW_ind),1));
[sorted, sortedInd] = sort(NW_exceed,'descend');
top3 = sorted(1);
NW_3 = find(seasonal_exceed >= top3 & stnlats > NW_lat(1) & stnlats < NW_lat(2) & stnlons > NW_lon(1) & stnlons < NW_lon(2) & NW_lats + NW_lons == 2);

% Southeastern India
SE_lon = [78 90];
SE_lat = [8 22];

% Find station indices  
SE_lats = ismember(stnlats, SA_locations(:,1));
SE_lons = ismember(stnlons, SA_locations(:,2));
SE_ind = find(SE_lats + SE_lons == 2 & stnlats > SE_lat(1) & stnlats < SE_lat(2) & stnlons > SE_lon(1) & stnlons < SE_lon(2));

% NEW: Only use three top stations
seasonal_exceed = squeeze(sum(OMI_thresh_1std(:,j,1,:),1));
SE_exceed = squeeze(sum(OMI_thresh_1std(:,j,1,SE_ind),1));
[sorted, sortedInd] = sort(SE_exceed,'descend');
top3 = sorted(1);
SE_3 = find(seasonal_exceed >= top3 & stnlats > SE_lat(1) & stnlats < SE_lat(2) & stnlons > SE_lon(1) & stnlons < SE_lon(2) & SE_lats + SE_lons == 2);

% NEW: Find matching indices for scrambled
stns = [coastal_ind; NW_ind; SE_ind];
[x, coastal_yes] = ismember(stns,coastal_3);
[x, NW_yes] = ismember(stns,NW_3);
[x, SE_yes] = ismember(stns,SE_3);

coastal_3scram = find(coastal_yes ~= 0);
NW_3scram = find(NW_yes ~= 0);
SE_3scram = find(SE_yes ~= 0);

%% Plotting
    
j = 5;
figure('Position', [10 10 800 900])
set(0,'defaultAxesFontSize',10)
thresh = 27+j;

% Calculate station exceedances
CME_num = sum(OMI_thresh_1std(:,j,2,coastal_3));

% Select PG data and calculate risk
CME_risk_alltemps = sum(OMI_thresh_1std(:,:,2,coastal_3),3)./sum(OMI_obs_len(:,2,coastal_3),2)*100;
scrambled_risk_CME = squeeze(sum(OMI_scrambled_thresh(:,:,2,coastal_3scram,:),3))./sum(OMI_MJOdays_scrambled_thresh(:,2,coastal_3scram,:),2)*100;

% Calculate 95th and 5th percentiles
scram_95_CME = mean(prctile(squeeze(scrambled_risk_CME(:,j,:)),95,2));
scram_5_CME = mean(prctile(squeeze(scrambled_risk_CME(:,j,:)),5,2));

%Create background shading
subplot(3,2,1)
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
a = scatter(phases,CME_risk_alltemps(:,j),1000,'k','.');
xlabel('Phase')
ylabel({'Likelihood Exceeds ',[num2str(thresh),' °C Wet Bulb (%)']})
xticks(phases)

uistack(a,'top')
hold off

% Plot station locations
addpath /home/ivanov/matlab/Code/m_map/

subplot(3,2,2);

m_proj('Equidistant Cylindrical','lon',CME_lon,'lat',CME_lat);
m_coast('patch',[.85 .85 .85],'edgecolor','black');
m_grid('xtick', 3, 'linestyle','none','tickdir','out');
hold on

m_scatter(stnlons(coastal_3),stnlats(coastal_3),100,'k','o','filled')
hold on

m_text(CME_lon(1)+0.5, CME_lat(2)-1.5,num2str(CME_num))
hold off

% Calculate station exceedances
NW_num = sum(OMI_thresh_1std(:,j,2,NW_3));

% Select SA data and calculate risk
NW_risk_alltemps = sum(OMI_thresh_1std(:,:,2,NW_3),3)./sum(OMI_obs_len(:,2,NW_3),2)*100;
scrambled_risk_NW = squeeze(sum(OMI_scrambled_thresh(:,:,2,NW_3scram+1,:),3))./sum(OMI_MJOdays_scrambled_thresh(:,2,NW_3scram+1,:),2)*100;

% Calculate 95th and 5th percentiles
scram_95_NW = mean(prctile(squeeze(scrambled_risk_NW(:,j,:)),95,2));
scram_5_NW = mean(prctile(squeeze(scrambled_risk_NW(:,j,:)),5,2));

%Create background shading
subplot(3,2,3)
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
a = scatter(phases,NW_risk_alltemps(:,j),1000,'k','.');
xlabel('Phase')
ylabel({'Likelihood Exceeds ',[num2str(thresh),' °C Wet Bulb (%)']})
xticks(phases)

uistack(a,'top')
hold off

% Plot station locations colored by whether they ever exceed and #
% of exceedences
subplot(3,2,4)
addpath /home/ivanov/matlab/Code/m_map/


m_proj('Equidistant Cylindrical','lon',NW_lon,'lat',NW_lat);
m_coast('patch',[.85 .85 .85],'edgecolor','black');
m_grid('xtick', 3, 'linestyle','none','tickdir','out');
hold on

m_scatter(stnlons(NW_3),stnlats(NW_3),100,'k','o','filled')
hold on

m_text(NW_lon(1)+0.5, NW_lat(2)-1,num2str(NW_num))
hold off

% Calculate station exceedances
SE_num = sum(OMI_thresh_1std(:,j,1,SE_3));

% Select SA data and calculate risk
SE_risk_alltemps = sum(OMI_thresh_1std(:,:,1,SE_3),3)./sum(OMI_obs_len(:,1,SE_3),2)*100;
scrambled_risk_SE = squeeze(sum(OMI_scrambled_thresh(:,:,1,SE_3scram +1,:),3))./sum(OMI_MJOdays_scrambled_thresh(:,1,SE_3scram+1,:),2)*100;

% Calculate 95th and 5th percentiles
scram_95_SE = mean(prctile(squeeze(scrambled_risk_SE(:,j,:)),95,2));
scram_5_SE = mean(prctile(squeeze(scrambled_risk_SE(:,j,:)),5,2));

%Create background shading
subplot(3,2,5)
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
a = scatter(phases,SE_risk_alltemps(:,j),1000,'k','.');
xlabel('Phase')
ylabel({'Likelihood Exceeds ',[num2str(thresh),' °C Wet Bulb (%)']})
xticks(phases)

uistack(a,'top')
hold off

% Plot station locations colored by whether they ever exceed and #
% of exceedences
subplot(3,2,6)
addpath /home/ivanov/matlab/Code/m_map/

m_proj('Equidistant Cylindrical','lon',SE_lon,'lat',SE_lat);
m_coast('patch',[.85 .85 .85],'edgecolor','black');
m_grid('xtick', 3, 'linestyle','none','tickdir','out');
hold on

m_scatter(stnlons(SE_3),stnlats(SE_3),100,'k','o','filled')
hold on

m_text(SE_lon(1)+0.5, SE_lat(2)-1,num2str(SE_num))
hold off