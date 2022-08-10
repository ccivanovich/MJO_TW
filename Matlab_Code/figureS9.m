%% Reliable Coastal ME
% Load data
addpath /home/ivanov/matlab/Code/MJO_TW

OMI_prep
OMI_stats_absthresh
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

% Southeastern India
SE_lon = [78 90];
SE_lat = [8 22];

% Find station indices  
SE_lats = ismember(stnlats, SA_locations(:,1));
SE_lons = ismember(stnlons, SA_locations(:,2));
SE_ind = find(SE_lats + SE_lons == 2 & stnlats > SE_lat(1) & stnlats < SE_lat(2) & stnlons > SE_lon(1) & stnlons < SE_lon(2));

%% Plotting
    
j = 1;
figure('Position', [10 10 800 900])
set(0,'defaultAxesFontSize',10)
thresh = 27+j;

% Select PG data and calculate risk
CME_risk_alltemps = sum(OMI_thresh_1std(:,:,2,coastal_ind),4)./sum(OMI_obs_len(:,2,coastal_ind),3)*100;
scrambled_risk_CME = squeeze(sum(OMI_scrambled_thresh(:,:,2,1:length(coastal_ind),:),4))./sum(OMI_MJOdays_scrambled_thresh(:,2,1:length(coastal_ind),:),3)*100;

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

JA_exceed = squeeze(sum(OMI_thresh_1std(:,j,2,:),1));
CME_exceed_ind = find(JA_exceed ~= 0 & stnlats > CME_lat(1) & stnlats < CME_lat(2) & stnlons > CME_lon(1) & stnlons < CME_lon(2) & coast_lats + coast_lons == 2);
CME_0_ind = find(JA_exceed == 0 & stnlats > CME_lat(1) & stnlats < CME_lat(2) & stnlons > CME_lon(1) & stnlons < CME_lon(2) & coast_lats + coast_lons == 2);

JA_colors = squeeze(sum(OMI_thresh_1std(:,j,2,:),1));
c = JA_colors(CME_exceed_ind);

subplot(3,2,2);

m_proj('Equidistant Cylindrical','lon',CME_lon,'lat',CME_lat);
m_coast('patch',[.85 .85 .85],'edgecolor','black');
m_grid('xtick', 3, 'linestyle','none','tickdir','out');
hold on

m_scatter(stnlons(CME_0_ind),stnlats(CME_0_ind),100,'k','o')
hold on

m_scatter(stnlons(CME_exceed_ind),stnlats(CME_exceed_ind),100,c,'o','filled')

colormap(flipud(hot))
h = colorbar;
ylabel(h, 'Total Exceedances')
hold off


% Select NW data and calculate risk
NW_risk_alltemps = sum(OMI_thresh_1std(:,:,2,NW_ind),4)./sum(OMI_obs_len(:,2,NW_ind),3)*100;
scrambled_risk_NW = squeeze(sum(OMI_scrambled_thresh(:,:,2,length(coastal_ind)+1:length(coastal_ind) + length(NW_ind),:),4))./sum(OMI_MJOdays_scrambled_thresh(:,2,length(coastal_ind)+1:length(coastal_ind) + length(NW_ind),:),3)*100;

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

% Plot station locations
addpath /home/ivanov/matlab/Code/m_map/

NW_exceed_ind = find(JA_exceed ~= 0 & stnlats > NW_lat(1) & stnlats < NW_lat(2) & stnlons > NW_lon(1) & stnlons < NW_lon(2) & NW_lats + NW_lons == 2);
NW_0_ind = find(JA_exceed == 0 & stnlats > NW_lat(1) & stnlats < NW_lat(2) & stnlons > NW_lon(1) & stnlons < NW_lon(2) & NW_lats + NW_lons == 2);

c = JA_colors(NW_exceed_ind);

subplot(3,2,4);

m_proj('Equidistant Cylindrical','lon',NW_lon,'lat',NW_lat);
m_coast('patch',[.85 .85 .85],'edgecolor','black');
m_grid('xtick', 3, 'linestyle','none','tickdir','out');
hold on

m_scatter(stnlons(NW_0_ind),stnlats(NW_0_ind),100,'k','o')
hold on

m_scatter(stnlons(NW_exceed_ind),stnlats(NW_exceed_ind),100,c,'o','filled')

colormap(flipud(hot))
h = colorbar;
ylabel(h, 'Total Exceedances')
hold off

% Select SE data and calculate risk
SE_risk_alltemps = sum(OMI_thresh_1std(:,:,2,SE_ind),4)./sum(OMI_obs_len(:,2,SE_ind),3)*100;
scrambled_risk_SE = squeeze(sum(OMI_scrambled_thresh(:,:,2,length(coastal_ind)+ length(NW_ind)+1:end,:),4))./sum(OMI_MJOdays_scrambled_thresh(:,2,length(coastal_ind)+ length(NW_ind)+1:end,:),3)*100;

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

% Plot station locations
addpath /home/ivanov/matlab/Code/m_map/

SE_exceed_ind = find(JA_exceed ~= 0 & stnlats > SE_lat(1) & stnlats < SE_lat(2) & stnlons > SE_lon(1) & stnlons < SE_lon(2) & SE_lats + SE_lons == 2);
SE_0_ind = find(JA_exceed == 0 & stnlats > SE_lat(1) & stnlats < SE_lat(2) & stnlons > SE_lon(1) & stnlons < SE_lon(2) & SE_lats + SE_lons == 2);

c = JA_colors(SE_exceed_ind);

subplot(3,2,6);

m_proj('Equidistant Cylindrical','lon',SE_lon,'lat',SE_lat);
m_coast('patch',[.85 .85 .85],'edgecolor','black');
m_grid('xtick', 3, 'linestyle','none','tickdir','out');
hold on

m_scatter(stnlons(SE_0_ind),stnlats(SE_0_ind),100,'k','o')
hold on

m_scatter(stnlons(SE_exceed_ind),stnlats(SE_exceed_ind),100,c,'o','filled')

colormap(flipud(hot))
h = colorbar;
ylabel(h, 'Total Exceedances')
hold off