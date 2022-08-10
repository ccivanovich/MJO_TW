%% Determine phase of each active BSISO day

addpath /home/ivanov/matlab/Code/MJO_TW

% Load data from file
fid = fopen('OMI_Index.txt','rt');
D = textscan(fid, '%f%f%f%f%f%f%f');
fclose(fid);
OMI_ind = cell2mat(D);

% Convert to polar coordinates and determine phase
[theta,rho] = cart2pol(OMI_ind(:,6),-OMI_ind(:,5));
angles = [-pi -3*pi/4 -2*pi/4 -pi/4 0 pi/4 2*pi/4 3*pi/4 pi];

phase = 0*theta;

for ii = 1:length(theta)
    if theta(ii) >= angles(1) && theta(ii) < angles(2)
        phase(ii) = 1;
    elseif theta(ii) >= angles(2) && theta(ii) < angles(3)
        phase(ii) = 2;
    elseif theta(ii) >= angles(3) && theta(ii) < angles(4)
        phase(ii) = 3;
    elseif theta(ii) >= angles(4) && theta(ii) < angles(5)
        phase(ii) = 4;
    elseif theta(ii) >= angles(5) && theta(ii) < angles(6)
        phase(ii) = 5;
    elseif theta(ii) >= angles(6) && theta(ii) < angles(7)
        phase(ii) = 6;
    elseif theta(ii) >= angles(7) && theta(ii) < angles(8)
        phase(ii) = 7;
    elseif theta(ii) >= angles(8) && theta(ii) < angles(9)
        phase(ii) = 8;
    end
end

% Add phase onto OMI index
OMI = horzcat(OMI_ind,phase);
%% Calculate number of days above 28C in each phase per year

addpath /home/ivanov/matlab/Code/MJO_TW
load('TW_allstations.mat')

load('stndata.mat')
%% Persian Gulf
fid = fopen('persiangulfstns.txt','rt');
P = textscan(fid, '%f%f');
fclose(fid);
PG_locations = cell2mat(P);
phases = 0:1:8;

% Identify coastal middle east stations
CME_lon = [45 60];
CME_lat = [20 36];
coastal = PG_locations(PG_locations(:,1) > CME_lat(1) & PG_locations(:,1) < CME_lat(2) & PG_locations(:,2) > CME_lon(1) & PG_locations(:,2) < CME_lon(2)-6,:);

% Find station indices  
coast_lats = ismember(stnlats, coastal(:,1));
coast_lons = ismember(stnlons, coastal(:,2));
coastal_ind = find(coast_lats + coast_lons == 2);

% Select PG station data and stitch together
PG_TW = [];

for i = 1:length(coastal_ind)
    PG_station = C{coastal_ind(i)};
    
    OMI_dates = datenum(double(OMI(:,1)),double(OMI(:,2)),double(OMI(:,3)));
    PG_station_dates = datenum(double(PG_station(:,2)), double(PG_station(:,3)), double(PG_station(:,4)));
    [inter,iO,iW] = intersect(OMI_dates, PG_station_dates);
    
    PG_station_phases = horzcat(PG_station(iW,:),OMI(iO,7:8));
    
    % Remove outliers
    outliers = find(PG_station_phases(:,5)>10);
    PG_station_phases(outliers,:) = [];
    PG_nomissing = rmmissing(PG_station_phases);
    
    % Assign inactive MJO/BSISO days to phase 0
    std_amp = std(PG_nomissing(:,5));
    inactive = find(PG_nomissing(:,5) < std_amp);
    PG_nomissing(inactive,6) = 0;
    
    PG_TW = vertcat(PG_TW,PG_station_phases);
    
end


%% Count number of days in each phase and number of days above 28C in each phase
count_28C = zeros(2019-1979,9);
phase_count = zeros(2019-1979,9);

for year = 1979:2019
    
    year_ind = year -1978;
    
    for phase = 0:8
        
        PG_phase28C = PG_TW(PG_TW(:,2) == year & PG_TW(:,6) == phase & PG_TW(:,1) >= 28,:);
        phase_days = PG_TW(PG_TW(:,2) == year & PG_TW(:,6) == phase,:);
        count_28C(year_ind,phase+1) = length(PG_phase28C);
        phase_count(year_ind,phase+1) = length(phase_days);
        
    end
end

%% Calculate trend in number of days in each phase and number of days in each phase above 28C
figure('Position', [10 10 900 700])
set(0,'defaultAxesFontSize',10)

for i = 1:8
    a(i) = subplot(4,2,i);
    years = [1979:2019];
    counts = phase_count(:,i+1)/length(coastal_ind);
    plot(years,counts,'Color', 'k')
    hold on
    
    % Linear regression
    f1 = fitlm(years, counts);
    counts_predict = f1.predict(years')';
    R2 = f1.Rsquared.Ordinary;
    coeff = f1.Coefficients.Estimate;
    
    plot(years, counts_predict);

    xlabel('Year')
    ylabel({'Daily Counts'})
    b(i) = annotation('textbox','String', sprintf('Slope = %g, R2 = %f', coeff(2), R2),'Position',a(i).Position,'Vert','top','EdgeColor','none','FontSize',8);
    xlim([1979,2019])
    ylim([0,100])
    title(num2str(i,'Phase %d'))
end

%%
% Extremes frequency
ext_freq = 100*count_28C./phase_count;

figure('Position', [10 10 900 700])
set(0,'defaultAxesFontSize',10)
titles = ['Phase 1','Phase 2','Phase 3','Phase 4','Phase 5','Phase 6','Phase 7','Phase 8'];

for i = 1:8
    a(i) = subplot(4,2,i);
    years = [1979:2019];
    freq = ext_freq(:,i+1);
    plot(years,freq,'Color', 'k')
    hold on
    
    % Linear regression
    f2 = fitlm(years, freq);
    counts_predict = f2.predict(years')';
    R2 = f2.Rsquared.Ordinary;
    coeff = f2.Coefficients.Estimate;
    
    plot(years, counts_predict);
    
    
    xlabel('Year')
    ylabel({'Frequency of'; 'Exceeding 28C'})
    b(i) = annotation('textbox','String', sprintf('Slope = %g, R2 = %f', coeff(2), R2),'Position',a(i).Position,'Vert','top','EdgeColor','none','FontSize',8);
    xlim([1979,2019])
    title(num2str(i,'Phase %d'))
end

%% NW South Asia
% Convert text file to matrix of lat/lon
fid = fopen('southasiastns.txt','rt');
S = textscan(fid, '%f%f');
fclose(fid);
SA_locations = cell2mat(S);

% Nothern South Asia
NW_lon = [68 78];
NW_lat = [22 32];

% Find station indices  
NW_lats = ismember(stnlats, SA_locations(:,1));
NW_lons = ismember(stnlons, SA_locations(:,2));
NW_ind = find(NW_lats + NW_lons == 2 & stnlats > NW_lat(1) & stnlats < NW_lat(2) & stnlons > NW_lon(1) & stnlons < NW_lon(2));

% Select PG station data and stitch together
NW_TW = [];

for i = 1:length(NW_ind)
    NW_station = C{NW_ind(i)};
    
    OMI_dates = datenum(double(OMI(:,1)),double(OMI(:,2)),double(OMI(:,3)));
    NW_station_dates = datenum(double(NW_station(:,2)), double(NW_station(:,3)), double(NW_station(:,4)));
    [inter,iO,iW] = intersect(OMI_dates, NW_station_dates);
    
    NW_station_phases = horzcat(NW_station(iW,:),OMI(iO,7:8));
    
    % Remove outliers
    outliers = find(NW_station_phases(:,5)>10);
    NW_station_phases(outliers,:) = [];
    NW_nomissing = rmmissing(NW_station_phases);
    
    % Assign inactive MJO/BSISO days to phase 0
    std_amp = std(NW_nomissing(:,5));
    inactive = find(NW_nomissing(:,5) < std_amp);
    NW_nomissing(inactive,6) = 0;
    
    NW_TW = vertcat(NW_TW,NW_station_phases);
    
end


% Count number of days in each phase and number of days above 28C in each phase
count_28C = zeros(2019-1979,9);
phase_count = zeros(2019-1979,9);

for year = 1979:2019
    
    year_ind = year -1978;
    
    for phase = 0:8
        
        NW_phase28C = NW_TW(NW_TW(:,2) == year & NW_TW(:,6) == phase & NW_TW(:,1) >= 28,:);
        phase_days = NW_TW(NW_TW(:,2) == year & NW_TW(:,6) == phase,:);
        count_28C(year_ind,phase+1) = length(NW_phase28C);
        phase_count(year_ind,phase+1) = length(phase_days);
        
    end
end

% Calculate trend in number of days in each phase and number of days in each phase above 28C
figure('Position', [10 10 900 700])
set(0,'defaultAxesFontSize',10)

for i = 1:8
    a(i) = subplot(4,2,i);
    years = [1979:2019];
    counts = phase_count(:,i+1)/length(NW_ind);
    plot(years,counts,'Color', 'k')
    hold on
    
    % Linear regression
    f1 = fitlm(years, counts);
    counts_predict = f1.predict(years')';
    R2 = f1.Rsquared.Ordinary;
    coeff = f1.Coefficients.Estimate;
    
    plot(years, counts_predict);

    xlabel('Year')
    ylabel({'Daily Counts'})
    b(i) = annotation('textbox','String', sprintf('Slope = %g, R2 = %f', coeff(2), R2),'Position',a(i).Position,'Vert','top','EdgeColor','none','FontSize',8);
    xlim([1979,2019])
    ylim([0,100])
    title(num2str(i,'Phase %d'))
end

% Extremes frequency
ext_freq = 100*count_28C./phase_count;

figure('Position', [10 10 900 700])
set(0,'defaultAxesFontSize',10)
titles = ['Phase 1','Phase 2','Phase 3','Phase 4','Phase 5','Phase 6','Phase 7','Phase 8'];

for i = 1:8
    a(i) = subplot(4,2,i);
    years = [1979:2019];
    freq = ext_freq(:,i+1);
    plot(years,freq,'Color', 'k')
    hold on
    
    % Linear regression
    f2 = fitlm(years, freq);
    counts_predict = f2.predict(years')';
    R2 = f2.Rsquared.Ordinary;
    coeff = f2.Coefficients.Estimate;
    
    plot(years, counts_predict);
    
    
    xlabel('Year')
    ylabel({'Frequency of'; 'Exceeding 28C'})
    b(i) = annotation('textbox','String', sprintf('Slope = %g, R2 = %f', coeff(2), R2),'Position',a(i).Position,'Vert','top','EdgeColor','none','FontSize',8);
    xlim([1979,2019])
    title(num2str(i,'Phase %d'))
end

%% SE India
% Convert text file to matrix of lat/lon
% SE South Asia
SE_lon = [78 90];
SE_lat = [8 22];

% Find station indices  
SE_lats = ismember(stnlats, SA_locations(:,1));
SE_lons = ismember(stnlons, SA_locations(:,2));
SE_ind = find(SE_lats + SE_lons == 2 & stnlats > SE_lat(1) & stnlats < SE_lat(2) & stnlons > SE_lon(1) & stnlons < SE_lon(2));

% Select PG station data and stitch together
SE_TW = [];

for i = 1:length(SE_ind)
    SE_station = C{SE_ind(i)};
    
    OMI_dates = datenum(double(OMI(:,1)),double(OMI(:,2)),double(OMI(:,3)));
    SE_station_dates = datenum(double(SE_station(:,2)), double(SE_station(:,3)), double(SE_station(:,4)));
    [inter,iO,iW] = intersect(OMI_dates, SE_station_dates);
    
    SE_station_phases = horzcat(SE_station(iW,:),OMI(iO,7:8));
    
    % Remove outliers
    outliers = find(SE_station_phases(:,5)>10);
    SE_station_phases(outliers,:) = [];
    SE_nomissing = rmmissing(SE_station_phases);
    
    % Assign inactive MJO/BSISO days to phase 0
    std_amp = std(SE_nomissing(:,5));
    inactive = find(SE_nomissing(:,5) < std_amp);
    SE_nomissing(inactive,6) = 0;
    
    SE_TW = vertcat(SE_TW,SE_station_phases);
    
end


% Count number of days in each phase and number of days above 28C in each phase
count_28C = zeros(2019-1979,9);
phase_count = zeros(2019-1979,9);

for year = 1979:2019
    
    year_ind = year -1978;
    
    for phase = 0:8
        
        SE_phase28C = SE_TW(SE_TW(:,2) == year & SE_TW(:,6) == phase & SE_TW(:,1) >= 28,:);
        phase_days = SE_TW(SE_TW(:,2) == year & SE_TW(:,6) == phase,:);
        count_28C(year_ind,phase+1) = length(SE_phase28C);
        phase_count(year_ind,phase+1) = length(phase_days);
        
    end
end

% Calculate trend in number of days in each phase and number of days in each phase above 28C
figure('Position', [10 10 900 700])
set(0,'defaultAxesFontSize',10)

for i = 1:8
    a(i) = subplot(4,2,i);
    years = [1979:2019];
    counts = phase_count(:,i+1)/length(SE_ind);
    plot(years,counts,'Color', 'k')
    hold on
    
    % Linear regression
    f1 = fitlm(years, counts);
    counts_predict = f1.predict(years')';
    R2 = f1.Rsquared.Ordinary;
    coeff = f1.Coefficients.Estimate;
    
    plot(years, counts_predict);

    xlabel('Year')
    ylabel({'Daily Counts'})
    b(i) = annotation('textbox','String', sprintf('Slope = %g, R2 = %f', coeff(2), R2),'Position',a(i).Position,'Vert','top','EdgeColor','none','FontSize',8);
    xlim([1979,2019])
    ylim([0,100])
    title(num2str(i,'Phase %d'))
end

% Extremes frequency
ext_freq = 100*count_28C./phase_count;

figure('Position', [10 10 900 700])
set(0,'defaultAxesFontSize',10)
titles = ['Phase 1','Phase 2','Phase 3','Phase 4','Phase 5','Phase 6','Phase 7','Phase 8'];

for i = 1:8
    a(i) = subplot(4,2,i);
    years = [1979:2019];
    freq = ext_freq(:,i+1);
    plot(years,freq,'Color', 'k')
    hold on
    
    % Linear regression
    f2 = fitlm(years, freq);
    counts_predict = f2.predict(years')';
    R2 = f2.Rsquared.Ordinary;
    coeff = f2.Coefficients.Estimate;
    
    plot(years, counts_predict);
    
    
    xlabel('Year')
    ylabel({'Frequency of'; 'Exceeding 28C'})
    b(i) = annotation('textbox','String', sprintf('Slope = %g, R2 = %f', coeff(2), R2),'Position',a(i).Position,'Vert','top','EdgeColor','none','FontSize',8);
    xlim([1979,2019])
    title(num2str(i,'Phase %d'))
end