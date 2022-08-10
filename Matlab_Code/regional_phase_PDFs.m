%% Create full wbt cell
addpath /home/ivanov/matlab/Code/MJO_TW
load('fullcell_updated_part1.mat')
load('fullcell_updated_part2.mat')

wbtmat = {fullcell_updated_part1, fullcell_updated_part2};
fullwbtcell = horzcat(wbtmat{:});

load('stndata.mat')

%% Determine phase of each active BSISO day

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

%% Doha

% Select CME data from cell
CME_ind = 3726;
CME_alldata = fullwbtcell{1,CME_ind};
    
% Find OMI where dates match that of station data   
OMI_dates = datenum(double(OMI(:,1)),double(OMI(:,2)),double(OMI(:,3)));
CME_dates = datenum(double(CME_alldata(:,1)), double(CME_alldata(:,2)), double(CME_alldata(:,3)));
[inter,iO,iW] = intersect(OMI_dates, CME_dates);

CME_table = horzcat(CME_alldata(iW,:),OMI(iO,7:8));

% Remove outliers
outliers = find(CME_table(:,8)>10);
CME_table(outliers,:) = [];
CME_nomissing = rmmissing(CME_table);

% Assign inactive MJO days to phase 0
std_amp = std(CME_nomissing(:,8));
inactive = find(CME_nomissing(:,8) < std_amp);
CME_nomissing(inactive,9) = 0;

% Select JA data
JA_CMEindex = find(CME_nomissing(:,2) == 7 | CME_nomissing(:,2) == 8);
CME_JA = CME_nomissing(JA_CMEindex,:);

%% Jacobabad

NW_ind = 3770;
NW_alldata = fullwbtcell{1,NW_ind};
    
% Find OMI where dates match that of station data   
OMI_dates = datenum(double(OMI(:,1)),double(OMI(:,2)),double(OMI(:,3)));
NW_dates = datenum(double(NW_alldata(:,1)), double(NW_alldata(:,2)), double(NW_alldata(:,3)));
[inter,iO,iW] = intersect(OMI_dates, NW_dates);

NW_table = horzcat(NW_alldata(iW,:),OMI(iO,7:8));

% Remove outliers
outliers = find(NW_table(:,8)>10);
NW_table(outliers,:) = [];
NW_nomissing = rmmissing(NW_table);

% Assign inactive MJO days to phase 0
std_amp = std(NW_nomissing(:,8));
inactive = find(NW_nomissing(:,8) < std_amp);
NW_nomissing(inactive,9) = 0;

% Select JA data
JA_NWindex = find(NW_nomissing(:,2) == 7 | NW_nomissing(:,2) == 8);
NW_JA = NW_nomissing(JA_NWindex,:);

%% Southeastern India

SE_ind = 3824;
SE_alldata = fullwbtcell{1,SE_ind};
    
% Find OMI where dates match that of station data   
OMI_dates = datenum(double(OMI(:,1)),double(OMI(:,2)),double(OMI(:,3)));
SE_dates = datenum(double(SE_alldata(:,1)), double(SE_alldata(:,2)), double(SE_alldata(:,3)));
[inter,iO,iW] = intersect(OMI_dates, SE_dates);

SE_table = horzcat(SE_alldata(iW,:),OMI(iO,7:8));

% Remove outliers
outliers = find(SE_table(:,8)>10);
SE_table(outliers,:) = [];
SE_nomissing = rmmissing(SE_table);

% Assign inactive MJO days to phase 0
std_amp = std(SE_nomissing(:,8));
inactive = find(SE_nomissing(:,8) < std_amp);
SE_nomissing(inactive,9) = 0;

% Select JJA data
MJ_SEindex = find(SE_nomissing(:,2) == 5 |SE_nomissing(:,2) == 6);
SE_MJ = SE_nomissing(MJ_SEindex,:);


%% Plot just most and least likely phases

%%%%%%% CME
figure('Position', [10 10 600 900])
subplot(3,1,1)

CME_high = CME_JA(CME_JA(:,9) == 3,5);
CME_high_mean = mean(CME_high);
h1 = histogram(CME_high,40, 'FaceColor','b','FaceAlpha', 0.5, 'EdgeColor', 'b');
xline(CME_high_mean,'--','Color','b', 'LineWidth', 4);

hold on

CME_low = CME_JA(CME_JA(:,9) == 8,5);
CME_low_mean = mean(CME_low);
h2 = histogram(CME_low,40, 'FaceColor',[0.3660 0.5740 0.0880],'FaceAlpha', 0.5, 'EdgeColor', [0.3660 0.5740 0.0880]);
xline(CME_low_mean, '--','Color',[0.3660 0.5740 0.0880], 'LineWidth',4);

legend([h1,h2],{'Phase 3','Phase 8'}, 'Location','northwest')
title('Doha')
ylabel('Frequency')
xlim([20 37])
ylim([0 40])
hold off

%%%%%%% NWSA
subplot(3,1,2)
NW_high = NW_JA(NW_JA(:,9) == 6,5);
NW_high_mean = mean(NW_high);
h3 = histogram(NW_high,40, 'FaceColor','b','FaceAlpha', 0.5, 'EdgeColor', 'b');
xline(NW_high_mean,'--','Color','b', 'LineWidth', 4);

hold on

NW_low = NW_JA(NW_JA(:,9) == 4,5);
NW_low_mean = mean(NW_low);
h4 = histogram(NW_low,40, 'FaceColor',[0.3660 0.5740 0.0880],'FaceAlpha', 0.5, 'EdgeColor', [0.3660 0.5740 0.0880]);
xline(NW_low_mean, '--','Color',[0.3660 0.5740 0.0880], 'LineWidth',4);

legend([h3,h4],{'Phase 6','Phase 4'}, 'Location','northwest')
title('Jacobabad')
ylabel('Frequency')
xlim([20 37])
ylim([0 40])

hold off

%%%%%%% SEI
subplot(3,1,3)
SE_high = SE_MJ(SE_MJ(:,9) == 1,5);
SE_high_mean = mean(SE_high);
h5 = histogram(SE_high,40, 'FaceColor','b','FaceAlpha', 0.5, 'EdgeColor', 'b');
xline(SE_high_mean,'--','Color','b', 'LineWidth', 4);

hold on

SE_low = SE_MJ(SE_MJ(:,9) == 5,5);
SE_low_mean = mean(NW_low);
h6 = histogram(SE_low,40, 'FaceColor',[0.3660 0.5740 0.0880],'FaceAlpha', 0.5, 'EdgeColor', [0.3660 0.5740 0.0880]);
xline(SE_low_mean, '--','Color',[0.3660 0.5740 0.0880], 'LineWidth',4);
    
legend([h5,h6],{'Phase 1','Phase 5'}, 'Location','northwest')
title('Bhubaneswar')
ylabel('Frequency')
xlim([20 37])
ylim([0 40])
xlabel('Wet Bulb Temperature (C)')

hold off










