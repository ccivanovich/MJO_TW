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

%% Prep for constant Tw lines

% Constant Tw lines prep
slp=1010;
CME_pressure = mean(pressurefromelev(stnelevs(coastal_ind)).*slp./1000.*100);
N = 1000;
temperature = linspace(20, 50, N);
spechum = flipud(linspace(0,0.03, N).');

wbt = [];

for i = 1:N
    for j = 1:N
    
        wbt(j,i) = calcwbt_daviesjones(temperature(i),CME_pressure,spechum(j));
    
    end
    
    display(i)

end

% Round to 1 decimal place
wbt_1dec = round(wbt,1);

%% Find "average grey squares" for all phases

% CME
% Select top 5 percent
topten = prctile(CME_JA(:,5),90);
topfifty = prctile(CME_JA(:,5),50);
CME_topten = CME_JA(CME_JA(:,5) >= topten,:);
CME_avg_topten = mean(CME_topten);

% Select bottom 50%
CME_bottom50 = CME_JA(CME_JA(:,5) < topfifty,:);
CME_avg_bottom50 = mean(CME_bottom50);

% NW
% Select top 5 percent
topten = prctile(NW_JA(:,5),90);
topfifty = prctile(NW_JA(:,5),50);
NW_topten = NW_JA(NW_JA(:,5) >= topten,:);
NW_avg_topten = mean(NW_topten);

% Select bottom 50%
NW_bottom50 = NW_JA(NW_JA(:,5) < topfifty,:);
NW_avg_bottom50 = mean(NW_bottom50);

% SE
% Select top 5 percent
topten = prctile(SE_MJ(:,5),90);
topfifty = prctile(SE_MJ(:,5),50);
SE_topten = SE_MJ(SE_MJ(:,5) >= topten,:);
SE_avg_topten = mean(SE_topten);

% Select bottom 50%
SE_bottom50 = SE_MJ(SE_MJ(:,5) < topfifty,:);
SE_avg_bottom50 = mean(SE_bottom50,1);


%% Plotted Together

% Create contours of full JJA dataset for CME
x_CME = [CME_JA(:,6), CME_JA(:,7)];   % Set x = two column matrix with temp and SH for full dataset
[n_CME, c_CME] = hist3(x_CME, 'Edges', {20:5:55 -0.0025:0.0025:0.0325});

% Create contours of full JJA dataset for NW
x_NW = [NW_JA(:,6), NW_JA(:,7)];   % Set x = two column matrix with temp and SH for full dataset
[n_NW, c_NW] = hist3(x_NW, 'Edges', {20:5:55 -0.0025:0.0025:0.0325});

% Create contours of full JJA dataset for SE
x_SE = [SE_MJ(:,6), SE_MJ(:,7)];   % Set x = two column matrix with temp and SH for full dataset
[n_SE, c_SE] = hist3(x_SE, 'Edges', {20:5:55 -0.0025:0.0025:0.0325});

figure('Position', [10 10 600 700])
set(0,'defaultAxesFontSize',8)
RGB = [0.0 0.6 1];
RGB2 = [0.5 0.5 0.5];

% Plot
levels = 20:5:35;
RGB3 = [0 0 0.5];

% CME
phasei = 3;
JJA_byphase = CME_JA(CME_JA(:,9) == phasei,:);

% Select top 5 percent
topten = prctile(JJA_byphase(:,5),90);
topfifty = prctile(JJA_byphase(:,5),50);
JJA_topten = JJA_byphase(JJA_byphase(:,5) >= topten,:);
avg_topten = mean(JJA_topten, 1);

% Select bottom 50%
JJA_bottom50 = JJA_byphase(JJA_byphase(:,5) < topfifty,:);
avg_bottom50 = mean(JJA_bottom50, 1);

a = subplot(3, 2, 1);

% Plot lines of constant Tw
contour(temperature, spechum, wbt_1dec, levels,'LineColor', RGB3,'LineWidth',2)
hold on

% Plot PDF contours
contour(c_CME{1},c_CME{2},n_CME.', 'LineColor', RGB2, 'LineWidth', 0.7)
hold on

% Scatter plots
scatter(JJA_topten(:,6), JJA_topten(:,7), 50,'r', '.')
hold on
scatter(JJA_bottom50(:,6), JJA_bottom50(:,7), 50,'k', '.')
hold on

% Grey and blue averages
scatter(CME_avg_topten(:,6), CME_avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
hold on
scatter(CME_avg_bottom50(:,6), CME_avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
hold on

scatter(avg_topten(:,6), avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
hold on
scatter(avg_bottom50(:,6), avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
hold on

b = annotation('textbox','String',num2str(phasei,'Phase %d'),'Position',a.Position,'Vert','bot','Horiz','left','FitBoxToText','off');

xlim([20 50])
ylim([0 0.03])
yticks([0 0.01 0.02 0.03])
ylabel('Specific Humidity (kg/kg)')
hold off

% CME
phasei = 8;
JJA_byphase = CME_JA(CME_JA(:,9) == phasei,:);

% Select top 5 percent
topten = prctile(JJA_byphase(:,5),90);
topfifty = prctile(JJA_byphase(:,5),50);
JJA_topten = JJA_byphase(JJA_byphase(:,5) >= topten,:);
avg_topten = mean(JJA_topten, 1);

% Select bottom 50%
JJA_bottom50 = JJA_byphase(JJA_byphase(:,5) < topfifty,:);
avg_bottom50 = mean(JJA_bottom50, 1);

a = subplot(3, 2, 2);

% Plot lines of constant Tw
contour(temperature, spechum, wbt_1dec, levels,'LineColor', RGB3,'LineWidth',2)
hold on

% Plot PDF contours
contour(c_CME{1},c_CME{2},n_CME.', 'LineColor', RGB2, 'LineWidth', 0.7)
hold on

% Scatter plots
scatter(JJA_topten(:,6), JJA_topten(:,7), 50,'r', '.')
hold on
scatter(JJA_bottom50(:,6), JJA_bottom50(:,7), 50,'k', '.')
hold on

% Grey and blue averages
scatter(CME_avg_topten(:,6), CME_avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
hold on
scatter(CME_avg_bottom50(:,6), CME_avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
hold on
scatter(avg_topten(:,6), avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
hold on
scatter(avg_bottom50(:,6), avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
hold on

b = annotation('textbox','String',num2str(phasei,'Phase %d'),'Position',a.Position,'Vert','bot','Horiz','left','FitBoxToText','off');

xlim([20 50])
ylim([0 0.03])
yticks([0 0.01 0.02 0.03])
hold off


% Pakistan
phasei = 6;
JJA_byphase = NW_JA(NW_JA(:,9) == phasei,:);

% Select top 5 percent
topten = prctile(JJA_byphase(:,5),90);
topfifty = prctile(JJA_byphase(:,5),50);
JJA_topten = JJA_byphase(JJA_byphase(:,5) >= topten,:);
avg_topten = mean(JJA_topten, 1);

% Select bottom 50%
JJA_bottom50 = JJA_byphase(JJA_byphase(:,5) < topfifty,:);
avg_bottom50 = mean(JJA_bottom50, 1);

a = subplot(3, 2, 3);

% Plot lines of constant Tw
contour(temperature, spechum, wbt_1dec, levels,'LineColor', RGB3,'LineWidth',2)
hold on

% Plot PDF contours
contour(c_NW{1},c_NW{2},n_NW.', 'LineColor', RGB2, 'LineWidth', 0.7)
hold on

% Scatter plots
scatter(JJA_topten(:,6), JJA_topten(:,7), 50,'r', '.')
hold on
scatter(JJA_bottom50(:,6), JJA_bottom50(:,7), 50,'k', '.')
hold on

% Grey and blue averages
scatter(NW_avg_topten(:,6), NW_avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
hold on
scatter(NW_avg_bottom50(:,6), NW_avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
hold on
scatter(avg_topten(:,6), avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
hold on
scatter(avg_bottom50(:,6), avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
hold on

b = annotation('textbox','String',num2str(phasei,'Phase %d'),'Position',a.Position,'Vert','bot','Horiz','left','FitBoxToText','off');

xlim([20 50])
ylim([0 0.03])
yticks([0 0.01 0.02 0.03])
ylabel('Specific Humidity (kg/kg)')
hold off

% Pakistan
phasei = 4;
JJA_byphase = NW_JA(NW_JA(:,9) == phasei,:);

% Select top 5 percent
topten = prctile(JJA_byphase(:,5),90);
topfifty = prctile(JJA_byphase(:,5),50);
JJA_topten = JJA_byphase(JJA_byphase(:,5) >= topten,:);
avg_topten = mean(JJA_topten, 1);

% Select bottom 50%
JJA_bottom50 = JJA_byphase(JJA_byphase(:,5) < topfifty,:);
avg_bottom50 = mean(JJA_bottom50, 1);

a = subplot(3, 2, 4);

% Plot lines of constant Tw
contour(temperature, spechum, wbt_1dec, levels,'LineColor', RGB3,'LineWidth',2)
hold on

% Plot PDF contours
contour(c_NW{1},c_NW{2},n_NW.', 'LineColor', RGB2, 'LineWidth', 0.7)
hold on

% Scatter plots
scatter(JJA_topten(:,6), JJA_topten(:,7), 50,'r', '.')
hold on
scatter(JJA_bottom50(:,6), JJA_bottom50(:,7), 50,'k', '.')
hold on

% Grey and blue averages
scatter(NW_avg_topten(:,6), NW_avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
hold on
scatter(NW_avg_bottom50(:,6), NW_avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
hold on
scatter(avg_topten(:,6), avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
hold on
scatter(avg_bottom50(:,6), avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
hold on

b = annotation('textbox','String',num2str(phasei,'Phase %d'),'Position',a.Position,'Vert','bot','Horiz','left','FitBoxToText','off');

xlim([20 50])
ylim([0 0.03])
yticks([0 0.01 0.02 0.03])
hold off

% India
phasei = 1;
JJA_byphase = SE_MJ(SE_MJ(:,9) == phasei,:);

% Select top 5 percent
topten = prctile(JJA_byphase(:,5),90);
topfifty = prctile(JJA_byphase(:,5),50);
JJA_topten = JJA_byphase(JJA_byphase(:,5) >= topten,:);
avg_topten = mean(JJA_topten, 1);

% Select bottom 50%
JJA_bottom50 = JJA_byphase(JJA_byphase(:,5) < topfifty,:);
avg_bottom50 = mean(JJA_bottom50, 1);

a = subplot(3, 2, 5);

% Plot lines of constant Tw
contour(temperature, spechum, wbt_1dec, levels,'LineColor', RGB3,'LineWidth',2)
hold on

% Plot PDF contours
contour(c_SE{1},c_SE{2},n_SE.', 'LineColor', RGB2, 'LineWidth', 0.7)
hold on

% Scatter plots
scatter(JJA_topten(:,6), JJA_topten(:,7), 50,'r', '.')
hold on
scatter(JJA_bottom50(:,6), JJA_bottom50(:,7), 50,'k', '.')
hold on

% Grey and blue averages
scatter(SE_avg_topten(:,6), SE_avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
hold on
scatter(SE_avg_bottom50(:,6), SE_avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
hold on
scatter(avg_topten(:,6), avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
hold on
scatter(avg_bottom50(:,6), avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
hold on

b = annotation('textbox','String',num2str(phasei,'Phase %d'),'Position',a.Position,'Vert','bot','Horiz','left','FitBoxToText','off');

xlim([20 50])
ylim([0 0.03])
yticks([0 0.01 0.02 0.03])
xlabel('Temperature (C)')
ylabel('Specific Humidity (kg/kg)')
hold off

% India
phasei = 5;
JJA_byphase = SE_MJ(SE_MJ(:,9) == phasei,:);

% Select top 5 percent
topten = prctile(JJA_byphase(:,5),90);
topfifty = prctile(JJA_byphase(:,5),50);
JJA_topten = JJA_byphase(JJA_byphase(:,5) >= topten,:);
avg_topten = mean(JJA_topten, 1);

% Select bottom 50%
JJA_bottom50 = JJA_byphase(JJA_byphase(:,5) < topfifty,:);
avg_bottom50 = mean(JJA_bottom50, 1);

a = subplot(3, 2, 6);

% Plot lines of constant Tw
contour(temperature, spechum, wbt_1dec, levels,'LineColor', RGB3,'LineWidth',2)
hold on

% Plot PDF contours
contour(c_SE{1},c_SE{2},n_SE.', 'LineColor', RGB2, 'LineWidth', 0.7)
hold on

% Scatter plots
scatter(JJA_topten(:,6), JJA_topten(:,7), 50,'r', '.')
hold on
scatter(JJA_bottom50(:,6), JJA_bottom50(:,7), 50,'k', '.')
hold on

% Grey and blue averages
scatter(SE_avg_topten(:,6), SE_avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
hold on
scatter(SE_avg_bottom50(:,6), SE_avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
hold on
scatter(avg_topten(:,6), avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
hold on
scatter(avg_bottom50(:,6), avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
hold on

b = annotation('textbox','String',num2str(phasei,'Phase %d'),'Position',a.Position,'Vert','bot','Horiz','left','FitBoxToText','off');

xlim([20 50])
ylim([0 0.03])
yticks([0 0.01 0.02 0.03])
xlabel('Temperature (C)')
hold off

%% CME: Plot with top 10% and bottom 50%

% Create contours of full JJA dataset for CME
x_CME = [CME_JA(:,6), CME_JA(:,7)];   % Set x = two column matrix with temp and SH for full dataset
[n_CME, c_CME] = hist3(x_CME, 'Edges', {20:5:55 -0.0025:0.0025:0.0325});

figure('Position', [10 10 600 900])
set(0,'defaultAxesFontSize',8)
RGB = [0.0 0.6 1];
RGB2 = [0.5 0.5 0.5];

% Plot
levels = 20:5:35;
RGB3 = [0 0 0.5];

for phasei = 1:8

    JJA_byphase = CME_JA(CME_JA(:,9) == phasei,:);
    
    % Select top 5 percent
    topten = prctile(JJA_byphase(:,5),90);
    topfifty = prctile(JJA_byphase(:,5),50);
    JJA_topten = JJA_byphase(JJA_byphase(:,5) >= topten,:);
    avg_topten = mean(JJA_topten, 1);
    
    % Select bottom 50%
    JJA_bottom50 = JJA_byphase(JJA_byphase(:,5) < topfifty,:);
    avg_bottom50 = mean(JJA_bottom50, 1);
    
    a = subplot(4, 2, phasei);
    
    % Plot lines of constant Tw
    contour(temperature, spechum, wbt_1dec, levels,'LineColor', RGB3,'LineWidth',2)
    hold on
    
    % Plot PDF contours
    contour(c_CME{1},c_CME{2},n_CME.', 'LineColor', RGB2, 'LineWidth', 0.7)
    hold on
    
    % Scatter plots
    scatter(JJA_topten(:,6), JJA_topten(:,7), 50,'r', '.')
    hold on
    scatter(JJA_bottom50(:,6), JJA_bottom50(:,7), 50,'k', '.')
    hold on
    
    % Grey and blue averages
    scatter(CME_avg_topten(:,6), CME_avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
    hold on
    scatter(CME_avg_bottom50(:,6), CME_avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
    hold on
    scatter(avg_topten(:,6), avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
    hold on
    scatter(avg_bottom50(:,6), avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
    hold on
    xlim([25 50])
    ylim([0 0.03])
    yticks([0 0.01 0.02 0.03])
    b = annotation('textbox','String',num2str(phasei,'Phase %d'),'Position',a.Position,'Vert','bot','Horiz','left','FitBoxToText','off');
    
    if rem(phasei,2) ~= 0
        ylabel('Specific Humidity (g/kg)')
    end
        
    if (phasei == 7 | phasei == 8)
        xlabel('Temperature (C)')
    end
    
    hold off

end

%% Pakistan: Plot with top 10% and bottom 50%

% Create contours of full JJA dataset for CME
x_NW = [NW_JA(:,6), NW_JA(:,7)];   % Set x = two column matrix with temp and SH for full dataset
[n_NW, c_NW] = hist3(x_NW, 'Edges', {20:5:55 -0.0025:0.0025:0.0325});

figure('Position', [10 10 600 900])
set(0,'defaultAxesFontSize',8)
RGB = [0.0 0.6 1];
RGB2 = [0.5 0.5 0.5];

% Plot
levels = 20:5:35;
RGB3 = [0 0 0.5];

for phasei = 1:8

    JJA_byphase = NW_JA(NW_JA(:,9) == phasei,:);
    
    % Select top 5 percent
    topten = prctile(JJA_byphase(:,5),90);
    topfifty = prctile(JJA_byphase(:,5),50);
    JJA_topten = JJA_byphase(JJA_byphase(:,5) >= topten,:);
    avg_topten = mean(JJA_topten, 1);
    
    % Select bottom 50%
    JJA_bottom50 = JJA_byphase(JJA_byphase(:,5) < topfifty,:);
    avg_bottom50 = mean(JJA_bottom50, 1);
    
    a = subplot(4, 2, phasei)
    
    % Plot lines of constant Tw
    contour(temperature, spechum, wbt_1dec, levels,'LineColor', RGB3,'LineWidth',2)
    hold on
    
    % Plot PDF contours
    contour(c_NW{1},c_NW{2},n_NW.', 'LineColor', RGB2, 'LineWidth', 0.7)
    hold on
    
    % Scatter plots
    scatter(JJA_topten(:,6), JJA_topten(:,7), 50,'r', '.')
    hold on
    scatter(JJA_bottom50(:,6), JJA_bottom50(:,7), 50,'k', '.')
    hold on
    
    % Grey and blue averages
    scatter(NW_avg_topten(:,6), NW_avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
    hold on
    scatter(NW_avg_bottom50(:,6), NW_avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
    hold on
    scatter(avg_topten(:,6), avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
    hold on
    scatter(avg_bottom50(:,6), avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
    hold on
    xlim([25 50])
    ylim([0 0.03])
    yticks([0 0.01 0.02 0.03])
    b = annotation('textbox','String',num2str(phasei,'Phase %d'),'Position',a.Position,'Vert','bot','Horiz','left','FitBoxToText','off');
    
    if rem(phasei,2) ~= 0
        ylabel('Specific Humidity (kg/kg)')
    end
        
    if (phasei == 7 | phasei == 8)
        xlabel('Temperature (C)')
    end
    
    hold off

end

%% Southeastern India: Plot with top 10% and bottom 50%

% Create contours of full JJA dataset for CME
x_SE = [SE_MJ(:,6), SE_MJ(:,7)];   % Set x = two column matrix with temp and SH for full dataset
[n_SE, c_SE] = hist3(x_SE, 'Edges', {20:5:55 -0.0025:0.0025:0.0325});

figure('Position', [10 10 600 900])
set(0,'defaultAxesFontSize',8)
RGB = [0.0 0.6 1];
RGB2 = [0.5 0.5 0.5];

% Plot
levels = 20:5:35;
RGB3 = [0 0 0.5];

for phasei = 1:8

    JJA_byphase = SE_MJ(SE_MJ(:,9) == phasei,:);
    
    % Select top 5 percent
    topten = prctile(JJA_byphase(:,5),90);
    topfifty = prctile(JJA_byphase(:,5),50);
    JJA_topten = JJA_byphase(JJA_byphase(:,5) >= topten,:);
    avg_topten = mean(JJA_topten, 1);
    
    % Select bottom 50%
    JJA_bottom50 = JJA_byphase(JJA_byphase(:,5) < topfifty,:);
    avg_bottom50 = mean(JJA_bottom50, 1);
    
    a = subplot(4, 2, phasei);
    
    % Plot lines of constant Tw
    contour(temperature, spechum, wbt_1dec, levels,'LineColor', RGB3,'LineWidth',2)
    hold on
    
    % Plot PDF contours
    contour(c_SE{1},c_SE{2},n_SE.', 'LineColor', RGB2, 'LineWidth', 0.7)
    hold on
    
    % Scatter plots
    scatter(JJA_topten(:,6), JJA_topten(:,7), 50,'r', '.')
    hold on
    scatter(JJA_bottom50(:,6), JJA_bottom50(:,7), 50,'k', '.')
    hold on
    
    % Grey and blue averages
    scatter(SE_avg_topten(:,6), SE_avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
    hold on
    scatter(SE_avg_bottom50(:,6), SE_avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB2, 'MarkerFaceColor', RGB2)
    hold on
    scatter(avg_topten(:,6), avg_topten(:,7),125, 's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
    hold on
    scatter(avg_bottom50(:,6), avg_bottom50(:,7),125,  's', 'filled', 'MarkerEdgeColor', RGB, 'MarkerFaceColor', RGB)
    hold on
    xlim([25 50])
    ylim([0 0.03])
    yticks([0 0.01 0.02 0.03])
    b = annotation('textbox','String',num2str(phasei,'Phase %d'),'Position',a.Position,'Vert','bot','Horiz','left','FitBoxToText','off');
    
    if rem(phasei,2) ~= 0
        ylabel('Specific Humidity (kg/kg)')
    end
        
    if (phasei == 7 | phasei == 8)
        xlabel('Temperature (C)')
    end
    
    hold off

end

