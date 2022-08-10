%% Stations
addpath /home/ivanov/matlab/Code/MJO_TW
load('stndata.mat')
load('TW_allstations.mat')

% Qatar
fid = fopen('persiangulfstns.txt','rt');
P = textscan(fid, '%f%f');
fclose(fid);
PG_locations = cell2mat(P);
phases = 0:1:8;
season = ["DJF", "MAM", "JJA", "SON"];

% Identify CME stations
CME_lon = [45 60];
CME_lat = [20 36];
CME = PG_locations(PG_locations(:,1) > CME_lat(1) & PG_locations(:,1) < CME_lat(2) & PG_locations(:,2) > CME_lon(1) & PG_locations(:,2) < CME_lon(2),:);

% Find station indices  
coastal_lats = ismember(stnlats, CME(:,1));
coastal_lons = ismember(stnlons, CME(:,2));
coastal_ind = find(coastal_lats + coastal_lons == 2);

%NW SA
% Convert text file to matrix of lat/lon
fid = fopen('southasiastns.txt','rt');
S = textscan(fid, '%f%f');
fclose(fid);
SA_locations = cell2mat(S);

% Identify station locations
NW_lon = [68 78];
NW_lat = [22 32];
NW = SA_locations(SA_locations(:,1) > NW_lat(1) & SA_locations(:,1) < NW_lat(2) & SA_locations(:,2) > NW_lon(1) & SA_locations(:,2) < NW_lon(2),:);

% Find station indices 
NW_lats = ismember(stnlats, NW(:,1));
NW_lons = ismember(stnlons, NW(:,2));
NW_ind = find(NW_lats + NW_lons == 2);

% SE
% Identify station locations
SE_lon = [78 90];
SE_lat = [8 22];
SE = SA_locations(SA_locations(:,1) > SE_lat(1) & SA_locations(:,1) < SE_lat(2) & SA_locations(:,2) > SE_lon(1) & SA_locations(:,2) < SE_lon(2),:);

% Find station indices 
SE_lats = ismember(stnlats, SE(:,1));
SE_lons = ismember(stnlons, SE(:,2));
SE_ind = find(SE_lats + SE_lons == 2);

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

%% Run bootstrapping

% Use only stations we're interested in
stns = [coastal_ind; NW_ind; SE_ind];

n = 500;
OMI_scrambled = zeros(9,2,length(stns),n);
OMI_MJOdays_scrambled = zeros(9,2,length(stns),n);

for ind = 1:n
    
    for stn=1:length(stns)
        
        if isempty(C{1,stns(stn)}) == 0

            wbt = C{1,stn};

            % Find OMI where dates match that of station data   
            OMI_dates = datenum(double(OMI(:,1)),double(OMI(:,2)),double(OMI(:,3)));
            TW_dates = datenum(double(wbt(:,2)), double(wbt(:,3)), double(wbt(:,4)));
            [inter,iO,iW] = intersect(OMI_dates, TW_dates);

            wbt_table = horzcat(wbt(iW,:),OMI(iO,7:8));

            % Remove outliers
            outliers = find(wbt_table(:,5)>10);
            wbt_table(outliers,:) = [];
            wbt_nomissing = rmmissing(wbt_table);

            % Assign inactive MJO/BSISO days to phase 0
            std_amp = std(wbt_nomissing(:,5));
            inactive = find(wbt_nomissing(:,5) < std_amp);
            wbt_nomissing(inactive,6) = 0;
            
            % Select monthly data
            JA = wbt_nomissing(wbt_nomissing(:,3) == 7 | wbt_nomissing(:,3) == 8,:);
            MJ = wbt_nomissing(wbt_nomissing(:,3) == 5 | wbt_nomissing(:,3) == 6,:);

            % Calculate seasonal TW anomalies
            JA_mean = mean(JA(:,1));
            MJ_mean = mean(MJ(:,1));

            JA_anom = [(JA(:,1) - JA_mean) JA(:,2:6)];
            MJ_anom = [(MJ(:,1) - MJ_mean) MJ(:,2:6)];            

            % Reassign phase to dates at random, adding a new column to "active"
            phase_orig_JA = JA_anom(:,6);
            phase_rand_JA = phase_orig_JA(randperm(length(phase_orig_JA)));
            JA_rand = horzcat(JA_anom,phase_rand_JA);
            JA_rand = rmmissing(JA_rand);
            
            phase_orig_MJ = MJ_anom(:,6);
            phase_rand_MJ = phase_orig_MJ(randperm(length(phase_orig_MJ)));
            MJ_rand = horzcat(MJ_anom,phase_rand_MJ);
            MJ_rand = rmmissing(MJ_rand);           
            
            % Calculate average anomaly in each phase
            for phasei = 1:9

                mjo_JA = find(JA_rand(:,7) == phasei-1);
                OMI_scrambled(phasei,2,stn,ind) = nanmean(JA_rand(mjo_JA),1);

                mjo_MJ = find(MJ_rand(:,7) == phasei-1);
                OMI_scrambled(phasei,1,stn,ind) = nanmean(MJ_rand(mjo_MJ),1);


            end
        
        end

    end
    
    display(['# = ' num2str(ind)])
    
end
