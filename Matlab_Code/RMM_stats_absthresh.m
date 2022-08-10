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

MJO_amplitude = netcdf.open('http://iridl.ldeo.columbia.edu/SOURCES/.BoM/.MJO/.RMM/.amplitude/dods','NC_NOWRITE');
MJO_phase = netcdf.open('http://iridl.ldeo.columbia.edu/SOURCES/.BoM/.MJO/.RMM/.phase/dods','NC_NOWRITE');
MJO_RMM1 = netcdf.open('http://iridl.ldeo.columbia.edu/SOURCES/.BoM/.MJO/.RMM/.RMM1/dods','NC_NOWRITE');
MJO_RMM2 = netcdf.open('http://iridl.ldeo.columbia.edu/SOURCES/.BoM/.MJO/.RMM/.RMM2/dods','NC_NOWRITE');

% Establish variables through December 31, 2019
time1 = netcdf.getVar(MJO_amplitude,0,0,16650);

amp = netcdf.getVar(MJO_amplitude,1,0,16650);
phase = netcdf.getVar(MJO_phase,1,0,16650);
rmm1 = netcdf.getVar(MJO_RMM1,1,0,16650);
rmm2 = netcdf.getVar(MJO_RMM2,1,0,16650);

% Determine monthly time series

time_serialnum = time1 - 1721059;
time_datetime = datetime(time_serialnum,'ConvertFrom','datenum');
time_datevec = datevec(time_datetime);

% Create matrix with MJO data
A = horzcat(rmm1, rmm2, amp, phase, time_datevec(:,1:3));
%% Run bootstrapping

% Use only stations we're interested in
stns = [coastal_ind; NW_ind; SE_ind];

n = 500;
RMM_scrambled_thresh = zeros(9,6,2,length(stns),n);
RMM_MJOdays_scrambled_thresh = zeros(9,2,length(stns),n);

for ind = 1:n
    
    for stn=1:length(stns)
        
        if isempty(C{1,stns(stn)}) == 0

        wbt = C{1,stns(stn)};

        % Find RMM where dates match that of station data   
        RMM_dates = datenum(double(A(:,5)),double(A(:,6)),double(A(:,7)));
        TW_dates = datenum(double(wbt(:,2)), double(wbt(:,3)), double(wbt(:,4)));
        [inter,iR,iW] = intersect(RMM_dates, TW_dates);

        wbt_table = horzcat(wbt(iW,:),A(iR,3:4));

        % Remove outliers
        outliers = find(wbt_table(:,5)>10);
        wbt_table(outliers,:) = [];
        wbt_nomissing = rmmissing(wbt_table);

        % Assign inactive MJO/BSISO days to phase 0
        std_amp = std(wbt_nomissing(:,5));
        inactive = find(wbt_nomissing(:,5) < std_amp);
        wbt_nomissing(inactive,6) = 0;
        
        for season = 1:2
            
            if season == 1
                
                % Select seasonal data
                season_index = find(wbt_nomissing(:,3) == 5 | wbt_nomissing(:,3) == 6);
                MJJA = wbt_nomissing(season_index,:);
    
                % Reassign phase to dates at random, adding a new column to "active"
                phase_orig = MJJA(:,6);
                phase_rand = phase_orig(randperm(length(phase_orig)));

                MJJA_rand = horzcat(MJJA,phase_rand);
        
                % Count number MJO days above threshold
                % in each month
                for thresh = 28:33
                    for phasei = 1:9

                        RMM_scrambled_thresh(phasei,thresh-27,season,stn,ind) = length(MJJA_rand(MJJA_rand(:,7) == phasei-1 & MJJA_rand(:,1) > thresh,1));

                    end
                end
        
                % Calculate number of MJO days in each phase for each station in
                % each month

                for phaseind = 1:9

                    RMM_MJOdays_scrambled_thresh(phaseind,season,stn,ind) = length(MJJA_rand(MJJA_rand(:,7) == phaseind-1,1));

                end
                
            end
            
            if season == 2
                
                % Select seasonal data
                season_index = find(wbt_nomissing(:,3) == 7 | wbt_nomissing(:,3) == 8);
                MJJA = wbt_nomissing(season_index,:);
    
                % Reassign phase to dates at random, adding a new column to "active"
                phase_orig = MJJA(:,6);
                phase_rand = phase_orig(randperm(length(phase_orig)));

                MJJA_rand = horzcat(MJJA,phase_rand);
        
                % Count number MJO days above threshold
                % in each month
                for thresh = 28:33
                    for phasei = 1:9

                        RMM_scrambled_thresh(phasei,thresh-27,season,stn,ind) = length(MJJA_rand(MJJA_rand(:,7) == phasei-1 & MJJA_rand(:,1) > thresh,1));

                    end
                end
        
                % Calculate number of MJO days in each phase for each station in
                % each month

                for phaseind = 1:9

                    RMM_MJOdays_scrambled_thresh(phaseind,season,stn,ind) = length(MJJA_rand(MJJA_rand(:,7) == phaseind-1,1));

                end
                
            end
            
        end
        


        end

    end
    
    display(['# = ' num2str(ind)])
    
end