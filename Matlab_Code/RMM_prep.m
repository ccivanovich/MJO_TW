%% Determine phase of each active BSISO day

addpath /home/ivanov/matlab/Code/MJO_TW

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

%% WBT Data Prep

load('stndata.mat')

[stncodes,stnlats,stnlons,stnelevs]=...
        textread('hadisd_station_metadata_v2019.txt','%12c %7f %8f %6f'); %from metadata file 

slp=1010;
datadir = '/Users/Casey/Columbia/Research/MJO_WetBulb/Data/station_data/all_stations/';
curnumstns = size(stncodes,1);

finalstnelev=NaN.*ones(curnumstns,1);

for stn=1:size(stncodes,1)

        if exist(strcat(datadir,'hadisd.3.1.0.2019f_19310101-20200101_',stncodes(stn,:),'.nc'),'file')==2

            finalstnelev(stn)=ncread(strcat(datadir,'hadisd.3.1.0.2019f_19310101-20200101_',stncodes(stn,:),'.nc'),'elevation');

        end
        
        display(['elevation # = ' num2str(stn)])
        
end

%% Calculate Wet Bulb Temperature

load('TW_allstations.mat')

% Compute WBT using Davies-Jones method
% Note that in each file, time is given as hours since 1931-01-01 00:00

RMM_wbt_MJJA = zeros(9,2,length(stnelevs));
RMM_obs_len = zeros(9,2,length(stnelevs));
RMM_thresh_1std = zeros(9,6,2,length(stnelevs));
RMM_MJOdays = zeros(9,2,length(stnelevs));

for stn=1:size(stncodes,1)
        
    if isempty(C{1,stn}) == 0

    wbt = C{1,stn};
    
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
    
    % Calculate average WBT for each MJO phase (seasonally)
    % Define indices for each season
    
     for season = 1:2
        
        if season == 1 % MJ season
    
            season_index = find(wbt_nomissing(:,3) == 5 | wbt_nomissing(:,3) == 6);
            MJJA = wbt_nomissing(season_index,:);

            % Count number MJO days and number MJO days above threshold
            % in each season
            for thresh = 28:33
                for phasei = 1:9

                    RMM_thresh_1std(phasei,thresh-27,season, stn) = length(MJJA(MJJA(:,6) == phasei-1 & MJJA(:,1) > thresh,1));

                end
            end


            % Calculate average TW
            for phasei = 1:9

                mjo_MJJA = find(MJJA(:,6) == phasei-1);
                RMM_wbt_MJJA(phasei,season,stn) = mean(MJJA(mjo_MJJA),1);
                RMM_obs_len(phasei,season,stn) = length(mjo_MJJA);
                RMM_MJOdays(phasei,season,stn) = length(MJJA(MJJA(:,6) == phasei-1,1));

            end


        end
        
        if season == 2
            
            season_index = find(wbt_nomissing(:,3) == 7 | wbt_nomissing(:,3) == 8);
            MJJA = wbt_nomissing(season_index,:);

            % Count number MJO days and number MJO days above threshold
            % in each season
            for thresh = 28:33
                for phasei = 1:9

                    RMM_thresh_1std(phasei,thresh-27,season, stn) = length(MJJA(MJJA(:,6) == phasei-1 & MJJA(:,1) > thresh,1));

                end
            end


            % Calculate average TW
            for phasei = 1:9

                mjo_MJJA = find(MJJA(:,6) == phasei-1);
                RMM_wbt_MJJA(phasei,season,stn) = mean(MJJA(mjo_MJJA),1);
                RMM_obs_len(phasei,season,stn) = length(mjo_MJJA);
                RMM_MJOdays(phasei,season,stn) = length(MJJA(MJJA(:,6) == phasei-1,1));


            end
            
        end
        
    end
    
    end
    

        display(['stn # = ' num2str(stn)])

end