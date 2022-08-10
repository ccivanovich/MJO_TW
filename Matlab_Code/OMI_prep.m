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

OMI_wbt_MJJA = zeros(9,2,length(stnelevs));
OMI_obs_len = zeros(9,2,length(stnelevs));
OMI_thresh_1std = zeros(9,6,2,length(stnelevs));
OMI_MJOdays = zeros(9,2,length(stnelevs));

for stn=1:size(stncodes,1)
        
    if isempty(C{1,stn}) == 0

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

                    OMI_thresh_1std(phasei,thresh-27,season, stn) = length(MJJA(MJJA(:,6) == phasei-1 & MJJA(:,1) > thresh,1));

                end
            end


            % Calculate average TW
            for phasei = 1:9

                mjo_MJJA = find(MJJA(:,6) == phasei-1);
                OMI_wbt_MJJA(phasei,season,stn) = mean(MJJA(mjo_MJJA),1);
                OMI_obs_len(phasei,season,stn) = length(mjo_MJJA);
                OMI_MJOdays(phasei,season,stn) = length(MJJA(MJJA(:,6) == phasei-1,1));

            end


        end
        
        if season == 2
            
            season_index = find(wbt_nomissing(:,3) == 7 | wbt_nomissing(:,3) == 8);
            MJJA = wbt_nomissing(season_index,:);

            % Count number MJO days and number MJO days above threshold
            % in each season
            for thresh = 28:33
                for phasei = 1:9

                    OMI_thresh_1std(phasei,thresh-27,season, stn) = length(MJJA(MJJA(:,6) == phasei-1 & MJJA(:,1) > thresh,1));

                end
            end


            % Calculate average TW
            for phasei = 1:9

                mjo_MJJA = find(MJJA(:,6) == phasei-1);
                OMI_wbt_MJJA(phasei,season,stn) = mean(MJJA(mjo_MJJA),1);
                OMI_obs_len(phasei,season,stn) = length(mjo_MJJA);
                OMI_MJOdays(phasei,season,stn) = length(MJJA(MJJA(:,6) == phasei-1,1));


            end
            
        end
        
    end
    
    end
    

        display(['stn # = ' num2str(stn)])

end
