%%% Loading MEOP data
load("SealData_All.mat");

%%% Loading bathymetry data
RTOPO.lat = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'lat'));
RTOPO.lon = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'lon'))';
RTOPO.bedrock_topography = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'bedrock_topography'))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUALITY CONTROL THRESHOLDS
qc.min_depth    = 350; % minimum depth to be considered a 'good' profile (350)
qc.min_profiles = 50;   % minimum number of 'good' profiles for a timeseries (100)
qc.max_time_gap = 5;   % max gap (days) between 'good' profiles before rejection (5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCE PROFILE SETTINGS
ref_settings.inner_window = 2;
ref_settings.outer_window = 11;
ref_settings.bathymetry_mean = 400;
ref_settings.bathymetry_std = 400;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filtering Out "Bad" Data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Identifying MEOP Profiles in the Southern Ocean (Below 60S)
idy = find([sealdata_all.LAT] < -60);
SO_sealdata = sealdata_all(idy);
clear idy

%%% Quality Controlling Southern Ocean MEOP Profiles (Identifying profiles
%%% with consistent temperature and salinity quality ratings of 1 - i.e. the best rating)

%%% Extracting rating data
time_qc = NaN(length(SO_sealdata),1);
pres_qc = cell(length(SO_sealdata),1);
salt_qc = cell(length(SO_sealdata),1);
temp_qc = cell(length(SO_sealdata),1);

for i = 1:length(SO_sealdata)
    time_qc(i) = str2double(SO_sealdata(i).TIME_QC);
    pres_qc{i,:} = str2num(SO_sealdata(i).PRES_QC);
    salt_qc{i,:} = str2num(SO_sealdata(i).SALT_QC);
    temp_qc{i,:} = str2num(SO_sealdata(i).TEMP_QC);
end

%%% Isolating profiles with the best data
for i = 1:length(SO_sealdata)
    if time_qc(i) == 1
        good_data{i,:} = find((pres_qc{i,1} == 1) & (salt_qc{i,1} == 1) & (temp_qc{i,1} == 1));
    else 
        good_data{i,:} = [];
    end
end

%%% Creating new structure with QC profiles
clear SO_sealdata_qc
u = 1;
for i = 1:length(SO_sealdata)
    if isempty(good_data{i,1}) == 0
        SO_sealdata_qc(u).lat = SO_sealdata(i).LAT;
        SO_sealdata_qc(u).lon = SO_sealdata(i).LON;
        SO_sealdata_qc(u).time = SO_sealdata(i).TIME;
        SO_sealdata_qc(u).pres = SO_sealdata(i).PRES(good_data{i,1});
        SO_sealdata_qc(u).salt = SO_sealdata(i).SALT(good_data{i,1});
        SO_sealdata_qc(u).temp = SO_sealdata(i).TEMP(good_data{i,1});
        SO_sealdata_qc(u).tag = str2double(SO_sealdata(i).TAG);
        u = u + 1;
    end
end

clear time_qc pres_qc salt_qc temp_qc u i x SO_sealdata good_data

%%% Extracting the max pressure that each profile reaches
for i = 1:length(SO_sealdata_qc)
    maxpres(i) = max(SO_sealdata_qc(i).pres);
end

%%% Removing profiles that don't go deep enough
flagit = (maxpres >= qc.min_depth);
SO_sealdata_qc = SO_sealdata_qc(flagit);

clear maxpres flagit i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolating to Uniform Pressure Grid %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear interp_sealdata_qc interp_salt interp_temp

depth_grid = 1:10:800;
depth_grid = depth_grid';

for i = 1:length(SO_sealdata_qc)
    interp_salt(i,:) = interp1(SO_sealdata_qc(i).pres, SO_sealdata_qc(i).salt, depth_grid);
    interp_temp(i,:) = interp1(SO_sealdata_qc(i).pres, SO_sealdata_qc(i).temp, depth_grid);
    interp_sealdata_qc(i) = struct('tag', SO_sealdata_qc(i).tag, 'lat',SO_sealdata_qc(i).lat,'lon',SO_sealdata_qc(i).lon, 'time', SO_sealdata_qc(i).time, 'salt', interp_salt(i,:), 'temp', interp_temp(i,:));
end

clear interp_salt interp_temp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Organizing Data into Time Series %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meop_data = interp_sealdata_qc;

%%% Grab MEOP meta data (tag num, date) and group them into time-series
for i = 1:length(meop_data)
    tagnum(i)   = meop_data(i).tag; % tag num
    profdate(i) = datenum(meop_data(i).time); % date
end

%%% Get unique tag numbers
taglist = unique(tagnum);

%%% Go through each tagnum and assign a cast number
meop_ts = [];
for i = 1:length(taglist)
	%%% Find profiles with that tagnum
	tagidx = find(tagnum == taglist(i));

	%%% Start meop time-series
	meop_ts(i).tag = taglist(i);

	%%% Sort the profiles by date
	[a,b] = sort(profdate(tagidx),'ascend');

	%%% Go through all instances of the same tag and fill in the data
	for j = 1:length(b)
		meop_ts(i).cast(j)   = j;						  % cast number
		meop_ts(i).lat(j)    = meop_data(tagidx(j)).lat;  % latitude
		meop_ts(i).lon(j)    = meop_data(tagidx(j)).lon;  % longitude
		meop_ts(i).time(j,:) = meop_data(tagidx(j)).time; % cast date
		meop_ts(i).salt(:,j) = meop_data(tagidx(j)).salt; % salinity
		meop_ts(i).temp(:,j) = meop_data(tagidx(j)).temp; % temperature
	end	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% QC FLAG (time gap) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Find where gap between casts is too large
%%% Set up a tag modifier to separate one T-S into many
tslist  = {'A','B','C','D','E','F','G','H',...
		   'I','J','K','L','M','N','O','P'};

%%% Set up a total time-series counter
tscount = 1;
for i = 1:length(meop_ts)
	% Clear loop variables
	ts_dates = []; didx = []; ind = []; tmpidx = [];
	% Find difference in dates between casts
	ts_dates = datenum(meop_ts(i).time);
	% Find where difference exceeds QC threshold (max time between casts)
	didx     = diff(ts_dates)<qc.max_time_gap;
	% Find the gaps greater than QC threshold 
	ind      = find(didx == 1); 
	% Clear start of time-series
	ts_start = []; ts_end = []; ts_cnt = 1;
	if isempty(ind);
		% No good groups found, go to next MEOP time-series
		continue
	elseif length(ind)==1
		% Only one group of consecutive casts exist
		tmpidx{ts_cnt} = [ind:ind+1];
	else
		% Many groups exist, separate them
		for j = 1:length(ind);
			if isempty(ts_start)
				ts_start = ind(j);
			end
			if j < length(ind) & ismember(ind(j)+1,ind);
				continue;
			else
				ts_end = ind(j)+1;
				tmpidx{ts_cnt} = [ts_start:ts_end];
				ts_start = []; ts_end = [];
				ts_cnt = ts_cnt + 1;
			end
		end
	end
	% Assign 'good' consecutive casts into separate time-series	
	for j = 1:length(tmpidx)
		qc_ts(tscount).tag  = [num2str(meop_ts(i).tag),tslist{j}];
		%%% Save as QC'd timeseries
		qc_ts(tscount).cast = meop_ts(i).cast(tmpidx{j});
		qc_ts(tscount).lat  = meop_ts(i).lat(tmpidx{j});
		qc_ts(tscount).lon  = meop_ts(i).lon(tmpidx{j});
		qc_ts(tscount).time = meop_ts(i).time(tmpidx{j},:);
		qc_ts(tscount).temp = meop_ts(i).temp(:,tmpidx{j});
		qc_ts(tscount).salt = meop_ts(i).salt(:,tmpidx{j});
		%%% Increase counter
		tscount = tscount + 1;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% QC FLAG (min profiles) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Find timeseries without enough profiles
Np = [];
for i = 1:length(qc_ts)
	Np(i) = length(qc_ts(i).cast);
end
ind = find(Np < qc.min_profiles);
qc_ts(ind) = [];

clear meop_data profdate tagidx taglist tagnum j i a b tmpdat tmpidx ts_cnt ...
    ts_end ts_start tscount tslist Np ind didx qc flagit ts_dates dbar_grid meop_ts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Additional Variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = 1:length(qc_ts)

    %%% Creating pressure matrix
    qc_ts(tag_no).pres = depth_grid .* ones(size(qc_ts(tag_no).salt));

    %%% Calculating absolute salinity and conservative temperature
    qc_ts(tag_no).salt_absolute = gsw_SA_from_SP(qc_ts(tag_no).salt, depth_grid, qc_ts(tag_no).lon, qc_ts(tag_no).lat);
    qc_ts(tag_no).temp_conservative = gsw_CT_from_t(qc_ts(tag_no).salt_absolute, qc_ts(tag_no).temp, qc_ts(tag_no).pres);

    %%% Calculating density
    qc_ts(tag_no).density = gsw_rho(qc_ts(tag_no).salt_absolute, qc_ts(tag_no).temp_conservative, qc_ts(tag_no).pres);

    %%% Calculating Potential Density Anomaly
    qc_ts(tag_no).sigma0 = gsw_sigma0(qc_ts(tag_no).salt_absolute, qc_ts(tag_no).temp_conservative);

    %%% Calculating N^2
    [qc_ts(tag_no).N2,~] = gsw_Nsquared(qc_ts(tag_no).salt_absolute, qc_ts(tag_no).temp_conservative, qc_ts(tag_no).pres, qc_ts(tag_no).lat .* ones(size(qc_ts(tag_no).salt)));

    %%% Calculating Spice
    qc_ts(tag_no).spice = gsw_spiciness0(qc_ts(tag_no).salt_absolute, qc_ts(tag_no).temp_conservative);

    %%% Calculating Bathymetry
    qc_ts(tag_no).bathymetry = interp2(RTOPO.lon, RTOPO.lat', RTOPO.bedrock_topography, qc_ts(tag_no).lon, qc_ts(tag_no).lat);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Building Reference Profiles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

for tag_no = 1:100
    
    %%% Finding indices to use for background profile calculation
    for i = 1:length(qc_ts(tag_no).cast)
        mean_ind(1,i) = {(i-ref_settings.outer_window):(i-ref_settings.inner_window)};
        mean_ind{1,i}(mean_ind{1,i} < 1) = []; %%% Making sure indices remain within ts boundaries
        mean_ind(2,i) = {(i+ref_settings.inner_window):(i+ref_settings.outer_window)};
        mean_ind{2,i}(mean_ind{2,i} > length(qc_ts(tag_no).cast)) = []; %%% Making sure indices remain within ts boundaries
    
%         %%% Only considering profiles that have a complete window on both
%         %%% sides
%         if (length(mean_ind{1,i}) + length(mean_ind{2,i})) < (2*ref_settings.outer_window - ref_settings.inner_window)
%             mean_ind{1,i} = [];
%             mean_ind{2,i} = [];
%         end

        %%% Creating special reference profiles for the start and end of
        %%% the time series
        if (length(mean_ind{1,i}) < ref_settings.outer_window - ref_settings.inner_window+1)
            mean_ind{2,i} = [mean_ind{2,i} max(mean_ind{2,i})+1:max(mean_ind{2,i})+(ref_settings.outer_window)-(ref_settings.inner_window-1)-length(mean_ind{1,i})];
        end

        if (length(mean_ind{2,i}) < ref_settings.outer_window - ref_settings.inner_window+1)
            mean_ind{1,i} = sort([mean_ind{1,i} min(mean_ind{1,i})-1:-1:min(mean_ind{1,i})-(ref_settings.outer_window)+(ref_settings.inner_window-1)+length(mean_ind{2,i})]);
        end

    end

    %%% Creating reference profiles for each of the time series profiles
    qc_ts(tag_no).ref_salt = cell(1,length(qc_ts(tag_no).cast));
    qc_ts(tag_no).ref_temp = cell(1,length(qc_ts(tag_no).cast));
    qc_ts(tag_no).ref_N2 = cell(1,length(qc_ts(tag_no).cast));
    qc_ts(tag_no).ref_spice = cell(1,length(qc_ts(tag_no).cast));

    for i = 1:length(qc_ts(tag_no).cast)

        %%% Removing NaNs from POI
        density_grid = qc_ts(tag_no).density(~isnan(qc_ts(tag_no).density(:,i)),i);
    
        %%% Extracting the assigned reference profiles
        tmp_density = qc_ts(tag_no).density(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_salt = qc_ts(tag_no).salt(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_temp = qc_ts(tag_no).temp(:,[mean_ind{1,i} mean_ind{2,i}]);
        %tmp_N2 = qc_ts(tag_no).N2(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_spice = qc_ts(tag_no).spice(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_bathymetry = qc_ts(tag_no).bathymetry(:,[mean_ind{1,i} mean_ind{2,i}]);

        if isempty(tmp_density)
            continue 
        end
        
%         %%%%%%%%%%%%%%%%%%%%%%%%
%         %%% Bathymetry Check %%%
%         %%%%%%%%%%%%%%%%%%%%%%%%
% 
%         %%% Checking each profile against assigned background profile
%         bathymetry_change = abs(qc_ts(tag_no).bathymetry(i) - mean(tmp_bathymetry));
% 
%         %%% Flagging profiles for removal based on bathymetry check
%         if ((bathymetry_change > ref_settings.bathymetry_mean) || (std(tmp_bathymetry) > ref_settings.bathymetry_mean))
%             qc_ts(tag_no).ref_salt{1,i} = [];
%             qc_ts(tag_no).ref_temp{1,i} = [];
%             qc_ts(tag_no).ref_N2{1,i} = [];
%             qc_ts(tag_no).ref_spice{1,i} = [];
%             continue
%         end

        %%% Interpolating the reference profiles to the POI's density grid
        clear interp_tmp_salt interp_tmp_temp interp_tmp_N2 interp_tmp_spice

        for j = 1:size(tmp_salt,2)

            %%% Removing NaNs
            tmp_density_prof = tmp_density(~isnan(tmp_density(:,j)),j);
            tmp_salt_prof = tmp_salt(~isnan(tmp_salt(:,j)),j);
            tmp_temp_prof = tmp_temp(~isnan(tmp_temp(:,j)),j);
            %tmp_N2_prof = tmp_N2(~isnan(tmp_N2(:,j)),j);
            tmp_spice_prof = tmp_spice(~isnan(tmp_spice(:,j)),j);

            %%% Interpolating data
            interp_tmp_salt(:,j) = interp1(tmp_density_prof, tmp_salt_prof, density_grid);
            interp_tmp_temp(:,j) = interp1(tmp_density_prof, tmp_temp_prof, density_grid);
            %interp_tmp_N2(:,j) = interp1(tmp_density_prof, tmp_N2_prof, density_grid);
            interp_tmp_spice(:,j) = interp1(tmp_density_prof, tmp_spice_prof, density_grid);
        end

        qc_ts(tag_no).ref_salt{1,i} = median(interp_tmp_salt,2);
        qc_ts(tag_no).ref_temp{1,i} = median(interp_tmp_temp,2);
        %qc_ts(tag_no).ref_N2{1,i} = median(interp_tmp_N2,2);
        qc_ts(tag_no).ref_spice{1,i} = median(interp_tmp_spice,2);

    end

end

clear tmp_density_prof tmp_density tmp_salt_prof tmp_salt tmp_temp_prof tmp_temp...
    tmp_N2_prof tmp_N2 tmp_spice_prof tmp_spice interp_tmp_salt interp_tmp_temp...
    interp_tmp_N2 interp_tmp_spice j i density_grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Anomalies %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = 1:100

    %%% Creating anomaly profiles for each of the time series profiles
    qc_ts(tag_no).salt_anom = cell(1,length(qc_ts(tag_no).cast));
    qc_ts(tag_no).temp_anom = cell(1,length(qc_ts(tag_no).cast));
    qc_ts(tag_no).N2_anom = cell(1,length(qc_ts(tag_no).cast));
    qc_ts(tag_no).spice_anom = cell(1,length(qc_ts(tag_no).cast));

    for i = 1:length(qc_ts(tag_no).cast)

        if isempty(qc_ts(tag_no).ref_salt{1,i})
            continue
        end

        %%% Removing NaNs from POI
        salt_prof = qc_ts(tag_no).salt(~isnan(qc_ts(tag_no).salt(:,i)),i);
        temp_prof = qc_ts(tag_no).temp(~isnan(qc_ts(tag_no).temp(:,i)),i);
        spice_prof = qc_ts(tag_no).spice(~isnan(qc_ts(tag_no).spice(:,i)),i);

        %%% Calculating Anomalies
        qc_ts(tag_no).salt_anom{1,i} = salt_prof - qc_ts(tag_no).ref_salt{1,i};
        qc_ts(tag_no).temp_anom{1,i} = temp_prof - qc_ts(tag_no).ref_temp{1,i};
        qc_ts(tag_no).spice_anom{1,i} = spice_prof - qc_ts(tag_no).ref_spice{1,i};
    end
end

clear salt_prof temp_prof spice_prof
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
%%% IQR Check %%%
%%%%%%%%%%%%%%%%%

for tag_no = 1

    %%% Finding indices to use for background profile calculation
    clear mean_ind
    for i = 1:length(qc_ts(tag_no).cast)
        mean_ind(1,i) = {(i-ref_settings.outer_window):(i-ref_settings.inner_window)};
        mean_ind{1,i}(mean_ind{1,i} < 1) = []; %%% Making sure indices remain within ts boundaries
        mean_ind(2,i) = {(i+ref_settings.inner_window):(i+ref_settings.outer_window)};
        mean_ind{2,i}(mean_ind{2,i} > length(qc_ts(tag_no).cast)) = []; %%% Making sure indices remain within ts boundaries

        %%% Creating special reference profiles for the start and end of
        %%% the time series
        if (length(mean_ind{1,i}) < ref_settings.outer_window - ref_settings.inner_window+1)
            mean_ind{2,i} = [mean_ind{2,i} max(mean_ind{2,i})+1:max(mean_ind{2,i})+(ref_settings.outer_window)-(ref_settings.inner_window-1)-length(mean_ind{1,i})];
        end

        if (length(mean_ind{2,i}) < ref_settings.outer_window - ref_settings.inner_window+1)
            mean_ind{1,i} = sort([mean_ind{1,i} min(mean_ind{1,i})-1:-1:min(mean_ind{1,i})-(ref_settings.outer_window)+(ref_settings.inner_window-1)+length(mean_ind{2,i})]);
        end

    end

    %%% Creating iqr profiles for each of the time series profiles
    qc_ts(tag_no).salt_iqr = cell(1,length(qc_ts(tag_no).cast));
    qc_ts(tag_no).temp_iqr = cell(1,length(qc_ts(tag_no).cast));
    qc_ts(tag_no).N2_iqr = cell(1,length(qc_ts(tag_no).cast));
    qc_ts(tag_no).spice_iqr = cell(1,length(qc_ts(tag_no).cast));

    for i = 1:length(qc_ts(tag_no).cast)

        %%% Removing NaNs from POI
        density_grid = qc_ts(tag_no).density(~isnan(qc_ts(tag_no).density(:,i)),i);

        %%% Extracting the assigned anomaly profiles
        tmp_density = qc_ts(tag_no).density(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_salt = qc_ts(tag_no).salt_anom(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_temp = qc_ts(tag_no).temp_anom(:,[mean_ind{1,i} mean_ind{2,i}]);
        %tmp_N2 = qc_ts(tag_no).N2_anom(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_spice = qc_ts(tag_no).spice_anom(:,[mean_ind{1,i} mean_ind{2,i}]);

        %%% Interpolating the anomaly profiles to the POI's density grid
        clear interp_tmp_salt interp_tmp_temp interp_tmp_N2 interp_tmp_spice

        for j = 1:size(tmp_salt,2)

            %%% Removing NaNs
            tmp_density_prof = tmp_density(~isnan(tmp_density(:,j)),j);
            tmp_salt_prof = tmp_salt{1,j};
            tmp_density_prof = tmp_density_prof(~isnan(tmp_salt_prof(:)));
            tmp_salt_prof = tmp_salt_prof(~isnan(tmp_salt_prof));
            tmp_temp_prof = tmp_temp{1,j};
            tmp_temp_prof = tmp_temp_prof(~isnan(tmp_temp_prof));
            %tmp_N2_prof = tmp_N2(~isnan(tmp_N2(:,j)),j);
            tmp_spice_prof = tmp_spice{1,j};
            tmp_spice_prof = tmp_spice_prof(~isnan(tmp_spice_prof));

            %%% Interpolating data
            interp_tmp_salt(:,j) = interp1(tmp_density_prof, tmp_salt_prof, density_grid);
            interp_tmp_temp(:,j) = interp1(tmp_density_prof, tmp_temp_prof, density_grid);
            %interp_tmp_N2(:,j) = interp1(tmp_density_prof, tmp_N2_prof, density_grid);
            interp_tmp_spice(:,j) = interp1(tmp_density_prof, tmp_spice_prof, density_grid);
        end

        %%% Calculating IQR
        qc_ts(tag_no).salt_iqr{1,i} = std(interp_tmp_salt,0,2);
        qc_ts(tag_no).temp_iqr{1,i} = std(interp_tmp_temp,0,2);
        %qc_ts(tag_no).N2_iqr{1,i} = iqr(tmp_N2,2);
        qc_ts(tag_no).spice_iqr{1,i} = std(interp_tmp_spice,0,2);
    end
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on
p = pcolor(datenum(qc_ts(tag_no).time),qc_ts(tag_no).pres, qc_ts(tag_no).spice);
set(p, 'EdgeColor', 'none');
set(gca, 'YDir','reverse');
set(gca, 'Layer','top');
center_date = qc_ts(tag_no).time(55,:);
xline(datenum(center_date), 'red', 'LineWidth', 2)
datetick('x', 'mm/dd');


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gaussian Fit Check %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = 1

    for i = 60
        ds_pres = qc_ts(tag_no).pres(~isnan(qc_ts(tag_no).density(:,i)),i);
        ds_spice_anom = qc_ts(tag_no).spice_anom{1,i};

        % Grab amplitude and depth of max spice anomaly
        spike.A  = max(ds_spice_anom);
        spike.P = ds_pres(find(ds_spice_anom == spike.A));

        % Get range of allowable parameters
        prng = [-0.2:0.05:0.2];  % allow pressure peak to vary between +- 20% of height
        arng  = [0.8:0.05:1.2];  % allow amplitude range of +- 20% of spice anomaly peak
        hrng  = [10:10:500];     % allow height to vary  between 50 and 850m

        % Set up matrix for least-squared error calculation
        lse = nan(length(prng),length(arng),length(hrng));

        % Go through all possible combinations
        hcnt = 0; % reset h counter
        for h = hrng
            hcnt = hcnt + 1; % increase 'h' counter
            acnt = 0;        % reset 'a' counter
            for a = arng
                acnt = acnt + 1; % increase 'a' counter
                pcnt = 0;        % reset 'p' counter
                for p = prng
                    pcnt = pcnt + 1; % increase 'p'

                    % Center Gaussian model around spike.P + p*h
                    zo = [];
                    zo = double(ds_pres - [spike.P + p*(4)*sqrt(h^2/2)]);
                    sa = double(ds_spice_anom);

                    % Reduce to where data exists
                    datcheck = sa + zo;
                    sa       = sa(~isnan(datcheck));
                    zo       = zo(~isnan(datcheck));

                    % Generate gaussian model using updated amplitude, center, and height
                    gauss = (spike.A*a)*exp((-(zo.^2))/(h.^2));

                    % Get gaussian limits for testing
                    pl = [spike.P + p*(4)*sqrt(h^2/2)] - 2*sqrt((h^2)/2); pl  = round(pl/10)*10;
                    ph = [spike.P + p*(4)*sqrt(h^2/2)] + 2*sqrt((h^2)/2); ph  = round(ph/10)*10;

                    % Grab results
                    zp     = [zo + spike.P + p*(4)*sqrt(h^2/2)];
                    dataX  = ds_spice_anom(pl <= ds_pres & ds_pres <= ph);
                    dataY  = ds_pres(pl <= ds_pres & ds_pres <= ph);
                    modelX = gauss(pl <= zp & zp <= ph);
                    modelY = zp(pl <= zp & zp <= ph);

                    % Check that depths of model and data intersect
                    if length(dataX) < length(modelX) | length(modelX) < length(dataX)
                        [c,~,~] = intersect(dataY,modelY);
                        ind     = find(min(c) <= dataY & dataY <= max(c));
                        dataX   = dataX(ind);   dataY = dataY(ind);
                        ind     = find(min(c) <= modelY & modelY <= max(c));
                        modelX  = modelX(ind); modelY = modelY(ind);
                    end

                    % Calculate R^2
                    R2(pcnt,acnt,hcnt) = corr2(dataX,modelX).^2;

                    % Save least-squared error results (ignore if bad R2 value (i.e. < 0.5))
                    if R2(pcnt,acnt,hcnt) < 0.5
                        lse(pcnt,acnt,hcnt) = NaN;
                    else
                        lse(pcnt,acnt,hcnt) = sum([dataX-modelX].^2);
                    end
                end
            end
        end

        % Find best zo,A,H combo according to lse
        [minlse,idxlse] = min(lse(:));
        [a,b,c] = ind2sub(size(lse),idxlse);

        % Update parameters
        results.A    = spike.A*arng(b);
        results.H    = hrng(c);
        results.P    = spike.P + prng(a)*(4)*sqrt(results.H^2/2);
        results.Plow = spike.P - 2*sqrt((results.H^2)/2);
        results.Phih = spike.P + 2*sqrt((results.H^2)/2);
        results.Plow = round(results.Plow/10)*10;
        results.Phih = round(results.Phih/10)*10;

        % Update zo,zp,gauss for final model
        zo    = double(ds_pres - [results.P]);
        zp    = zo + results.P;
        gauss = results.A*exp((-(zo.^2))/(results.H.^2));

        % Save final model
        results.X = gauss;
        results.Y = zp;

        % Plot
        figure()
        plot(ds_spice_anom,ds_pres,'k','linewidth',2)
        hold on; grid on; set(gca,'YDir','Reverse')
        plot(results.X,results.Y,'Color','Blue','LineWidth',3,'LineStyle','-.')
        xlabel('kg/m^3')
        ylabel('dbar');
        set(gca,'fontsize',10,'fontname','Helvetica')
    end
end

