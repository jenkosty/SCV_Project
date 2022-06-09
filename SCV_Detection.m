%%% Loading MEOP data
%load("SealData_All.mat");

%%% Loading bathymetry data
RTOPO.lat = double(ncread('/Users/jenkosty/Research/detectSCV-main/RTOPO2.nc', 'lat'));
RTOPO.lon = double(ncread('/Users/jenkosty/Research/detectSCV-main/RTOPO2.nc', 'lon'))';
RTOPO.bedrock_topography = double(ncread('/Users/jenkosty/Research/detectSCV-main/RTOPO2.nc', 'bedrock_topography'))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME SERIES QUALITY CONTROL THRESHOLDS
qc.min_depth    = 350; % minimum depth to be considered a 'good' profile 
qc.min_profiles = 50;   % minimum number of 'good' profiles for a timeseries 
qc.max_time_gap = 5;   % max gap (days) between 'good' profiles before rejection

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE CALCULATION SETTINGS
var.isopycnal_sep = 0.005;
var.den_grid_sep = 0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCE PROFILE SETTINGS
ref_settings.inner_window = 2;
ref_settings.outer_window = 15; %%% number of days away from profile of interest
ref_settings.bathymetry_mean = 400;
ref_settings.bathymetry_std = 400;
ref_settings.no_profiles = 15;
ref_settings.bottom_depth = 350; %%% depth at which there must be at least x # profiles (bottom-up approach)
ref_settings.top_depth = 100; %%% depth at which there must still be at least x # profiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IQR CHECK SETTINGS
iqr_settings.density_levels = 1;
iqr_settings.min_pres = 100;
iqr_settings.max_pres = 500;
iqr_settings.no_profiles = 15;

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

test_prof = 60:70;

for tag_no = test_prof

    %%% Creating pressure matrix
    qc_ts(tag_no).pres = depth_grid .* ones(size(qc_ts(tag_no).salt));

    %%% Calculating absolute salinity and conservative temperature
    qc_ts(tag_no).salt_absolute = gsw_SA_from_SP(qc_ts(tag_no).salt, depth_grid, qc_ts(tag_no).lon, qc_ts(tag_no).lat);
    qc_ts(tag_no).temp_conservative = gsw_CT_from_t(qc_ts(tag_no).salt_absolute, qc_ts(tag_no).temp, qc_ts(tag_no).pres);

    %%% Calculating density
    qc_ts(tag_no).density = gsw_rho(qc_ts(tag_no).salt_absolute, qc_ts(tag_no).temp_conservative, zeros(size(qc_ts(tag_no).salt)));

    %%% Calculating Potential Density Anomaly
    qc_ts(tag_no).sigma0 = gsw_sigma0(qc_ts(tag_no).salt_absolute, qc_ts(tag_no).temp_conservative);

    %%% Calculating N^2
    [N2, mid_pres] = gsw_Nsquared(qc_ts(tag_no).salt_absolute, qc_ts(tag_no).temp_conservative, qc_ts(tag_no).pres, qc_ts(tag_no).lat .* ones(size(qc_ts(tag_no).salt)));
    for i = 1:length(qc_ts(tag_no).cast)
        qc_ts(tag_no).N2(:,i) = interp1(mid_pres(:,i), N2(:,i), qc_ts(tag_no).pres(:,i));
    end

    %%% Calculating Spice
    qc_ts(tag_no).spice = gsw_spiciness0(qc_ts(tag_no).salt_absolute, qc_ts(tag_no).temp_conservative);

    %%% Calculating isopycnal separation
    for k = 1:length(qc_ts(tag_no).cast)
        %%% Creating y-axis for isopycnal separation calculation
        isopycnal_sep_y_axis = (1026:var.den_grid_sep:1028); %min(qc_ts(tag_no).density(:,k)):var.den_grid_sep:max(qc_ts(tag_no).density(:,k));

        %%% Density-space isopycnal separation calculation
        isopycnal_sep = NaN(length(isopycnal_sep_y_axis), 1);
        u = 1;
        for j = isopycnal_sep_y_axis
            isopycnal_sep(u) = length(find(qc_ts(tag_no).density(:,k) <= (j+var.isopycnal_sep) & qc_ts(tag_no).density(:,k) >= (j-var.isopycnal_sep)));
            isopycnal_sep(u) = isopycnal_sep(u) * mean(diff(depth_grid));
            u = u + 1;
        end

        %%% Calculating pressure along isopycnal separation y-axis
        pres = NaN(length(isopycnal_sep_y_axis), 1);
        density = qc_ts(tag_no).density(~isnan(qc_ts(tag_no).density(:,k)),k);
        depths = depth_grid(~isnan(qc_ts(tag_no).density(:,k)));
        [~,b] = unique(density);
        if length(b) > 2
            pres = interp1(density(b), depths(b), isopycnal_sep_y_axis');
        end

        %%% Calculating isopycnal separation in pressure-space
        isopycnal_sep_final = NaN(length(depth_grid), 1);
        pres_ds = pres(~isnan(pres));
        isopycnal_sep_ds = isopycnal_sep(~isnan(pres));
        [~,b] = unique(pres_ds);
        if length(b) > 2
            isopycnal_sep_final = interp1(pres_ds(b), isopycnal_sep_ds(b), depth_grid);
        end

        qc_ts(tag_no).isopycnal_sep(:,k) = isopycnal_sep_final;
    end

        clear ts_isopycnal_sep ts_pres u pres_final pres j isopycnal_sep isopycnal_sep_ds isopycnal_sep_y_axis...
            i a b density depths density_final k isopycnal_sep_final pres_ds ts_density

    %%% Calculating Bathymetry
    qc_ts(tag_no).bathymetry = interp2(RTOPO.lon, RTOPO.lat', RTOPO.bedrock_topography, qc_ts(tag_no).lon, qc_ts(tag_no).lat);

end

clear N2 mid_pres

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Building Reference Profiles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = test_prof
    
     mean_ind = cell(2,length(qc_ts(tag_no).cast));
%     prof_dates = datenum(qc_ts(tag_no).time);
%     u = 1;

    %%% Finding indices to use for background profile calculation
%     for i = prof_dates'
% 
%         mean_ind(1,u) = {find(prof_dates < i & prof_dates > i-ref_settings.outer_window)'};
%         mean_ind(2,u) = {find(prof_dates > i & prof_dates < i+ref_settings.outer_window)'};
% 
%         u = u + 1;
% 
%     end

    %%% Finding indices to use for background profile calculation
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

    %%% Creating reference profiles for each of the time series profiles
    qc_ts(tag_no).ref_salt = NaN(size(qc_ts(tag_no).salt));
    qc_ts(tag_no).ref_temp = NaN(size(qc_ts(tag_no).temp));
    qc_ts(tag_no).ref_N2 = NaN(size(qc_ts(tag_no).N2));
    qc_ts(tag_no).ref_spice = NaN(size(qc_ts(tag_no).spice));
    qc_ts(tag_no).ref_isopycnal_sep = NaN(size(qc_ts(tag_no).isopycnal_sep));

    for i = 1:length(qc_ts(tag_no).cast)
    
        %%% Extracting the profiles to build the reference profile
        tmp_density = qc_ts(tag_no).density(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_salt = qc_ts(tag_no).salt(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_temp = qc_ts(tag_no).temp(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_N2 = qc_ts(tag_no).N2(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_spice = qc_ts(tag_no).spice(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_isopycnal_sep = qc_ts(tag_no).isopycnal_sep(:,[mean_ind{1,i} mean_ind{2,i}]);

        %%% Interpolating the reference profiles to the POI's density grid
        clear interp_tmp_salt interp_tmp_temp interp_tmp_N2 interp_tmp_spice interp_tmp_isopycnal_sep

        for j = 1:size(tmp_salt,2)
            %%% Removing NaNs
            tmp_density_prof = tmp_density(~isnan(tmp_density(:,j)),j);
            tmp_salt_prof = tmp_salt(~isnan(tmp_salt(:,j)),j);
            tmp_temp_prof = tmp_temp(~isnan(tmp_temp(:,j)),j);
            tmp_N2_prof = tmp_N2(~isnan(tmp_density(:,j)) & ~isnan(tmp_N2(:,j)),j);
            tmp_density_prof_N2 = tmp_density(~isnan(tmp_density(:,j)) & ~isnan(tmp_N2(:,j)),j);
            tmp_spice_prof = tmp_spice(~isnan(tmp_spice(:,j)),j);
            tmp_isopycnal_sep_prof = tmp_isopycnal_sep(~isnan(tmp_isopycnal_sep(:,j)),j);
            tmp_density_prof_isopycnal_sep = tmp_density(~isnan(tmp_density(:,j)) & ~isnan(tmp_isopycnal_sep(:,j)),j);
            
            if isempty(tmp_isopycnal_sep_prof)
                continue
            end

            %%% Interpolating data
            interp_tmp_salt(:,j) = interp1(tmp_density_prof, tmp_salt_prof, qc_ts(tag_no).density(:,i));
            interp_tmp_temp(:,j) = interp1(tmp_density_prof, tmp_temp_prof, qc_ts(tag_no).density(:,i));
            interp_tmp_N2(:,j) = interp1(tmp_density_prof_N2, tmp_N2_prof, qc_ts(tag_no).density(:,i));
            interp_tmp_spice(:,j) = interp1(tmp_density_prof, tmp_spice_prof, qc_ts(tag_no).density(:,i));
            interp_tmp_isopycnal_sep(:,j) = interp1(tmp_density_prof_isopycnal_sep, tmp_isopycnal_sep_prof, qc_ts(tag_no).density(:,i));
        end

        %%% Checking data coverage at each depth level
        for j = length(depth_grid):-1:1
            good_prof(j,:) = {find(~isnan(interp_tmp_spice(j,:)) & ~isnan(interp_tmp_N2(j,:)) & ~isnan(interp_tmp_isopycnal_sep(j,:)))};
            length_check(j) = length(good_prof{j,1}) > ref_settings.no_profiles;
        end
        
        %%% Only building reference profiles if there is enough coverage
        if min(depth_grid(length_check)) < ref_settings.top_depth & max(depth_grid(length_check)) > ref_settings.bottom_depth %%% Checking min/max depths
            data = good_prof(length_check,:);
            depth_data = depth_grid(length_check);
            last_ind = find(depth_data < ref_settings.top_depth, 1, 'last');
            overlap = intersect(data{last_ind,1}, data{end,1}); %%% Checking number of consistent profiles over depth range
            if length(overlap) > ref_settings.no_profiles %%% Generating reference profiles for those that meet the requirements
                qc_ts(tag_no).ref_salt(:,i) = median(interp_tmp_salt(:,overlap),2);
                qc_ts(tag_no).ref_temp(:,i) = median(interp_tmp_temp(:,overlap),2);
                qc_ts(tag_no).ref_N2(:,i) = median(interp_tmp_N2(:,overlap),2);
                qc_ts(tag_no).ref_spice(:,i) = median(interp_tmp_spice(:,overlap),2);
                qc_ts(tag_no).ref_isopycnal_sep(:,i) = median(interp_tmp_isopycnal_sep(:,overlap),2);
            end
        end
    end

end

clear tmp_density_prof tmp_density tmp_salt_prof tmp_salt tmp_temp_prof tmp_temp...
    tmp_N2_prof tmp_N2 tmp_density_prof_N2 tmp_spice_prof tmp_spice interp_tmp_salt interp_tmp_temp...
    interp_tmp_N2 interp_tmp_spice j i density_grid tmp_bathymetry u prof_dates good_prof top_data bottom_data...
    depth_check  depth_data last_ind length_check overlap data tmp_isopycnal_sep_prof tmp_isopycnal_sep ...
    tmp_density_prof_isopycnal_sep interp_tmp_isopycnal_sep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Anomalies %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = test_prof

    %%% Creating anomaly profiles for each of the time series profiles
    qc_ts(tag_no).salt_anom = NaN(size(qc_ts(tag_no).salt));
    qc_ts(tag_no).temp_anom = NaN(size(qc_ts(tag_no).temp));
    qc_ts(tag_no).N2_anom = NaN(size(qc_ts(tag_no).N2));
    qc_ts(tag_no).spice_anom = NaN(size(qc_ts(tag_no).spice));
    qc_ts(tag_no).isopycnal_sep_anom = NaN(size(qc_ts(tag_no).isopycnal_sep));

    for i = 1:length(qc_ts(tag_no).cast)

        %%% Calculating Anomalies
        qc_ts(tag_no).salt_anom(:,i) = qc_ts(tag_no).salt(:,i) - qc_ts(tag_no).ref_salt(:,i);
        qc_ts(tag_no).temp_anom(:,i) = qc_ts(tag_no).temp(:,i) - qc_ts(tag_no).ref_temp(:,i);
        qc_ts(tag_no).N2_anom(:,i) = qc_ts(tag_no).N2(:,i) - qc_ts(tag_no).ref_N2(:,i);
        qc_ts(tag_no).spice_anom(:,i) = qc_ts(tag_no).spice(:,i) - qc_ts(tag_no).ref_spice(:,i);
        qc_ts(tag_no).isopycnal_sep_anom(:,i) = qc_ts(tag_no).isopycnal_sep(:,i) - qc_ts(tag_no).ref_isopycnal_sep(:,i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating IQR %%%
%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = test_prof

    mean_ind = cell(2,length(qc_ts(tag_no).cast));
% %     prof_dates = datenum(qc_ts(tag_no).time);
% %     u = 1;
% % 
% %     %%% Finding indices to use for background profile calculation
% %     for i = prof_dates'
% % 
% %         mean_ind(1,u) = {find(prof_dates < i & prof_dates > i-ref_settings.outer_window)'};
% %         mean_ind(2,u) = {find(prof_dates > i & prof_dates < i+ref_settings.outer_window)'};
% % 
% %         u = u + 1;
% % 
% %     end

    %%% Finding indices to use for background profile calculation
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
    qc_ts(tag_no).salt_anom_iqr = NaN(size(qc_ts(tag_no).salt));
    qc_ts(tag_no).temp_anom_iqr = NaN(size(qc_ts(tag_no).temp));
    qc_ts(tag_no).N2_anom_iqr = NaN(size(qc_ts(tag_no).N2));
    qc_ts(tag_no).spice_anom_iqr = NaN(size(qc_ts(tag_no).spice));
    qc_ts(tag_no).isopycnal_sep_anom_iqr = NaN(size(qc_ts(tag_no).isopycnal_sep));

    %%% Creating lower iqr limit profiles for each of the time series profiles
    qc_ts(tag_no).salt_anom_lim_lo = NaN(size(qc_ts(tag_no).salt));
    qc_ts(tag_no).temp_anom_lim_lo = NaN(size(qc_ts(tag_no).temp));
    qc_ts(tag_no).N2_anom_lim_lo = NaN(size(qc_ts(tag_no).N2));
    qc_ts(tag_no).spice_anom_lim_lo = NaN(size(qc_ts(tag_no).spice));
    qc_ts(tag_no).isopycnal_sep_anom_lim_lo = NaN(size(qc_ts(tag_no).isopycnal_sep));

    %%% Creating upper iqr limit profiles for each of the time series profiles
    qc_ts(tag_no).salt_anom_lim_hi = NaN(size(qc_ts(tag_no).salt));
    qc_ts(tag_no).temp_anom_lim_hi = NaN(size(qc_ts(tag_no).temp));
    qc_ts(tag_no).N2_anom_lim_hi = NaN(size(qc_ts(tag_no).N2));
    qc_ts(tag_no).spice_anom_lim_hi = NaN(size(qc_ts(tag_no).spice));
    qc_ts(tag_no).isopycnal_sep_anom_lim_hi = NaN(size(qc_ts(tag_no).isopycnal_sep));

    for i = 1:length(qc_ts(tag_no).cast)

        %%% Extracting the assigned anomaly profiles
        tmp_density = qc_ts(tag_no).density(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_salt_anom = qc_ts(tag_no).salt_anom(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_temp_anom = qc_ts(tag_no).temp_anom(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_N2_anom = qc_ts(tag_no).N2_anom(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_spice_anom = qc_ts(tag_no).spice_anom(:,[mean_ind{1,i} mean_ind{2,i}]);
        tmp_isopycnal_sep_anom = qc_ts(tag_no).isopycnal_sep_anom(:,[mean_ind{1,i} mean_ind{2,i}]);

        %%% Interpolating the anomaly profiles to the POI's density grid
        clear interp_tmp_salt_anom interp_tmp_temp_anom interp_tmp_N2_anom interp_tmp_spice_anom interp_tmp_isopycnal_sep_anom

        for j = 1:size(tmp_salt_anom,2)

            if length(tmp_salt_anom(~isnan(tmp_density(:,j)) & ~isnan(tmp_salt_anom(:,j)),j)) < 2 || length(tmp_N2_anom(~isnan(tmp_density(:,j)) & ~isnan(tmp_N2_anom(:,j)),j)) < 2 ...
                    || length(tmp_isopycnal_sep_anom(~isnan(tmp_density(:,j)) & ~isnan(tmp_isopycnal_sep_anom(:,j)),j)) < 2 
                interp_tmp_salt_anom(:,j) = NaN(size(qc_ts(tag_no).density(:,i)));
                interp_tmp_temp_anom(:,j) = NaN(size(qc_ts(tag_no).density(:,i)));
                interp_tmp_N2_anom(:,j) = NaN(size(qc_ts(tag_no).density(:,i)));
                interp_tmp_spice_anom(:,j) = NaN(size(qc_ts(tag_no).density(:,i)));
                interp_tmp_isopycnal_sep_anom(:,j) = NaN(size(qc_ts(tag_no).density(:,i)));

                continue
            end

            %%% Interpolating data
            interp_tmp_salt_anom(:,j) = interp1(tmp_density(~isnan(tmp_density(:,j)) & ~isnan(tmp_salt_anom(:,j)),j), tmp_salt_anom(~isnan(tmp_density(:,j)) & ~isnan(tmp_salt_anom(:,j)),j), qc_ts(tag_no).density(:,i));
            interp_tmp_temp_anom(:,j) = interp1(tmp_density(~isnan(tmp_density(:,j)) & ~isnan(tmp_temp_anom(:,j)),j), tmp_temp_anom(~isnan(tmp_density(:,j)) & ~isnan(tmp_temp_anom(:,j)),j), qc_ts(tag_no).density(:,i));
            interp_tmp_N2_anom(:,j) = interp1(tmp_density(~isnan(tmp_density(:,j)) & ~isnan(tmp_N2_anom(:,j)),j), tmp_N2_anom(~isnan(tmp_density(:,j)) & ~isnan(tmp_N2_anom(:,j)),j), qc_ts(tag_no).density(:,i));
            interp_tmp_spice_anom(:,j) = interp1(tmp_density(~isnan(tmp_density(:,j)) & ~isnan(tmp_spice_anom(:,j)),j), tmp_spice_anom(~isnan(tmp_density(:,j)) & ~isnan(tmp_spice_anom(:,j)),j), qc_ts(tag_no).density(:,i));
            interp_tmp_isopycnal_sep_anom(:,j) = interp1(tmp_density(~isnan(tmp_density(:,j)) & ~isnan(tmp_isopycnal_sep_anom(:,j)),j), tmp_isopycnal_sep_anom(~isnan(tmp_density(:,j)) & ~isnan(tmp_isopycnal_sep_anom(:,j)),j), qc_ts(tag_no).density(:,i));
        end

        %%% Checking data coverage at each depth level
        for j = length(depth_grid):-1:1
            good_prof(j,:) = {find(~isnan(interp_tmp_spice_anom(j,:)) & ~isnan(interp_tmp_N2_anom(j,:)) & ~isnan(interp_tmp_isopycnal_sep_anom(j,:)))};
            length_check(j) = length(good_prof{j,1}) > ref_settings.no_profiles;
        end
        
        %%% Only building reference profiles if there is enough coverage
        if min(depth_grid(length_check)) < ref_settings.top_depth & max(depth_grid(length_check)) > ref_settings.bottom_depth %%% Checking min/max depths
            data = good_prof(length_check,:);
            depth_data = depth_grid(length_check);
            last_ind = find(depth_data < ref_settings.top_depth, 1, 'last');
            overlap = intersect(data{last_ind,1}, data{end,1}); %%% Checking number of consistent profiles over depth range
            if length(overlap) > ref_settings.no_profiles %%% Generating reference profiles for those that meet the requirements
                %%% Calculating IQR
                qc_ts(tag_no).salt_anom_iqr(:,i) = prctile(interp_tmp_salt_anom(:,overlap), 75, 2) - prctile(interp_tmp_salt_anom(:,overlap), 25, 2);
                qc_ts(tag_no).temp_anom_iqr(:,i) = prctile(interp_tmp_temp_anom(:,overlap), 75, 2) - prctile(interp_tmp_temp_anom(:,overlap), 25, 2);
                qc_ts(tag_no).N2_anom_iqr(:,i) = prctile(interp_tmp_N2_anom(:,overlap), 75, 2) - prctile(interp_tmp_N2_anom(:,overlap), 25, 2);
                qc_ts(tag_no).spice_anom_iqr(:,i) = prctile(interp_tmp_spice_anom(:,overlap), 75, 2) - prctile(interp_tmp_spice_anom(:,overlap), 25, 2);
                qc_ts(tag_no).isopycnal_sep_anom_iqr(:,i) = prctile(interp_tmp_isopycnal_sep_anom(:,overlap), 75, 2) - prctile(interp_tmp_isopycnal_sep_anom(:,overlap), 25, 2);

                %%% Calculating upper thresholds
                qc_ts(tag_no).salt_anom_lim_hi(:,i) = prctile(interp_tmp_salt_anom(:,overlap), 75, 2) + 1.5*qc_ts(tag_no).salt_anom_iqr(:,i);
                qc_ts(tag_no).temp_anom_lim_hi(:,i) = prctile(interp_tmp_temp_anom(:,overlap), 75, 2) + 1.5*qc_ts(tag_no).temp_anom_iqr(:,i);
                qc_ts(tag_no).N2_anom_lim_hi(:,i) = prctile(interp_tmp_N2_anom(:,overlap), 75, 2) + 1.5*qc_ts(tag_no).N2_anom_iqr(:,i);
                qc_ts(tag_no).spice_anom_lim_hi(:,i) = prctile(interp_tmp_spice_anom(:,overlap), 75, 2) + 1.5*qc_ts(tag_no).spice_anom_iqr(:,i);
                qc_ts(tag_no).isopycnal_sep_anom_lim_hi(:,i) = prctile(interp_tmp_isopycnal_sep_anom(:,overlap), 75, 2) + 1.5*qc_ts(tag_no).isopycnal_sep_anom_iqr(:,i);

                %%% Calculating lower thresholds
                qc_ts(tag_no).salt_anom_lim_lo(:,i) = prctile(interp_tmp_salt_anom(:,overlap), 25, 2) - 1.5*qc_ts(tag_no).salt_anom_iqr(:,i);
                qc_ts(tag_no).temp_anom_lim_lo(:,i) = prctile(interp_tmp_temp_anom(:,overlap), 25, 2) - 1.5*qc_ts(tag_no).temp_anom_iqr(:,i);
                qc_ts(tag_no).N2_anom_lim_lo(:,i) = prctile(interp_tmp_N2_anom(:,overlap), 25, 2) - 1.5*qc_ts(tag_no).N2_anom_iqr(:,i);
                qc_ts(tag_no).spice_anom_lim_lo(:,i) = prctile(interp_tmp_spice_anom(:,overlap), 25, 2) - 1.5*qc_ts(tag_no).spice_anom_iqr(:,i);
                qc_ts(tag_no).isopycnal_sep_anom_lim_lo(:,i) = prctile(interp_tmp_isopycnal_sep_anom(:,overlap), 25, 2) - 1.5*qc_ts(tag_no).isopycnal_sep_anom_iqr(:,i);
            end
        end
    end
end

clear density_grid interp_tmp_spice_anom interp_tmp_salt_anom interp_tmp_temp_anom interp_tmp_N2_anom j i tmp_density tmp_density_prof tmp_salt_anom tmp_salt_anom_prof ...
    tmp_spice_anom tmp_spice_anom_prof tmp_temp_anom tmp_temp_anom_prof tmp_N2_anom tmp_N2_anom_prof interp_tmp_isopycnal_sep_anom tmp_isopycnal_sep_anom_prof...
    depth_data data good_prof last_ind mean_ind length_check overlap prof_dates tmp_isopycnal_sep_anom 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
%%% IQR Check %%%
%%%%%%%%%%%%%%%%%

for tag_no = test_prof
    z = 1;
    u = 1;
    zz = 1;
    uu = 1;
    qc_ts(tag_no).pot_cyclones_N2 = [];
    qc_ts(tag_no).pot_cyclones_isopycnal_sep = [];
    qc_ts(tag_no).pot_anticyclones_N2 = [];
    qc_ts(tag_no).pot_anticyclones_isopcyanal_sep = [];

    for i = 1:length(qc_ts(tag_no).cast)
        pres_levels = qc_ts(tag_no).pres(:,i) > iqr_settings.min_pres & qc_ts(tag_no).pres(:,i) < iqr_settings.max_pres;
        N2_lt = sum(double(qc_ts(tag_no).N2_anom(pres_levels,i) < qc_ts(tag_no).N2_anom_lim_lo(pres_levels,i)));
        N2_gt = sum(double(qc_ts(tag_no).N2_anom(pres_levels,i) > qc_ts(tag_no).N2_anom_lim_hi(pres_levels,i)));
        isopycnal_sep_lt = sum(double(qc_ts(tag_no).isopycnal_sep_anom(pres_levels,i) < qc_ts(tag_no).isopycnal_sep_anom_lim_lo(pres_levels,i)));
        isopycnal_sep_gt = sum(double(qc_ts(tag_no).isopycnal_sep_anom(pres_levels,i) > qc_ts(tag_no).isopycnal_sep_anom_lim_hi(pres_levels,i)));
        
        if N2_lt > iqr_settings.density_levels
            qc_ts(tag_no).pot_anticyclones_N2(z) = i;
            z = z + 1;
        elseif isopycnal_sep_gt > iqr_settings.density_levels
            qc_ts(tag_no).pot_anticyclones_isopycnal_sep(zz) = i;
            zz = zz + 1;
        elseif N2_gt > iqr_settings.density_levels 
            qc_ts(tag_no).pot_cyclones_N2(u) = i;
            u = u + 1;
        elseif isopycnal_sep_lt > iqr_settings.density_levels
            qc_ts(tag_no).pot_cyclones_isopycnal_sep(uu) = i;
            uu = uu + 1;
        end

    end
end

clear u uu z zz pres_levels N2_lt N2_gt isopycnal_sep_lt isopycnal_sep_gt i 

%%
isopycnals = 0.01;

tag_no = 63;

figure('Renderer', 'painters', 'Position', [0 0 1000 950])
sgtitle('MEOP Seal ' + string(qc_ts(tag_no).tag), 'FontSize', 18, 'FontWeight', 'bold')

ts_density = gsw_rho(qc_ts(tag_no).salt, qc_ts(tag_no).temp, zeros(size(qc_ts(tag_no).salt,1), size(qc_ts(tag_no).salt,2)));

%%% Bathymetry Subplot
ax1 = subplot(4,1,1);
hold on
plot(datenum(qc_ts(tag_no).time), qc_ts(tag_no).bathymetry)
for i = qc_ts(tag_no).pot_anticyclones_N2
    xline(datenum(qc_ts(tag_no).time(i,:)), 'r')
end
for i = qc_ts(tag_no).pot_anticyclones_isopycnal_sep
    xline(datenum(qc_ts(tag_no).time(i,:)), 'r--')
end
for i = qc_ts(tag_no).pot_cyclones_N2
    xline(datenum(qc_ts(tag_no).time(i,:)), 'b')
end
for i = qc_ts(tag_no).pot_cyclones_isopycnal_sep
    xline(datenum(qc_ts(tag_no).time(i,:)), 'b--')
end
hold off
xticks(linspace(datenum(qc_ts(tag_no).time(1,:)), datenum(qc_ts(tag_no).time(end,:)), (datenum(qc_ts(tag_no).time(end,:)) - datenum(qc_ts(tag_no).time(1,:))) / 5))
datetick('x', 'mm/dd/yy', 'keepticks');
ylabel('Depth (m)','FontSize',12);
title('Bedrock Topography','FontSize', 12);

%%% Temperature Subplot
ax2 = subplot(4,1,2);
hold on
[~,IB] = unique(datenum(qc_ts(tag_no).time));
pp = pcolor(unique(datenum(qc_ts(tag_no).time)),qc_ts(tag_no).pres(:,1), qc_ts(tag_no).temp(:,IB));
set(pp, 'EdgeColor', 'none');
[C,h] = contour(ax2, unique(datenum(qc_ts(tag_no).time)), depth_grid, ts_density(:,IB), round(min(min(qc_ts(tag_no).density)):isopycnals:max(max(qc_ts(tag_no).density)), 2), 'k');
clabel(C,h,'LabelSpacing',500);
for i = qc_ts(tag_no).pot_anticyclones_N2
    xline(datenum(qc_ts(tag_no).time(i,:)), 'r')
end
for i = qc_ts(tag_no).pot_anticyclones_isopycnal_sep
    xline(datenum(qc_ts(tag_no).time(i,:)), 'r--')
end
for i = qc_ts(tag_no).pot_cyclones_N2
    xline(datenum(qc_ts(tag_no).time(i,:)), 'b')
end
for i = qc_ts(tag_no).pot_cyclones_isopycnal_sep
    xline(datenum(qc_ts(tag_no).time(i,:)), 'b--')
end
hold off
cmap = cmocean('thermal');
colormap(ax2, cmap);
colorbar;
caxis([min(min(qc_ts(tag_no).temp)) max(max(qc_ts(tag_no).temp))])
set(gca, 'YDir','reverse');
set(gca, 'Layer','top');
xticks(linspace(datenum(qc_ts(tag_no).time(1,:)), datenum(qc_ts(tag_no).time(end,:)), (datenum(qc_ts(tag_no).time(end,:)) - datenum(qc_ts(tag_no).time(1,:))) / 5))
datetick('x', 'mm/dd/yy', 'keepticks');
ylabel('Pressure (dbar)', 'FontSize', 12);
ylim([0 500])
title('Temperature', 'FontSize', 12);

%%% Salinity Subplot
ax3 = subplot(4,1,3);
hold on
[~,IB] = unique(datenum(qc_ts(tag_no).time));
pp = pcolor(unique(datenum(qc_ts(tag_no).time)),qc_ts(tag_no).pres(:,1), qc_ts(tag_no).salt(:,IB));
set(pp, 'EdgeColor', 'none');
[C,h] = contour(ax3, unique(datenum(qc_ts(tag_no).time)), depth_grid, qc_ts(tag_no).density(:,IB), round(min(min(qc_ts(tag_no).density)):isopycnals:max(max(qc_ts(tag_no).density)), 2), 'k');
clabel(C,h,'LabelSpacing',500);
for i = qc_ts(tag_no).pot_anticyclones_N2
    xline(datenum(qc_ts(tag_no).time(i,:)), 'r')
end
for i = qc_ts(tag_no).pot_anticyclones_isopycnal_sep
    xline(datenum(qc_ts(tag_no).time(i,:)), 'r--')
end
for i = qc_ts(tag_no).pot_cyclones_N2
    xline(datenum(qc_ts(tag_no).time(i,:)), 'b')
end
for i = qc_ts(tag_no).pot_cyclones_isopycnal_sep
    xline(datenum(qc_ts(tag_no).time(i,:)), 'b--')
end
hold off
cmap = cmocean('haline');
colormap(ax3, cmap);
colorbar;
caxis([min(min(qc_ts(tag_no).salt)) max(max(qc_ts(tag_no).salt))])
set(gca, 'YDir','reverse');
set(gca, 'Layer','top');
xticks(linspace(datenum(qc_ts(tag_no).time(1,:)), datenum(qc_ts(tag_no).time(end,:)), (datenum(qc_ts(tag_no).time(end,:)) - datenum(qc_ts(tag_no).time(1,:))) / 5))
datetick('x', 'mm/dd/yy', 'keepticks');
ylabel('Pressure (dbar)', 'FontSize', 12);
ylim([0 500])
title('Salinity', 'FontSize', 12);

%%% Isopycnal Separation Subplot
ax4 = subplot(4,1,4);
hold on
[~,IB] = unique(datenum(qc_ts(tag_no).time));
pp = pcolor(unique(datenum(qc_ts(tag_no).time)),qc_ts(tag_no).pres(:,1), qc_ts(tag_no).isopycnal_sep(:,IB));
set(pp, 'EdgeColor', 'none');
[C,h] = contour(ax4, unique(datenum(qc_ts(tag_no).time)), depth_grid, qc_ts(tag_no).density(:,IB), round(min(min(qc_ts(tag_no).density)):isopycnals:max(max(qc_ts(tag_no).density)), 2), 'k');
clabel(C,h,'LabelSpacing',500);
for i = qc_ts(tag_no).pot_anticyclones_N2
    xline(datenum(qc_ts(tag_no).time(i,:)), 'r')
end
for i = qc_ts(tag_no).pot_anticyclones_isopycnal_sep
    xline(datenum(qc_ts(tag_no).time(i,:)), 'r--')
end
for i = qc_ts(tag_no).pot_cyclones_N2
    xline(datenum(qc_ts(tag_no).time(i,:)), 'b')
end
for i = qc_ts(tag_no).pot_cyclones_isopycnal_sep
    xline(datenum(qc_ts(tag_no).time(i,:)), 'b--')
end
xline(datenum(qc_ts(tag_no).time(84,:)), 'w')
hold off
colorbar;
caxis([min(min(qc_ts(tag_no).isopycnal_sep)) max(max(qc_ts(tag_no).isopycnal_sep))])
set(gca, 'YDir','reverse');
set(gca, 'Layer','top');
xticks(linspace(datenum(qc_ts(tag_no).time(1,:)), datenum(qc_ts(tag_no).time(end,:)), (datenum(qc_ts(tag_no).time(end,:)) - datenum(qc_ts(tag_no).time(1,:))) / 5))
datetick('x', 'mm/dd/yy', 'keepticks');
ylabel('Pressure (dbar)', 'FontSize', 12);
ylim([0 500])
title('Isopycnal Separation', 'FontSize', 12);

% %% N^2 Subplot
% ax5 = subplot(5,1,5);
% hold on
% [~,IB] = unique(datenum(qc_ts(tag_no).time));
% pp = pcolor(unique(datenum(qc_ts(tag_no).time)),qc_ts(tag_no).pres(:,1), qc_ts(tag_no).N2(:,IB));
% set(pp, 'EdgeColor', 'none');
% [C,h] = contour(ax5, unique(datenum(qc_ts(tag_no).time)), depth_grid, qc_ts(tag_no).density(:,IB), round(min(min(qc_ts(tag_no).density)):isopycnals:max(max(qc_ts(tag_no).density)), 2), 'k');
% clabel(C,h,'LabelSpacing',500);
% for i = qc_ts(tag_no).pot_anticyclones_N2
%     xline(datenum(qc_ts(tag_no).time(i,:)), 'r')
% end
% for i = qc_ts(tag_no).pot_anticyclones_isopycnal_sep
%     xline(datenum(qc_ts(tag_no).time(i,:)), 'r--')
% end
% for i = qc_ts(tag_no).pot_cyclones_N2
%     xline(datenum(qc_ts(tag_no).time(i,:)), 'b')
% end
% for i = qc_ts(tag_no).pot_cyclones_isopycnal_sep
%     xline(datenum(qc_ts(tag_no).time(i,:)), 'b--')
% end
% hold off
% cmap = cmocean('thermal');
% colormap(ax5, cmap);
% colorbar;
% caxis([min(min(qc_ts(tag_no).N2)) max(max(qc_ts(tag_no).N2))])
% set(gca, 'YDir','reverse');
% set(gca, 'Layer','top');
% xticks(linspace(datenum(qc_ts(tag_no).time(1,:)), datenum(qc_ts(tag_no).time(end,:)), (datenum(qc_ts(tag_no).time(end,:)) - datenum(qc_ts(tag_no).time(1,:))) / 5))
% datetick('x', 'mm/dd/yy', 'keepticks');
% ylabel('Pressure (dbar)', 'FontSize', 12);
% ylim([0 500])
% title('N^2', 'FontSize', 12);

p1 = get(ax1, 'Position');
p2 = get(ax2, 'Position');
p3 = get(ax3, 'Position');
p4 = get(ax4, 'Position');
%p5 = get(ax5, 'Position');

p2(3) = p1(3);
p3(3) = p1(3);
p4(3) = p1(3);
%p5(3) = p1(3);

set(ax2, 'Position', p2);
set(ax3, 'Position', p3);
set(ax4, 'Position', p4);
%set(ax5, 'Position', p5);

clear ax1 ax2 ax3 ax4 ax5 C h cmap fig i h IB isopycnals p1 p2 p3 p4 p5 pp
%%
for i = qc_ts(tag_no).pot_anticyclones_isopycnal_sep

    fig = figure();
    subplot(121)
    hold on
    plot(qc_ts(tag_no).isopycnal_sep(:,i), qc_ts(tag_no).pres(:,i), 'k', 'DisplayName', 'Profile')
    plot(qc_ts(tag_no).ref_isopycnal_sep(:,i), qc_ts(tag_no).pres(:,i), 'r', 'DisplayName', 'Reference');
    hold off
    xlabel('Isopycnal Separation');
    ylabel('Pressure')
    ylim([0 max(qc_ts(tag_no).pres(~isnan(qc_ts(tag_no).isopycnal_sep_anom(:,i)),i))])
    set(gca, 'YDir', 'reverse')
    legend();


    subplot(122)
    hold on
    plot(qc_ts(tag_no).isopycnal_sep_anom(:,i), qc_ts(tag_no).pres(:,i), 'k', 'DisplayName', 'Profile Anomaly');
    plot(qc_ts(tag_no).isopycnal_sep_anom_lim_lo(:,i), qc_ts(tag_no).pres(:,i), 'b', 'DisplayName', 'Lower Threshold');
    plot(qc_ts(tag_no).isopycnal_sep_anom_lim_hi(:,i), qc_ts(tag_no).pres(:,i), 'r', 'DisplayName', 'UpperThreshold');
    hold off
    xlabel('Isopycnal Separation Anomaly');
    ylabel('Pressure')
    ylim([0 max(qc_ts(tag_no).pres(~isnan(qc_ts(tag_no).isopycnal_sep_anom(:,i)),i))])
    set(gca, 'YDir', 'reverse')
    legend();

    pause

    close(fig)

end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gaussian Fit Check %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = test_prof

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


%%
%%

