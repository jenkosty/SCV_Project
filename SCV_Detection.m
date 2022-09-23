 %%% Loading MEOP data
% load("SealData_All.mat");

%%% Loading bathymetry data
RTOPO.lat = double(ncread('/Users/jenkosty/Research/detectSCV-main/RTOPO2.nc', 'lat'));
RTOPO.lon = double(ncread('/Users/jenkosty/Research/detectSCV-main/RTOPO2.nc', 'lon'))';
RTOPO.bedrock_topography = double(ncread('/Users/jenkosty/Research/detectSCV-main/RTOPO2.nc', 'bedrock_topography'))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME SERIES QUALITY CONTROL THRESHOLDS
qc.min_depth    = 350; % minimum depth to be considered a 'good' profile 
qc.min_profiles = 50;   % minimum number of 'good' profiles for a time series 
qc.max_time_gap = 5;   % max gap (days) between 'good' profiles before rejection

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCE PROFILE SETTINGS
ref_settings.inner_window = 2;
ref_settings.outer_window = 12; %%% number of days away from profile of interest
ref_settings.bathymetry_mean = 400;
ref_settings.bathymetry_std = 400;
ref_settings.no_profiles = 10;
ref_settings.bottom_depth = 350; %%% depth at which there must be at least x # profiles (bottom-up approach)
ref_settings.top_depth = 100; %%% depth at which there must still be at least x # profiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IQR CHECK SETTINGS
iqr_settings.density_levels = 2;
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
        raw_data.pres = SO_sealdata(i).PRES;
        raw_data.salt = SO_sealdata(i).SALT;
        raw_data.temp = SO_sealdata(i).TEMP;
        raw_data.pres_qc = SO_sealdata(i).PRES_QC;
        raw_data.salt_qc = SO_sealdata(i).SALT_QC;
        raw_data.temp_qc = SO_sealdata(i).TEMP_QC;
        SO_sealdata_qc(u).raw_data = raw_data;
        u = u + 1;
    end
end

clear time_qc pres_qc salt_qc temp_qc u SO_sealdata i x good_data

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
    interp_sealdata_qc(i) = struct('tag', SO_sealdata_qc(i).tag, 'lat',SO_sealdata_qc(i).lat,...
        'lon',SO_sealdata_qc(i).lon, 'time', SO_sealdata_qc(i).time, 'salt', interp_salt(i,:), ...
        'temp', interp_temp(i,:), 'raw_data', SO_sealdata_qc(i).raw_data);
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
        meop_ts(i).raw_data(j) = meop_data(tagidx(j)).raw_data;
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
        qc_ts(tscount).raw_data = meop_ts(i).raw_data(tmpidx{j});
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Variables (Pressure Space) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_prof = [63, 99, 199, 270, 318, 345, 418, 436, 82, 91, 116, 201];

for tag_no = test_prof
    
    %%% Calculating Bathymetry
    qc_ts(tag_no).bathymetry = interp2(RTOPO.lon, RTOPO.lat', RTOPO.bedrock_topography, qc_ts(tag_no).lon, qc_ts(tag_no).lat);

    %%% Creating pressure space structure
    tmp_pres_space.pres = depth_grid .* ones(size(qc_ts(tag_no).salt));
    
    %%% Assigning temperature and salinity data to new pressure space
    %%% structure
    tmp_pres_space.salt = qc_ts(tag_no).salt;
    tmp_pres_space.temp = qc_ts(tag_no).temp;
    
    %%% Calculating bottom pressure
    for i = 1:length(qc_ts(tag_no).cast)
        tmp_pres_space.prof_bot_pres(i) = max(tmp_pres_space.pres(~isnan(tmp_pres_space.salt(:,i)),i));     
    end

    %%% Calculating absolute salinity and conservative temperature
    tmp_pres_space.salt_absolute = gsw_SA_from_SP(tmp_pres_space.salt, depth_grid, qc_ts(tag_no).lon, qc_ts(tag_no).lat);
    tmp_pres_space.temp_conservative = gsw_CT_from_t(tmp_pres_space.salt_absolute, tmp_pres_space.temp, tmp_pres_space.pres);

    %%% Calculating density
    tmp_pres_space.density = gsw_rho_CT_exact(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative, ones(size(tmp_pres_space.salt)).*400);

    %%% Calculating potential density anomaly
    tmp_pres_space.sigma0 = gsw_sigma0(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative);
    
    %%% Calculating neutral density
    [tmp_pres_space.gamma_n,~,~] = eos80_legacy_gamma_n(tmp_pres_space.salt, tmp_pres_space.temp, tmp_pres_space.pres, qc_ts(tag_no).lon', qc_ts(tag_no).lat');

    %%% Calculating N^2
    [N2, mid_pres] = gsw_Nsquared(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative, tmp_pres_space.pres, qc_ts(tag_no).lat .* ones(size(tmp_pres_space.salt)));
    for i = 1:length(qc_ts(tag_no).cast)
        tmp_pres_space.N2(:,i) = interp1(mid_pres(:,i), N2(:,i), tmp_pres_space.pres(:,i));
    end

    %%% Calculating Spice
    tmp_pres_space.spice = gsw_spiciness0(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative);
    
    %%% Assigning variables to structure
    qc_ts(tag_no).ps = tmp_pres_space;

    clear ts_isopycnal_sep ts_pres u pres_final pres j isopycnal_sep isopycnal_sep_ds isopycnal_sep_y_axis...
            i a b density depths density_final k isopycnal_sep_final pres_ds ts_density tmp_pres_space 

end

%%% Removing previous salinity and temperature fields for
%%% organizational purposes
% qc_ts = rmfield(qc_ts, 'temp');
% qc_ts = rmfield(qc_ts, 'salt');

clear N2 mid_pres

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Checking for Density Inversions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = test_prof
    
    qc_ts(tag_no).rejected = strings(size(qc_ts(tag_no).cast));
    qc_ts(tag_no).ps.max_sigma0_inversion = NaN(size(qc_ts(tag_no).cast));
    qc_ts(tag_no).ps.max_gamma_n_inversion = NaN(size(qc_ts(tag_no).cast));
    
    for i =  1:length(qc_ts(tag_no).cast)
        
        %%% Checking to see if sigma0 profile is ascending with depth
        if issorted(qc_ts(tag_no).ps.sigma0(~isnan(qc_ts(tag_no).ps.sigma0(:,i)),i), 'ascend') == 0
            
            %%% If not, sorting the profile
            sigma0_orig = qc_ts(tag_no).ps.sigma0(:,i);
            idx = ~isnan(sigma0_orig);
            sigma0_sort_1 = sort(sigma0_orig, 'ascend');
            sigma0_sort_1(isnan(sigma0_sort_1)) = [];
            sigma0_sort = sigma0_orig;
            sigma0_sort(idx) = sigma0_sort_1;

            %%% Checking max difference between sorted and original density
            %%% profile. If the max difference it too large, the profile
            %%% will be rejected. If the max difference is small, the
            %%% profile will be replaced with the sorted version. 
         
            qc_ts(tag_no).ps.max_sigma0_inversion(:,i) = max(abs(sigma0_orig-sigma0_sort));
       
            if max(abs(sigma0_orig-sigma0_sort)) > 0.1
                qc_ts(tag_no).ps.sigma0(:,i) = NaN;
                qc_ts(tag_no).rejected(i) = "Density Inversion";
            else
                qc_ts(tag_no).ps.sigma0(:,i) = sigma0_sort;
            end
           
        end
        
        %%% Checking to see if gamma_n profile is ascending with depth
        if issorted(qc_ts(tag_no).ps.gamma_n(~isnan(qc_ts(tag_no).ps.gamma_n(:,i)),i), 'ascend') == 0
            
            %%% If not, sorting the profile
            gamma_n_orig = qc_ts(tag_no).ps.gamma_n(:,i);
            idx = ~isnan(gamma_n_orig);
            gamma_n_sort_1 = sort(gamma_n_orig, 'ascend');
            gamma_n_sort_1(isnan(gamma_n_sort_1)) = [];
            gamma_n_sort = gamma_n_orig;
            gamma_n_sort(idx) = gamma_n_sort_1;

            %%% Checking max difference between sorted and original density
            %%% profile. If the max difference it too large, the profile
            %%% will be rejected. If the max difference is small, the
            %%% profile will be replaced with the sorted version. 
         
            qc_ts(tag_no).ps.max_gamma_n_inversion(:,i) = max(abs(gamma_n_orig-gamma_n_sort));
       
            if max(abs(gamma_n_orig-gamma_n_sort)) > 0.1
                qc_ts(tag_no).ps.gamma_n(:,i) = NaN;
                qc_ts(tag_no).rejected(i) = "Density Inversion";
            else
                qc_ts(tag_no).ps.gamma_n(:,i) = gamma_n_sort;
            end
            
        end
    end
end

clear sigma0_orig sigma0_sort idx sigma0_sort_1 i gamma_n_orig gamma_n_sort idx gamma_n_sort_1 i

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolating to Density Grid %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Creating density grid
density_grid = (26.6:0.01:28)';

for tag_no = test_prof
    
    %%% Creating matrices for density-interpolated data
    interp_salt = NaN(length(density_grid), length(qc_ts(tag_no).cast));
    interp_temp = NaN(length(density_grid), length(qc_ts(tag_no).cast));
    interp_N2 = NaN(length(density_grid), length(qc_ts(tag_no).cast));
    interp_spice = NaN(length(density_grid), length(qc_ts(tag_no).cast));
    interp_pres = NaN(length(density_grid), length(qc_ts(tag_no).cast));
    
    for i = 1:length(qc_ts(tag_no).cast)
        
        %%% Removing NaNs
        tmp_sigma0 = qc_ts(tag_no).ps.sigma0(~isnan(qc_ts(tag_no).ps.sigma0(:,i)),i);
        tmp_salt = qc_ts(tag_no).ps.salt(~isnan(qc_ts(tag_no).ps.salt(:,i)),i);
        tmp_temp = qc_ts(tag_no).ps.temp(~isnan(qc_ts(tag_no).ps.temp(:,i)),i);
        tmp_N2 = qc_ts(tag_no).ps.N2(~isnan(qc_ts(tag_no).ps.N2(:,i)),i);
        tmp_sigma0_N2 = qc_ts(tag_no).ps.sigma0(~isnan(qc_ts(tag_no).ps.sigma0(:,i)) & ~isnan(qc_ts(tag_no).ps.N2(:,i)),i);
        tmp_spice = qc_ts(tag_no).ps.spice(~isnan(qc_ts(tag_no).ps.spice(:,i)),i);
        tmp_pres = depth_grid(~isnan(qc_ts(tag_no).ps.sigma0(:,i)));
        
        %%% Interpolating data
        interp_salt(:,i) = interp1(tmp_sigma0, tmp_salt, density_grid);
        interp_temp(:,i) = interp1(tmp_sigma0, tmp_temp, density_grid);
        interp_N2(:,i) = interp1(tmp_sigma0_N2, tmp_N2, density_grid);
        interp_spice(:,i) = interp1(tmp_sigma0, tmp_spice, density_grid);
        interp_pres(:,i) = interp1(tmp_sigma0, tmp_pres, density_grid);
    end
    
    %%% Saving density-interpolated data
    tmp_sigma0_space.salt = interp_salt;
    tmp_sigma0_space.temp = interp_temp;
    tmp_sigma0_space.N2 = interp_N2;
    tmp_sigma0_space.spice = interp_spice;
    tmp_sigma0_space.pres = interp_pres;
    
    %%% Saving structure
    qc_ts(tag_no).ds = tmp_sigma0_space;
    
    clear tmp_density tmp_salt tmp_temp tmp_N2 tmp_density_N2 tmp_spice tmp_pres ...
        interp_salt interp_temp interp_N2 interp_spice interp_pres tmp_density_space ...
        tmp_sigma0 tmp_sigma0_N2 tmp_sigma0_space i 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Isopycnal Separation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = test_prof
    
    isopycnal_separation = NaN(length(density_grid), length(qc_ts(tag_no).cast));
    
    for i = 1:length(qc_ts(tag_no).cast)
        for j = 1:length(density_grid)
            isopycnal_separation(j,i) = qc_ts(tag_no).ds.pres(min(j+1, length(density_grid)), i) - qc_ts(tag_no).ds.pres(max(j-1, 1), i);
        end
    end
    
    qc_ts(tag_no).ds.isopycnal_separation = isopycnal_separation;
    
    clear isopycnal_separation i j 
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Finding Indices to Build Reference Profiles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = test_prof
    
    mean_ind = cell(2,length(qc_ts(tag_no).cast));
    
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
    
    %%% Saving indices
    qc_ts(tag_no).ref_ind = mean_ind;
    
end

clear mean_ind i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Building Reference Profiles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = test_prof
    
    %%% Creating reference profiles for each of the time series profiles
    qc_ts(tag_no).ds.ref_salt = NaN(size(qc_ts(tag_no).ds.salt));
    qc_ts(tag_no).ds.ref_temp = NaN(size(qc_ts(tag_no).ds.temp));
    qc_ts(tag_no).ds.ref_N2 = NaN(size(qc_ts(tag_no).ds.N2));
    qc_ts(tag_no).ds.ref_spice = NaN(size(qc_ts(tag_no).ds.spice));
    qc_ts(tag_no).ds.ref_isopycnal_separation = NaN(size(qc_ts(tag_no).ds.isopycnal_separation));

    for i = 1:length(qc_ts(tag_no).cast)
    
        %%% Extracting the profiles to build the reference profile
        tmp_salt_ds = qc_ts(tag_no).ds.salt(:,[qc_ts(tag_no).ref_ind{1,i} qc_ts(tag_no).ref_ind{2,i}]);
        tmp_temp_ds = qc_ts(tag_no).ds.temp(:,[qc_ts(tag_no).ref_ind{1,i} qc_ts(tag_no).ref_ind{2,i}]);
        tmp_N2_ds = qc_ts(tag_no).ds.N2(:,[qc_ts(tag_no).ref_ind{1,i} qc_ts(tag_no).ref_ind{2,i}]);
        tmp_spice_ds = qc_ts(tag_no).ds.spice(:,[qc_ts(tag_no).ref_ind{1,i} qc_ts(tag_no).ref_ind{2,i}]);
        tmp_isopycnal_separation_ds = qc_ts(tag_no).ds.isopycnal_separation(:,[qc_ts(tag_no).ref_ind{1,i} qc_ts(tag_no).ref_ind{2,i}]);

        %%% Calculating a climatological value for each density level if at least
        %%% 75% of the potential reference profiles have data at said
        %%% level
        for j = 1:size(tmp_salt_ds, 1)
            tmp_salt_ds_level = tmp_salt_ds(j,~isnan(tmp_salt_ds(j,:)));
            tmp_temp_ds_level = tmp_temp_ds(j,~isnan(tmp_temp_ds(j,:)));
            tmp_spice_ds_level = tmp_spice_ds(j,~isnan(tmp_spice_ds(j,:)));
            if length(tmp_salt_ds_level) > 0.75*size(tmp_salt_ds,2)
                qc_ts(tag_no).ds.ref_salt(j,i) = median(tmp_salt_ds_level);
                qc_ts(tag_no).ds.ref_temp(j,i) = median(tmp_temp_ds_level);
                qc_ts(tag_no).ds.ref_spice(j,i) = median(tmp_spice_ds_level); 
            end
            
            tmp_N2_ds_level = tmp_N2_ds(j,~isnan(tmp_N2_ds(j,:)));
            if length(tmp_N2_ds_level) > 0.75*size(tmp_N2_ds,2)
                qc_ts(tag_no).ds.ref_N2(j,i) = median(tmp_N2_ds_level);  
            end
            
            tmp_isopycnal_separation_ds_level = tmp_isopycnal_separation_ds(j,~isnan(tmp_isopycnal_separation_ds(j,:)));
            if length(tmp_isopycnal_separation_ds_level) > 0.75*size(tmp_isopycnal_separation_ds,2)
                qc_ts(tag_no).ds.ref_isopycnal_separation(j,i) = median(tmp_isopycnal_separation_ds_level);  
            end
        end
        
        
    end

end

clear tmp_salt_ds tmp_temp_ds tmp_N2_ds tmp_spice_ds tmp_isopycnal_separation_ds tmp_salt_ds_level ...
    tmp_temp_ds_level tmp_spice_ds_level tmp_N2_ds_level tmp_isopycnal_separation_ds_level i j

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Anomalies %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = test_prof

    %%% Creating anomaly profiles for each of the time series profiles
    qc_ts(tag_no).ds.salt_anom = NaN(size(qc_ts(tag_no).ds.salt));
    qc_ts(tag_no).ds.temp_anom = NaN(size(qc_ts(tag_no).ds.temp));
    qc_ts(tag_no).ds.N2_anom = NaN(size(qc_ts(tag_no).ds.N2));
    qc_ts(tag_no).ds.spice_anom = NaN(size(qc_ts(tag_no).ds.spice));
    qc_ts(tag_no).ds.isopycnal_separation_anom = NaN(size(qc_ts(tag_no).ds.isopycnal_separation));

    for i = 1:length(qc_ts(tag_no).cast)

        %%% Calculating Anomalies
        qc_ts(tag_no).ds.salt_anom(:,i) = qc_ts(tag_no).ds.salt(:,i) - qc_ts(tag_no).ds.ref_salt(:,i);
        qc_ts(tag_no).ds.temp_anom(:,i) = qc_ts(tag_no).ds.temp(:,i) - qc_ts(tag_no).ds.ref_temp(:,i);
        qc_ts(tag_no).ds.N2_anom(:,i) = qc_ts(tag_no).ds.N2(:,i) - qc_ts(tag_no).ds.ref_N2(:,i);
        qc_ts(tag_no).ds.spice_anom(:,i) = qc_ts(tag_no).ds.spice(:,i) - qc_ts(tag_no).ds.ref_spice(:,i);
        qc_ts(tag_no).ds.isopycnal_separation_anom(:,i) = qc_ts(tag_no).ds.isopycnal_separation(:,i) - qc_ts(tag_no).ds.ref_isopycnal_separation(:,i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating IQR %%%
%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = test_prof

    %%% Creating iqr profiles for each of the time series profiles
    qc_ts(tag_no).ds.salt_anom_iqr = NaN(size(qc_ts(tag_no).ds.salt));
    qc_ts(tag_no).ds.temp_anom_iqr = NaN(size(qc_ts(tag_no).ds.temp));
    qc_ts(tag_no).ds.N2_anom_iqr = NaN(size(qc_ts(tag_no).ds.N2));
    qc_ts(tag_no).ds.spice_anom_iqr = NaN(size(qc_ts(tag_no).ds.spice));
    qc_ts(tag_no).ds.isopycnal_separation_anom_iqr = NaN(size(qc_ts(tag_no).ds.isopycnal_separation));

    %%% Creating lower iqr limit profiles for each of the time series profiles
    qc_ts(tag_no).ds.salt_anom_lim_lo = NaN(size(qc_ts(tag_no).ds.salt));
    qc_ts(tag_no).ds.temp_anom_lim_lo = NaN(size(qc_ts(tag_no).ds.temp));
    qc_ts(tag_no).ds.N2_anom_lim_lo = NaN(size(qc_ts(tag_no).ds.N2));
    qc_ts(tag_no).ds.spice_anom_lim_lo = NaN(size(qc_ts(tag_no).ds.spice));
    qc_ts(tag_no).ds.isopycnal_separation_anom_lim_lo = NaN(size(qc_ts(tag_no).ds.isopycnal_separation));

    %%% Creating upper iqr limit profiles for each of the time series profiles
    qc_ts(tag_no).ds.salt_anom_lim_hi = NaN(size(qc_ts(tag_no).ds.salt));
    qc_ts(tag_no).ds.temp_anom_lim_hi = NaN(size(qc_ts(tag_no).ds.temp));
    qc_ts(tag_no).ds.N2_anom_lim_hi = NaN(size(qc_ts(tag_no).ds.N2));
    qc_ts(tag_no).ds.spice_anom_lim_hi = NaN(size(qc_ts(tag_no).ds.spice));
    qc_ts(tag_no).ds.isopycnal_separation_anom_lim_hi = NaN(size(qc_ts(tag_no).ds.isopycnal_separation));

    for i = 1:length(qc_ts(tag_no).cast)

        %%% Extracting the assigned anomaly profiles
        tmp_salt_anom_ds = qc_ts(tag_no).ds.salt_anom(:,[qc_ts(tag_no).ref_ind{1,i} qc_ts(tag_no).ref_ind{2,i}]);
        tmp_salt_anom_ds(sum(double(~isnan(tmp_salt_anom_ds)), 2) < 0.75*size(tmp_salt_anom_ds,2),:) = NaN;
        tmp_temp_anom_ds = qc_ts(tag_no).ds.temp_anom(:,[qc_ts(tag_no).ref_ind{1,i} qc_ts(tag_no).ref_ind{2,i}]);
        tmp_temp_anom_ds(sum(double(~isnan(tmp_temp_anom_ds)), 2) < 0.75*size(tmp_temp_anom_ds,2),:) = NaN;
        tmp_spice_anom_ds = qc_ts(tag_no).ds.spice_anom(:,[qc_ts(tag_no).ref_ind{1,i} qc_ts(tag_no).ref_ind{2,i}]);
        tmp_spice_anom_ds(sum(double(~isnan(tmp_spice_anom_ds)), 2) < 0.75*size(tmp_spice_anom_ds,2),:) = NaN;
        tmp_N2_anom_ds = qc_ts(tag_no).ds.N2_anom(:,[qc_ts(tag_no).ref_ind{1,i} qc_ts(tag_no).ref_ind{2,i}]);
        tmp_N2_anom_ds(sum(double(~isnan(tmp_N2_anom_ds)), 2) < 0.75*size(tmp_N2_anom_ds,2),:) = NaN;
        tmp_isopycnal_separation_anom_ds = qc_ts(tag_no).ds.isopycnal_separation_anom(:,[qc_ts(tag_no).ref_ind{1,i} qc_ts(tag_no).ref_ind{2,i}]);
        tmp_isopycnal_separation_anom_ds(sum(double(~isnan(tmp_isopycnal_separation_anom_ds)), 2) < 0.75*size(tmp_isopycnal_separation_anom_ds,2),:) = NaN;
              
        %%% Calculating IQR
        qc_ts(tag_no).ds.salt_anom_iqr(:,i) = prctile(tmp_salt_anom_ds, 75, 2) - prctile(tmp_salt_anom_ds, 25, 2);
        qc_ts(tag_no).ds.temp_anom_iqr(:,i) = prctile(tmp_temp_anom_ds, 75, 2) - prctile(tmp_temp_anom_ds, 25, 2);
        qc_ts(tag_no).ds.spice_anom_iqr(:,i) = prctile(tmp_spice_anom_ds, 75, 2) - prctile(tmp_spice_anom_ds, 25, 2);
        qc_ts(tag_no).ds.N2_anom_iqr(:,i) = prctile(tmp_N2_anom_ds, 75, 2) - prctile(tmp_N2_anom_ds, 25, 2);
        qc_ts(tag_no).ds.isopycnal_separation_anom_iqr(:,i) = prctile(tmp_isopycnal_separation_anom_ds, 75, 2) - prctile(tmp_isopycnal_separation_anom_ds, 25, 2);
        
        %%% Calculating upper thresholds
        qc_ts(tag_no).ds.salt_anom_lim_hi(:,i) = prctile(tmp_salt_anom_ds, 75, 2) + 1.5*qc_ts(tag_no).ds.salt_anom_iqr(:,i);
        qc_ts(tag_no).ds.temp_anom_lim_hi(:,i) = prctile(tmp_temp_anom_ds, 75, 2) + 1.5*qc_ts(tag_no).ds.temp_anom_iqr(:,i);
        qc_ts(tag_no).ds.spice_anom_lim_hi(:,i) = prctile(tmp_spice_anom_ds, 75, 2) + 1.5*qc_ts(tag_no).ds.spice_anom_iqr(:,i);
        qc_ts(tag_no).ds.N2_anom_lim_hi(:,i) = prctile(tmp_N2_anom_ds, 75, 2) + 1.5*qc_ts(tag_no).ds.N2_anom_iqr(:,i);
        qc_ts(tag_no).ds.isopycnal_separation_anom_lim_hi(:,i) = prctile(tmp_isopycnal_separation_anom_ds, 75, 2) + 1.5*qc_ts(tag_no).ds.isopycnal_separation_anom_iqr(:,i);
        
        %%% Calculating lower thresholds
        qc_ts(tag_no).ds.salt_anom_lim_lo(:,i) = prctile(tmp_salt_anom_ds, 25, 2) - 1.5*qc_ts(tag_no).ds.salt_anom_iqr(:,i);
        qc_ts(tag_no).ds.temp_anom_lim_lo(:,i) = prctile(tmp_temp_anom_ds, 25, 2) - 1.5*qc_ts(tag_no).ds.temp_anom_iqr(:,i);
        qc_ts(tag_no).ds.spice_anom_lim_lo(:,i) = prctile(tmp_spice_anom_ds, 25, 2) - 1.5*qc_ts(tag_no).ds.spice_anom_iqr(:,i);
        qc_ts(tag_no).ds.N2_anom_lim_lo(:,i) = prctile(tmp_N2_anom_ds, 25, 2) - 1.5*qc_ts(tag_no).ds.N2_anom_iqr(:,i);
        qc_ts(tag_no).ds.isopycnal_separation_anom_lim_lo(:,i) = prctile(tmp_isopycnal_separation_anom_ds, 25, 2) - 1.5*qc_ts(tag_no).ds.isopycnal_separation_anom_iqr(:,i);
    end
end

clear tmp_salt_anom_ds tmp_temp_anom_ds tmp_spice_anom_ds tmp_N2_anom_ds tmp_isopycnal_separation_anom_ds ...
    i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%% IQR Check %%%
%%%%%%%%%%%%%%%%%

for tag_no = test_prof
    u = 1;
    uu = 1;
    cyclones.N2 = [];
    cyclones.isopycnal_separation = [];
    anticyclones.N2 = [];
    anticyclones.isopycnal_separation = [];
    
    for i = 1:length(qc_ts(tag_no).cast)

        %%% Checking isopycnal separation
        isopycnal_separation_check = qc_ts(tag_no).ds.isopycnal_separation_anom(:,i) > qc_ts(tag_no).ds.isopycnal_separation_anom_lim_hi(:,i);
        if sum(double(isopycnal_separation_check)) > iqr_settings.density_levels
            
            %%% Checking pressure levels
            pres_levels = qc_ts(tag_no).ds.pres(isopycnal_separation_check,i);
            
            %%% Checking number of continuous anomalies
            y = diff(find([0 double(isopycnal_separation_check(qc_ts(tag_no).ds.pres(:,i) > 50))' 0]==0))-1;
            y(y==0) = [];
            
            if max(pres_levels) < iqr_settings.min_pres
                qc_ts(tag_no).rejected(i) = "IQR Anomaly Too Shallow";
                continue 
            elseif min(pres_levels) > iqr_settings.max_pres
                qc_ts(tag_no).rejected(i) = "IQR Anomaly Too Deep";
                continue
            elseif max(pres_levels) - min(pres_levels) < 100
                qc_ts(tag_no).rejected(i) = "IQR Anomaly Not Tall Enough";
                continue
            elseif max(y) < 4
                qc_ts(tag_no).rejected(i) = "IQR Anomaly Does Not Cross Enough Density Surfaces";
                continue
            end
            
            anticyclones.isopycnal_separation(u) = i;
            u = u + 1; 
        end
        
        %%% Checking N^2
        N2_check = qc_ts(tag_no).ds.N2_anom(:,i) > qc_ts(tag_no).ds.N2_anom_lim_hi(:,i);
        if sum(double(N2_check)) > iqr_settings.density_levels
            
            %%% Checking pressure levels
            pres_levels = qc_ts(tag_no).ds.pres(N2_check,i);
            
            %%% Checking number of continuous anomalies
            y = diff(find([0 double(N2_check(qc_ts(tag_no).ds.pres(:,i) > 100))' 0]==0))-1;
            y(y==0) = [];
            
            if max(pres_levels) < iqr_settings.min_pres
                qc_ts(tag_no).rejected(i) = "IQR Anomaly Too Shallow";
                continue 
            elseif min(pres_levels) > iqr_settings.max_pres
                qc_ts(tag_no).rejected(i) = "IQR Anomaly Too Deep";
                continue
            elseif max(pres_levels) - min(pres_levels) < 100
                qc_ts(tag_no).rejected(i) = "IQR Anomaly Not Tall Enough";
                continue
            elseif max(y) < 4
                qc_ts(tag_no).rejected(i) = "IQR Anomaly Does Not Cross Enough Density Surfaces";
                continue
            end
            
            cyclones.N2(uu) = i;
            uu = uu + 1; 
        end
        
    end
    
    qc_ts(tag_no).cyclones = cyclones;
    qc_ts(tag_no).anticyclones = anticyclones;
end

clear u uu z zz pres_levels N2_lt N2_gt isopycnal_sep_lt isopycnal_sep_gt i cyclones anticyclones ...
    y isopycnal_separation_check N2_check i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gaussian Fit Check %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = test_prof
    
    u = 1;
    qc_ts(tag_no).anticyclones.gaussian = [];
    
    for i = qc_ts(tag_no).anticyclones.isopycnal_separation
        
        % Grab amplitude and depth of max spice anomaly
        spike.A  = max(qc_ts(tag_no).ds.spice_anom(:,i));
        spike.P = qc_ts(tag_no).ds.pres(find(qc_ts(tag_no).ds.spice_anom(:,i) == spike.A),i);
        
        % Get range of allowable parameters
        prng = [-0.2:0.05:0.2];  % allow pressure peak to vary between +- 20% of height
        arng  = [0.8:0.05:1.2];  % allow amplitude range of +- 20% of spice anomaly peak
        hrng  = [50:10:500];     % allow height to vary  between 50 and 850m
        
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
                    zo = double(qc_ts(tag_no).ds.pres(:,i) - [spike.P + p*(4)*sqrt(h^2/2)]);
                    sa = double(qc_ts(tag_no).ds.spice_anom(:,i));
                    
                    % Reduce to where data exists
                    dat      = sa + zo;
                    sa       = sa(~isnan(dat));
                    zo       = zo(~isnan(dat));
                    
                    % Generate gaussian model using updated amplitude, center, and height
                    gauss = (spike.A*a)*exp((-(zo.^2))/(h.^2));
                    
                    % Get gaussian limits for testing
                    pl = [spike.P + p*(4)*sqrt(h^2/2)] - 2*sqrt((h^2)/2); pl  = round(pl/10)*10;
                    ph = [spike.P + p*(4)*sqrt(h^2/2)] + 2*sqrt((h^2)/2); ph  = round(ph/10)*10;
                    
                    % Grab results
                    zp     = [zo + spike.P + p*(4)*sqrt(h^2/2)];
                    dataX  = qc_ts(tag_no).ds.spice_anom(pl <= qc_ts(tag_no).ds.pres(:,i) & qc_ts(tag_no).ds.pres(:,i) <= ph, i);
                    dataY  = qc_ts(tag_no).ds.pres(pl <= qc_ts(tag_no).ds.pres(:,i) & qc_ts(tag_no).ds.pres(:,i) <= ph, i);
                    dataY  = round(dataY, 6);
                    modelX = gauss(pl <= zp & zp <= ph);
                    modelY = zp(pl <= zp & zp <= ph);
                    modelY = round(modelY, 6);
                    
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
        zo    = double(qc_ts(tag_no).ds.pres(:,i) - [results.P]);
        zp    = zo + results.P;
        gauss = results.A*exp((-(zo.^2))/(results.H.^2));
        
        % Save final model
        results.X = gauss;
        results.Y = zp;
        
        if minlse < 0.05 && results.A > 0.1
            qc_ts(tag_no).anticyclones.gaussian(u) = i;
            u = u +1;
            
            qc_ts(tag_no).gauss_fit(i) = results;
        else 
            qc_ts(tag_no).rejected(i) = "Failed Gaussian Fit Test";
        end
        
    end
end

% figure()
% plot(qc_ts(tag_no).ds.spice_anom(:,i),qc_ts(tag_no).ds.pres(:,i),'k','linewidth',2)
% hold on; grid on; set(gca,'YDir','Reverse')
% plot(results.X,results.Y,'Color','r','LineWidth',3,'LineStyle','-.')
% xlabel('kg/m^3')
% ylabel('dbar');
% set(gca,'fontsize',10,'fontname','Helvetica')

clear a acnt arng b c dat dataX dataY gauss hcnt hrng idxlse ind lse minlse modelX ...
    modelY p pcnt ph pl prng sa zo zp u spike R2 i results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Dynamic Height Anomaly %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Dynamic Height Anomaly calculation
for tag_no = test_prof
    
    qc_ts(tag_no).ps.dha_bot_pres = NaN(size(qc_ts(tag_no).cast));
    
    for i = 1:length(qc_ts(tag_no).cast)
        ref_ind = [qc_ts(tag_no).ref_ind{1,i} qc_ts(tag_no).ref_ind{2,i}];
        
        max_pres = NaN(size(ref_ind));
        k = 1;
        
        for j = ref_ind
            max_pres(k) = max(qc_ts(tag_no).ps.pres(~isnan(qc_ts(tag_no).ps.salt(:,j)),j));
            k = k + 1;
        end
        
        qc_ts(tag_no).ps.dha_bot_pres(i) = min([min(max_pres), qc_ts(tag_no).ps.prof_bot_pres]);
        
        qc_ts(tag_no).ps.dyn_height_anom(:,i) = gsw_geo_strf_dyn_height(qc_ts(tag_no).ps.salt_absolute(:,i), qc_ts(tag_no).ps.temp_conservative(:,i), qc_ts(tag_no).ps.pres(:,i), 0); %qc_ts(tag_no).bot_pres(i));
    end
    
    clear max_pres ref_ind k i j
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Reference Profiles (Pressure Space) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = test_prof
    
    %%% Creating reference profiles for each of the time series profiles
    qc_ts(tag_no).ps.ref_N2 = NaN(size(qc_ts(tag_no).ps.N2));
    
    for i = 1:length(qc_ts(tag_no).cast)
        
        %%% Extracting the profiles to build the reference profile
        tmp_N2_ps = qc_ts(tag_no).ps.N2(:,[qc_ts(tag_no).ref_ind{1,i} qc_ts(tag_no).ref_ind{2,i}]);
        
        %%%
        for j = 1:size(tmp_N2_ps, 1)
            
            tmp_N2_ps_level = tmp_N2_ps(j,:);
            if length(tmp_N2_ps_level) > 0.75*size(tmp_N2_ps,2)
                qc_ts(tag_no).ps.ref_N2(j,i) = median(tmp_N2_ps_level);
            end
            
        end
    end
    
    clear tmp_N2_ps tmp_N2_ps_level
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dynamic Height Anomaly Check %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script to apply dynmodes routine to recover the vertical velocity and horizontal structure modes
opts1 = optimset('display','off','UseParallel',false);

for tag_no = test_prof
    
    qc_ts(tag_no).anticyclones.dha = [];
    u = 1;
    
    for i = qc_ts(tag_no).anticyclones.gaussian
        
        %%% Extract vertical velocity and horizontal structure modes of climatology
        [meop_profile(i).ref.wmodes, meop_profile(i).ref.pmodes, ~, ~] = dynmodes(qc_ts(tag_no).ps.ref_N2(~isnan(qc_ts(tag_no).ps.ref_N2(:,i)),i), qc_ts(tag_no).ps.pres(~isnan(qc_ts(tag_no).ps.ref_N2(:,i)),i),1);

        %%% Grab pressure levels of mode decomposition, add zero level (surface)
        meop_profile(i).ref.mode_pres = qc_ts(tag_no).ps.pres(~isnan(qc_ts(tag_no).ps.ref_N2(:,i)),i);
        meop_profile(i).ref.mode_pres = [0;meop_profile(i).ref.mode_pres];

        % Interpolate 1st baroclinic mode to pressure of SCV cast
        meop_profile(i).dyn_pres = qc_ts(tag_no).ps.pres(~isnan(qc_ts(tag_no).ps.dyn_height_anom(:,i)),i);
        meop_profile(i).dyn_height_anom = qc_ts(tag_no).ps.dyn_height_anom(~isnan(qc_ts(tag_no).ps.dyn_height_anom(:,i)),i);
        meop_profile(i).ref.BC1_data = meop_profile(i).ref.pmodes(:,1);
        meop_profile(i).ref.BC1_pres = meop_profile(i).ref.mode_pres;
        meop_profile(i).ref.BC1_data = interp1(meop_profile(i).ref.mode_pres(~isnan(meop_profile(i).ref.BC1_data)),meop_profile(i).ref.BC1_data(~isnan(meop_profile(i).ref.BC1_data)),meop_profile(i).dyn_pres);
        meop_profile(i).ref.BC1_pres = meop_profile(i).dyn_pres;

        % Create function that describes residuals between projected BC1 and dyn_height_anom
        % Exclude data inbetween SCV limits for better fit to first mode
        dat = [];
        dat = [meop_profile(i).ref.BC1_data + meop_profile(i).dyn_height_anom];
        x_o = [];
        x_o = meop_profile(i).ref.BC1_data(~isnan(dat));
        x_p = [];
        x_p = meop_profile(i).ref.BC1_pres(~isnan(dat));
        x_f = [];
        x_f = meop_profile(i).dyn_height_anom(~isnan(dat));

        % Get limits
        pl = qc_ts(tag_no).gauss_fit(i).Plow;
        ph = qc_ts(tag_no).gauss_fit(i).Phih;
        
        % Remove values between upper/lower limits of SCV to avoid bad fit
        ind = [];
        ind = find(pl < x_p & x_p < ph);
        if ind(end) == length(x_p)
            ind = ind(1:end-1);
        end
        x_o(ind) = [];
        x_f(ind) = [];
        %x_p(ind) = [];

        % Remove mixed layer depths (Lynne Talley method, first density greater than 0.03 from sfc value
        ind      = [];
        mld_dens = qc_ts(tag_no).ps.sigma0(~isnan(qc_ts(tag_no).ps.sigma0(:,i)),i);
        mld_pres = qc_ts(tag_no).ps.pres(~isnan(qc_ts(tag_no).ps.sigma0(:,i)),i);
        ind      = find(mld_dens > mld_dens(1)+0.03);
        mld_pres = mld_pres(ind(1));
        ind      = find(x_p < mld_pres);
        if length(ind) > length(x_o)
            continue
        end
        x_o(ind) = [];
        x_f(ind) = [];

        % f simply evaluates a given alpha (modal amplitude) and returns the
        % difference between the input DHanom profile and the projected 1st mode
        % We want to restrict our solutions such that the bottom of the projected
        % profile is equal to the bottom of the DHanom profile
        % SO let alpha2 = DHanom(end) - alpha*BT1(end)
        f = [];
        f = @(alpha) (alpha*x_o - x_f + (x_f(end) - alpha*x_o(end)));
        x0  = 0.05; % First guess

        % Solve for best modal amplitude
        alpha = [];
        alpha = lsqnonlin(f,x0,[-1],[1],opts1);

        % Redfine x_o and x_f with full profile
        x_o = meop_profile(i).ref.BC1_data(~isnan(dat));
        x_p = meop_profile(i).ref.BC1_pres(~isnan(dat)); %%% NOTE: Changed from meop_profile(i).ref.mode_pres - need to check with Danny
        x_f = meop_profile(i).dyn_height_anom(~isnan(dat));

        % Fix dynamic height anomaly by removing projected 1st mode, add back in barotopic mode
        meop_profile(i).dyn_height_anom_BC1 = [x_f] - [x_o*alpha + (x_f(end) - alpha*x_o(end))];
        meop_profile(i).dyn_height_pres_BC1 = meop_profile(i).dyn_pres(~isnan(dat));

        % Save VMD results
        meop_profile(i).ref.VMD.x_f      = x_f;
        meop_profile(i).ref.VMD.x_o      = x_o;
        meop_profile(i).ref.VMD.x_p      = x_p;
        meop_profile(i).ref.VMD.alpha    = alpha;

        % Get mode decomposition results
        BC1 = meop_profile(i).ref.VMD.x_o*meop_profile(i).ref.VMD.alpha;
        BC1 = BC1 - BC1(end); %// Set bottom to zero
        
        %%% Checking if positive maximum is achieved within anomaly height
        local_max = islocalmax(meop_profile(i).dyn_height_anom_BC1);
        
        for j = find(local_max)'
            if meop_profile(i).dyn_height_anom_BC1(j) > 0 && meop_profile(i).dyn_height_pres_BC1(j) > qc_ts(tag_no).gauss_fit(i).Plow ...
                    && meop_profile(i).dyn_height_pres_BC1(j) < qc_ts(tag_no).gauss_fit(i).Phih
                qc_ts(tag_no).anticyclones.dha(u) = i;
                u = u + 1;
                break
            else
                qc_ts(tag_no).rejected(i) = "Failed DHA Test";
            end
        end
           
        % Plot results
%         figure();
%         subplot(121)
%         plot(-meop_profile(i).ref.pmodes(:,1),meop_profile(i).ref.mode_pres,'r','linewidth',2)
%         hold on; grid on; set(gca,'YDir','Reverse')
%         plot(meop_profile(i).ref.pmodes(:,2),meop_profile(i).ref.mode_pres,'b','linewidth',2)
%         plot(meop_profile(i).ref.pmodes(:,3),meop_profile(i).ref.mode_pres,'g','linewidth',2)
%         plot(meop_profile(i).ref.pmodes(:,4),meop_profile(i).ref.mode_pres,'y','linewidth',2)
%         plot(meop_profile(i).ref.pmodes(:,5),meop_profile(i).ref.mode_pres,'color',[0.5 0 0.5],'linewidth',2)
%         title({'\it\bf\fontsize{8}\fontname{Helvetica}Horizontal Velocity','Modes'})
%         set(gca,'XTick',[0])
%         ylabel('Pressure (dbar)')
%         [l,~] = legend('Mode-1','Mode-2','Mode-3','Mode-4','Mode-5','location','southeast');
%         l.Box = 'off';
%         ylim([0 400]);
%         
%         subplot(122)
%         plot(meop_profile(i).dyn_height_anom,meop_profile(i).dyn_pres,'k','linewidth',2)
%         hold on; grid on; set(gca,'YDir','Reverse')
%         plot(BC1,meop_profile(i).ref.VMD.x_p,':r','linewidth',2)
%         plot(meop_profile(i).dyn_height_anom_BC1,meop_profile(i).dyn_height_pres_BC1,':k','linewidth',2)
%         legend('DH''_{orig}','BC1_{fit}','DH''_{adj}','location','southeast');
%         xlabel('m^2/s^2')
%         title({'\it\bf\fontsize{8}\fontname{Helvetica}Dynamic Height','Anomaly'})
%         set(gca,'YTickLabel',[])
%         ylim([0 400]);
 
    end
    
    %%% Saving BC1 fitting step
    qc_ts(tag_no).DHA_check = meop_profile;
   
end

    clear x0 x_f x_o x_p ref_orig pl ph mld_dens mld_pres l ind icons f dat BC1 argo alpha flaggedprofs_anticyclones...
        flaggedprofs_cyclones anticyclones pot_cyclones local_max j opts1 u

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting time series %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isopycnals = 0.05;

for tag_no = test_prof
    
    figure('Renderer', 'painters', 'Position', [0 0 1000 950])
    sgtitle('MEOP Seal ' + string(qc_ts(tag_no).tag), 'FontSize', 18, 'FontWeight', 'bold')
    
    %%% Bathymetry Subplot
    ax1 = subplot(4,1,1);
    hold on
    plot(datenum(qc_ts(tag_no).time), qc_ts(tag_no).bathymetry)
    for i = qc_ts(tag_no).pot_anticyclones.dha
        xline(datenum(qc_ts(tag_no).time(i,:)), 'r')
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
    pp = pcolor(unique(datenum(qc_ts(tag_no).time)),qc_ts(tag_no).ps.pres(:,1), qc_ts(tag_no).ps.temp(:,IB));
    set(pp, 'EdgeColor', 'none');
    [C,h] = contour(ax2, unique(datenum(qc_ts(tag_no).time)), depth_grid, qc_ts(tag_no).ps.density(:,IB), round(min(min(qc_ts(tag_no).ps.density)):isopycnals:max(max(qc_ts(tag_no).ps.density)), 2), 'k');
    clabel(C,h,'LabelSpacing',500);
    for i = qc_ts(tag_no).pot_anticyclones.dha
        xline(datenum(qc_ts(tag_no).time(i,:)), 'r')
    end
    hold off
    cmap = cmocean('thermal');
    colormap(ax2, cmap);
    colorbar;
    caxis([min(min(qc_ts(tag_no).ps.temp)) max(max(qc_ts(tag_no).ps.temp))])
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
    pp = pcolor(unique(datenum(qc_ts(tag_no).time)),qc_ts(tag_no).ps.pres(:,1), qc_ts(tag_no).ps.salt(:,IB));
    set(pp, 'EdgeColor', 'none');
    [C,h] = contour(ax3, unique(datenum(qc_ts(tag_no).time)), depth_grid, qc_ts(tag_no).ps.density(:,IB), round(min(min(qc_ts(tag_no).ps.density)):isopycnals:max(max(qc_ts(tag_no).ps.density)), 2), 'k');
    clabel(C,h,'LabelSpacing',500);
    for i = qc_ts(tag_no).pot_anticyclones.dha
        xline(datenum(qc_ts(tag_no).time(i,:)), 'r')
    end
    hold off
    cmap = cmocean('haline');
    colormap(ax3, cmap);
    colorbar;
    caxis([min(min(qc_ts(tag_no).ps.salt)) max(max(qc_ts(tag_no).ps.salt))])
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
    pp = pcolor(unique(datenum(qc_ts(tag_no).time)),density_grid, qc_ts(tag_no).ds.isopycnal_separation(:,IB));
    set(pp, 'EdgeColor', 'none');
    for i = qc_ts(tag_no).pot_anticyclones.dha
        xline(datenum(qc_ts(tag_no).time(i,:)), 'r')
    end
    hold off
    colorbar;
    caxis([min(min(qc_ts(tag_no).ds.isopycnal_separation)) max(max(qc_ts(tag_no).ds.isopycnal_separation))])
    set(gca, 'YDir','reverse');
    set(gca, 'Layer','top');
    xticks(linspace(datenum(qc_ts(tag_no).time(1,:)), datenum(qc_ts(tag_no).time(end,:)), (datenum(qc_ts(tag_no).time(end,:)) - datenum(qc_ts(tag_no).time(1,:))) / 5))
    datetick('x', 'mm/dd/yy', 'keepticks');
    ylabel('Density (kg/m^3)', 'FontSize', 12);
    ylim([27.2 27.9])
    title('Isopycnal Separation', 'FontSize', 12);
    
    p1 = get(ax1, 'Position');
    p2 = get(ax2, 'Position');
    p3 = get(ax3, 'Position');
    p4 = get(ax4, 'Position');
    
    p2(3) = p1(3);
    p3(3) = p1(3);
    p4(3) = p1(3);
    
    set(ax2, 'Position', p2);
    set(ax3, 'Position', p3);
    set(ax4, 'Position', p4);
end


clear ax1 ax2 ax3 ax4 ax5 C h cmap fig i h IB isopycnals p1 p2 p3 p4 p5 pp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting time series %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isopycnals = 0.05;

tag_no = 116;

figure('Renderer', 'painters', 'Position', [0 0 1000 950])
sgtitle('MEOP Seal ' + string(qc_ts(tag_no).tag), 'FontSize', 18, 'FontWeight', 'bold')

%%% Bathymetry Subplot
ax1 = subplot(4,1,1);
hold on
plot(datenum(qc_ts(tag_no).time), qc_ts(tag_no).bathymetry)
for i = qc_ts(tag_no).pot_anticyclones.N2
    xline(datenum(qc_ts(tag_no).time(i,:)), 'r')
end
for i = qc_ts(tag_no).pot_anticyclones.isopycnal_separation
    xline(datenum(qc_ts(tag_no).time(i,:)), 'r--')
end
for i = qc_ts(tag_no).pot_cyclones.N2
    xline(datenum(qc_ts(tag_no).time(i,:)), 'b')
end
for i = qc_ts(tag_no).pot_cyclones.isopycnal_separation
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
pp = pcolor(unique(datenum(qc_ts(tag_no).time)),qc_ts(tag_no).ps.pres(:,1), qc_ts(tag_no).ps.temp(:,IB));
set(pp, 'EdgeColor', 'none');
[C,h] = contour(ax2, unique(datenum(qc_ts(tag_no).time)), depth_grid, qc_ts(tag_no).ps.density(:,IB), round(min(min(qc_ts(tag_no).ps.density)):isopycnals:max(max(qc_ts(tag_no).ps.density)), 2), 'k');
clabel(C,h,'LabelSpacing',500);
for i = qc_ts(tag_no).pot_anticyclones.N2
    xline(datenum(qc_ts(tag_no).time(i,:)), 'r')
end
for i = qc_ts(tag_no).pot_anticyclones.isopycnal_separation
    xline(datenum(qc_ts(tag_no).time(i,:)), 'r--')
end
for i = qc_ts(tag_no).pot_cyclones.N2
    xline(datenum(qc_ts(tag_no).time(i,:)), 'b')
end
for i = qc_ts(tag_no).pot_cyclones.isopycnal_separation
    xline(datenum(qc_ts(tag_no).time(i,:)), 'b--')
end
hold off
cmap = cmocean('thermal');
colormap(ax2, cmap);
colorbar;
caxis([min(min(qc_ts(tag_no).ps.temp)) max(max(qc_ts(tag_no).ps.temp))])
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
pp = pcolor(unique(datenum(qc_ts(tag_no).time)),qc_ts(tag_no).ps.pres(:,1), qc_ts(tag_no).ps.salt(:,IB));
set(pp, 'EdgeColor', 'none');
[C,h] = contour(ax3, unique(datenum(qc_ts(tag_no).time)), depth_grid, qc_ts(tag_no).ps.density(:,IB), round(min(min(qc_ts(tag_no).ps.density)):isopycnals:max(max(qc_ts(tag_no).ps.density)), 2), 'k');
clabel(C,h,'LabelSpacing',500);
for i = qc_ts(tag_no).pot_anticyclones.N2
    xline(datenum(qc_ts(tag_no).time(i,:)), 'r')
end
for i = qc_ts(tag_no).pot_anticyclones.isopycnal_separation
    xline(datenum(qc_ts(tag_no).time(i,:)), 'r--')
end
for i = qc_ts(tag_no).pot_cyclones.N2
    xline(datenum(qc_ts(tag_no).time(i,:)), 'b')
end
for i = qc_ts(tag_no).pot_cyclones.isopycnal_separation
    xline(datenum(qc_ts(tag_no).time(i,:)), 'b--')
end
hold off
cmap = cmocean('haline');
colormap(ax3, cmap);
colorbar;
caxis([min(min(qc_ts(tag_no).ps.salt)) max(max(qc_ts(tag_no).ps.salt))])
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
pp = pcolor(unique(datenum(qc_ts(tag_no).time)),density_grid, qc_ts(tag_no).ds.isopycnal_separation(:,IB));
set(pp, 'EdgeColor', 'none');
for i = qc_ts(tag_no).pot_anticyclones.N2
    xline(datenum(qc_ts(tag_no).time(i,:)), 'r')
end
for i = qc_ts(tag_no).pot_anticyclones.isopycnal_separation
    xline(datenum(qc_ts(tag_no).time(i,:)), 'r--')
end
for i = qc_ts(tag_no).pot_cyclones.N2
    xline(datenum(qc_ts(tag_no).time(i,:)), 'b')
end
for i = qc_ts(tag_no).pot_cyclones.isopycnal_separation
    xline(datenum(qc_ts(tag_no).time(i,:)), 'b--')
end
hold off
colorbar;
caxis([min(min(qc_ts(tag_no).ds.isopycnal_separation)) max(max(qc_ts(tag_no).ds.isopycnal_separation))])
set(gca, 'YDir','reverse');
set(gca, 'Layer','top');
xticks(linspace(datenum(qc_ts(tag_no).time(1,:)), datenum(qc_ts(tag_no).time(end,:)), (datenum(qc_ts(tag_no).time(end,:)) - datenum(qc_ts(tag_no).time(1,:))) / 5))
datetick('x', 'mm/dd/yy', 'keepticks');
ylabel('Pressure (dbar)', 'FontSize', 12);
ylim([27.2 27.9])
title('Isopycnal Separation', 'FontSize', 12);

p1 = get(ax1, 'Position');
p2 = get(ax2, 'Position');
p3 = get(ax3, 'Position');
p4 = get(ax4, 'Position');

p2(3) = p1(3);
p3(3) = p1(3);
p4(3) = p1(3);

set(ax2, 'Position', p2);
set(ax3, 'Position', p3);
set(ax4, 'Position', p4);


clear ax1 ax2 ax3 ax4 ax5 C h cmap fig i h IB isopycnals p1 p2 p3 p4 p5 pp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for i = qc_ts(tag_no).pot_anticyclones.isopycnal_separation
    figure()
    hold on
    plot(qc_ts(tag_no).ds.isopycnal_separation_anom(:,i), qc_ts(tag_no).ds.pres(:,i), 'k')
    plot(qc_ts(tag_no).ds.isopycnal_separation_anom_lim_lo(:,i), qc_ts(tag_no).ds.pres(:,i), 'r')
    plot(qc_ts(tag_no).ds.isopycnal_separation_anom_lim_hi(:,i), qc_ts(tag_no).ds.pres(:,i), 'g')
    set(gca, 'YDir', 'reverse');
    
    pause 
   
    
    close all
    
end

%%
for i = qc_ts(tag_no).pot_cyclones.N2
    figure()
    hold on
    plot(qc_ts(tag_no).ds.N2_anom(:,i), qc_ts(tag_no).ds.pres(:,i), 'k')
    plot(qc_ts(tag_no).ds.N2_anom_lim_lo(:,i), qc_ts(tag_no).ds.pres(:,i), 'r')
    plot(qc_ts(tag_no).ds.N2_anom_lim_hi(:,i), qc_ts(tag_no).ds.pres(:,i), 'g')
    set(gca, 'YDir', 'reverse');
    
    pause 
    
    close all
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
tag_no = 63;
prof_no = 84;
isopycnals = 0.01;

fig = figure('Renderer', 'painters', 'Position', [0 0 1000 950]);
sgtitle('MEOP Seal ' + string(qc_ts(tag_no).tag), 'FontSize', 18, 'FontWeight', 'bold')

%%% Temperature Subplot
ax1 = subplot(4,3,1:3);
hold on
[~,IB] = unique(datenum(qc_ts(tag_no).time));
pp = pcolor(unique(datenum(qc_ts(tag_no).time)),qc_ts(tag_no).ps.pres(:,1), qc_ts(tag_no).temp(:,IB));
set(pp, 'EdgeColor', 'none');
[C,h] = contour(ax1, unique(datenum(qc_ts(tag_no).time)), depth_grid, qc_ts(tag_no).ps.sigma0(:,IB), round(min(min(qc_ts(tag_no).ps.sigma0)):isopycnals:max(max(qc_ts(tag_no).ps.sigma0)), 2), 'k');
clabel(C,h,'LabelSpacing',500);
xline(datenum(qc_ts(tag_no).time(prof_no,:)), 'r', 'LineWidth', 1.5)
hold off
cmap = cmocean('thermal');
colormap(ax1, cmap);
colorbar;
caxis([min(min(qc_ts(tag_no).temp)) max(max(qc_ts(tag_no).temp))])
set(gca, 'YDir','reverse');
set(gca, 'Layer','top');
xticks(linspace(datenum(qc_ts(tag_no).time(1,:)), datenum(qc_ts(tag_no).time(end,:)), (datenum(qc_ts(tag_no).time(end,:)) - datenum(qc_ts(tag_no).time(1,:))) / 5))
datetick('x', 'mm/dd/yy', 'keepticks');
ylabel('Pressure (dbar)', 'FontSize', 12);
ylim([0 500])
xlim([datenum(qc_ts(tag_no).time(max(prof_no-12, 1),:)) datenum(qc_ts(tag_no).time(prof_no+12,:))])
title('Temperature', 'FontSize', 12);

%%% Salinity Subplot
ax2 = subplot(4,3,4:6);
hold on
[~,IB] = unique(datenum(qc_ts(tag_no).time));
pp = pcolor(unique(datenum(qc_ts(tag_no).time)),qc_ts(tag_no).ps.pres(:,1), qc_ts(tag_no).salt(:,IB));
set(pp, 'EdgeColor', 'none');
[C,h] = contour(ax2, unique(datenum(qc_ts(tag_no).time)), depth_grid, qc_ts(tag_no).ps.sigma0(:,IB), round(min(min(qc_ts(tag_no).ps.sigma0)):isopycnals:max(max(qc_ts(tag_no).ps.sigma0)), 2), 'k');
clabel(C,h,'LabelSpacing',500);
xline(datenum(qc_ts(tag_no).time(prof_no,:)), 'r', 'LineWidth', 1.5)
hold off
cmap = cmocean('haline');
colormap(ax2, cmap);
colorbar;
caxis([min(min(qc_ts(tag_no).salt)) max(max(qc_ts(tag_no).salt))])
set(gca, 'YDir','reverse');
set(gca, 'Layer','top');
xticks(linspace(datenum(qc_ts(tag_no).time(1,:)), datenum(qc_ts(tag_no).time(end,:)), (datenum(qc_ts(tag_no).time(end,:)) - datenum(qc_ts(tag_no).time(1,:))) / 5))
datetick('x', 'mm/dd/yy', 'keepticks');
ylabel('Pressure (dbar)', 'FontSize', 12);
ylim([0 500])
title('Salinity', 'FontSize', 12);
xlim([datenum(qc_ts(tag_no).time(max(prof_no-12, 1),:)) datenum(qc_ts(tag_no).time(prof_no+12,:))])

%linkaxes([ax1 ax2]);

%%% Spice Profile
ax3 = subplot(4,3,7);
hold on
plot(qc_ts(tag_no).ds.spice(:,prof_no), qc_ts(tag_no).ds.pres(:,prof_no), 'r', 'DisplayName', 'Profile', 'LineWidth', 1.5)
plot(qc_ts(tag_no).ds.ref_spice(:,prof_no),qc_ts(tag_no).ds.pres(:,prof_no), 'k', 'DisplayName', 'Reference', 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse');
xlabel('Spice');
hold off
legend('Location', 'best')

%%% Isopycnal Separation Profile
ax4 = subplot(4,3,8);
hold on
plot(qc_ts(tag_no).ds.isopycnal_separation(:,prof_no), qc_ts(tag_no).ds.pres(:,prof_no), 'r', 'DisplayName', 'Profile', 'LineWidth', 1.5)
plot(qc_ts(tag_no).ds.ref_isopycnal_separation(:,prof_no), qc_ts(tag_no).ds.pres(:,prof_no), 'k','DisplayName', 'Reference', 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse');
xlabel('Isopycnal Separation');
hold off
legend('Location', 'best')

%%% N2 Profile
ax5 = subplot(4,3,9);
hold on
plot(qc_ts(tag_no).ds.N2(:,prof_no), qc_ts(tag_no).ds.pres(:,prof_no), 'r', 'DisplayName', 'Profile', 'LineWidth', 1.5)
plot(qc_ts(tag_no).ds.ref_N2(:,prof_no), qc_ts(tag_no).ds.pres(:,prof_no), 'k','DisplayName', 'Reference', 'LineWidth', 1.5)
set(gca, 'YDir', 'reverse');
xlabel('N^2');
hold off
legend('Location', 'best')

%%% Spice Anomaly Profile
ax6 = subplot(4,3,10);
hold on
plot(qc_ts(tag_no).ds.spice_anom(:,prof_no), qc_ts(tag_no).ds.pres(:,prof_no), 'b', 'DisplayName', 'Profile', 'LineWidth', 1.5)
xline(0, '--k', 'LineWidth', 1);
set(gca, 'YDir', 'reverse');
xlabel('Spice Anomaly');
hold off

%%% Isopycnal Separation Anomaly Profile
ax7 = subplot(4,3,11);
hold on
plot(qc_ts(tag_no).ds.isopycnal_separation_anom(:,prof_no), qc_ts(tag_no).ds.pres(:,prof_no), 'b', 'DisplayName', 'Profile', 'LineWidth', 1.5)
xline(0, '--k', 'LineWidth', 1);
set(gca, 'YDir', 'reverse');
xlabel('Isopycnal Separation Anomaly');
hold off

%%% N2 Anomaly Profile
ax8 = subplot(4,3,12);
hold on
plot(qc_ts(tag_no).ds.N2_anom(:,prof_no), qc_ts(tag_no).ds.pres(:,prof_no), 'b', 'DisplayName', 'Profile', 'LineWidth', 1.5)
xline(0, '--k', 'LineWidth', 1);
set(gca, 'YDir', 'reverse');
xlabel('N^2 Anomaly');
hold off

saveas(fig, '/Users/jenkosty/Research/SCV_Project/Figures/15Sept2022/' + string(qc_ts(tag_no).tag), 'png')
%linkaxes([ax3 ax4 ax5 ax6 ax7 ax8], 'y')

clear ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 C h cmap fig i h IB isopycnals p1 p2 p3 p4 p5 pp

