
RTOPO_lat = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'lat'));
RTOPO_lon = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'lon'))';
RTOPO_bedrock_topography = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'bedrock_topography'))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUALITY CONTROL THRESHOLDS
qc.min_depth    = [350]; % minimum depth to be considered a 'good' profile
qc.max_time_gap = [5];   % max gap (days) between 'good' profiles before rejection
qc.bathymetry_change = [1000]; % max bathymatry change between 'good' profiles before rejection
qc.min_profiles = [30];   % minimum number of 'good' profiles for a timeseries

dbar_grid = 1:500;
qc_ts_initial = [];
qc_ts = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA LOADING AND MANIPULATION

meop_data = interp_meop_data;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIRST QC FLAG (max depth) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Flag all casts that don't go deep enough
flagit = [];
for i = 1:length(meop_ts)
	for j = 1:length(meop_ts(i).cast)
		tmpdat = [meop_ts(i).temp(:,j) + meop_ts(i).salt(:,j)];
		ind = find(isnan(tmpdat)==0);
		if dbar_grid(ind(end)) < qc.min_depth
			flagit = [flagit; i j];
		end
	end
end

%%% Remove
for i = size(flagit,1):-1:1 % (NOTE: Flipped direction of deletion (issue with deleting correct profile)
	meop_ts(flagit(i,1)).cast(flagit(i,2))   = [];
	meop_ts(flagit(i,1)).lat(flagit(i,2))    = [];
	meop_ts(flagit(i,1)).lon(flagit(i,2))    = [];
	meop_ts(flagit(i,1)).time(flagit(i,2),:) = [];
	meop_ts(flagit(i,1)).salt(:,flagit(i,2)) = [];
	meop_ts(flagit(i,1)).temp(:,flagit(i,2)) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SECOND QC FLAG (time gap) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
		qc_ts_initial(tscount).tag  = [num2str(meop_ts(i).tag),tslist{j}];
		%%% Save as QC'd timeseries
		qc_ts_initial(tscount).cast = meop_ts(i).cast(tmpidx{j});
		qc_ts_initial(tscount).lat  = meop_ts(i).lat(tmpidx{j});
		qc_ts_initial(tscount).lon  = meop_ts(i).lon(tmpidx{j});
		qc_ts_initial(tscount).time = meop_ts(i).time(tmpidx{j},:);
		qc_ts_initial(tscount).temp = meop_ts(i).temp(:,tmpidx{j});
		qc_ts_initial(tscount).salt = meop_ts(i).salt(:,tmpidx{j});
		%%% Increase counter
		tscount = tscount + 1;
	end
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIRD QC FLAG (min profiles) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Find timeseries without enough profiles
Np = [];
for i = 1:length(qc_ts_initial)
	Np(i) = length(qc_ts_initial(i).cast);
end
ind = find(Np < qc.min_profiles);
qc_ts_initial(ind) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FOURTH QC FLAG (bathymetry change) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Find where gap between casts is too large
%%% Set up a tag modifier to separate one T-S into many
tslist  = {'A','B','C','D','E','F','G','H',...
		   'I','J','K','L','M','N','O','P',...
           'Q','R','S','T','U','V','W','X',...
           'Y','Z'};

%%% Find bathymetry depth at each profile
for i = 1:50
    qc_ts_initial(i).bathymetry = interp2(RTOPO_lon, RTOPO_lat', RTOPO_bedrock_topography, qc_ts_initial(i).lon, qc_ts_initial(i).lat);
end

tscount = 1;
for i = 1:50
    %%% Clear loop variables
    bathymetry_buffer = []; ts_bathymetry_diff = []; didx = []; ind = [];
    %%% Find difference in bathymetry between profiles
    ts_bathymetry_diff = diff(qc_ts_initial(i).bathymetry);
    %%% Find where difference greater than QC threshold
    didx = abs(ts_bathymetry_diff) < qc.bathymetry_change;
    %%% Find the gaps greater than QC threshold 
	ind = find(didx == 1); 
    %%% Clear start of time-series
	ts_start = []; ts_end = []; ts_cnt = 1; tmpidx = [];
	if isempty(ind);
		%%% No good groups found, go to next MEOP time-series
		continue
	elseif length(ind)==1
		%%% Only one group of consecutive casts exist
		tmpidx{ts_cnt} = [ind:ind+1];
	else
		%%% Many groups exist, separate them
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

    %%% Assign 'good' consecutive casts into separate time-series	
    for j = 1:length(tmpidx)
        qc_ts(tscount).tag  = [qc_ts_initial(i).tag, tslist{j}];
        %%% Save as QC'd timeseries
        qc_ts(tscount).cast = qc_ts_initial(i).cast(tmpidx{j});
        qc_ts(tscount).lat  = qc_ts_initial(i).lat(tmpidx{j});
        qc_ts(tscount).lon  = qc_ts_initial(i).lon(tmpidx{j});
        qc_ts(tscount).time = qc_ts_initial(i).time(tmpidx{j},:);
        qc_ts(tscount).temp = qc_ts_initial(i).temp(:,tmpidx{j});
        qc_ts(tscount).salt = qc_ts_initial(i).salt(:,tmpidx{j});
        qc_ts(tscount).bathymetry = qc_ts_initial(i).bathymetry(:,tmpidx{j});
        %%% Increase counter
        tscount = tscount + 1;
    end
end





