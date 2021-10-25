%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Algorithm to Detect MEOP SCVs - Following Methods of Zhang et al %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUALITY CONTROL THRESHOLDS
qc.min_length    = [100]; % minimum density to be considered a 'good' profile
qc.min_profiles = [100];   % minimum number of 'good' profiles for a timeseries
qc.max_time_gap = [5];   % max gap (days) between 'good' profiles before rejection

dbar_grid = depth_grid;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIRD QC FLAG (min profiles) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Find timeseries without enough profiles
Np = [];
for i = 1:length(qc_ts)
	Np(i) = length(qc_ts(i).cast);
end
ind = find(Np < qc.min_profiles);
qc_ts(ind) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETECTION ALGORTHIM START

for tag_no = 63
        
    %%% Converting time data into datenum
    ts_time = datenum(qc_ts(tag_no).time);
    
    %%% Density calculation
    ts_density = gsw_rho(qc_ts(tag_no).salt, qc_ts(tag_no).temp, zeros(size(qc_ts(tag_no).salt,1), size(qc_ts(tag_no).salt,2)));

    %%% Interpolating into density space
    for i = 1:length(ts_density)
        [a,b] = unique(ts_density(:,i));
        if length(b) > 2
            ds_qc_ts = struct('tag', qc_ts(tag_no).tag, 'cast', qc_ts(tag_no).cast, 'lat', qc_ts(tag_no).lat, 'lon', qc_ts(tag_no).lon, 'time', qc_ts(tag_no).time, 'salt', interp1(ts_density(:,i), qc_ts(tag_no).salt(:,i),density_grid), 'temp',  interp1(ts_density(:,i), qc_ts(tag_no).temp(:,i),density_grid), 'pres', interp1(ts_density(:,i), depth_grid, density_grid));
        end
    end
    
    
    
    %%% Calculating PV
    ts_dz_dp = diff(-qc_ts(tag_no).pres,1) ./ mean(diff(density_grid));
    ts_dp_dz = 1 ./ ts_dz_dp;
    ts_pv = -gsw_f(qc_ts(tag_no).lat) .* ts_dp_dz;
    
    %%% Identifying outliers (using a moving mean across 12 days)
    ts_dates_tda = {};
    j = 1;
    for i = ts_time'
        ts_dates_tda(1,j) = {find(ts_time < i-1.5 & ts_time > i-10)};
        ts_dates_tda(2,j) = {find(ts_time > i+1.5 & ts_time < i+10)};
        j = j + 1;
    end
    
    %%% Calculating number of profiles included in the running mean
    no_profiles = [];
    for i = 1:size(ts_dates_tda,2)
        no_profiles(i) = length(cat(1, ts_dates_tda{1,i}, ts_dates_tda{2,i}));
    end
    
    "Mean number of profiles included in isoutlier detection: " + string(mean(no_profiles));
    "STD of number of profiles included in isoutlier detection: " + string(std(no_profiles));
   
    %%% Identifying PV outliers
    outliers_pv = outlier_detection_date_version(ts_pv, ts_dates_tda);
    
    
    
    figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);
    ax1 = subplot(2,1,1);
    pp = pcolor(ts_time, density_grid(1:size(ts_pv,1)), ts_pv);
    set(pp, 'EdgeColor', 'none');
    set(gca, 'YDir','reverse');
    set(gca, 'Layer','top');
    xlim([ts_time(1) ts_time(end)]);
    xticks(linspace(ts_time(1), ts_time(end), (ts_time(end) - ts_time(1)) / 5));
    set(ax1,'Xticklabel',[])
    ylim([27 27.7]);
    cmap = cmocean('thermal');
    colormap(ax1, cmap);
    cb1 = colorbar('eastoutside');
    title('PV');
    
    ax2 = subplot(2,1,2);
    pp = pcolor(ts_time, density_grid(1:size(outliers_pv,1)), outliers_pv);
    set(pp, 'EdgeColor', 'none');
    set(gca, 'YDir','reverse');
    set(gca, 'Layer','top');
    set(gca,'XColor','w');
    xlim([ts_time(1) ts_time(end)]);
    xticks(linspace(ts_time(1), ts_time(end), (ts_time(end) - ts_time(1)) / 5));
    datetick('x', 'mm/dd', 'keepticks');
    ylim([27 27.7]);
    cmap = cmocean('thermal');
    colormap(ax2, cmap);
    cb2 = colorbar;
    title('PV Outliers');
    
    p1 = get(ax1, 'Position');
    p2 = get(ax2, 'Position');
    p1(2) = p2(2) + p2(4) + 0.03;
    p1(3) = p2(3);
    set(ax1, 'pos', p1);
    
    c1 = get(cb1, 'Position');
    c2 = get(cb2, 'Position');
    c1(1) = c2(1);
    set(cb1, 'pos', c1); 
   
end

function o = outlier_detection_date_version(colorbar_variable, date_indexes)
        %%% Initializing outlier matrix
        o = zeros(size(colorbar_variable));
        
        %%% Loop through data to find outliers over the course of a given
        %%% time window
        for i = 1:size(colorbar_variable,1)
            for j = 1:size(colorbar_variable,2)
                
                movMat_before = colorbar_variable(i, date_indexes{1,j});
                
                movMat_after = colorbar_variable(i, date_indexes{2,j});
                
                no_nans = double(isnan([movMat_before movMat_after]));
                
                if sum(double(no_nans == 0)) < 30
                    
                    o(i,j) = 0;
                    
                else
                    
                    outliers_gen = isoutlier([movMat_before colorbar_variable(i,j) movMat_after]);
                    
                    o(i,j) = double(outliers_gen(length(movMat_before) + 1));
                    
                end
               
            end
        end
end