
%%% Loading MEOP data
% load("SealData_All.mat");

%%% Loading bathymetry data
RTOPO.lat = double(ncread('/Volumes/Elements/Raw Data/RTOPO2.nc', 'lat'));
RTOPO.lon = double(ncread('/Volumes/Elements/Raw Data/RTOPO2.nc', 'lon'))';
RTOPO.bedrock_topography = double(ncread('/Volumes/Elements/Raw Data/RTOPO2.nc', 'bedrock_topography'))';

%%% Loading BedMachineAntarctica data
BedMachineAntarctica.x = double(ncread('/Volumes/Elements/Raw Data/BedMachineAntarctica-v3.nc', 'x'));
BedMachineAntarctica.y = double(ncread('/Volumes/Elements/Raw Data/BedMachineAntarctica-v3.nc', 'y'));
BedMachineAntarctica.bedrock_topography = double(ncread('/Volumes/Elements/Raw Data/BedMachineAntarctica-v3.nc', 'bed'));

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
    pres_qc{i,:} = str2num(SO_sealdata(i).PRES_RAW_QC);
    salt_qc{i,:} = str2num(SO_sealdata(i).SALT_RAW_QC);
    temp_qc{i,:} = str2num(SO_sealdata(i).TEMP_RAW_QC);
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
        SO_sealdata_qc(u).pres = SO_sealdata(i).PRES_RAW(good_data{i,1});
        SO_sealdata_qc(u).salt = SO_sealdata(i).SALT_RAW(good_data{i,1});
        SO_sealdata_qc(u).temp = SO_sealdata(i).TEMP_RAW(good_data{i,1});
        SO_sealdata_qc(u).tag = str2double(SO_sealdata(i).TAG);
        raw_data.pres = SO_sealdata(i).PRES;
        raw_data.salt = SO_sealdata(i).SALT;
        raw_data.temp = SO_sealdata(i).TEMP;
        raw_data.pres_qc = SO_sealdata(i).PRES_QC;
        raw_data.salt_qc = SO_sealdata(i).SALT_QC;
        raw_data.temp_qc = SO_sealdata(i).TEMP_QC;
        raw_data.pres_raw = SO_sealdata(i).PRES_RAW;
        raw_data.salt_raw = SO_sealdata(i).SALT_RAW;
        raw_data.temp_raw = SO_sealdata(i).TEMP_RAW;
        raw_data.pres_raw_qc = SO_sealdata(i).PRES_RAW_QC;
        raw_data.salt_raw_qc = SO_sealdata(i).SALT_RAW_QC;
        raw_data.temp_raw_qc = SO_sealdata(i).TEMP_RAW_QC;
        SO_sealdata_qc(u).raw_data = raw_data;
        u = u + 1;
    end
end

clear time_qc pres_qc salt_qc temp_qc u SO_sealdata i x good_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Quality Controlling Individual Profiles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Imposing max pressure difference of qc.max_pres_diff
for i = 1:length(SO_sealdata_qc)
    pres_diff = [0; diff(SO_sealdata_qc(i).pres)];
    ind = (pres_diff) <= prms.qc.max_pres_diff_shallow;
    SO_sealdata_qc(i).pres = SO_sealdata_qc(i).pres(ind);
    SO_sealdata_qc(i).salt = SO_sealdata_qc(i).salt(ind);
    SO_sealdata_qc(i).temp = SO_sealdata_qc(i).temp(ind);
end

%%% Extracting max pressure achieved by each profile
max_pres = NaN(1,length(SO_sealdata_qc));
for i = 1:length(SO_sealdata_qc)
    max_pres(i) = max(SO_sealdata_qc(i).pres);
end

%%% Only keeping profiles that go deep enough
keepit = (max_pres >= prms.qc.min_pres);
SO_sealdata_qc = SO_sealdata_qc(keepit);

clear keepit pres_diff ind max_pres

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolating to Uniform Pressure Grid %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear interp_sealdata_qc interp_salt interp_temp

depth_grid = 0:1:800;
depth_grid = depth_grid';

parfor i = 1:length(SO_sealdata_qc)
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
	[~,b] = sort(profdate(tagidx),'ascend');

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
	didx     = diff(ts_dates)<prms.qc.max_time_gap;
	% Find the gaps greater than QC threshold 
	ind      = find(didx == 1); 
	% Clear start of time-series
	ts_start = []; ts_end = []; ts_cnt = 1;
	if isempty(ind)
		% No good groups found, go to next MEOP time-series
		continue
	elseif length(ind)==1
		% Only one group of consecutive casts exist
		tmpidx{ts_cnt} = [ind:ind+1];
	else
		% Many groups exist, separate them
		for j = 1:length(ind)
			if isempty(ts_start)
				ts_start = ind(j);
			end
			if j < length(ind) & ismember(ind(j)+1,ind)
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
ind = find(Np < prms.qc.min_profiles);
qc_ts(ind) = [];

clear meop_data profdate tagidx taglist tagnum j i a b tmpdat tmpidx ts_cnt ...
    ts_end ts_start tscount tslist Np ind didx qc flagit ts_dates dbar_grid meop_ts

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Start of Detection Algorithm %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Keep an eye on tag 197, prof 163
test_prof = 115; %1:length(qc_ts); %[115 197 62 63 65 77 81 92 25 100 101 129 122 118 199 52 205 207 208 257 376 441];

%%% Loading algorithm settings
run("MEOPseals_algorithm_settings.m")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------------')
disp('Calculating Variables')
disp('---------------------')

for tag_no = test_prof

    %%% Creating arrays to hold rejection and justification data
    rejected.anticyclones = zeros(size(qc_ts(tag_no).cast)); % Anticyclones
    reason.anticyclones = strings(size(qc_ts(tag_no).cast));
    rejected.cyclones = zeros(size(qc_ts(tag_no).cast)); % Cyclones
    reason.cyclones = strings(size(qc_ts(tag_no).cast));
    qc_ts(tag_no).rejected = rejected;
    qc_ts(tag_no).reason = reason;
    
    %%% Calculating bathymetry
    qc_ts(tag_no).bathymetry = interp2(RTOPO.lon, RTOPO.lat', RTOPO.bedrock_topography, qc_ts(tag_no).lon, qc_ts(tag_no).lat);

    %%% Calculating BedMachineAntarctica bathymetry
    [x,y] = ll2xy(qc_ts(tag_no).lat, qc_ts(tag_no).lon, -1);
    qc_ts(tag_no).BedMachine = interpBedmachineAntarctica(x, y, 'bed');

    %%% Calculating pressure space variables
    qc_ts(tag_no).ps = calc_pres_space_vars(qc_ts, tag_no, depth_grid, 0);

end

clear rejected reason

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Checking for Density Inversions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------------')
disp('Checking for Density Inversions')
disp('-------------------------------')

for tag_no = test_prof

    %%% Checking for density inversions. Small inversions are corrected.
    %%% Profiles with large inversions will be excluded.
    [qc_ts(tag_no).ps.sigma0, qc_ts(tag_no).ps.max_sigma0_inversion,...
        qc_ts(tag_no).rejected, qc_ts(tag_no).reason]...
        = checking_for_density_inversions(qc_ts, tag_no);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolating to Density Grid %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('------------------------------')
disp('Interpolating to Density Space')
disp('------------------------------')

%%% Creating density grid
density_grid = (26.0:0.001:28.5)';

for tag_no = test_prof
    
    %%% Interpolating pressure space variables to density grid
    qc_ts(tag_no).ds = interp_to_sigma0_space(qc_ts, tag_no, density_grid);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Isopycnal Separation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--------------------------------')
disp('Calculating Isopycnal Separation')
disp('--------------------------------')

for tag_no = test_prof
    
    %%% Calculating isopycnal separation
    qc_ts(tag_no).ds.isopycnal_separation = calc_isopycnal_separation(qc_ts, tag_no);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Finding Indices to Build Reference Profiles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--------------------------------------')
disp('Getting Indices for Reference Profiles')
disp('--------------------------------------')
 
for tag_no = test_prof
        
    %%% Finding indices to build reference profiles
    qc_ts(tag_no).ref_ind = indices_for_ref_profiles(qc_ts,tag_no, prms.refprof);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Building Reference Profiles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------------------')
disp('Building Reference Profiles')
disp('---------------------------')

for tag_no = test_prof

    %%% Building density-space reference profiles
    [qc_ts(tag_no).ds.ref_salt, qc_ts(tag_no).ds.ref_temp,...
        qc_ts(tag_no).ds.ref_N2, qc_ts(tag_no).ds.ref_spice,...
        qc_ts(tag_no).ds.ref_isopycnal_separation] = build_density_space_ref_profiles(qc_ts, tag_no);

    %%% Building pressure-space reference profiles
    [qc_ts(tag_no).ps.ref_dyn_height_anom, qc_ts(tag_no).ps.ref_N2]...
        = build_pres_space_ref_profiles(qc_ts, tag_no);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Anomalies %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------------')
disp('Calculating Anomalies')
disp('---------------------')

for tag_no = test_prof

    %%% Calculating density-space anomalies
    qc_ts(tag_no).ds.anoms = calc_density_space_anomalies(qc_ts, tag_no);

    %%% Calculating pressure-space anomalies
    qc_ts(tag_no).ps.anoms = calc_pres_space_anomalies(qc_ts, tag_no);

end

%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating IQR %%%
%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------')
disp('Calculating IQR')
disp('---------------')

for tag_no = test_prof

    %%% Calculating iqr, lower-threshold, and upper threshold
    qc_ts(tag_no).ds.iqrs = calc_iqr(qc_ts, tag_no);
    
end

%save('/Volumes/Elements/MEOPData.mat', 'qc_ts', 'depth_grid');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating Arrays for Detection Results %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------')
disp('Creating Detection Arrays')
disp('-------------------------')

for tag_no = test_prof

    %%% Creating array to hold anticyclonic detections
    anticyclones.isa_iqr = [];
    anticyclones.isa_gaussian = [];
    anticyclones.dha = [];
    anticyclones.isopycnal_stability = [];
    anticyclones.bathymetric_stability = [];
    anticyclones.MLD = [];
    anticyclones.final = [];
    qc_ts(tag_no).anticyclones = anticyclones;

    %%% Creating array to hold cyclonic detections
    cyclones.spice_iqr = [];
    cyclones.spice_gaussian = [];
    cyclones.dha = [];
    cyclones.isopycnal_stability = [];
    cyclones.bathymetric_stability = [];
    cyclones.MLD = [];
    cyclones.final = [];
    qc_ts(tag_no).cyclones = cyclones;

    %%% Creating arrays to hold profile rejections and reasons
    rejected.anticyclones = zeros(size(qc_ts(tag_no).cast));
    reason.anticyclones = strings(size(qc_ts(tag_no).cast));
    rejected.cyclones = zeros(size(qc_ts(tag_no).cast));
    reason.cyclones = strings(size(qc_ts(tag_no).cast));
    qc_ts(tag_no).reason = reason;
    qc_ts(tag_no).rejected = rejected;

end

clear anticyclones cyclones
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Anticyclonic Detections %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------------------------')
disp('Beginning Anticyclonic Detections')
disp('---------------------------------')

qc_ts = anticyclonic_checks(qc_ts, test_prof, prms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cyclonic Detection Checks %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-----------------------------')
disp('Beginning Cyclonic Detections')
disp('-----------------------------')

qc_ts = cyclonic_checks(qc_ts, test_prof, prms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Isopycnal Stability Check (anticyclones AND cyclones) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------')
disp('Isopycnal Stability Check')
disp('-------------------------')

for tag_no = test_prof

    %%% Creating array to hold isopycnal separation variance data
    qc_ts(tag_no).isopycnal_stability = NaN(size(qc_ts(tag_no).cast));

    %%% Isopycnal stability check (anticyclones AND cyclones)
    [qc_ts(tag_no).anticyclones, qc_ts(tag_no).cyclones, ...
        qc_ts(tag_no).isopycnal_stability, qc_ts(tag_no).rejected, ...
        qc_ts(tag_no).reason] = isopycnal_stability_check(qc_ts, tag_no, prms);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bathymetric Stability Check (anticyclones AND cyclones) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------------------');
disp('Bathymetric Stability Check');
disp('---------------------------');

for tag_no = test_prof

    %%% Creating array to hold bathymetric variance data
    qc_ts(tag_no).bathymetric_var = NaN(size(qc_ts(tag_no).cast));

    %%% Bathymetric stability check (anticyclones AND cyclones)
    [qc_ts(tag_no).anticyclones, qc_ts(tag_no).cyclones,...
        qc_ts(tag_no).bathymetric_var, qc_ts(tag_no).rejected,...
        qc_ts(tag_no).reason] = bathymetry_check(qc_ts, tag_no, prms);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mixed Layer Depth Check %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-----------------------')
disp('Mixed Layer Depth Check')
disp('-----------------------')

for tag_no = test_prof
    
    qc_ts(tag_no).MLD = NaN(size(qc_ts(tag_no).cast));

    [qc_ts(tag_no).anticyclones, qc_ts(tag_no).cyclones,...
        qc_ts(tag_no).MLD, qc_ts(tag_no).rejected,...
        qc_ts(tag_no).reason] = MLD_check(qc_ts, tag_no, prms);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Removing Edge Cases %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------');
disp('Removing Edge Cases');
disp('-------------------');

for tag_no = test_prof

    %%% Removing edge cases (anticyclones AND cyclones)
    [qc_ts(tag_no).anticyclones, qc_ts(tag_no).cyclones,...
       qc_ts(tag_no).rejected, qc_ts(tag_no).reason]...
       = removing_edge_cases(qc_ts, tag_no);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Summarizing Detected SCVs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('----------------------')
disp('Summarizing Detections')
disp('----------------------')

MEOPscvs = [];
u = 1;
for tag_no = test_prof
    for i = unique([qc_ts(tag_no).anticyclones.final qc_ts(tag_no).cyclones.final])
        if isempty(unique([qc_ts(tag_no).anticyclones.final qc_ts(tag_no).cyclones.final]))
            break
        end
        MEOPscvs(u).tag = qc_ts(tag_no).tag;
        MEOPscvs(u).tag_no = tag_no;
        MEOPscvs(u).cast = i;
        MEOPscvs(u).time = qc_ts(tag_no).time(i,:);
        MEOPscvs(u).lat = qc_ts(tag_no).lat(i);
        MEOPscvs(u).lon = qc_ts(tag_no).lon(i);
        if ismember(i, qc_ts(tag_no).anticyclones.final)
            MEOPscvs(u).type = "Anticyclonic";
        else
            MEOPscvs(u).type = "Cyclonic";
        end
        MEOPscvs(u).isopycnal_stability = qc_ts(tag_no).isopycnal_stability(i);
        MEOPscvs(u).bathymetric_variance = qc_ts(tag_no).bathymetric_var(i);

        u = u + 1;

    end
end

%%% Noting region of detection
for i = 1:length(MEOPscvs)
    if MEOPscvs(i).lon > -60 && MEOPscvs(i).lon < 0
        MEOPscvs(i).region = "Weddell";
    elseif MEOPscvs(i).lon <= -60 && MEOPscvs(i).lon > -120
        MEOPscvs(i).region = "WAP";
    elseif MEOPscvs(i).lon <= -120 || MEOPscvs(i).lon > 170
        MEOPscvs(i).region = "Ross";
    else
        MEOPscvs(i).region = "East";
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Summarizing Figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%MEOPanticyclones = MEOPscvs(strcmp([MEOPscvs.type], "Anticyclonic") & [MEOPscvs.isopycnal_stability] > 75);

isopycnals = 0.03;
lw = 2;

j = 0;
for tag_no = 199 %vertcat(MEOPscvs.tag_no)'
        j = j + 1;
        i = 179 %MEOPscvs(j).cast;

        casts = [qc_ts(tag_no).ref_ind{1,i} qc_ts(tag_no).ref_ind{2,i}];
        ind = min(casts):max(casts);
        x_datenum = datenum(qc_ts(tag_no).time(ind,:));

        fig = figure('Position', [0 0 1300 950]);
        %sgtitle('MEOP Seal Track, Tag = ' + string(qc_ts(tag_no).tag) + ', Tag # = ' + string(tag_no) + ', Cast = ' + string(i) + ', Type = ' + string(MEOPscvs(j).type), 'FontSize', 18, 'FontWeight', 'bold')
        sgtitle('MEOP Seal Track, Tag: ' + string(qc_ts(tag_no).tag) + ', Cast: ' + string(i) + ', Isopycnal Stability: ' + string(qc_ts(tag_no).isopycnal_stability(i)), 'FontSize', 18, 'FontWeight', 'bold')

        %%% Bathymetry Subplot
        ax0 = subplot(5,3, [1 2]);
        hold on
        a = plot(x_datenum, qc_ts(tag_no).bathymetry(ind), 'k-o', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'DisplayName', 'RTOPO2');
        b = plot(x_datenum, qc_ts(tag_no).BedMachine(ind), 'r-o', 'MarkerFaceColor', 'r', 'LineWidth', 2, 'DisplayName', 'BedMachine');
        xline(datenum(qc_ts(tag_no).time(i,:)), 'r', 'LineWidth', 1.5)
        hold off
        xlim([qc_ts(tag_no).cast(ind(1)) qc_ts(tag_no).cast(ind(end))])
        ylabel('Bedrock Topography (m)', 'FontSize', 12)
        title('Bathymetry', 'FontSize', 12)
        datetick('x', 'mm/dd')
        xlim([x_datenum(1) x_datenum(end)])
        legend([a b])

        %%% Temperature Subplot
        ax1 = subplot(5,3, [4 5]);
        hold on
        pp = pcolor(x_datenum, qc_ts(tag_no).ps.pres(:,1), qc_ts(tag_no).temp(:,ind));
        set(pp, 'EdgeColor', 'none');
        [~, dd] = unique(x_datenum);
        [C,h] = contour(ax1, x_datenum(dd), depth_grid, qc_ts(tag_no).ps.sigma0(:,ind(dd)), round(min(min(qc_ts(tag_no).ps.sigma0)):isopycnals:max(max(qc_ts(tag_no).ps.sigma0)), 2), 'k');
        clabel(C,h,'LabelSpacing',500);
        xline(datenum(qc_ts(tag_no).time(i,:)), 'r', 'LineWidth', 1.5)
        hold off
        cmap = cmocean('thermal'); colormap(ax1, cmap); colorbar;
        clim([min(min(qc_ts(tag_no).temp(:,ind))) max(max(qc_ts(tag_no).temp(:,ind)))])
        set(gca, 'YDir','reverse'); set(gca, 'Layer','top');
        ylabel('Pressure (dbar)', 'FontSize', 12); ylim([0 500])
        title('Temperature', 'FontSize', 12);
        datetick('x', 'mm/dd')
        xlim([x_datenum(1) x_datenum(end)])

        %%% Salinity Subplot
        ax2 = subplot(5,3,[7 8]);
        hold on
        pp = pcolor(x_datenum, qc_ts(tag_no).ps.pres(:,1), qc_ts(tag_no).salt(:,ind));
        set(pp, 'EdgeColor', 'none');
        [~, dd] = unique(x_datenum);
        [C,h] = contour(ax2, x_datenum(dd), depth_grid, qc_ts(tag_no).ps.sigma0(:,ind(dd)), round(min(min(qc_ts(tag_no).ps.sigma0)):isopycnals:max(max(qc_ts(tag_no).ps.sigma0)), 2), 'k');
        clabel(C,h,'LabelSpacing',500);
        xline(datenum(qc_ts(tag_no).time(i,:)), 'r', 'LineWidth', 1.5)
        hold off
        cmap = cmocean('haline'); colormap(ax2, cmap); colorbar;
        clim([min(min(qc_ts(tag_no).salt(:,ind))) max(max(qc_ts(tag_no).salt(:,ind)))])
        set(gca, 'YDir','reverse'); set(gca, 'Layer','top');
        ylabel('Pressure (dbar)', 'FontSize', 12); ylim([0 500]);
        title('Salinity', 'FontSize', 12);
        datetick('x', 'mm/dd')
        xlim([x_datenum(1) x_datenum(end)])

        %%% Realigning Plots
        p0 = get(ax0, 'Position');
        p1 = get(ax1, 'Position');
        p2 = get(ax2, 'Position');
        p1(3) = p0(3);
        p2(3) = p0(3);
        set(ax1, 'Position', p1);
        set(ax2, 'Position', p2);
        clear p0 p1 p2

        %%% Map Subplot
        subplot(5,3,[3 6 9])
        load coastlines
        axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -57]);
        axis off; framem on; 
        geoshow(coastlat, coastlon, 'Color', 'k')
        scatterm(qc_ts(tag_no).lat(i), qc_ts(tag_no).lon(i), 20, 'ko', 'MarkerFaceColor', 'g')

%         %%% Zoomed Map Subplot
%         ax3 = subplot(5,3,[6 9]);
% 
% %         lat_min = min(qc_ts(tag_no).lat(ind)) - 0.1;
% %         lat_max = max(qc_ts(tag_no).lat(ind)) + 0.1;
% %         lon_min = min(qc_ts(tag_no).lon(ind)) - 0.5;
% %         lon_max = max(qc_ts(tag_no).lon(ind)) + 0.5;
% 
%         lat_min = qc_ts(tag_no).lat(i) - 0.5;
%         lat_max = qc_ts(tag_no).lat(i) + 0.5;
%         lon_min = qc_ts(tag_no).lon(i) - 1;
%         lon_max = qc_ts(tag_no).lon(i) + 1;
% 
%         if lon_max - lon_min > 340
%             lon_sec = qc_ts(tag_no).lon(ind);
%             pos_min = min(lon_sec(lon_sec > 0));
%             lon_min = -180 - (180-pos_min) - 0.5;
%             lon_max = max(lon_sec(lon_sec < 0)) + 0.5;
%         end
% 
%         %%% Figure settings
%         load coastlines
%         hold on
% 
%         %%% Grabbing sector of LLC time series
%         plot(qc_ts(tag_no).lat(ind), qc_ts(tag_no).lon(ind), '-s', 'Color', [0.5 0.5 0.5],...
%             'MarkerSize', 4, 'MarkerFaceColor', [0.5 0.5 0.5], 'LineWidth', 2)
%         scatter(qc_ts(tag_no).lat(i), qc_ts(tag_no).lon(i), 6, 'gs', 'MarkerFaceColor', 'g')
%         xlabel('Latitude'); xlim([lat_min lat_max]);
%         ylabel('Longitude'); ylim([lon_min lon_max]);

        %%% Spice Profile
        subplot(5,4,13);
        hold on
        plot(qc_ts(tag_no).ds.spice(:,i), qc_ts(tag_no).ds.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', lw)
        plot(qc_ts(tag_no).ds.ref_spice(:,i), qc_ts(tag_no).ds.pres(:,i), 'k', 'DisplayName', 'Reference', 'LineWidth', lw)
        hold off
        xlabel('Spice');
        set(gca, 'YDir', 'reverse');   
        ylim([0 500]);
        legend('Location', 'best')

        %%% Isopycnal Separation Profile
        subplot(5,4,14);
        hold on
        plot(qc_ts(tag_no).ds.isopycnal_separation(:,i), qc_ts(tag_no).ds.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', lw)
        plot(qc_ts(tag_no).ds.ref_isopycnal_separation(:,i), qc_ts(tag_no).ds.pres(:,i), 'k','DisplayName', 'Reference', 'LineWidth', lw)
        hold off
        xlabel('Isopycnal Separation');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);
        legend('Location', 'best')
        
        %%% N2 Profile
        subplot(5,4,15);
        hold on
        plot(qc_ts(tag_no).ds.N2(:,i), qc_ts(tag_no).ds.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', lw)
        plot(qc_ts(tag_no).ds.ref_N2(:,i), qc_ts(tag_no).ds.pres(:,i), 'k','DisplayName', 'Reference', 'LineWidth', lw)
        hold off
        xlabel('N^2');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);
        legend('Location', 'best')

        %%% DHA Profile
        subplot(5,4,16)
        hold on
        plot(qc_ts(tag_no).ps.dyn_height_anom(:,i), qc_ts(tag_no).ps.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', lw)
        plot(qc_ts(tag_no).ps.ref_dyn_height_anom(:,i), qc_ts(tag_no).ps.pres(:,i), 'k','DisplayName', 'Reference', 'LineWidth', lw)
        hold off
        xlabel('Dynamic Height Anomaly');
        set(gca, 'YDir', 'reverse');
        legend('Location', 'best')
        ylim([0 500])

        %%% Spice Anomaly Profile
        subplot(5,4,17);
        hold on
        plot(qc_ts(tag_no).ds.anoms.spice(:,i), qc_ts(tag_no).ds.pres(:,i), 'b', 'DisplayName', 'Profile', 'LineWidth', lw)
        if ~isempty(qc_ts(tag_no).spice_gauss) && (length(qc_ts(tag_no).spice_gauss) >= i) && ~isempty(qc_ts(tag_no).spice_gauss(i).rejected)
            plot(qc_ts(tag_no).spice_gauss(i).X,qc_ts(tag_no).spice_gauss(i).Y,'--','Color', 'm', 'LineWidth',lw);
        end
        xline(0, '--k', 'LineWidth', 1);
        hold off
        xlabel('Spice Anomaly');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);

        %%% Isopycnal Separation Anomaly Profile
        subplot(5,4,18);
        hold on
        plot(qc_ts(tag_no).ds.anoms.isopycnal_separation_normalized(:,i), qc_ts(tag_no).ds.pres(:,i), 'b', 'DisplayName', 'Profile', 'LineWidth', lw)
        if ~isempty(qc_ts(tag_no).isa_gauss) && (length(qc_ts(tag_no).isa_gauss) >= i) && ~isempty(qc_ts(tag_no).isa_gauss(i).rejected)
            plot(qc_ts(tag_no).isa_gauss(i).X,qc_ts(tag_no).isa_gauss(i).Y,'--','Color', 'm', 'LineWidth',lw);
        end
        xline(0, '--k', 'LineWidth', 1);
        hold off
        xlabel('Isopycnal Separation Anomaly');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);

        %%% N2 Anomaly Profile
        subplot(5,4,19);
        hold on
        plot(qc_ts(tag_no).ds.anoms.N2(:,i), qc_ts(tag_no).ds.pres(:,i), 'b', 'DisplayName', 'Profile', 'LineWidth', lw)
        xline(0, '--k', 'LineWidth', 1);
        hold off
        xlabel('N^2 Anomaly');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);

        %%% Dynamic Height Anomaly Profile
        subplot(5,4,20);
        hold on
        xline(0, '--k', 'LineWidth', 1);
        a = plot(qc_ts(tag_no).ps.anoms.dyn_height_anom(:,i), qc_ts(tag_no).ps.pres(:,i),'b', 'DisplayName', 'PS','LineWidth',lw);
        if ~isempty(qc_ts(tag_no).dha_gauss) && (length(qc_ts(tag_no).dha_gauss) >= i) && ~isempty(qc_ts(tag_no).dha_gauss(i).rejected)
            b = plot(qc_ts(tag_no).dha_gauss(i).dataX,qc_ts(tag_no).dha_gauss(i).dataY,'Color', [0.5 0.5 0.5],'LineWidth',lw, 'DisplayName', 'Adjusted');
            plot(qc_ts(tag_no).dha_gauss(i).X,qc_ts(tag_no).dha_gauss(i).Y,'--','Color', 'm', 'LineWidth',lw);
            legend([a b], 'Location', 'best')
        end
        hold off
        xlabel('Dynamic Height Anomaly Anomaly');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);
        
        %saveas(fig, '/Volumes/Elements/MEOP SCVs/' + string(qc_ts(tag_no).tag) + '_' + string(i), 'png')
end
        clear ax0 ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 C h cmap fig h IB isopycnals p1 p2 p3 p4 p5 pp casts lon_max lon_min lat_max lat_min coastlat coastlon i j a b

%%
for ii  = 1:length(MEOPscvs)
    tag_no = MEOPscvs(ii).tag_no;
    i = MEOPscvs(ii).cast;
    MEOPscvs(ii).MLD = MLD_calc(qc_ts, tag_no, i);
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Summarizing SCV Statistics %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
sgtitle('MEOP Anticyclonic Detections', 'FontSize', 15, 'FontWeight', 'bold')
for j = 1:length(MEOPscvs)
    color_matrix(j,:) = rand(1,3);
end


new_pres = 0:10:500;

subplot(4,2,1)
hold on
xline(0, '--k', 'LineWidth', 1);
for j = 1:length(MEOPscvs)
    tag_no = MEOPscvs(j).tag_no;
    i = MEOPscvs(j).cast;
    plot(qc_ts(tag_no).ds.anoms.salt(:,i), qc_ts(tag_no).ds.pres(:,i), 'Color', color_matrix(j,:));
    ind = ~isnan(qc_ts(tag_no).ds.anoms.salt(:,i));
    salt_all(:,j) = interp1(qc_ts(tag_no).ds.pres(ind,i), qc_ts(tag_no).ds.anoms.salt(ind,i), new_pres);
end
plot(mean(salt_all, 2, 'omitnan'), new_pres, 'k', 'LineWidth', 2)
xlabel('Salinity Anomaly (psu)');
ylabel('Pressure (dbar)')
ylim([0 600])
set(gca, 'YDir', 'reverse');

subplot(4,2,3)
hold on
xline(0, '--k', 'LineWidth', 1);
for j = 1:length(MEOPscvs)
    tag_no = MEOPscvs(j).tag_no;
    i = MEOPscvs(j).cast;
    plot(qc_ts(tag_no).ds.anoms.temp(:,i), qc_ts(tag_no).ds.pres(:,i), 'Color', color_matrix(j,:));
    ind = ~isnan(qc_ts(tag_no).ds.anoms.temp(:,i));
    temp_all(:,j) = interp1(qc_ts(tag_no).ds.pres(ind,i), qc_ts(tag_no).ds.anoms.temp(ind,i), new_pres);
end
plot(mean(temp_all, 2, 'omitnan'), new_pres, 'k', 'LineWidth', 2)
xlabel('Temperature Anomaly (degC)');
ylabel('Pressure (dbar)')
ylim([0 600])
set(gca, 'YDir', 'reverse');

subplot(4,2,2)
hold on
xline(0, '--k', 'LineWidth', 1);
for j = 1:length(MEOPscvs)
    tag_no = MEOPscvs(j).tag_no;
    i = MEOPscvs(j).cast;
    plot(qc_ts(tag_no).ds.anoms.spice(:,i), qc_ts(tag_no).ds.pres(:,i), 'Color', color_matrix(j,:));
    ind = ~isnan(qc_ts(tag_no).ds.anoms.spice(:,i));
    spice_all(:,j) = interp1(qc_ts(tag_no).ds.pres(ind,i), qc_ts(tag_no).ds.anoms.spice(ind,i), new_pres);
end
plot(mean(spice_all, 2, 'omitnan'), new_pres, 'k', 'LineWidth', 2)
xlabel('Spice Anomaly (kg/m^3)');
ylabel('Pressure (dbar)')
ylim([0 600])
set(gca, 'YDir', 'reverse');

subplot(4,2,4)
hold on
xline(0, '--k', 'LineWidth', 1);
for j = 1:length(MEOPscvs)
    tag_no = MEOPscvs(j).tag_no;
    i = MEOPscvs(j).cast;
    plot(qc_ts(tag_no).ds.anoms.isopycnal_separation(:,i), qc_ts(tag_no).ds.pres(:,i), 'Color', color_matrix(j,:));
    ind = ~isnan(qc_ts(tag_no).ds.anoms.isopycnal_separation(:,i));
    isa_all(:,j) = interp1(qc_ts(tag_no).ds.pres(ind,i), qc_ts(tag_no).ds.anoms.isopycnal_separation(ind,i), new_pres);
end
plot(mean(isa_all, 2, 'omitnan'), new_pres, 'k', 'LineWidth', 2)
xlabel('Isopycnal Separation Anomaly');
ylabel('Pressure')
ylim([0 600])
set(gca, 'YDir', 'reverse');

subplot(4,2,5)
hold on
xline(0, '--k', 'LineWidth', 1);
for j = 1:length(MEOPscvs)
    tag_no = MEOPscvs(j).tag_no;
    i = MEOPscvs(j).cast;
    plot(qc_ts(tag_no).ps.anoms.dyn_height_anom(:,i), qc_ts(tag_no).ps.pres(:,i), 'Color', color_matrix(j,:));
    ind = ~isnan(qc_ts(tag_no).ps.anoms.dyn_height_anom(:,i));
    dha_all(:,j) = interp1(qc_ts(tag_no).ps.pres(ind,i), qc_ts(tag_no).ps.anoms.dyn_height_anom(ind,i), new_pres);
end
plot(mean(dha_all, 2, 'omitnan'), new_pres, 'k', 'LineWidth', 2)
xlabel('Dynamic Height Anomaly Anomaly (m^2/s^2)');
ylabel('Pressure (dbar)')
ylim([0 600])
set(gca, 'YDir', 'reverse');

subplot(4,2,7)
hold on
xline(0, '--k', 'LineWidth', 1);
for j = 1:length(MEOPscvs)
    tag_no = MEOPscvs(j).tag_no;
    i = MEOPscvs(j).cast;
    plot(qc_ts(tag_no).ds.anoms.N2(:,i), qc_ts(tag_no).ds.pres(:,i), 'Color', color_matrix(j,:));
    ind = ~isnan(qc_ts(tag_no).ds.anoms.N2(:,i));
    N2_all(:,j) = interp1(qc_ts(tag_no).ds.pres(ind,i), qc_ts(tag_no).ds.anoms.N2(ind,i), new_pres);
end
plot(mean(N2_all, 2, 'omitnan'), new_pres, 'k', 'LineWidth', 2)
xlabel('N^2 Anomaly (1/min^2)');
ylabel('Pressure (dbar)')
ylim([0 600])
set(gca, 'YDir', 'reverse');

subplot(4,2,6)
for j = 1:length(MEOPscvs)
    tag_no = MEOPscvs(j).tag_no;
    i = MEOPscvs(j).cast;
    prof_datetime = datetime(qc_ts(tag_no).time(i,:));
    doy(j) = day(prof_datetime, 'dayofyear');
end
histogram(doy, 365);
xlabel('Day of Year')
ylabel('Number of Detections')

subplot(4,2,8)
load coastlines
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -57]);
axis off; framem on;
geoshow(coastlat, coastlon, 'Color', 'k')
for j = 1:length(MEOPscvs)
    tag_no = MEOPscvs(j).tag_no;
    i = MEOPscvs(j).cast;
    scatterm(qc_ts(tag_no).lat(i), qc_ts(tag_no).lon(i), 50, 'ko', 'MarkerFaceColor', color_matrix(j,:))
end

