
RTOPO_lat = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'lat'));
RTOPO_lon = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'lon'))';
RTOPO_bedrock_topography = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'bedrock_topography'))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inner_window = 2; %%% setting inner "exclusion" window for background calculation
outer_window = 7; %%% setting outer "inclusion" window for background calculation

for tag_no = 1
    
    %%% Creating new structure for seal of interest
    meop_ts = qc_ts(tag_no);

    %%% Finding indices to use for background profile calculation 
    for i = 1:length(qc_ts(tag_no).cast)
        mean_ind(1,i) = {(i-outer_window):(i-inner_window)};
        mean_ind{1,i}(mean_ind{1,i} < 1) = []; %%% Making sure indices remain within ts boundaries
        mean_ind(2,i) = {(i+inner_window):(i+outer_window)};
        mean_ind{2,i}(mean_ind{2,i} > length(qc_ts(tag_no).cast)) = []; %%% Making sure indices remain within ts boundaries
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Calculating Variables for Time Series %%%%%%%%%%%%%%%%%%

    %%% Calculating N^2, Dynamic Heigh Anomaly

    depth_grid = 1:500;

    %%% Calculating absolute salinity and conservative temperature
    meop_ts.salt_absolute = gsw_SA_from_SP(meop_ts.salt, depth_grid', meop_ts.lon, meop_ts.lat);
    meop_ts.temp_conservative = gsw_CT_from_t(meop_ts.salt_absolute, meop_ts.temp, depth_grid');

    %%% Calculating N^2
    [meop_ts.N2,~] = gsw_Nsquared(meop_ts.salt_absolute, meop_ts.temp_conservative, depth_grid' .* ones(size(meop_ts.salt)), meop_ts.lat .* ones(size(meop_ts.salt)));

    for i = 1:length(meop_ts.cast)
        %%% Finding maximum depth
        meop_ts.max_depth(i) = find(~isnan(meop_ts.temp(:,i)),1,'last');

        %%% calculating dynamic height anomaly
        meop_ts.dyn_height_anom(:,i) = gsw_geo_strf_dyn_height(meop_ts.salt_absolute(:,i), meop_ts.temp_conservative(:,i), depth_grid', meop_ts.max_depth(i));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Calculating Reference Profiles %%%%%%%%%%%%%%%%%%%%%%%%

    %%% Creating reference profiles against which the anomalies can be
    %%% calculated
    for i = 1:length(meop_ts.cast)
        ref_meop_ts.salt(:,i) = mean(meop_ts.salt(:,[mean_ind{1,i} mean_ind{2,i}]), 2, 'omitnan');
        ref_meop_ts.temp(:,i) = mean(meop_ts.temp(:,[mean_ind{1,i} mean_ind{2,i}]), 2, 'omitnan');
        ref_meop_ts.N2(:,i) = mean(meop_ts.N2(:,[mean_ind{1,i} mean_ind{2,i}]), 2, 'omitnan');
        ref_meop_ts.dyn_height_anom(:,i) = mean(meop_ts.dyn_height_anom(:,[mean_ind{1,i} mean_ind{2,i}]), 2, 'omitnan');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Bathymetry Check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Creates new time series with profiles that may be checked for SCVs

    %%% Finding bathymetry along seal time series
    ts_bathymetry = interp2(RTOPO_lon, RTOPO_lat', RTOPO_bedrock_topography, qc_ts(tag_no).lon,qc_ts(tag_no).lat);
    
    %%% Calculating bathymetry of background profile
    ref_bathymetry_mean = NaN(size(ts_bathymetry));
    ref_bathymetry_std = NaN(size(ts_bathymetry));
    for i = 1:length(qc_ts(tag_no).cast)
        ref_bathymetry_mean(i) = mean(ts_bathymetry([mean_ind{1,i}, mean_ind{2,i}]));
        ref_bathymetry_std(i) = std(ts_bathymetry([mean_ind{1,i}, mean_ind{2,i}]));
    end

    %%% Checking each profile against assigned background profile
    bathymetry_change = abs(ts_bathymetry - ref_bathymetry_mean);

    %%% Flagging profiles for removal based on bathymetry check
    flagit = find((bathymetry_change > 500) | (ref_bathymetry_std > 500));

    %%% Removing flagged profiles from seal time series
    meop_ts.cast(flagit) = [];
    meop_ts.lat(flagit) = [];
    meop_ts.lon(flagit) = [];
    meop_ts.time(flagit,:) = [];
    meop_ts.salt(:,flagit) = [];
    meop_ts.temp(:,flagit) = [];
    meop_ts.salt_absolute(:,flagit) = [];
    meop_ts.temp_conservative(:,flagit) = [];
    meop_ts.N2(:,flagit) = [];
    meop_ts.dyn_height_anom(:,flagit) = [];
    meop_ts.max_depth(:,flagit) = [];
    
    %%% Removing corresponding reference profiles
    ref_meop_ts.salt(:,flagit) = [];
    ref_meop_ts.temp(:,flagit) = [];
    ref_meop_ts.N2(:,flagit) = [];
    ref_meop_ts.dyn_height_anom(:,flagit) = [];

    clear ts_bathymetry ref_bathymetry_std ref_bathymetry_mean bathymetry_change flagit i

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Figure to see Profile vs Reference %%%%%%%%%%%%%%%%%%%%%
    i = 1;

    figure()
    subplot(1,2,1);
    hold on
    plot(meop_ts.N2(:,i), 1:499, 'k', 'LineWidth', 2)
    plot(ref_meop_ts.N2(:,i), 1:499, 'r', 'LineWidth', 2)
    hold off
    set(gca,'YDir','Reverse');
    title('N^2')

    subplot(1,2,2);
    hold on
    plot(meop_ts.dyn_height_anom(:,i), 1:500, 'k', 'LineWidth', 2)
    plot(ref_meop_ts.dyn_height_anom(:,i), 1:500, 'r', 'LineWidth', 2)
    hold off
    set(gca,'YDir','Reverse');
    title('Dynamic Height Anomaly')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% Dynamic Height Anomaly Check %%%%%%%%%%%%%%%%%%%%

    for i = 1:length(meop_ts.cast)

        %%% Creating 
        ref_meop_profile(i).salt = ref_meop_ts.salt(:,i);
        ref_meop_profile(i).temp = ref_meop_ts.temp(:,i);
        ref_meop_profile(i).N2 = ref_meop_ts.N2(:,i);
        ref_meop_profile(i).dyn_height_anom = ref_meop_ts.dyn_height_anom(:,i);

        %%% Calculating dynamic modes
        [ref_meop_profile(i).wmodes, ref_meop_profile(i).pmodes, ~, ~] = dynmodes(ref_meop_profile(i).N2(~isnan(ref_meop_profile(i).N2)), depth_grid(~isnan(ref_meop_profile(i).N2)),1);
    end

end
