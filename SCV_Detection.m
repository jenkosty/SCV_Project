
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
    %%%%%%%%%%%%%%%%%%%%%%% Bathymetry Check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Finding bathymetry along seal time series
    ts_bathymetry = interp2(RTOPO_lon, RTOPO_lat', RTOPO_bedrock_topography, qc_ts(tag_no).lon,qc_ts(tag_no).lat);
    
    %%% Calculating bathymetry of background profile
    background_bathymetry_mean = NaN(size(ts_bathymetry));
    background_bathymetry_std = NaN(size(ts_bathymetry));
    for i = 1:length(qc_ts(tag_no).cast)
        background_bathymetry_mean(i) = mean(ts_bathymetry([mean_ind{1,i}, mean_ind{2,i}]));
        background_bathymetry_std(i) = std(ts_bathymetry([mean_ind{1,i}, mean_ind{2,i}]));
    end

    %%% Checking each profile against assigned background profile
    bathymetry_change = abs(ts_bathymetry - background_bathymetry_mean);

    %%% Flagging profiles for removal based on bathymetry check
    flagit = find((bathymetry_change > 500) | (background_bathymetry_std > 500));

    %%% Removing flagged profiles 
    meop_ts.cast(flagit) = [];
    meop_ts.lat(flagit) = [];
    meop_ts.lon(flagit) = [];
    meop_ts.time(flagit,:) = [];
    meop_ts.salt(:,flagit) = [];
    meop_ts.temp(:,flagit) = [];

    clear ts_bathymetry background_bathymetry_std background_bathymetry_mean bathymetry_change flagit i

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% Dynamic Height Anomaly Check %%%%%%%%%%%%%%%%%%%%

    depth_grid = 1:500;

    %%% Calculating absolute salinity and conservative temperature
    ts_salt_absolute = gsw_SA_from_SP(meop_ts.salt, depth_grid', meop_ts.lon, meop_ts.lat);
    ts_temp_conservative = gsw_CT_from_t(ts_salt_absolute, meop_ts.temp, depth_grid');
    

    for i = 1:length(meop_ts.cast)
        %%% Finding maximum depth
        ts_max_depth(i) = find(~isnan(meop_ts.temp(:,i)),1,'last');

        %%% calculating dynamic height anomaly
        ts_dynamic_height(:,i) = gsw_geo_strf_dyn_height(ts_salt_absolute(:,i), ts_temp_conservative(:,i), depth_grid', ts_max_depth(i));
    end

end
