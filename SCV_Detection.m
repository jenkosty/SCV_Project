
RTOPO_lat = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'lat'));
RTOPO_lon = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'lon'))';
RTOPO_bedrock_topography = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'bedrock_topography'))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inner_window = 2; %%% setting inner "exclusion" window for background calculation
outer_window = 7; %%% setting outer "inclusion" window for background calculation

for tag_no = 87

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

    %%% Calculating Potential Density Anomaly
    meop_ts.sigma0 = gsw_sigma0(meop_ts.salt_absolute, meop_ts.temp_conservative);

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
    meop_ts.sigma0(:,flagit) = [];

    %%% Removing corresponding reference profiles
    ref_meop_ts.salt(:,flagit) = [];
    ref_meop_ts.temp(:,flagit) = [];
    ref_meop_ts.N2(:,flagit) = [];
    ref_meop_ts.dyn_height_anom(:,flagit) = [];

    clear ts_bathymetry ref_bathymetry_std ref_bathymetry_mean bathymetry_change i

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Figure to see Profile vs Reference %%%%%%%%%%%%%%%%%%%%%
%     i = 1;
% 
%     figure()
%     subplot(1,2,1);
%     hold on
%     plot(meop_ts.N2(:,i), 1:499, 'k', 'LineWidth', 2)
%     plot(ref_meop_ts.N2(:,i), 1:499, 'r', 'LineWidth', 2)
%     hold off
%     set(gca,'YDir','Reverse');
%     title('N^2')
% 
%     subplot(1,2,2);
%     hold on
%     plot(meop_ts.dyn_height_anom(:,i), 1:500, 'k', 'LineWidth', 2)
%     plot(ref_meop_ts.dyn_height_anom(:,i), 1:500, 'r', 'LineWidth', 2)
%     hold off
%     set(gca,'YDir','Reverse');
%     title('Dynamic Height Anomaly')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% Dynamic Height Anomaly Check %%%%%%%%%%%%%%%%%%%%
    clear ref_meop_profile meop_profile

    for i = 1:10
        %%% Creating structure of reference profiles
        ref_meop_profile(i).salt = ref_meop_ts.salt(:,i);
        ref_meop_profile(i).temp = ref_meop_ts.temp(:,i);
        ref_meop_profile(i).N2 = ref_meop_ts.N2(:,i);
        ref_meop_profile(i).dyn_height_anom = ref_meop_ts.dyn_height_anom(:,i);

        %%% Extract vertical velocity and horizontal structure modes of climatology
        [ref_meop_profile(i).wmodes, ref_meop_profile(i).pmodes, ~, ~] = dynmodes(ref_meop_profile(i).N2(~isnan(ref_meop_profile(i).N2)), depth_grid(~isnan(ref_meop_profile(i).N2)),1);

        %%% Grab pressure levels of mode decomposition, add zero level (surface)
        ref_meop_profile(i).mode_pres = depth_grid(~isnan(ref_meop_profile(i).N2));
        ref_meop_profile(i).mode_pres = [0 ref_meop_profile(i).mode_pres]';

        %%% Creating structure of MEOP profiles
        meop_profile(i).salt = meop_ts.salt(:,i);
        meop_profile(i).temp = meop_ts.temp(:,i);
        meop_profile(i).N2 = meop_ts.N2(:,i);
        meop_profile(i).dyn_height_anom = meop_ts.dyn_height_anom(:,i);
        meop_profile(i).sigma0 = meop_ts.sigma0(:,i);
        meop_profile(i).pres = 1:500;

        % Interpolate 1st baroclinic mode to pressure of SCV cast
        meop_profile(i).dyn_pres = depth_grid(~isnan(meop_profile(i).dyn_height_anom));
        meop_profile(i).dyn_height_anom = meop_profile(i).dyn_height_anom(~isnan(meop_profile(i).dyn_height_anom));
        ref_meop_profile(i).BC1_data = ref_meop_profile(i).pmodes(:,1);
        ref_meop_profile(i).BC1_pres = ref_meop_profile(i).mode_pres;
        ref_meop_profile(i).BC1_data = interp1(ref_meop_profile(i).mode_pres(~isnan(ref_meop_profile(i).BC1_data)),ref_meop_profile(i).BC1_data(~isnan(ref_meop_profile(i).BC1_data)),meop_profile(i).dyn_pres)';
        ref_meop_profile(i).BC1_pres = meop_profile(i).dyn_pres';
    end

    for i = 9

    % Create function that describes residuals between projected BC1 and dyn_height_anom
    % Exclude data inbetween SCV limits for better fit to first mode
    dat = [];
    dat = [ref_meop_profile(i).BC1_data + meop_profile(i).dyn_height_anom];
    x_o = [];
    x_o = ref_meop_profile(i).BC1_data(~isnan(dat));
    x_p = [];
    x_p = ref_meop_profile(i).BC1_pres(~isnan(dat));
    x_f = [];
    x_f = meop_profile(i).dyn_height_anom(~isnan(dat));

%     % Get limits
%     pl = meop_profile(i).limits.shallow_pres;
%     ph = meop_profile(i).limits.deep_pres;
% 
%     % Remove values between upper/lower limits of SCV to avoid bad fit
%     ind = [];
%     ind = find(pl < x_p & x_p < ph);
%     if ind(end) == length(x_p)
%         ind = ind(1:end-1);
%     end
%     x_o(ind) = [];
%     x_f(ind) = [];

    % Remove mixed layer depths (Lynne Talley method, first density greater than 0.03 from sfc value
    ind      = [];
    mld_dens = meop_profile(i).sigma0(~isnan(meop_profile(i).sigma0));
    mld_pres = meop_profile(i).pres(~isnan(meop_profile(i).sigma0));
    ind      = find(mld_dens > mld_dens(1)+0.03);
    mld_pres = mld_pres(ind(1));
    ind      = find(x_p < mld_pres);
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
    x_o = ref_meop_profile(i).BC1_data(~isnan(dat));
    x_p = ref_meop_profile(i).mode_pres(~isnan(dat));
    x_f = meop_profile(i).dyn_height_anom(~isnan(dat));

    % Fix dynamic height anomaly by removing projected 1st mode, add back in barotopic mode
    meop_profile(i).dyn_height_anom_BC1 = [x_f] - [x_o*alpha + (x_f(end) - alpha*x_o(end))];
    meop_profile(i).dyn_height_pres_BC1 = meop_profile(i).dyn_pres(~isnan(dat));

    % Save VMD results
    ref_meop_profile(i).VMD.x_f      = x_f;
    ref_meop_profile(i).VMD.x_o      = x_o;
    ref_meop_profile(i).VMD.x_p      = x_p;
    ref_meop_profile(i).VMD.alpha    = alpha;

    % Get mode decomposition results
    BC1 = ref_meop_profile(i).VMD.x_o*ref_meop_profile(i).VMD.alpha;
    BC1 = BC1 - BC1(end); %// Set bottom to zero

    % Plot results
    	figure();
    	subplot(131)
    	plot(-ref_meop_profile(i).pmodes(:,1),ref_meop_profile(i).mode_pres,'r','linewidth',2)
    	hold on; grid on; set(gca,'YDir','Reverse')
    	plot(ref_meop_profile(i).pmodes(:,2),ref_meop_profile(i).mode_pres,'b','linewidth',2)
    	plot(ref_meop_profile(i).pmodes(:,3),ref_meop_profile(i).mode_pres,'g','linewidth',2)
    	plot(ref_meop_profile(i).pmodes(:,4),ref_meop_profile(i).mode_pres,'y','linewidth',2)
    	plot(ref_meop_profile(i).pmodes(:,5),ref_meop_profile(i).mode_pres,'color',[0.5 0 0.5],'linewidth',2)
    	title({'\it\bf\fontsize{8}\fontname{Helvetica}Horizontal Velocity','Modes'})
    	set(gca,'XTick',[0])
    	ylabel('Pressure (dbar)')
    	[l,icons] = legend('Mode-1','Mode-2','Mode-3','Mode-4','Mode-5','location','southeast');
    	l.Box = 'off';

    	subplot(132)
    	plot(meop_profile(i).dyn_height_anom,meop_profile(i).dyn_pres,'k','linewidth',2)
    	hold on; grid on; set(gca,'YDir','Reverse')
    	plot(BC1,ref_meop_profile(i).VMD.x_p,':r','linewidth',2)
    	plot(meop_profile(i).dyn_height_anom_BC1,meop_profile(i).dyn_pres(1:494),':k','linewidth',2)
    	[l,icons] = legend('DH''_{orig}','BC1_{fit}','DH''_{adj}','location','southeast');
    	xlabel('m^2/s^2')
    	title({'\it\bf\fontsize{8}\fontname{Helvetica}Dynamic Height','Anomaly'})
    	set(gca,'YTickLabel',[])

        subplot(133)
        plot(meop_profile(i).temp - ref_meop_profile(i).temp,meop_profile(i).pres,'k','linewidth',2)
        set(gca,'YDir','Reverse')
    end


end

