
RTOPO_lat = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'lat'));
RTOPO_lon = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'lon'))';
RTOPO_bedrock_topography = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'bedrock_topography'))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inner_window = 2; %%% setting inner "exclusion" window for background calculation
outer_window = 7; %%% setting outer "inclusion" window for background calculation

for tag_no = 362
    clear meop_ts

    %%% Creating new structure for current seal
    meop_ts = qc_ts(tag_no);

    %%% Finding indices to use for background profile calculation
    for i = 1:length(qc_ts(tag_no).cast)
        mean_ind(1,i) = {(i-outer_window):(i-inner_window)};
        mean_ind{1,i}(mean_ind{1,i} < 1) = []; %%% Making sure indices remain within ts boundaries
        mean_ind(2,i) = {(i+inner_window):(i+outer_window)};
        mean_ind{2,i}(mean_ind{2,i} > length(qc_ts(tag_no).cast)) = []; %%% Making sure indices remain within ts boundaries
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Calculating Reference Profiles %%%%%%%%%%%%%%%%%%%%%%%%

    %%% Creating reference profiles against which the anomalies can be
    %%% calculated
    for i = 1:length(meop_ts.cast)
        meop_ts.ref.salt(:,i) = mean(meop_ts.salt(:,[mean_ind{1,i} mean_ind{2,i}]), 2, 'omitnan');
        meop_ts.ref.temp(:,i) = mean(meop_ts.temp(:,[mean_ind{1,i} mean_ind{2,i}]), 2, 'omitnan');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Bathymetry Check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Creates new time series with profiles that may be checked for SCVs

    %%% Finding bathymetry along seal time series
    ts_bathymetry = interp2(RTOPO.lon, RTOPO.lat', RTOPO.bedrock_topography, meop_ts.lon, meop_ts.lat);

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

    %%% Removing corresponding reference profiles
    meop_ts.ref.salt(:,flagit) = [];
    meop_ts.ref.temp(:,flagit) = [];

    clear ts_bathymetry ref_bathymetry_std ref_bathymetry_mean bathymetry_change i
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Calculating Pressure Space Variables %%%%%%%%%%%%%%%%%%%
    
    %%% Creating pressure time series
    depth_grid = 1:500;
    depth_grid = depth_grid';
    meop_ts.pres = depth_grid .* ones(size(meop_ts.salt));

    %%% Calculating absolute salinity and conservative temperature
    meop_ts.salt_absolute = gsw_SA_from_SP(meop_ts.salt, meop_ts.pres, meop_ts.lon, meop_ts.lat);
    meop_ts.temp_conservative = gsw_CT_from_t(meop_ts.salt_absolute, meop_ts.temp, meop_ts.pres);

    for i = 1:length(meop_ts.cast)
        %%% Finding maximum depth
        meop_ts.max_depth(i) = find(~isnan(meop_ts.temp(:,i)),1,'last');

        %%% calculating dynamic height anomaly
        meop_ts.dyn_height_anom(:,i) = gsw_geo_strf_dyn_height(meop_ts.salt_absolute(:,i), meop_ts.temp_conservative(:,i), meop_ts.pres(:,i), meop_ts.max_depth(i));
    end

    %%% Calculating density
    meop_ts.density = gsw_rho(meop_ts.salt_absolute, meop_ts.temp_conservative, meop_ts.pres);

    %%% Calculating Potential Density Anomaly
    meop_ts.sigma0 = gsw_sigma0(meop_ts.salt_absolute, meop_ts.temp_conservative);

    %%% Calculating N^2
    [meop_ts.N2,~] = gsw_Nsquared(meop_ts.salt_absolute, meop_ts.temp_conservative, meop_ts.pres, meop_ts.lat .* ones(size(meop_ts.salt)));
    
    %%% Creating density grid for interpolation
    density_grid(:,1) = linspace(1026.5, 1030, 500);

    %%% Interpolating into density space
    for i = 1:length(meop_ts.cast)
        meop_ts.ds.salt(:,i) = interp1(meop_ts.density(~isnan(meop_ts.density(:,i)),i), meop_ts.salt(~isnan(meop_ts.density(:,i)),i), density_grid);
        meop_ts.ds.temp(:,i) = interp1(meop_ts.density(~isnan(meop_ts.density(:,i)),i), meop_ts.temp(~isnan(meop_ts.density(:,i)),i), density_grid);
        meop_ts.ds.pres(:,i) = interp1(meop_ts.density(~isnan(meop_ts.density(:,i)),i), meop_ts.pres(~isnan(meop_ts.density(:,i)),i), density_grid);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Calculating Reference Pressure Space Variables %%%%%%%%%%%%%

    %%% Creating pressure time series
    meop_ts.ref.pres = depth_grid .* ones(size(meop_ts.salt));

    %%% Calculating absolute salinity and conservative temperature
    meop_ts.ref.salt_absolute = gsw_SA_from_SP(meop_ts.ref.salt, meop_ts.ref.pres, meop_ts.lon, meop_ts.lat);
    meop_ts.ref.temp_conservative = gsw_CT_from_t(meop_ts.ref.salt_absolute, meop_ts.ref.temp, meop_ts.ref.pres);

    for i = 1:length(meop_ts.cast)
        %%% Finding maximum depth
        meop_ts.ref.max_depth(i) = find(~isnan(meop_ts.ref.temp(:,i)),1,'last');

        %%% calculating dynamic height anomaly
        meop_ts.ref.dyn_height_anom(:,i) = gsw_geo_strf_dyn_height(meop_ts.ref.salt_absolute(:,i), meop_ts.ref.temp_conservative(:,i), meop_ts.ref.pres(:,i), meop_ts.ref.max_depth(i));
    end

    %%% Calculating density
    meop_ts.ref.density = gsw_rho(meop_ts.ref.salt_absolute, meop_ts.ref.temp_conservative, meop_ts.ref.pres);

    %%% Calculating Potential Density Anomaly
    meop_ts.ref.sigma0 = gsw_sigma0(meop_ts.ref.salt_absolute, meop_ts.ref.temp_conservative);

    %%% Calculating N^2
    [meop_ts.ref.N2,~] = gsw_Nsquared(meop_ts.ref.salt_absolute, meop_ts.ref.temp_conservative, meop_ts.ref.pres, meop_ts.lat .* ones(size(meop_ts.ref.salt)));

    %%% Interpolating into density space
    for i = 1:length(meop_ts.cast)
        meop_ts.ds.ref.salt(:,i) = interp1(meop_ts.ref.density(~isnan(meop_ts.ref.density(:,i)),i), meop_ts.ref.salt(~isnan(meop_ts.ref.density(:,i)),i), density_grid);
        meop_ts.ds.ref.temp(:,i) = interp1(meop_ts.ref.density(~isnan(meop_ts.ref.density(:,i)),i), meop_ts.ref.temp(~isnan(meop_ts.ref.density(:,i)),i), density_grid);
        meop_ts.ds.ref.pres(:,i) = interp1(meop_ts.ref.density(~isnan(meop_ts.ref.density(:,i)),i), meop_ts.ref.pres(~isnan(meop_ts.ref.density(:,i)),i), density_grid);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Calculating Density Space Variables %%%%%%%%%%%%%%%%%%

    %%% Calculating absolute salinity and conservative temperature
    meop_ts.ds.salt_absolute = gsw_SA_from_SP(meop_ts.ds.salt, meop_ts.ds.pres, meop_ts.lon, meop_ts.lat);
    meop_ts.ds.temp_conservative = gsw_CT_from_t(meop_ts.ds.salt_absolute, meop_ts.ds.temp, meop_ts.ds.pres);

    %%% Calculating Potential Density Anomaly
    %meop_ts.ds.sigma0 = gsw_sigma0(meop_ts.ds.salt_absolute, meop_ts.ds.temp_conservative);

    %%% Calculating Spice
    meop_ts.ds.spice = gsw_spiciness0(meop_ts.ds.salt_absolute, meop_ts.ds.temp_conservative);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Calculating Reference Density Space Variables %%%%%%%%%%%%%%%
    
    %%% Calculating absolute salinity and conservative temperature
    meop_ts.ds.ref.salt_absolute = gsw_SA_from_SP(meop_ts.ds.ref.salt, meop_ts.ds.ref.pres, meop_ts.lon, meop_ts.lat);
    meop_ts.ds.ref.temp_conservative = gsw_CT_from_t(meop_ts.ds.ref.salt_absolute, meop_ts.ds.ref.temp, meop_ts.ds.ref.pres);

    %%% Calculating Spice
    meop_ts.ds.ref.spice = gsw_spiciness0(meop_ts.ds.ref.salt_absolute, meop_ts.ds.ref.temp_conservative);

 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%% Figure to see Profile vs Reference %%%%%%%%%%%%%%%%%%%%%
%     
%         figure()
%         subplot(211);
%         pcolor(datenum(meop_ts.time), depth_grid, meop_ts.temp)
%         datetick('x', 'dd/mm');
%         set(gca,'YDir','Reverse');
%         shading interp
%     
%         subplot(212);
%         pcolor(datenum(meop_ts.time), depth_grid, ref_meop_ts.temp)
%         datetick('x', 'dd/mm');
%         set(gca,'YDir','Reverse');
%         shading interp

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

    clear ref_meop_profile meop_profile

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Isolating profiles to be checked for SCVs %%%%%%%%%%%%%%%%%%%

    for i = 1:length(meop_ts.cast)
        %%% Creating structure of individual MEOP profiles
        meop_profile(i).salt = meop_ts.salt(:,i);
        meop_profile(i).temp = meop_ts.temp(:,i);
        meop_profile(i).dyn_height_anom = meop_ts.dyn_height_anom(:,i);
        meop_profile(i).sigma0 = meop_ts.sigma0(:,i);
        meop_profile(i).N2 = meop_ts.N2(:,i);
        meop_profile(i).pres = meop_ts.pres(:,i);

        meop_profile(i).ds.salt = meop_ts.ds.salt(:,i);
        meop_profile(i).ds.temp = meop_ts.ds.temp(:,i);
        meop_profile(i).ds.pres = meop_ts.ds.pres(:,i);
        meop_profile(i).ds.spice = meop_ts.ds.spice(:,i);

        %%% Creating structure of individual reference profiles
        meop_profile(i).ref.salt = meop_ts.ref.salt(:,i);
        meop_profile(i).ref.temp = meop_ts.ref.temp(:,i);
        meop_profile(i).ref.dyn_height_anom = meop_ts.ref.dyn_height_anom(:,i);
        meop_profile(i).ref.sigma0 = meop_ts.ref.sigma0(:,i);
        meop_profile(i).ref.N2 = meop_ts.ref.N2(:,i);

        meop_profile(i).ds.ref.salt = meop_ts.ds.ref.salt(:,i);
        meop_profile(i).ds.ref.temp = meop_ts.ds.ref.temp(:,i);
        meop_profile(i).ds.ref.pres = meop_ts.ds.ref.pres(:,i);
        meop_profile(i).ds.ref.spice = meop_ts.ds.ref.spice(:,i);

        %%% Calculating anomalies
        meop_profile(i).salt_anom = meop_profile(i).salt - meop_profile(i).ref.salt;
        meop_profile(i).temp_anom = meop_profile(i).temp - meop_profile(i).ref.temp;
        meop_profile(i).sigma0_anom = meop_profile(i).sigma0 - meop_profile(i).ref.sigma0;

        %%% Calculating anomalies
        meop_profile(i).ds.salt_anom = meop_profile(i).ds.salt - meop_profile(i).ds.ref.salt;
        meop_profile(i).ds.temp_anom = meop_profile(i).ds.temp - meop_profile(i).ds.ref.temp;
        meop_profile(i).ds.spice_anom = meop_profile(i).ds.spice - meop_profile(i).ds.ref.spice;
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% Gaussian Fit Check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 27

        % Grab amplitude and depth of max spice anomaly
        spike.A  = max(meop_profile(i).ds.spice_anom);
        spike.P = meop_profile(i).ds.pres(find(meop_profile(i).ds.spice_anom == spike.A));

        % Get range of allowable parameters
        prng = [-0.2:0.05:0.2];  % allow pressure peak to vary between +- 20% of height
        arng  = [0.8:0.05:1.2];  % allow amplitude range of +- 20% of spice anomaly peak
        hrng  = [0:10:500];     % allow height to vary  between 50 and 850m

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
        			zo = double(meop_profile(i).ds.pres - [spike.P + p*(4)*sqrt(h^2/2)]);
        			sa = double(meop_profile(i).ds.spice_anom);

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
        			dataX  = meop_profile(i).ds.spice_anom(pl <= meop_profile(i).ds.pres & meop_profile(i).ds.pres <= ph);
        			dataY  = meop_profile(i).ds.pres(pl <= meop_profile(i).ds.pres & meop_profile(i).ds.pres <= ph);
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
        zo    = double(meop_profile(i).ds.pres - [results.P]);
        zp    = zo + results.P;
        gauss = results.A*exp((-(zo.^2))/(results.H.^2));

        % Save final model
        results.X = gauss;
        results.Y = zp;

        % Plot
        figure()
        plot(meop_profile(i).ds.spice_anom,meop_profile(i).ds.pres,'k','linewidth',2)
        hold on; grid on; set(gca,'YDir','Reverse')
        plot(results.X,results.Y,'Color','Blue','LineWidth',3,'LineStyle','-.')
        xlabel('kg/m^3')
        ylabel('dbar');
        set(gca,'fontsize',10,'fontname','Helvetica')

    end
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% Dynamic Height Anomaly Check %%%%%%%%%%%%%%%%%%%%

    for i = 25

        %%% Extract vertical velocity and horizontal structure modes of climatology
        [ref_meop_profile(i).wmodes, ref_meop_profile(i).pmodes, ~, ~] = dynmodes(ref_meop_profile(i).N2(~isnan(ref_meop_profile(i).N2)), depth_grid(~isnan(ref_meop_profile(i).N2)),1);

        %%% Grab pressure levels of mode decomposition, add zero level (surface)
        ref_meop_profile(i).mode_pres = depth_grid(~isnan(ref_meop_profile(i).N2));
        ref_meop_profile(i).mode_pres = [0 ref_meop_profile(i).mode_pres]';

        % Interpolate 1st baroclinic mode to pressure of SCV cast
        meop_profile(i).dyn_pres = depth_grid(~isnan(meop_profile(i).dyn_height_anom));
        meop_profile(i).dyn_height_anom = meop_profile(i).dyn_height_anom(~isnan(meop_profile(i).dyn_height_anom));
        ref_meop_profile(i).BC1_data = ref_meop_profile(i).pmodes(:,1);
        ref_meop_profile(i).BC1_pres = ref_meop_profile(i).mode_pres;
        ref_meop_profile(i).BC1_data = interp1(ref_meop_profile(i).mode_pres(~isnan(ref_meop_profile(i).BC1_data)),ref_meop_profile(i).BC1_data(~isnan(ref_meop_profile(i).BC1_data)),meop_profile(i).dyn_pres)';
        ref_meop_profile(i).BC1_pres = meop_profile(i).dyn_pres';

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
        plot(meop_profile(i).dyn_height_anom_BC1,meop_profile(i).dyn_pres,':k','linewidth',2)
        [l,icons] = legend('DH''_{orig}','BC1_{fit}','DH''_{adj}','location','southeast');
        xlabel('m^2/s^2')
        title({'\it\bf\fontsize{8}\fontname{Helvetica}Dynamic Height','Anomaly'})
        set(gca,'YTickLabel',[])

        subplot(133)
        plot(meop_profile(i).sigma0 - ref_meop_profile(i).sigma0,meop_profile(i).pres,'k','linewidth',2)
        set(gca,'YDir','Reverse')
        grid on
        xlim([-0.2 0.2]);
    end


end


