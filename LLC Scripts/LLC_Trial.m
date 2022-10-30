clear; close all; clc

%%% Loading MEOP seal data
load("qc_ts.mat")
test_prof = 4;

%%
%%% Loading LLC data

LLC_1.mat = matfile('/Volumes/Elements/LLC_Snapshots/SO_snapshots_NSF_LLC4320_k1-86_face1_01-Dec-2011.mat');
LLC_2.mat = matfile('/Volumes/Elements/LLC_Snapshots/SO_snapshots_NSF_LLC4320_k1-86_face2_01-Dec-2011.mat');
LLC_4.mat = matfile('/Volumes/Elements/LLC_Snapshots/SO_snapshots_NSF_LLC4320_k1-86_face4_01-Dec-2011.mat');
LLC_5.mat = matfile('/Volumes/Elements/LLC_Snapshots/SO_snapshots_NSF_LLC4320_k1-86_face5_01-Dec-2011.mat');

LLC_1.edge_lats = [double(LLC_1.mat.yc(1,:)), double(LLC_1.mat.yc(:,end))', flip(double(LLC_1.mat.yc(end,:))), flip(double(LLC_1.mat.yc(:,1)))'];
LLC_1.edge_lons = [double(LLC_1.mat.xc(1,:)), double(LLC_1.mat.xc(:,end))', flip(double(LLC_1.mat.xc(end,:))), flip(double(LLC_1.mat.xc(:,1)))'];
LLC_1.polygon = geopolyshape(LLC_1.edge_lats, LLC_1.edge_lons);
LLC_1.lats = double(LLC_1.mat.yc);
LLC_1.lons = double(LLC_1.mat.xc);

LLC_2.edge_lats = [double(LLC_2.mat.yc(1,:)), double(LLC_2.mat.yc(:,end))', flip(double(LLC_2.mat.yc(end,:))), flip(double(LLC_2.mat.yc(:,1)))'];
LLC_2.edge_lons = [double(LLC_2.mat.xc(1,:)), double(LLC_2.mat.xc(:,end))', flip(double(LLC_2.mat.xc(end,:))), flip(double(LLC_2.mat.xc(:,1)))'];
LLC_2.polygon = geopolyshape(LLC_2.edge_lats, LLC_2.edge_lons);
LLC_2.lats = double(LLC_2.mat.yc);
LLC_2.lons = double(LLC_2.mat.xc);

LLC_4.edge_lats = [double(LLC_4.mat.yc(1,:)), double(LLC_4.mat.yc(:,end))', flip(double(LLC_4.mat.yc(end,:))), flip(double(LLC_4.mat.yc(:,1)))'];
LLC_4.edge_lons = [double(LLC_4.mat.xc(1,:)), double(LLC_4.mat.xc(:,end))', flip(double(LLC_4.mat.xc(end,:))), flip(double(LLC_4.mat.xc(:,1)))'];
LLC_4.polygon = geopolyshape(LLC_4.edge_lats, LLC_4.edge_lons);
LLC_4.lats = double(LLC_4.mat.yc);
LLC_4.lons = double(LLC_4.mat.xc);

LLC_5.edge_lats = [double(LLC_5.mat.yc(1,:)), double(LLC_5.mat.yc(:,end))', flip(double(LLC_5.mat.yc(end,:))), flip(double(LLC_5.mat.yc(:,1)))'];
LLC_5.edge_lons = [double(LLC_5.mat.xc(1,:)), double(LLC_5.mat.xc(:,end))', flip(double(LLC_5.mat.xc(end,:))), flip(double(LLC_5.mat.xc(:,1)))'];
LLC_5.polygon = geopolyshape(LLC_5.edge_lats, LLC_5.edge_lons);
LLC_5.lats = double(LLC_5.mat.yc);
LLC_5.lons = double(LLC_5.mat.xc);

%%% Calculating the Okubo Weiss Parameter
ind = [27 32];
LLC_1.OW = Okubo_Weiss(LLC_1.mat.dxc(:,:), LLC_1.mat.dyc(:,:), LLC_1.mat.u(:,:,ind), LLC_1.mat.v(:,:,ind));
LLC_2.OW = Okubo_Weiss(LLC_2.mat.dxc(:,:), LLC_2.mat.dyc(:,:), LLC_2.mat.u(:,:,ind), LLC_2.mat.v(:,:,ind));
LLC_4.OW = Okubo_Weiss(LLC_4.mat.dxc(:,:), LLC_4.mat.dyc(:,:), LLC_4.mat.u(:,:,ind), LLC_4.mat.v(:,:,ind));
LLC_5.OW = Okubo_Weiss(LLC_5.mat.dxc(:,:), LLC_5.mat.dyc(:,:), LLC_5.mat.u(:,:,ind), LLC_5.mat.v(:,:,ind));

%%
test_prof = 1:length(qc_ts);

for tag_no = test_prof
    
    clear LLCseal

    %%% Grabbing data on MEOP seal time series
    LLCseal.tag = qc_ts(tag_no).tag;
    LLCseal.cast = qc_ts(tag_no).cast;
    LLCseal.lat = qc_ts(tag_no).lat;
    LLCseal.lon = qc_ts(tag_no).lon;
    LLCseal.time = qc_ts(tag_no).time;

    sector = NaN(1,length(LLCseal.cast));
    %%% Identifying the LLC sector of each MEOP profile
    for i = 1:length(LLCseal.cast)
        if isinterior(LLC_1.polygon, geopointshape(LLCseal.lat(i), LLCseal.lon(i)))
            sector(i) = 1;
        elseif isinterior(LLC_2.polygon, geopointshape(LLCseal.lat(i), LLCseal.lon(i)))
            sector(i) = 2;
        elseif isinterior(LLC_5.polygon, geopointshape(LLCseal.lat(i), LLCseal.lon(i)))
            sector(i) = 5;
        else
            sector(i) = 4;
        end
    end

    LLCseal.sector = sector;

    %%% Initializing matrices for interpolated data
    LLCseal_salt = NaN(size(qc_ts(tag_no).salt));
    LLCseal_temp = NaN(size(qc_ts(tag_no).temp));
    LLCseal_vort = NaN(size(qc_ts(tag_no).temp));
    LLCseal_OW = NaN(size(LLC_1.OW, 3), length(qc_ts(tag_no).cast));

    for i = 1:length(qc_ts(tag_no).cast)

        %%% Grabbing LLC data for profile sector
        if sector(i) == 1
            LLC = LLC_1;
        elseif sector(i) == 2
            LLC = LLC_2;
        elseif sector(i) == 4
            LLC = LLC_4;
        elseif sector(i) == 5
            LLC = LLC_5;
        end

        %%% Finding LLC points close to MEOP profile
        LLC_lats = LLC.lats;
        LLC_lons = LLC.lons;
        nan_ind = LLC_lats > LLCseal.lat(i) + 0.025 | LLC_lats < LLCseal.lat(i) - 0.025 | LLC_lons > LLCseal.lon(i) + 0.04 | LLC_lons < LLCseal.lon(i) - 0.04;
        LLC_lats(nan_ind) = NaN;
        close_ind = find(~isnan(LLC_lats));
        [rows, cols] = ind2sub(size(LLC_lats), close_ind);

        lats_unfmt = NaN(1,length(rows));
        lons_unfmt = NaN(1,length(rows));
        salt = NaN(86, length(rows));
        temp = NaN(86, length(rows));
        vort = NaN(86, length(rows));
        okubo_weiss = NaN(size(LLC.OW, 3), length(rows));

        %%% Extracting LLC data close to MEOP profile
        for j = 1:length(rows)
            lats_unfmt(j) = double(LLC.mat.yc(rows(j), cols(j)));
            lons_unfmt(j) = double(LLC.mat.xc(rows(j), cols(j)));
            salt(:,j) = squeeze(LLC.mat.s(rows(j), cols(j), :));
            temp(:,j) = squeeze(LLC.mat.t(rows(j), cols(j), :));
            vort(:,j) = squeeze(LLC.mat.Ro(rows(j), cols(j), :));
            okubo_weiss(:,j) = squeeze(LLC.OW(rows(j), cols(j),:));
        end
        temp(salt == 0) = NaN;
        vort(salt == 0) = NaN;
        salt(salt == 0) = NaN;

        lats = double(lats_unfmt .* ones(size(salt)));
        lons = double(lons_unfmt .* ones(size(salt)));
        depths = LLC.mat.rc(1:size(salt, 1), 1 ) .* ones(size(salt));
        
        %%% Interpolating LLC data
        S = scatteredInterpolant(lats(:), lons(:), depths(:), salt(:));
        LLCseal_salt(:,i) = S(qc_ts(tag_no).lat(i).*ones(size(depth_grid)), qc_ts(tag_no).lon(i).*ones(size(depth_grid)), -depth_grid);

        T = scatteredInterpolant(lats(:), lons(:), depths(:), temp(:));
        LLCseal_temp(:,i) = T(qc_ts(tag_no).lat(i).*ones(size(depth_grid)), qc_ts(tag_no).lon(i).*ones(size(depth_grid)), -depth_grid);

        V = scatteredInterpolant(lats(:), lons(:), depths(:), vort(:));
        LLCseal_vort(:,i) = V(qc_ts(tag_no).lat(i).*ones(size(depth_grid)), qc_ts(tag_no).lon(i).*ones(size(depth_grid)), -depth_grid);

        for j = 1:size(LLC.OW, 3)
            OW = scatteredInterpolant(lats_unfmt', lons_unfmt', okubo_weiss(j,:)');
            LLCseal_OW(j,i) = OW(qc_ts(tag_no).lat(i), qc_ts(tag_no).lon(i));
        end

        clear lats_unfmt lons_unfmt salt vort temp okubo_weiss cols rows nan_ind lats lons depths S T V OW j close_ind 
    end

    LLCseal.salt = LLCseal_salt;
    LLCseal.temp = LLCseal_temp;
    LLCseal.vort = LLCseal_vort;
    LLCseal.OW = LLCseal_OW;

    %%% Using Okubo Weiss to flag profiles as eddies
    for i = 1:length(LLCseal.cast)
        if abs(LLCseal.OW(1,i)) > 0.4e-8 && abs(LLCseal.OW(2,i)) > 0.4e-8
            LLCseal.scv(i) = 1;
        else
            LLCseal.scv(i) = 0;
        end
    end

    %%% Saving interpolated LLC data to a structure
    LLCsealdata(tag_no) = LLCseal;

    clear sector i LLCseal LLCseal_salt LLCseal_temp LLCseal_vort LLCseal_OW LLC_lats LLC_lons
end

save("LLCsealdata", "LLCsealdata")

%%

for tag_no = test_prof
    
    %%% Calculating Bathymetry
%     LLC_data(tag_no).bathymetry = interp2(RTOPO.lon, RTOPO.lat', RTOPO.bedrock_topography, qc_ts(tag_no).lon, qc_ts(tag_no).lat);

    %%% Creating pressure space structure
    tmp_pres_space.pres = depth_grid .* ones(size(LLCsealdata(tag_no).salt));
    
    %%% Assigning temperature and salinity data to new pressure space
    %%% structure
    tmp_pres_space.salt = LLCsealdata(tag_no).salt;
    tmp_pres_space.temp = LLCsealdata(tag_no).temp;
    tmp_pres_space.vort = LLCsealdata(tag_no).vort;
    
    %%% Calculating bottom pressure
    for i = 1:length(LLCsealdata(tag_no).cast)
        prof_pres = tmp_pres_space.pres(~isnan(tmp_pres_space.salt(:,i)),i);
        if isempty(prof_pres)
            tmp_pres_space.bot_prof_pres = NaN;
        else
            tmp_pres_space.prof_bot_pres(i) = max(prof_pres); 
        end
    end

    %%% Calculating absolute salinity and conservative temperature
    tmp_pres_space.salt_absolute = gsw_SA_from_SP(tmp_pres_space.salt, depth_grid, LLCsealdata(tag_no).lon, LLCsealdata(tag_no).lat);
    tmp_pres_space.temp_conservative = gsw_CT_from_t(tmp_pres_space.salt_absolute, tmp_pres_space.temp, tmp_pres_space.pres);

    %%% Calculating density
    tmp_pres_space.density = gsw_rho_CT_exact(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative, ones(size(tmp_pres_space.salt)).*400);

    %%% Calculating potential density anomaly
    tmp_pres_space.sigma0 = gsw_sigma0(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative);
    
    %%% Calculating neutral density
    %[tmp_pres_space.gamma_n,~,~] = eos80_legacy_gamma_n(tmp_pres_space.salt, tmp_pres_space.temp, tmp_pres_space.pres, qc_ts(tag_no).lon, qc_ts(tag_no).lat);

    %%% Calculating N^2
    [N2, mid_pres] = gsw_Nsquared(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative, tmp_pres_space.pres, LLCsealdata(tag_no).lat .* ones(size(tmp_pres_space.salt)));
    for i = 1:length(qc_ts(tag_no).cast)
        tmp_pres_space.N2(:,i) = interp1(mid_pres(:,i), N2(:,i), tmp_pres_space.pres(:,i));
    end

    %%% Calculating Spice
    tmp_pres_space.spice = gsw_spiciness0(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative);
    
    %%% Assigning variables to structure
    LLCsealdata(tag_no).ps = tmp_pres_space;

    clear ts_isopycnal_sep ts_pres u pres_final pres j isopycnal_sep isopycnal_sep_ds isopycnal_sep_y_axis...
            i a b density depths density_final k isopycnal_sep_final pres_ds ts_density tmp_pres_space prof_pres

end

clear N2 mid_pres

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolating to Density Grid %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Creating density grid
density_grid = (26.6:0.001:28.3)';

for tag_no = test_prof
    
    %%% Creating matrices for density-interpolated data
    interp_salt = NaN(length(density_grid), length(LLCsealdata(tag_no).cast));
    interp_temp = NaN(length(density_grid), length(LLCsealdata(tag_no).cast));
    interp_N2 = NaN(length(density_grid), length(LLCsealdata(tag_no).cast));
    interp_spice = NaN(length(density_grid), length(LLCsealdata(tag_no).cast));
    interp_pres = NaN(length(density_grid), length(LLCsealdata(tag_no).cast));
    
    for i = 1:length(LLCsealdata(tag_no).cast)
        
        %%% Removing NaNs
        tmp_sigma0 = LLCsealdata(tag_no).ps.sigma0(~isnan(LLCsealdata(tag_no).ps.sigma0(:,i)),i);
        tmp_salt = LLCsealdata(tag_no).ps.salt(~isnan(LLCsealdata(tag_no).ps.salt(:,i)),i);
        tmp_temp = LLCsealdata(tag_no).ps.temp(~isnan(LLCsealdata(tag_no).ps.temp(:,i)),i);
        tmp_N2 = LLCsealdata(tag_no).ps.N2(~isnan(LLCsealdata(tag_no).ps.N2(:,i)),i);
        tmp_sigma0_N2 = LLCsealdata(tag_no).ps.sigma0(~isnan(LLCsealdata(tag_no).ps.sigma0(:,i)) & ~isnan(LLCsealdata(tag_no).ps.N2(:,i)),i);
        tmp_spice = LLCsealdata(tag_no).ps.spice(~isnan(LLCsealdata(tag_no).ps.spice(:,i)),i);
        tmp_pres = depth_grid(~isnan(LLCsealdata(tag_no).ps.sigma0(:,i)));
        
        if isempty(tmp_salt)
            interp_salt(:,i) = NaN(length(density_grid), 1);
            interp_temp(:,i) = NaN(length(density_grid), 1);
            interp_N2(:,i) = NaN(length(density_grid), 1);
            interp_spice(:,i) = NaN(length(density_grid), 1);
            interp_pres(:,i) = NaN(length(density_grid), 1);
        else
            %%% Interpolating data
            interp_salt(:,i) = interp1(tmp_sigma0, tmp_salt, density_grid);
            interp_temp(:,i) = interp1(tmp_sigma0, tmp_temp, density_grid);
            interp_N2(:,i) = interp1(tmp_sigma0_N2, tmp_N2, density_grid);
            interp_spice(:,i) = interp1(tmp_sigma0, tmp_spice, density_grid);
            interp_pres(:,i) = interp1(tmp_sigma0, tmp_pres, density_grid);
        end
    end
    
    %%% Saving density-interpolated data
    tmp_sigma0_space.salt = interp_salt;
    tmp_sigma0_space.temp = interp_temp;
    tmp_sigma0_space.N2 = interp_N2;
    tmp_sigma0_space.spice = interp_spice;
    tmp_sigma0_space.pres = interp_pres;
    
    %%% Saving structure
    LLCsealdata(tag_no).ds = tmp_sigma0_space;
    
    clear tmp_density tmp_salt tmp_temp tmp_N2 tmp_density_N2 tmp_spice tmp_pres ...
        interp_salt interp_temp interp_N2 interp_spice interp_pres tmp_density_space ...
        tmp_sigma0 tmp_sigma0_N2 tmp_sigma0_space i
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Isopycnal Separation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = 1:length(LLCsealdata)
    
    isopycnal_separation = NaN(length(density_grid), length(LLCsealdata(tag_no).cast));
    
    for i = 1:length(LLCsealdata(tag_no).cast)
        for j = 1:length(density_grid)
            isopycnal_separation(j,i) = LLCsealdata(tag_no).ds.pres(min(j+1, length(density_grid)), i) - LLCsealdata(tag_no).ds.pres(max(j-1, 1), i);
        end
    end
    
    LLCsealdata(tag_no).ds.isopycnal_separation = isopycnal_separation;
    
    clear isopycnal_separation i j 
end

%%
%%%%%%%%%%%%%%%%%%
%%% Map figure %%%
%%%%%%%%%%%%%%%%%%

for tag_no = test_prof

    lat_min = min(LLCsealdata(tag_no).lat) - 0.1;
    lat_max = max(LLCsealdata(tag_no).lat) + 0.1;
    lon_min = min(LLCsealdata(tag_no).lon) - 0.5;
    lon_max = max(LLCsealdata(tag_no).lon) + 0.5;

    if lon_max - lon_min > 360
        pos_min = min(LLCsealdata(tag_no).lon(LLCsealdata(tag_no).lon > 0));
        lon_min = -180 - (180-pos_min) - 0.5;
        lon_max = max(LLCsealdata(tag_no).lon(LLCsealdata(tag_no).lon < 0)) + 0.5;
    end

    %%% Figure settings
    load coastlines
    figure('Position', [500 100 1000 850])
    axesm('lambertstd', 'MapParallels',[-75 -15],'MapLatLimit',[lat_min lat_max], ...
        'MapLonLimit', [lon_min lon_max], 'MLineLocation', 30, 'PLineLocation', 10, 'FontSize',10);
    axis off; framem on; gridm on; mlabel on; plabel on;
    colormap(cmocean('balance')); colorbar; clim([-1e-8 1e-8]);

    %%% Plotting data
    hold on

    %%% Grabbing sector of LLC time series
    for j = unique(LLCsealdata(tag_no).sector)
        if j == 1
            LLC = LLC_1;
            pcolorm(LLC.lats(:,:), LLC.lons(:,:), LLC.OW(:,:, 2));
        elseif j == 2
            LLC = LLC_2;
            pcolorm(LLC.lats(:,:), LLC.lons(:,:), LLC.OW(:,:, 2));
        elseif j == 4
            LLC = LLC_4;
            pcolorm(LLC.lats(:,:), LLC.lons(:,:), LLC.OW(:,:, 2));
        elseif j == 5
            LLC = LLC_5;
            pcolorm(LLC.lats(:,:), LLC.lons(:,:), LLC.OW(:,:, 2));
        end
    end
    geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor','white');
    plotm(LLCsealdata(tag_no).lat, LLCsealdata(tag_no).lon, '-s', 'Color', [0.5 0.5 0.5],...
        'MarkerSize', 4, 'MarkerFaceColor', [0.5 0.5 0.5], 'LineWidth', 2)
    scatterm(LLCsealdata(tag_no).lat(LLCsealdata(tag_no).scv == 1), LLCsealdata(tag_no).lon(LLCsealdata(tag_no).scv == 1), 6, 'gs', 'MarkerFaceColor', 'g')

end

clear lat_min lat_max lon_min lon_max pos_min coastlat coastlon
