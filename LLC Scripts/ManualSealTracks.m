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

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generating fake seal data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plotting vorticity data
load coastlines
figure('Position', [500 100 1000 850])
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -57],'MLineLocation', 30, 'PLineLocation', 10, 'FontSize',10);
axis off; framem on; gridm on; mlabel on; plabel on;
colormap(cmocean('balance')); colorbar; clim([-0.3 0.3]);
hold on
pcolorm(LLC_5.lats(:,:), LLC_5.lons(:,:), LLC_5.mat.Ro(:,:,32));
geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor','white');
hold off

%%

%%% Selecting lat/lon to samply eddy
ind = 4;
points = inputm(25);

%%% Saving data
fakesealdata(ind).tag = ind;
fakesealdata(ind).lat = points(:,1);
fakesealdata(ind).lon = points(:,2);
fakesealdata(ind).cast = 1:length(points);
region = input("Approximate Region? ");
fakesealdata(ind).region = region;
eddysamples = input("How many data points were selected in the eddy center? ");
fakesealdata(ind).eddysamples = eddysamples;
centerindex = input("Eddy center index: ");
fakesealdata(ind).centerindex = centerindex;
rotationtype = input("Rossby number color? ");
fakesealdata(ind).rotationtype = rotationtype;

save("FakeSealData.mat", "fakesealdata")

clear region eddysamples centerindex rotationtype

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolating LLC data to create "time series" of selected lat/lons %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear LLCsealdata

for tag_no = 1:length(fakesealdata)
    
    clear LLCseal

    %%% Grabbing data on MEOP seal time series
    LLCseal.tag = fakesealdata(tag_no).tag;
    LLCseal.cast = fakesealdata(tag_no).cast;
    LLCseal.region = fakesealdata(tag_no).region;
    LLCseal.eddysamples = fakesealdata(tag_no).eddysamples;
    LLCseal.centerindex = fakesealdata(tag_no).centerindex;
    LLCseal.rotationtype = fakesealdata(tag_no).rotationtype;
    LLCseal.lat = fakesealdata(tag_no).lat;
    LLCseal.lon = fakesealdata(tag_no).lon;

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
    LLCseal_salt = NaN(length(depth_grid), size(fakesealdata(tag_no).lat, 1));
    LLCseal_temp = NaN(length(depth_grid), size(fakesealdata(tag_no).lat, 1));
    LLCseal_vort = NaN(length(depth_grid), size(fakesealdata(tag_no).lat, 1));

    for i = 1:length(LLCseal.cast)

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
        nan_ind = LLC_lats > LLCseal.lat(i) + 0.025 | LLC_lats < LLCseal.lat(i) - 0.025 | LLC_lons > LLCseal.lon(i) + 0.025 | LLC_lons < LLCseal.lon(i) - 0.025;
        LLC_lats(nan_ind) = NaN;
        LLC_lons(nan_ind) = NaN;
        close_ind = find(~isnan(LLC_lats));
        [rows, cols] = ind2sub(size(LLC_lats), close_ind);

        lats_unfmt = NaN(1,length(rows));
        lons_unfmt = NaN(1,length(rows));
        salt = NaN(86, length(rows));
        temp = NaN(86, length(rows));
        vort = NaN(86, length(rows));

        %%% Extracting LLC data close to MEOP profile
        for j = 1:length(rows)
            lats_unfmt(j) = double(LLC.mat.yc(rows(j), cols(j)));
            lons_unfmt(j) = double(LLC.mat.xc(rows(j), cols(j)));
            salt(:,j) = squeeze(LLC.mat.s(rows(j), cols(j), :));
            temp(:,j) = squeeze(LLC.mat.t(rows(j), cols(j), :));
            vort(:,j) = squeeze(LLC.mat.Ro(rows(j), cols(j), :));
        end
        temp(salt == 0) = NaN;
        vort(salt == 0) = NaN;
        salt(salt == 0) = NaN;

        lats = double(lats_unfmt .* ones(size(salt)));
        lons = double(lons_unfmt .* ones(size(salt)));
        depths = LLC.mat.rc(1:size(salt, 1), 1 ) .* ones(size(salt));
        
        %%% Interpolating LLC data
        S = scatteredInterpolant(lats(:), lons(:), depths(:), salt(:));
        LLCseal_salt(:,i) = S(LLCseal.lat(i).*ones(size(depth_grid)), LLCseal.lon(i).*ones(size(depth_grid)), -depth_grid);

        T = scatteredInterpolant(lats(:), lons(:), depths(:), temp(:));
        LLCseal_temp(:,i) = T(LLCseal.lat(i).*ones(size(depth_grid)), LLCseal.lon(i).*ones(size(depth_grid)), -depth_grid);

        V = scatteredInterpolant(lats(:), lons(:), depths(:), vort(:));
        LLCseal_vort(:,i) = V(LLCseal.lat(i).*ones(size(depth_grid)), LLCseal.lon(i).*ones(size(depth_grid)), -depth_grid);

        clear lats_unfmt lons_unfmt salt vort temp cols rows nan_ind lats lons depths S T V j close_ind 
    end

    LLCseal.salt = LLCseal_salt;
    LLCseal.temp = LLCseal_temp;
    LLCseal.vort = LLCseal_vort;

    %%% Saving interpolated LLC data to a structure
    LLCsealdata(tag_no) = LLCseal;

    clear sector i LLCseal
end

%%
for tag_no = 1:length(LLCsealdata)
    
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
        tmp_pres_space.prof_bot_pres(i) = max(tmp_pres_space.pres(~isnan(tmp_pres_space.salt(:,i)),i));     
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
    [N2, mid_pres] = gsw_Nsquared(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative, tmp_pres_space.pres, LLCsealdata(tag_no).lat' .* ones(size(tmp_pres_space.salt)));
    for i = 1:length(LLCsealdata(tag_no).cast)
        tmp_pres_space.N2(:,i) = interp1(mid_pres(:,i), N2(:,i), tmp_pres_space.pres(:,i));
    end

    %%% Calculating Dynamic Height Anomaly 
    tmp_pres_space.dyn_height_anom = gsw_geo_strf_dyn_height(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative, tmp_pres_space.pres, 0);

    %%% Calculating Spice
    tmp_pres_space.spice = gsw_spiciness0(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative);
    
    %%% Assigning variables to structure
    LLCsealdata(tag_no).ps = tmp_pres_space;

    clear ts_isopycnal_sep ts_pres u pres_final pres j isopycnal_sep isopycnal_sep_ds isopycnal_sep_y_axis...
            i a b density depths density_final k isopycnal_sep_final pres_ds ts_density tmp_pres_space 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolating to Density Grid %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Creating density grid
density_grid = (26.6:0.01:28.3)';

for tag_no = 1:length(LLCsealdata)
    
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Finding Indices to Build Reference Profiles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ref_settings.inner_window = 2;
ref_settings.outer_window = 12;

for tag_no = 1:length(LLCsealdata)

    LLCsealdata(tag_no).ref_ind = [];
    mean_ind = cell(2,1);
    
    %%% Finding indices to use for background profile calculation
    i = LLCsealdata(tag_no).centerindex;

    mean_ind(1,:) = {(i-ref_settings.outer_window):(i-ref_settings.inner_window)};
    mean_ind(2,:) = {(i+ref_settings.inner_window):(i+ref_settings.outer_window)};

    %%% Saving indices
    LLCsealdata(tag_no).ref_ind = mean_ind;
    
end

clear mean_ind i

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Building Reference Profiles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = 1:length(LLCsealdata)
    
    %%% Creating reference profiles for each of the time series profiles
    LLCsealdata(tag_no).ds.ref_salt = NaN(size(LLCsealdata(tag_no).ds.salt,1),1);
    LLCsealdata(tag_no).ds.ref_temp = NaN(size(LLCsealdata(tag_no).ds.temp,1),1);
    LLCsealdata(tag_no).ds.ref_N2 = NaN(size(LLCsealdata(tag_no).ds.N2,1),1);
    LLCsealdata(tag_no).ds.ref_spice = NaN(size(LLCsealdata(tag_no).ds.spice,1),1);
    LLCsealdata(tag_no).ds.ref_isopycnal_separation = NaN(size(LLCsealdata(tag_no).ds.isopycnal_separation,1),1);

    for i = LLCsealdata(tag_no).centerindex
    
        %%% Extracting the profiles to build the reference profile
        tmp_salt_ds = LLCsealdata(tag_no).ds.salt(:,[LLCsealdata(tag_no).ref_ind{1,1} LLCsealdata(tag_no).ref_ind{2,1}]);
        tmp_temp_ds = LLCsealdata(tag_no).ds.temp(:,[LLCsealdata(tag_no).ref_ind{1,1} LLCsealdata(tag_no).ref_ind{2,1}]);
        tmp_N2_ds = LLCsealdata(tag_no).ds.N2(:,[LLCsealdata(tag_no).ref_ind{1,1} LLCsealdata(tag_no).ref_ind{2,1}]);
        tmp_spice_ds = LLCsealdata(tag_no).ds.spice(:,[LLCsealdata(tag_no).ref_ind{1,1} LLCsealdata(tag_no).ref_ind{2,1}]);
        tmp_isopycnal_separation_ds = LLCsealdata(tag_no).ds.isopycnal_separation(:,[LLCsealdata(tag_no).ref_ind{1,1} LLCsealdata(tag_no).ref_ind{2,1}]);

        %%% Calculating a climatological value for each density level if at least
        %%% 75% of the potential reference profiles have data at said
        %%% level
        for j = 1:size(tmp_salt_ds, 1)
            tmp_salt_ds_level = tmp_salt_ds(j,~isnan(tmp_salt_ds(j,:)));
            tmp_temp_ds_level = tmp_temp_ds(j,~isnan(tmp_temp_ds(j,:)));
            tmp_spice_ds_level = tmp_spice_ds(j,~isnan(tmp_spice_ds(j,:)));
            if length(tmp_salt_ds_level) > 0.75*size(tmp_salt_ds,2)
                LLCsealdata(tag_no).ds.ref_salt(j,1) = median(tmp_salt_ds_level);
                LLCsealdata(tag_no).ds.ref_temp(j,1) = median(tmp_temp_ds_level);
                LLCsealdata(tag_no).ds.ref_spice(j,1) = median(tmp_spice_ds_level); 
            end
            
            tmp_N2_ds_level = tmp_N2_ds(j,~isnan(tmp_N2_ds(j,:)));
            if length(tmp_N2_ds_level) > 0.75*size(tmp_N2_ds,2)
                LLCsealdata(tag_no).ds.ref_N2(j,1) = median(tmp_N2_ds_level);  
            end
            
            tmp_isopycnal_separation_ds_level = tmp_isopycnal_separation_ds(j,~isnan(tmp_isopycnal_separation_ds(j,:)));
            if length(tmp_isopycnal_separation_ds_level) > 0.75*size(tmp_isopycnal_separation_ds,2)
                LLCsealdata(tag_no).ds.ref_isopycnal_separation(j,1) = median(tmp_isopycnal_separation_ds_level);  
            end
        end
        
        
    end

end

clear tmp_pres_ds tmp_salt_ds tmp_temp_ds tmp_N2_ds tmp_spice_ds tmp_isopycnal_separation_ds tmp_pres_ds_level tmp_salt_ds_level ...
    tmp_temp_ds_level tmp_spice_ds_level tmp_N2_ds_level tmp_isopycnal_separation_ds_level i j

for tag_no = 1:length(LLCsealdata)

    %%% Creating reference profiles for each of the time series profiles
    LLCsealdata(tag_no).ps.ref_dyn_height_anom = NaN(size(LLCsealdata(tag_no).ps.dyn_height_anom,1),1);

    for i = LLCsealdata(tag_no).centerindex
    
        %%% Extracting the profiles to build the reference profile
        tmp_dyn_height_anom_ps = LLCsealdata(tag_no).ps.dyn_height_anom(:,[LLCsealdata(tag_no).ref_ind{1,1} LLCsealdata(tag_no).ref_ind{2,1}]);

        %%% Calculating a climatological value for each density level if at least
        %%% 75% of the potential reference profiles have data at said
        %%% level
        for j = 1:size(tmp_dyn_height_anom_ps, 1)
            tmp_dyn_height_anom_ps_level = tmp_dyn_height_anom_ps(j,~isnan(tmp_dyn_height_anom_ps(j,:)));
            if length(tmp_dyn_height_anom_ps_level) > 0.75*size(tmp_dyn_height_anom_ps,2)
                LLCsealdata(tag_no).ps.ref_dyn_height_anom(j,1) = median(tmp_dyn_height_anom_ps_level);  
            end
        end
        
    end

end

clear tmp_dyn_height_anom_ps tmp_dyn_height_anom_ps_level

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Anomalies %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = 1:length(LLCsealdata)

    %%% Creating anomaly profiles for each of the time series profiles
    LLCsealdata(tag_no).ds.salt_anom = NaN(size(LLCsealdata(tag_no).ds.salt,1),1);
    LLCsealdata(tag_no).ds.temp_anom = NaN(size(LLCsealdata(tag_no).ds.temp,1),1);
    LLCsealdata(tag_no).ds.N2_anom = NaN(size(LLCsealdata(tag_no).ds.N2,1),1);
    LLCsealdata(tag_no).ds.spice_anom = NaN(size(LLCsealdata(tag_no).ds.spice,1),1);
    LLCsealdata(tag_no).ds.isopycnal_separation_anom = NaN(size(LLCsealdata(tag_no).ds.isopycnal_separation,1),1);

    for i = LLCsealdata(tag_no).centerindex

        %%% Calculating Anomalies
        LLCsealdata(tag_no).ds.salt_anom(:,1) = LLCsealdata(tag_no).ds.salt(:,i) - LLCsealdata(tag_no).ds.ref_salt(:,1);
        LLCsealdata(tag_no).ds.temp_anom(:,1) = LLCsealdata(tag_no).ds.temp(:,i) - LLCsealdata(tag_no).ds.ref_temp(:,1);
        LLCsealdata(tag_no).ds.N2_anom(:,1) = LLCsealdata(tag_no).ds.N2(:,i) - LLCsealdata(tag_no).ds.ref_N2(:,1);
        LLCsealdata(tag_no).ds.spice_anom(:,1) = LLCsealdata(tag_no).ds.spice(:,i) - LLCsealdata(tag_no).ds.ref_spice(:,1);
        LLCsealdata(tag_no).ds.isopycnal_separation_anom(:,1) = LLCsealdata(tag_no).ds.isopycnal_separation(:,i) - LLCsealdata(tag_no).ds.ref_isopycnal_separation(:,1);
    end
end

for tag_no = 1:length(LLCsealdata)
    LLCsealdata(tag_no).ps.dyn_height_anom_anom = NaN(size(LLCsealdata(tag_no).ps.dyn_height_anom,1),1);

    for i = LLCsealdata(tag_no).centerindex
        LLCsealdata(tag_no).ps.dyn_height_anom_anom(:,1) = LLCsealdata(tag_no).ps.dyn_height_anom(:,i) - LLCsealdata(tag_no).ps.ref_dyn_height_anom(:,1);
    end
end

clear i

%%
isopycnals = 0.01;

for tag_no = 1:length(LLCsealdata)
        i = LLCsealdata(tag_no).centerindex;
        
        fig = figure('Position', [0 0 1000 850]);
        sgtitle('Manually Created "Seal" Track', 'FontSize', 18, 'FontWeight', 'bold')
        
        %%% Temperature Subplot
        ax1 = subplot(4,4,1:3);
        hold on
        pp = pcolor(LLCsealdata(tag_no).cast, LLCsealdata(tag_no).ps.pres(:,1), LLCsealdata(tag_no).temp);
        set(pp, 'EdgeColor', 'none');
        [C,h] = contour(ax1, LLCsealdata(tag_no).cast, depth_grid, LLCsealdata(tag_no).ps.sigma0, round(min(min(LLCsealdata(tag_no).ps.sigma0)):isopycnals:max(max(LLCsealdata(tag_no).ps.sigma0)), 2), 'k');
        clabel(C,h,'LabelSpacing',500);
        xline(LLCsealdata(tag_no).cast(i), 'r', 'LineWidth', 1.5)
        hold off
        cmap = cmocean('thermal'); colormap(ax1, cmap); colorbar;
        clim([min(min(LLCsealdata(tag_no).temp)) max(max(LLCsealdata(tag_no).temp))])
        set(gca, 'YDir','reverse'); set(gca, 'Layer','top');
        ylabel('Pressure (dbar)', 'FontSize', 12);
        ylim([0 800])
        title('Temperature', 'FontSize', 12);
        
        %%% Salinity Subplot
        ax2 = subplot(4,4,5:7);
        hold on
        pp = pcolor(LLCsealdata(tag_no).cast, LLCsealdata(tag_no).ps.pres(:,1), LLCsealdata(tag_no).salt);
        set(pp, 'EdgeColor', 'none');
        [C,h] = contour(ax2, LLCsealdata(tag_no).cast, depth_grid, LLCsealdata(tag_no).ps.sigma0, round(min(min(LLCsealdata(tag_no).ps.sigma0)):isopycnals:max(max(LLCsealdata(tag_no).ps.sigma0)), 2), 'k');
        clabel(C,h,'LabelSpacing',500);
        xline(LLCsealdata(tag_no).cast(i), 'r', 'LineWidth', 1.5)
        hold off
        cmap = cmocean('haline'); colormap(ax2, cmap); colorbar;
        clim([min(min(LLCsealdata(tag_no).salt)) max(max(LLCsealdata(tag_no).salt))])
        set(gca, 'YDir','reverse'); set(gca, 'Layer','top');
        ylabel('Pressure (dbar)', 'FontSize', 12);
        ylim([0 800])
        title('Salinity', 'FontSize', 12);
        
        %%% Map Subplot
        ax3 = subplot(4,4,4);
        axesm('lambertstd','MapParallels',[-75 -15],'MapLatLimit',[min(LLCsealdata(tag_no).lat)-0.1 max(LLCsealdata(tag_no).lat)+0.1],'MapLonLimit',[min(LLCsealdata(tag_no).lon)-0.2 max(LLCsealdata(tag_no).lon)+0.2], 'MLineLocation', 2, 'PLineLocation', 1, 'FontSize',13);
        axis off; framem on; gridm on; mlabel on; plabel on;
        hold on
        colormap(ax3, cmocean('balance')); colorbar; clim([-0.3 0.3]);
        if LLCsealdata(tag_no).region == "WAP"
            pcolorm(LLC_5.lats(:,:), LLC_5.lons(:,:), LLC_5.mat.Ro(:,:,32));
        elseif LLCsealdata(tag_no).region == "East"
            pcolorm(LLC_1.lats(:,:), LLC_1.lons(:,:), LLC_1.mat.Ro(:,:,32));
        end
        geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor','white');
        plotm(LLCsealdata(tag_no).lat, LLCsealdata(tag_no).lon, 'linestyle', '-', 'Color','g', 'Marker','o', 'MarkerFaceColor','w', 'MarkerSize', 1);

        %%% Map Subplot
        subplot(4,4,8)
        axesm('lambertstd','MapParallels',[-75 -15],'MapLatLimit',[min(LLCsealdata(tag_no).lat)-10 max(LLCsealdata(tag_no).lat)+10],'MapLonLimit',[min(LLCsealdata(tag_no).lon)-20 max(LLCsealdata(tag_no).lon)+20], 'MLineLocation', 20, 'PLineLocation', 5, 'FontSize',13);
        axis off; framem on; gridm on; mlabel on; plabel on;
        hold on
        geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor','white');
        plotm(LLCsealdata(tag_no).lat, LLCsealdata(tag_no).lon, 'linestyle', '-', 'Color','g', 'Marker','o', 'MarkerFaceColor','w', 'MarkerSize', 1);

        %%% Spice Profile
        subplot(4,4,9);
        hold on
        plot(LLCsealdata(tag_no).ds.spice(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        plot(LLCsealdata(tag_no).ds.ref_spice(:,1), LLCsealdata(tag_no).ds.pres(:,1), 'k', 'DisplayName', 'Reference', 'LineWidth', 1.5)
        set(gca, 'YDir', 'reverse');
        xlabel('Spice');
        hold off
        legend('Location', 'best')
        ylim([0 800]);
        
        %%% Isopycnal Separation Profile
        subplot(4,4,10);
        hold on
        plot(LLCsealdata(tag_no).ds.isopycnal_separation(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        plot(LLCsealdata(tag_no).ds.ref_isopycnal_separation(:,1), LLCsealdata(tag_no).ds.pres(:,1), 'k','DisplayName', 'Reference', 'LineWidth', 1.5)
        set(gca, 'YDir', 'reverse');
        xlabel('Isopycnal Separation');
        hold off
        legend('Location', 'best')
        ylim([0 800]);
        
        %%% N2 Profile
        subplot(4,4,11);
        hold on
        plot(LLCsealdata(tag_no).ds.N2(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        plot(LLCsealdata(tag_no).ds.ref_N2(:,1), LLCsealdata(tag_no).ds.pres(:,1), 'k','DisplayName', 'Reference', 'LineWidth', 1.5)
        set(gca, 'YDir', 'reverse');
        xlabel('N^2');
        hold off
        legend('Location', 'best')
        ylim([0 800]);

        %%% DHA Profile
        subplot(4,4,12)
        hold on
        plot(LLCsealdata(tag_no).ps.dyn_height_anom(:,i), LLCsealdata(tag_no).ps.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        plot(LLCsealdata(tag_no).ps.ref_dyn_height_anom(:,1), LLCsealdata(tag_no).ps.pres(:,1), 'k','DisplayName', 'Reference', 'LineWidth', 1.5)
        set(gca, 'YDir', 'reverse');
        xlabel('Dynamic Height Anomaly');
        hold off
        legend('Location', 'best')
        ylim([0 800]);
        
        %%% Spice Anomaly Profile
        subplot(4,4,13);
        hold on
        plot(LLCsealdata(tag_no).ds.spice_anom(:,1), LLCsealdata(tag_no).ds.pres(:,1), 'b', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        xline(0, '--k', 'LineWidth', 1);
        set(gca, 'YDir', 'reverse');
        xlabel('Spice Anomaly');
        hold off
        ylim([0 800]);
        
        %%% Isopycnal Separation Anomaly Profile
        subplot(4,4,14);
        hold on
        plot(LLCsealdata(tag_no).ds.isopycnal_separation_anom(:,1), LLCsealdata(tag_no).ds.pres(:,1), 'b', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        xline(0, '--k', 'LineWidth', 1);
        set(gca, 'YDir', 'reverse');
        xlabel('Isopycnal Separation Anomaly');
        hold off
        ylim([0 800]);
        
        %%% N2 Anomaly Profile
        subplot(4,4,15);
        hold on
        plot(LLCsealdata(tag_no).ds.N2_anom(:,1), LLCsealdata(tag_no).ds.pres(:,1), 'b', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        xline(0, '--k', 'LineWidth', 1);
        set(gca, 'YDir', 'reverse');
        xlabel('N^2 Anomaly');
        hold off
        ylim([0 800]);
             
        %%% Dynamic Height Anomaly Profile
        subplot(4,4,16);
        hold on
        plot(LLCsealdata(tag_no).ps.dyn_height_anom_anom(:,1), LLCsealdata(tag_no).ps.pres(:,1),'b', 'DisplayName', 'Adjusted Profile','LineWidth',1.5)
        xline(0, '--k', 'LineWidth', 1);
        set(gca, 'YDir', 'reverse');
        xlabel('Dynamic Height Anomaly Anomaly');
        hold off
        ylim([0 800]);

        saveas(fig, '/Users/jenkosty/Documents/Research/SCV_Project/Figures/25Oct2022/' + string(LLCsealdata(tag_no).tag), 'png')
        
end

        clear ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 C h cmap fig h IB isopycnals p1 p2 p3 p4 p5 pp

