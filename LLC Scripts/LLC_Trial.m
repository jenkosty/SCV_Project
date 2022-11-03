% clear; close all; clc

%%% Loading MEOP seal data
load("qc_ts.mat")

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
LLC_1.salt = double(LLC_1.mat.s);
LLC_1.temp = double(LLC_1.mat.t);
LLC_1.vort = double(LLC_1.mat.Ro);
LLC_1.depth = LLC_1.mat.rc;

LLC_2.edge_lats = [double(LLC_2.mat.yc(1,:)), double(LLC_2.mat.yc(:,end))', flip(double(LLC_2.mat.yc(end,:))), flip(double(LLC_2.mat.yc(:,1)))'];
LLC_2.edge_lons = [double(LLC_2.mat.xc(1,:)), double(LLC_2.mat.xc(:,end))', flip(double(LLC_2.mat.xc(end,:))), flip(double(LLC_2.mat.xc(:,1)))'];
LLC_2.polygon = geopolyshape(LLC_2.edge_lats, LLC_2.edge_lons);
LLC_2.lats = double(LLC_2.mat.yc);
LLC_2.lons = double(LLC_2.mat.xc);
LLC_2.salt = double(LLC_2.mat.s);
LLC_2.temp = double(LLC_2.mat.t);
LLC_2.vort = double(LLC_2.mat.Ro);
LLC_2.depth = LLC_2.mat.rc;

LLC_4.edge_lats = [double(LLC_4.mat.yc(1,:)), double(LLC_4.mat.yc(:,end))', flip(double(LLC_4.mat.yc(end,:))), flip(double(LLC_4.mat.yc(:,1)))'];
LLC_4.edge_lons = [double(LLC_4.mat.xc(1,:)), double(LLC_4.mat.xc(:,end))', flip(double(LLC_4.mat.xc(end,:))), flip(double(LLC_4.mat.xc(:,1)))'];
LLC_4.polygon = geopolyshape(LLC_4.edge_lats, LLC_4.edge_lons);
LLC_4.lats = double(LLC_4.mat.yc);
LLC_4.lons = double(LLC_4.mat.xc);
LLC_4.salt = double(LLC_4.mat.s);
LLC_4.temp = double(LLC_4.mat.t);
LLC_4.vort = double(LLC_4.mat.Ro);
LLC_4.depth = LLC_4.mat.rc;

LLC_5.edge_lats = [double(LLC_5.mat.yc(1,:)), double(LLC_5.mat.yc(:,end))', flip(double(LLC_5.mat.yc(end,:))), flip(double(LLC_5.mat.yc(:,1)))'];
LLC_5.edge_lons = [double(LLC_5.mat.xc(1,:)), double(LLC_5.mat.xc(:,end))', flip(double(LLC_5.mat.xc(end,:))), flip(double(LLC_5.mat.xc(:,1)))'];
LLC_5.polygon = geopolyshape(LLC_5.edge_lats, LLC_5.edge_lons);
LLC_5.lats = double(LLC_5.mat.yc);
LLC_5.lons = double(LLC_5.mat.xc);
LLC_5.salt = double(LLC_5.mat.s);
LLC_5.temp = double(LLC_5.mat.t);
LLC_5.vort = double(LLC_5.mat.Ro);
LLC_5.depth = LLC_5.mat.rc;

%%% Calculating the Okubo Weiss Parameter
ind = [27 32];
LLC_1.OW = Okubo_Weiss(LLC_1.mat.dxc(:,:), LLC_1.mat.dyc(:,:), LLC_1.mat.u(:,:,ind), LLC_1.mat.v(:,:,ind));
LLC_2.OW = Okubo_Weiss(LLC_2.mat.dxc(:,:), LLC_2.mat.dyc(:,:), LLC_2.mat.u(:,:,ind), LLC_2.mat.v(:,:,ind));
LLC_4.OW = Okubo_Weiss(LLC_4.mat.dxc(:,:), LLC_4.mat.dyc(:,:), LLC_4.mat.u(:,:,ind), LLC_4.mat.v(:,:,ind));
LLC_5.OW = Okubo_Weiss(LLC_5.mat.dxc(:,:), LLC_5.mat.dyc(:,:), LLC_5.mat.u(:,:,ind), LLC_5.mat.v(:,:,ind));

%%
test_prof = 1:length(qc_ts);

depth_grid = 1:800;
depth_grid = depth_grid';

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
    LLCseal_salt = NaN(length(depth_grid), length(qc_ts(tag_no).cast));
    LLCseal_temp = NaN(length(depth_grid), length(qc_ts(tag_no).cast));
    LLCseal_vort = NaN(length(depth_grid), length(qc_ts(tag_no).cast));
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
            lats_unfmt(j) = double(LLC.lats(rows(j), cols(j)));
            lons_unfmt(j) = double(LLC.lons(rows(j), cols(j)));
            salt(:,j) = squeeze(LLC.salt(rows(j), cols(j), :));
            temp(:,j) = squeeze(LLC.temp(rows(j), cols(j), :));
            vort(:,j) = squeeze(LLC.vort(rows(j), cols(j), :));
            okubo_weiss(:,j) = squeeze(LLC.OW(rows(j), cols(j),:));
        end
        temp(salt == 0) = NaN;
        vort(salt == 0) = NaN;
        salt(salt == 0) = NaN;

        lats = double(lats_unfmt .* ones(size(salt)));
        lons = double(lons_unfmt .* ones(size(salt)));
        depths = LLC.depth(1:size(salt, 1), 1 ) .* ones(size(salt));
        
        %%% Interpolating LLC data
        S = scatteredInterpolant(lats(:), lons(:), depths(:), salt(:));
        S_prof = S(qc_ts(tag_no).lat(i).*ones(size(depth_grid)), qc_ts(tag_no).lon(i).*ones(size(depth_grid)), -depth_grid);
        if isempty(S_prof)
            LLCseal_salt(:,i) = NaN(size(depth_grid));
        else
            LLCseal_salt(:,i) = S_prof;
        end

        T = scatteredInterpolant(lats(:), lons(:), depths(:), temp(:));
        T_prof = T(qc_ts(tag_no).lat(i).*ones(size(depth_grid)), qc_ts(tag_no).lon(i).*ones(size(depth_grid)), -depth_grid);
        if isempty(T_prof)
            LLCseal_temp(:,i) = NaN(size(depth_grid));
        else
            LLCseal_temp(:,i) = T_prof;
        end

        V = scatteredInterpolant(lats(:), lons(:), depths(:), vort(:));
        V_prof = V(qc_ts(tag_no).lat(i).*ones(size(depth_grid)), qc_ts(tag_no).lon(i).*ones(size(depth_grid)), -depth_grid);
        if isempty(V_prof)
            LLCseal_vort(:,i) = NaN(size(depth_grid));
        else
            LLCseal_vort(:,i) = V_prof;
        end

        for j = 1:size(LLC.OW, 3)
            OW = scatteredInterpolant(lats_unfmt', lons_unfmt', okubo_weiss(j,:)');
            OW_prof = OW(qc_ts(tag_no).lat(i), qc_ts(tag_no).lon(i));
            if isempty(OW_prof)
                LLCseal_OW(j,i) = NaN;
            else
                LLCseal_OW(j,i) = OW_prof;
            end
        end

        clear lats_unfmt lons_unfmt salt vort temp okubo_weiss cols rows nan_ind lats lons depths S T V OW j close_ind S_prof T_prof V_prof 
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

    %%% Calculating Dynamic Height Anomaly 
    tmp_pres_space.dyn_height_anom = gsw_geo_strf_dyn_height(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative, tmp_pres_space.pres, 0);

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

for tag_no = test_prof
    
    isopycnal_separation = NaN(length(density_grid), length(LLCsealdata(tag_no).cast));
    
    for i = 1:length(LLCsealdata(tag_no).cast)
        for j = 1:length(density_grid)
            isopycnal_separation(j,i) = LLCsealdata(tag_no).ds.pres(min(j+10, length(density_grid)), i) - LLCsealdata(tag_no).ds.pres(max(j-10, 1), i);
        end
    end
    
    LLCsealdata(tag_no).ds.isopycnal_separation = isopycnal_separation;
    
    clear isopycnal_separation i j 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Finding Indices to Build Reference Profiles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ref_settings.inner_window = 2;
ref_settings.outer_window = 12;
 
for tag_no = test_prof
    
    mean_ind = cell(2,length(LLCsealdata(tag_no).cast));
    
    %%% Finding indices to use for background profile calculation
    for i = 1:length(LLCsealdata(tag_no).cast)
        
        mean_ind(1,i) = {(i-ref_settings.outer_window):(i-ref_settings.inner_window)};
        mean_ind{1,i}(mean_ind{1,i} < 1) = []; %%% Making sure indices remain within ts boundaries
        mean_ind(2,i) = {(i+ref_settings.inner_window):(i+ref_settings.outer_window)};
        mean_ind{2,i}(mean_ind{2,i} > length(LLCsealdata(tag_no).cast)) = []; %%% Making sure indices remain within ts boundaries
        
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
    LLCsealdata(tag_no).ref_ind = mean_ind;
    
end

clear mean_ind i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Building Reference Profiles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = test_prof
    
    %%% Creating reference profiles for each of the time series profiles
    LLCsealdata(tag_no).ds.ref_salt = NaN(size(LLCsealdata(tag_no).ds.salt));
    LLCsealdata(tag_no).ds.ref_temp = NaN(size(LLCsealdata(tag_no).ds.temp));
    LLCsealdata(tag_no).ds.ref_N2 = NaN(size(LLCsealdata(tag_no).ds.N2));
    LLCsealdata(tag_no).ds.ref_spice = NaN(size(LLCsealdata(tag_no).ds.spice));
    LLCsealdata(tag_no).ds.ref_isopycnal_separation = NaN(size(LLCsealdata(tag_no).ds.isopycnal_separation));

    for i = 1:length(LLCsealdata(tag_no).cast)
    
        %%% Extracting the profiles to build the reference profile
        tmp_salt_ds = LLCsealdata(tag_no).ds.salt(:,[LLCsealdata(tag_no).ref_ind{1,i} LLCsealdata(tag_no).ref_ind{2,i}]);
        tmp_temp_ds = LLCsealdata(tag_no).ds.temp(:,[LLCsealdata(tag_no).ref_ind{1,i} LLCsealdata(tag_no).ref_ind{2,i}]);
        tmp_N2_ds = LLCsealdata(tag_no).ds.N2(:,[LLCsealdata(tag_no).ref_ind{1,i} LLCsealdata(tag_no).ref_ind{2,i}]);
        tmp_spice_ds = LLCsealdata(tag_no).ds.spice(:,[LLCsealdata(tag_no).ref_ind{1,i} LLCsealdata(tag_no).ref_ind{2,i}]);
        tmp_isopycnal_separation_ds = LLCsealdata(tag_no).ds.isopycnal_separation(:,[LLCsealdata(tag_no).ref_ind{1,i} LLCsealdata(tag_no).ref_ind{2,i}]);

        %%% Calculating a climatological value for each density level if at least
        %%% 75% of the potential reference profiles have data at said
        %%% level
        for j = 1:size(tmp_salt_ds, 1)
            tmp_salt_ds_level = tmp_salt_ds(j,~isnan(tmp_salt_ds(j,:)));
            tmp_temp_ds_level = tmp_temp_ds(j,~isnan(tmp_temp_ds(j,:)));
            tmp_spice_ds_level = tmp_spice_ds(j,~isnan(tmp_spice_ds(j,:)));
            if length(tmp_salt_ds_level) > 0.75*size(tmp_salt_ds,2)
                LLCsealdata(tag_no).ds.ref_salt(j,i) = median(tmp_salt_ds_level);
                LLCsealdata(tag_no).ds.ref_temp(j,i) = median(tmp_temp_ds_level);
                LLCsealdata(tag_no).ds.ref_spice(j,i) = median(tmp_spice_ds_level); 
            end
            
            tmp_N2_ds_level = tmp_N2_ds(j,~isnan(tmp_N2_ds(j,:)));
            if length(tmp_N2_ds_level) > 0.75*size(tmp_N2_ds,2)
                LLCsealdata(tag_no).ds.ref_N2(j,i) = median(tmp_N2_ds_level);  
            end
            
            tmp_isopycnal_separation_ds_level = tmp_isopycnal_separation_ds(j,~isnan(tmp_isopycnal_separation_ds(j,:)));
            if length(tmp_isopycnal_separation_ds_level) > 0.75*size(tmp_isopycnal_separation_ds,2)
                LLCsealdata(tag_no).ds.ref_isopycnal_separation(j,i) = median(tmp_isopycnal_separation_ds_level);  
            end
        end
        
        
    end

end

clear tmp_pres_ds tmp_salt_ds tmp_temp_ds tmp_N2_ds tmp_spice_ds tmp_isopycnal_separation_ds tmp_pres_ds_level tmp_salt_ds_level ...
    tmp_temp_ds_level tmp_spice_ds_level tmp_N2_ds_level tmp_isopycnal_separation_ds_level i j

for tag_no = test_prof

    %%% Creating reference profiles for each of the time series profiles
    LLCsealdata(tag_no).ps.ref_dyn_height_anom = NaN(size(LLCsealdata(tag_no).ps.dyn_height_anom));

    for i = 1:length(LLCsealdata(tag_no).cast)
    
        %%% Extracting the profiles to build the reference profile
        tmp_dyn_height_anom_ps = LLCsealdata(tag_no).ps.dyn_height_anom(:,[LLCsealdata(tag_no).ref_ind{1,i} LLCsealdata(tag_no).ref_ind{2,i}]);

        %%% Calculating a climatological value for each density level if at least
        %%% 75% of the potential reference profiles have data at said
        %%% level
        for j = 1:size(tmp_dyn_height_anom_ps, 1)
            tmp_dyn_height_anom_ps_level = tmp_dyn_height_anom_ps(j,~isnan(tmp_dyn_height_anom_ps(j,:)));
            if length(tmp_dyn_height_anom_ps_level) > 0.75*size(tmp_dyn_height_anom_ps,2)
                LLCsealdata(tag_no).ps.ref_dyn_height_anom(j,i) = median(tmp_dyn_height_anom_ps_level);  
            end
        end
        
    end

end

clear tmp_dyn_height_anom_ps tmp_dyn_height_anom_ps_level

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Anomalies %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = test_prof

    %%% Creating anomaly profiles for each of the time series profiles
    LLCsealdata(tag_no).ds.salt_anom = NaN(size(LLCsealdata(tag_no).ds.salt));
    LLCsealdata(tag_no).ds.temp_anom = NaN(size(LLCsealdata(tag_no).ds.temp));
    LLCsealdata(tag_no).ds.N2_anom = NaN(size(LLCsealdata(tag_no).ds.N2));
    LLCsealdata(tag_no).ds.spice_anom = NaN(size(LLCsealdata(tag_no).ds.spice));
    LLCsealdata(tag_no).ds.isopycnal_separation_anom = NaN(size(LLCsealdata(tag_no).ds.isopycnal_separation));

    for i = 1:length(LLCsealdata(tag_no).cast)

        %%% Calculating Anomalies
        LLCsealdata(tag_no).ds.salt_anom(:,i) = LLCsealdata(tag_no).ds.salt(:,i) - LLCsealdata(tag_no).ds.ref_salt(:,i);
        LLCsealdata(tag_no).ds.temp_anom(:,i) = LLCsealdata(tag_no).ds.temp(:,i) - LLCsealdata(tag_no).ds.ref_temp(:,i);
        LLCsealdata(tag_no).ds.N2_anom(:,i) = LLCsealdata(tag_no).ds.N2(:,i) - LLCsealdata(tag_no).ds.ref_N2(:,i);
        LLCsealdata(tag_no).ds.spice_anom(:,i) = LLCsealdata(tag_no).ds.spice(:,i) - LLCsealdata(tag_no).ds.ref_spice(:,i);
        LLCsealdata(tag_no).ds.isopycnal_separation_anom(:,i) = LLCsealdata(tag_no).ds.isopycnal_separation(:,i) - LLCsealdata(tag_no).ds.ref_isopycnal_separation(:,i);
    end
end

for tag_no = test_prof
    LLCsealdata(tag_no).ps.dyn_height_anom_anom = NaN(size(LLCsealdata(tag_no).ps.dyn_height_anom));

    for i = 1:length(LLCsealdata(tag_no).cast)
        LLCsealdata(tag_no).ps.dyn_height_anom_anom(:,i) = LLCsealdata(tag_no).ps.dyn_height_anom(:,i) - LLCsealdata(tag_no).ps.ref_dyn_height_anom(:,i);
    end
end

clear i

% 
%%
isopycnals = 0.01;

for tag_no = test_prof

    for i =  [57] %find(LLCsealdata(tag_no).scv)

        ind = min(LLCsealdata(tag_no).ref_ind{1,i}):max(LLCsealdata(tag_no).ref_ind{2,i});

        fig = figure('Position', [0 0 1000 850]);
        sgtitle('LLC Seal Track, Tag = ' + string(LLCsealdata(tag_no).tag) + ', Cast = ' + string(i), 'FontSize', 18, 'FontWeight', 'bold')

        %%% Temperature Subplot
        ax1 = subplot(4,3, [1 2]);
        hold on
        pp = pcolor(LLCsealdata(tag_no).cast(ind), LLCsealdata(tag_no).ps.pres(:,1), LLCsealdata(tag_no).temp(:,ind));
        set(pp, 'EdgeColor', 'none');
        [C,h] = contour(ax1, LLCsealdata(tag_no).cast(ind), depth_grid, LLCsealdata(tag_no).ps.sigma0(:,ind), round(min(min(LLCsealdata(tag_no).ps.sigma0)):isopycnals:max(max(LLCsealdata(tag_no).ps.sigma0)), 2), 'k');
        clabel(C,h,'LabelSpacing',500);
        xline(LLCsealdata(tag_no).cast(i), 'r', 'LineWidth', 1.5)
        hold off
        cmap = cmocean('thermal'); colormap(ax1, cmap); colorbar;
        clim([min(min(LLCsealdata(tag_no).temp(ind))) max(max(LLCsealdata(tag_no).temp(ind)))])
        set(gca, 'YDir','reverse'); set(gca, 'Layer','top');
        ylabel('Pressure (dbar)', 'FontSize', 12);
        ylim([0 500])
        title('Temperature', 'FontSize', 12);

        %%% Salinity Subplot
        ax2 = subplot(4,3,[4 5]);
        hold on
        pp = pcolor(LLCsealdata(tag_no).cast(ind), LLCsealdata(tag_no).ps.pres(:,1), LLCsealdata(tag_no).salt(:,ind));
        set(pp, 'EdgeColor', 'none');
        [C,h] = contour(ax2, LLCsealdata(tag_no).cast(ind), depth_grid, LLCsealdata(tag_no).ps.sigma0(:,ind), round(min(min(LLCsealdata(tag_no).ps.sigma0)):isopycnals:max(max(LLCsealdata(tag_no).ps.sigma0)), 2), 'k');
        clabel(C,h,'LabelSpacing',500);
        xline(LLCsealdata(tag_no).cast(i), 'r', 'LineWidth', 1.5)
        hold off
        cmap = cmocean('haline'); colormap(ax2, cmap); colorbar;
        clim([min(min(LLCsealdata(tag_no).salt(ind))) max(max(LLCsealdata(tag_no).salt(ind())))])
        set(gca, 'YDir','reverse'); set(gca, 'Layer','top');
        ylabel('Pressure (dbar)', 'FontSize', 12);
        ylim([0 500]);
        title('Salinity', 'FontSize', 12);

        %%% Map Subplot
        ax3 = subplot(4,3,[3 6]);
        lat_min = min(LLCsealdata(tag_no).lat(ind)) - 0.1;
        lat_max = max(LLCsealdata(tag_no).lat(ind)) + 0.1;
        lon_min = min(LLCsealdata(tag_no).lon(ind)) - 0.5;
        lon_max = max(LLCsealdata(tag_no).lon(ind)) + 0.5;

        if lon_max - lon_min > 340
            lon_sec = LLCsealdata(tag_no).lon(ind);
            pos_min = min(lon_sec(lon_sec > 0));
            lon_min = -180 - (180-pos_min) - 0.5;
            lon_max = max(lon_sec(lon_sec < 0)) + 0.5;
        end

        %%% Figure settings
        load coastlines

        %%% Plotting data
        hold on

        %%% Grabbing sector of LLC time series
        for j = unique(LLCsealdata(tag_no).sector(ind))
            if j == 1
                LLC = LLC_1;
                pp = pcolor(LLC_1.lats(:,:), LLC_1.lons(:,:), LLC_1.vort(:,:, 36));
                set(pp, 'EdgeColor', 'none');
            elseif j == 2
                LLC = LLC_2;
                pp = pcolor(LLC_2.lats(:,:), LLC_2.lons(:,:), LLC_2.vort(:,:, 36));
                set(pp, 'EdgeColor', 'none');
            elseif j == 4
                LLC = LLC_4;
                pp = pcolor(LLC_4.lats(:,:), LLC_4.lons(:,:), LLC_4.vort(:,:, 36));
                set(pp, 'EdgeColor', 'none');
            elseif j == 5
                LLC = LLC_5;
                pp = pcolor(LLC_5.lats(:,:), LLC_5.lons(:,:), LLC_5.vort(:,:, 36));
                set(pp, 'EdgeColor', 'none');
            end
        end
        colormap(ax3, cmocean('balance')); colorbar; clim([-0.4 0.4]);
        plot(coastlat, coastlon);
        plot(LLCsealdata(tag_no).lat(ind), LLCsealdata(tag_no).lon(ind), '-s', 'Color', [0.5 0.5 0.5],...
            'MarkerSize', 4, 'MarkerFaceColor', [0.5 0.5 0.5], 'LineWidth', 2)
        scatter(LLCsealdata(tag_no).lat(i), LLCsealdata(tag_no).lon(i), 6, 'gs', 'MarkerFaceColor', 'g')
        xlim([lat_min lat_max]);
        xlabel('Latitude');
        ylim([lon_min lon_max]);
        ylabel('Longitude');

        %%% Spice Profile
        subplot(4,4,9);
        hold on
        plot(LLCsealdata(tag_no).ds.spice(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        plot(LLCsealdata(tag_no).ds.ref_spice(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'k', 'DisplayName', 'Reference', 'LineWidth', 1.5)
        set(gca, 'YDir', 'reverse');
        xlabel('Spice');
        hold off
        legend('Location', 'best')
        ylim([0 500]);

        %%% Isopycnal Separation Profile
        subplot(4,4,10);
        hold on
        plot(LLCsealdata(tag_no).ds.isopycnal_separation(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        plot(LLCsealdata(tag_no).ds.ref_isopycnal_separation(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'k','DisplayName', 'Reference', 'LineWidth', 1.5)
        set(gca, 'YDir', 'reverse');
        xlabel('Isopycnal Separation');
        hold off
        legend('Location', 'best')
        ylim([0 500]);

        %%% N2 Profile
        subplot(4,4,11);
        hold on
        plot(LLCsealdata(tag_no).ds.N2(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        plot(LLCsealdata(tag_no).ds.ref_N2(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'k','DisplayName', 'Reference', 'LineWidth', 1.5)
        set(gca, 'YDir', 'reverse');
        xlabel('N^2');
        hold off
        legend('Location', 'best')
        ylim([0 500]);

        %%% DHA Profile
        subplot(4,4,12)
        hold on
        plot(LLCsealdata(tag_no).ps.dyn_height_anom(:,i), LLCsealdata(tag_no).ps.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        plot(LLCsealdata(tag_no).ps.ref_dyn_height_anom(:,i), LLCsealdata(tag_no).ps.pres(:,i), 'k','DisplayName', 'Reference', 'LineWidth', 1.5)
        set(gca, 'YDir', 'reverse');
        xlabel('Dynamic Height Anomaly');
        hold off
        legend('Location', 'best')
        ylim([0 500]);

        %%% Spice Anomaly Profile
        subplot(4,4,13);
        hold on
        plot(LLCsealdata(tag_no).ds.spice_anom(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'b', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        xline(0, '--k', 'LineWidth', 1);
        set(gca, 'YDir', 'reverse');
        xlabel('Spice Anomaly');
        hold off
        ylim([0 500]);

        %%% Isopycnal Separation Anomaly Profile
        subplot(4,4,14);
        hold on
        plot(LLCsealdata(tag_no).ds.isopycnal_separation_anom(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'b', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        xline(0, '--k', 'LineWidth', 1);
        set(gca, 'YDir', 'reverse');
        xlabel('Isopycnal Separation Anomaly');
        hold off
        ylim([0 500]);

        %%% N2 Anomaly Profile
        subplot(4,4,15);
        hold on
        plot(LLCsealdata(tag_no).ds.N2_anom(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'b', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        xline(0, '--k', 'LineWidth', 1);
        set(gca, 'YDir', 'reverse');
        xlabel('N^2 Anomaly');
        hold off
        ylim([0 500]);

        %%% Dynamic Height Anomaly Profile
        subplot(4,4,16);
        hold on
        plot(LLCsealdata(tag_no).ps.dyn_height_anom_anom(:,i), LLCsealdata(tag_no).ps.pres(:,i),'b', 'DisplayName', 'Adjusted Profile','LineWidth',1.5)
        xline(0, '--k', 'LineWidth', 1);
        set(gca, 'YDir', 'reverse');
        xlabel('Dynamic Height Anomaly Anomaly');
        hold off
        ylim([0 500]);

        %saveas(fig, '/Users/jenkosty/Documents/Research/SCV_Project/Figures/25Oct2022/' + string(LLCsealdata(tag_no).tag), 'png')
    end
end

        clear ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 C h cmap fig h IB isopycnals p1 p2 p3 p4 p5 pp
% 
% 
%%
%%%%%%%%%%%%%%%%%%
%%% Map figure %%%
%%%%%%%%%%%%%%%%%%

for tag_no = 221

    lat_min = min(LLCsealdata(tag_no).lat) - 0.1;
    lat_max = max(LLCsealdata(tag_no).lat) + 0.1;
    lon_min = min(LLCsealdata(tag_no).lon) - 0.5;
    lon_max = max(LLCsealdata(tag_no).lon) + 0.5;

    if lon_max - lon_min > 180
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
            pcolorm(LLC_1.lats(:,:), LLC_1.lons(:,:), LLC_1.OW(:,:, 1));
        elseif j == 2
            LLC = LLC_2;
            pcolorm(LLC_2.lats(:,:), LLC_2.lons(:,:), LLC_2.OW(:,:, 1));
        elseif j == 4
            LLC = LLC_4;
            pcolorm(LLC_4.lats(:,:), LLC_4.lons(:,:), LLC_4.OW(:,:, 1));
        elseif j == 5
            LLC = LLC_5;
            pcolorm(LLC_5.lats(:,:), LLC_5.lons(:,:), LLC_5.OW(:,:, 1));
        end
    end
    geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor','white');
    plotm(LLCsealdata(tag_no).lat, LLCsealdata(tag_no).lon, '-s', 'Color', [0.5 0.5 0.5],...
        'MarkerSize', 4, 'MarkerFaceColor', [0.5 0.5 0.5], 'LineWidth', 2)
    scatterm(LLCsealdata(tag_no).lat(LLCsealdata(tag_no).scv == 1), LLCsealdata(tag_no).lon(LLCsealdata(tag_no).scv == 1), 6, 'gs', 'MarkerFaceColor', 'g')

end

clear lat_min lat_max lon_min lon_max pos_min coastlat coastlon


%%

%%% Figure settings
load coastlines
figure('Position', [500 100 1000 850])
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -57],'MLineLocation', 30, 'PLineLocation', 10, 'FontSize',10);
axis off; framem on; gridm on; mlabel on; plabel on;
colormap(cmocean('balance')); colorbar; clim([-1e-8 1e-8]);

%%% Plotting data
hold on
pcolorm(LLC_1.lats(:,:), LLC_1.lons(:,:), LLC_1.OW(:,:, 1));
pcolorm(LLC_2.lats(:,:), LLC_2.lons(:,:), LLC_2.OW(:,:, 1));
pcolorm(LLC_4.lats(:,:), LLC_4.lons(:,:), LLC_4.OW(:,:, 1));
pcolorm(LLC_5.lats(:,:), LLC_5.lons(:,:), LLC_5.OW(:,:, 1));
geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor','white');

scv_lat = NaN;
scv_lon = NaN;
scv_tag = string;
scv_cast = string;
for tag_no = test_prof
    scv_lat_i = LLCsealdata(tag_no).lat(find(LLCsealdata(tag_no).scv));
    scv_lon_i = LLCsealdata(tag_no).lon(find(LLCsealdata(tag_no).scv));
    scv_tag_i = repmat(string(tag_no), 1, length(find(LLCsealdata(tag_no).scv)));
    scv_cast_i = string(find(LLCsealdata(tag_no).scv));
    scv_lat = horzcat(scv_lat, scv_lat_i);
    scv_lon = horzcat(scv_lon, scv_lon_i);
    scv_tag = horzcat(scv_tag, scv_tag_i);
    scv_cast = horzcat(scv_cast, scv_cast_i);
end
scv_tag = scv_tag(~isnan(scv_lon));
scv_cast = scv_cast(~isnan(scv_lon));
scv_lat = scv_lat(~isnan(scv_lat));
scv_lon = scv_lon(~isnan(scv_lon));

scv_label = string;
for i = 1:length(scv_tag)
    scv_label(i) = scv_tag(i) + ', ' + scv_cast(i);
end

scatterm(scv_lat, scv_lon, 6, 'gs', 'filled')
textm(scv_lat, scv_lon, scv_label)




