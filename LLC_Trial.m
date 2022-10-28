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

%%
%%% Calculating the Okubo Weiss Parameter
ind = [27 32];
LLC_1.OW = Okubo_Weiss(LLC_1.mat.dxc(:,:), LLC_1.mat.dyc(:,:), LLC_1.mat.u(:,:,ind), LLC_1.mat.v(:,:,ind));
LLC_2.OW = Okubo_Weiss(LLC_2.mat.dxc(:,:), LLC_2.mat.dyc(:,:), LLC_2.mat.u(:,:,ind), LLC_2.mat.v(:,:,ind));
LLC_4.OW = Okubo_Weiss(LLC_4.mat.dxc(:,:), LLC_4.mat.dyc(:,:), LLC_4.mat.u(:,:,ind), LLC_4.mat.v(:,:,ind));
LLC_5.OW = Okubo_Weiss(LLC_5.mat.dxc(:,:), LLC_5.mat.dyc(:,:), LLC_5.mat.u(:,:,ind), LLC_5.mat.v(:,:,ind));

%%
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

    %%% Saving interpolated LLC data to a structure
    LLCsealdata(tag_no) = LLCseal;

    clear sector i LLCseal LLCseal_salt LLCseal_temp LLCseal_vort LLCseal_OW LLC_lats LLC_lons
end

%%
%%% Using Okubo Weiss to flag profiles as eddies
for tag_no = test_prof
    for i = 1:length(LLCsealdata(tag_no).cast)
        if abs(LLCsealdata(tag_no).OW(1,i)) > 0.4e-8 && abs(LLCsealdata(tag_no).OW(2,i)) > 0.4e-8
            LLCsealdata(tag_no).scv(i) = 1;
        else
            LLCsealdata(tag_no).scv(i) = 0;
        end
    end
end



%%

for tag_no = test_prof
    
    %%% Calculating Bathymetry
%     LLC_data(tag_no).bathymetry = interp2(RTOPO.lon, RTOPO.lat', RTOPO.bedrock_topography, qc_ts(tag_no).lon, qc_ts(tag_no).lat);

    %%% Creating pressure space structure
    tmp_pres_space.pres = depth_grid .* ones(size(LLC_data(tag_no).salt));
    
    %%% Assigning temperature and salinity data to new pressure space
    %%% structure
    tmp_pres_space.salt = LLC_data(tag_no).salt;
    tmp_pres_space.temp = LLC_data(tag_no).temp;
    tmp_pres_space.vort = LLC_data(tag_no).vort;
    
    %%% Calculating bottom pressure
    for i = 1:length(qc_ts(tag_no).cast)
        tmp_pres_space.prof_bot_pres(i) = max(tmp_pres_space.pres(~isnan(tmp_pres_space.salt(:,i)),i));     
    end

    %%% Calculating absolute salinity and conservative temperature
    tmp_pres_space.salt_absolute = gsw_SA_from_SP(tmp_pres_space.salt, depth_grid, LLC_data(tag_no).lon, LLC_data(tag_no).lat);
    tmp_pres_space.temp_conservative = gsw_CT_from_t(tmp_pres_space.salt_absolute, tmp_pres_space.temp, tmp_pres_space.pres);

    %%% Calculating density
    tmp_pres_space.density = gsw_rho_CT_exact(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative, ones(size(tmp_pres_space.salt)).*400);

    %%% Calculating potential density anomaly
    tmp_pres_space.sigma0 = gsw_sigma0(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative);
    
    %%% Calculating neutral density
    %[tmp_pres_space.gamma_n,~,~] = eos80_legacy_gamma_n(tmp_pres_space.salt, tmp_pres_space.temp, tmp_pres_space.pres, qc_ts(tag_no).lon, qc_ts(tag_no).lat);

    %%% Calculating N^2
    [N2, mid_pres] = gsw_Nsquared(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative, tmp_pres_space.pres, LLC_data(tag_no).lat .* ones(size(tmp_pres_space.salt)));
    for i = 1:length(qc_ts(tag_no).cast)
        tmp_pres_space.N2(:,i) = interp1(mid_pres(:,i), N2(:,i), tmp_pres_space.pres(:,i));
    end

    %%% Calculating Spice
    tmp_pres_space.spice = gsw_spiciness0(tmp_pres_space.salt_absolute, tmp_pres_space.temp_conservative);
    
    %%% Assigning variables to structure
    LLC_data(tag_no).ps = tmp_pres_space;

    clear ts_isopycnal_sep ts_pres u pres_final pres j isopycnal_sep isopycnal_sep_ds isopycnal_sep_y_axis...
            i a b density depths density_final k isopycnal_sep_final pres_ds ts_density tmp_pres_space 

end

clear N2 mid_pres

%%
tag_no = test_prof;
i = 154;

isopycnals = 0.02;

figure('Renderer', 'painters', 'Position', [0 0 1000 850])

ax1 = subplot(311);
hold on
[~,IB] = unique(datenum(LLC_data(tag_no).time));
pp = pcolor(unique(datenum(LLC_data(tag_no).time)),LLC_data(tag_no).ps.pres(:,1), LLC_data(tag_no).salt(:,IB));
set(pp, 'EdgeColor', 'none');
[C,h] = contour(ax1, unique(datenum(LLC_data(tag_no).time)), depth_grid, LLC_data(tag_no).ps.sigma0(:,IB), round(min(min(LLC_data(tag_no).ps.sigma0)):isopycnals:max(max(LLC_data(tag_no).ps.sigma0)), 2), 'k');
clabel(C,h,'LabelSpacing',500);
xline(datenum(LLC_data(tag_no).time(i,:)), 'r', 'LineWidth', 1.5)
hold off
cmap = cmocean('haline'); colormap(ax1, cmap); colorbar; caxis([min(min(LLC_data(tag_no).salt)) max(max(LLC_data(tag_no).salt))])
set(gca, 'YDir','reverse'); set(gca, 'Layer','top');
xticks(linspace(datenum(LLC_data(tag_no).time(1,:)), datenum(LLC_data(tag_no).time(end,:)), (datenum(LLC_data(tag_no).time(end,:)) - datenum(LLC_data(tag_no).time(1,:))) / 5))
datetick('x', 'mm/dd/yy', 'keepticks');
ylim([0 500])
ylabel('Pressure (dbar)', 'FontSize', 12);
title('Salinity', 'FontSize', 12);

ax2 = subplot(312);
hold on
[~,IB] = unique(datenum(LLC_data(tag_no).time));
pp = pcolor(unique(datenum(LLC_data(tag_no).time)),LLC_data(tag_no).ps.pres(:,1), LLC_data(tag_no).temp(:,IB));
set(pp, 'EdgeColor', 'none');
[C,h] = contour(ax2, unique(datenum(LLC_data(tag_no).time)), depth_grid, LLC_data(tag_no).ps.sigma0(:,IB), round(min(min(LLC_data(tag_no).ps.sigma0)):isopycnals:max(max(LLC_data(tag_no).ps.sigma0)), 2), 'k');
clabel(C,h,'LabelSpacing',500);
xline(datenum(LLC_data(tag_no).time(i,:)), 'r', 'LineWidth', 1.5)
hold off
cmap = cmocean('thermal'); colormap(ax2, cmap); colorbar; caxis([min(min(LLC_data(tag_no).temp)) max(max(LLC_data(tag_no).temp))])
set(gca, 'YDir','reverse'); set(gca, 'Layer','top');
xticks(linspace(datenum(LLC_data(tag_no).time(1,:)), datenum(LLC_data(tag_no).time(end,:)), (datenum(LLC_data(tag_no).time(end,:)) - datenum(LLC_data(tag_no).time(1,:))) / 5))
datetick('x', 'mm/dd/yy', 'keepticks');
ylim([0 500])
ylabel('Pressure (dbar)', 'FontSize', 12);
title('Temperature', 'FontSize', 12);

ax3 = subplot(313);
hold on
[~,IB] = unique(datenum(LLC_data(tag_no).time));
pp = pcolor(unique(datenum(LLC_data(tag_no).time)),LLC_data(tag_no).ps.pres(:,1), LLC_data(tag_no).vort(:,IB));
set(pp, 'EdgeColor', 'none');
[C,h] = contour(ax3, unique(datenum(LLC_data(tag_no).time)), depth_grid, LLC_data(tag_no).ps.sigma0(:,IB), round(min(min(LLC_data(tag_no).ps.sigma0)):isopycnals:max(max(LLC_data(tag_no).ps.sigma0)), 2), 'k');
clabel(C,h,'LabelSpacing',500);
xline(datenum(LLC_data(tag_no).time(i,:)), 'r', 'LineWidth', 1.5)
hold off
cmap = cmocean('balance'); colormap(ax3, cmap); colorbar; caxis([-max(abs(LLC_data(tag_no).vort(:))) max(abs(LLC_data(tag_no).vort(:)))])
set(gca, 'YDir','reverse'); set(gca, 'Layer','top');
xticks(linspace(datenum(LLC_data(tag_no).time(1,:)), datenum(LLC_data(tag_no).time(end,:)), (datenum(LLC_data(tag_no).time(end,:)) - datenum(LLC_data(tag_no).time(1,:))) / 5))
datetick('x', 'mm/dd/yy', 'keepticks');
ylim([0 500])
ylabel('Pressure (dbar)', 'FontSize', 12);
title('Rossby Number', 'FontSize', 12);




%%

%%%
sea_ice = shaperead('/Users/jenkosty/Documents/Research/Sea_Ice_Extent/median_extent_S_12_1981-2010_polyline_v3.0/median_extent_S_12_1981-2010_polyline_v3.0');
proj = shapeinfo('/Users/jenkosty/Documents/Research/Sea_Ice_Extent/median_extent_S_12_1981-2010_polyline_v3.0/median_extent_S_12_1981-2010_polyline_v3.0.shp').CoordinateReferenceSystem;
[sea_ice_lat, sea_ice_lon] = projinv(proj, [sea_ice.X], [sea_ice.Y]);

tag_no = 92;
%%% Plotting vorticity data
load coastlines
figure('Position', [500 100 1000 850])
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -57],'MLineLocation', 30, 'PLineLocation', 10, 'FontSize',10);
axis off; framem on; gridm on; mlabel on; plabel on;
colormap(cmocean('balance')); colorbar; clim([-1e-8 1e-8]);
hold on
pcolorm(LLC_1.lats(:,:), LLC_1.lons(:,:), LLC_1.OW(:,:, 1));
plotm(sea_ice_lat, sea_ice_lon, '-g')
geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor','white');
%plotm(qc_ts(tag_no).lat, qc_ts(tag_no).lon, 'g-','Marker', '.','MarkerSize', 10, 'LineWidth', 2)
hold off
%title('LLC Vorticity Data: ' + string(-LLC.depth(depth_ind)) + ' m')



