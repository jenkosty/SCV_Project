clear; close all; clc

%%% Loading LLC data
LLC = loadLLCdata(5, 1:46);

%%% Loading seal data
load("qc_ts.mat")
test_prof = 164;

%%
%%% Generating synthetic seal track
for tag_no = test_prof

    %%% Cutting LLC Data by Latitude
    [~,col_lat] = find(LLC.lat < max(qc_ts(tag_no).lat+0.1) & LLC.lat > min(qc_ts(tag_no).lat-0.1));
    LON = LLC.lon(:,min(col_lat):max(col_lat));
    LAT = LLC.lat(:,min(col_lat):max(col_lat));
    SALT = LLC.salt(:,min(col_lat):max(col_lat),:);
    TEMP = LLC.temp(:,min(col_lat):max(col_lat),:);
    VORT = LLC.vort(:,min(col_lat):max(col_lat),:);

    %%% Cutting LLC Data by Longitude
    [row_lon,~] = find(LON < max(qc_ts(tag_no).lon+0.1) & LON > min(qc_ts(tag_no).lon-0.1));
    LON_final = LON(min(row_lon):max(row_lon),:);
    LAT_final = LAT(min(row_lon):max(row_lon),:);
    SALT_final = SALT(min(row_lon):max(row_lon),:,:);
    TEMP_final = TEMP(min(row_lon):max(row_lon),:,:);
    VORT_final = VORT(min(row_lon):max(row_lon),:,:);
    clear LON LAT SALT TEMP col_lat row_lon

    %%% Formatting LLC Data for Interpolation
    LON_final = LON_final .* ones(size(LON_final,1), size(LON_final,2), size(LLC.depth,1));
    LAT_final = LAT_final .* ones(size(LAT_final,1), size(LAT_final,2), size(LLC.depth,1));
    DEPTH(1,1,:) = -1 .* LLC.depth;
    DEPTH_final = DEPTH .* ones(size(LAT_final,1), size(LAT_final,2), size(LLC.depth,1));
    clear DEPTH

    salt_interp = cell(size(qc_ts(tag_no).cast));
    temp_interp = cell(size(qc_ts(tag_no).cast));
    vort_interp = cell(size(qc_ts(tag_no).cast));
    parfor i = 1:length(qc_ts(tag_no).cast)
        salt_interp{:,i} = squeeze(interp3(LAT_final, LON_final, DEPTH_final, SALT_final,qc_ts(tag_no).lat(i), qc_ts(tag_no).lon(i),qc_ts(tag_no).raw_data(i).pres(~isnan(qc_ts(tag_no).raw_data(i).pres))));
        temp_interp{:,i} = squeeze(interp3(LAT_final, LON_final, DEPTH_final, TEMP_final,qc_ts(tag_no).lat(i), qc_ts(tag_no).lon(i),qc_ts(tag_no).raw_data(i).pres(~isnan(qc_ts(tag_no).raw_data(i).pres))));
        vort_interp{:,i} = squeeze(interp3(LAT_final, LON_final, DEPTH_final, VORT_final,qc_ts(tag_no).lat(i), qc_ts(tag_no).lon(i),qc_ts(tag_no).raw_data(i).pres(~isnan(qc_ts(tag_no).raw_data(i).pres))));
    end

    salt_unfmt = NaN(length(depth_grid), size(qc_ts(tag_no).cast,2));
    temp_unfmt = NaN(length(depth_grid), size(qc_ts(tag_no).cast,2));
    vort_unfmt = NaN(length(depth_grid), size(qc_ts(tag_no).cast,2));
    for i = 1:length(qc_ts(tag_no).cast)
        salt_unfmt(:,i) = interp1(qc_ts(tag_no).raw_data(i).pres(~isnan(qc_ts(tag_no).raw_data(i).pres)), salt_interp{:,i}, depth_grid);
        temp_unfmt(:,i) = interp1(qc_ts(tag_no).raw_data(i).pres(~isnan(qc_ts(tag_no).raw_data(i).pres)), temp_interp{:,i}, depth_grid);
        vort_unfmt(:,i) = interp1(qc_ts(tag_no).raw_data(i).pres(~isnan(qc_ts(tag_no).raw_data(i).pres)), vort_interp{:,i}, depth_grid);
    end

    LLC_data(tag_no) = struct("tag", qc_ts(tag_no).tag, "cast", qc_ts(tag_no).cast,...
        "lat", qc_ts(tag_no).lat, "lon", qc_ts(tag_no).lon, "time", qc_ts(tag_no).time, "salt", salt_unfmt, "temp", temp_unfmt, ...
        "vort", vort_unfmt);
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
tag_no = 92;
depth_ind = 30;
%%% Plotting vorticity data
load coastlines
figure()
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -57],'MLineLocation', 30, 'PLineLocation', 10, 'FontSize',10);
axis off; framem on; gridm on; mlabel on; plabel on;
colormap(cmocean('balance'));
caxis([-0.3 0.3]);
colorbar;
hold on
pcolorm(LLC.lat(1:1:end,1:1:end), LLC.lon(1:1:end,1:1:end), LLC.vort(1:1:end,1:1:end,depth_ind));
geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor','white');
plotm(qc_ts(tag_no).lat, qc_ts(tag_no).lon, 'g-','Marker', '.','MarkerSize', 10, 'LineWidth', 2)
plotm(qc_ts(tag_no).lat(154), qc_ts(tag_no).lon(154),'c', 'Marker', '*','MarkerSize', 10, 'LineWidth', 2)
hold off
title('LLC Vorticity Data: ' + string(-LLC.depth(depth_ind)) + ' m')
