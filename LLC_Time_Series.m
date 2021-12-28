%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LLC Data Near Seal Track %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tag_no = 131;
LLC = LLC_4;

%%% Cutting LLC Data by Latitude
[~,col_lat] = find(LLC.lat < max(qc_ts(tag_no).lat) & LLC.lat > min(qc_ts(tag_no).lat));
LON = LLC.lon(:,min(col_lat):max(col_lat));
LAT = LLC.lat(:,min(col_lat):max(col_lat));
SALT = LLC.salt(:,min(col_lat):max(col_lat),:);
TEMP = LLC.temp(:,min(col_lat):max(col_lat),:);

%%% Cutting LLC Data by Longitude
[row_lon,~] = find(LON < max(qc_ts(tag_no).lon) & LON > min(qc_ts(tag_no).lon));
LON_final = LON(min(row_lon):max(row_lon),:);
LAT_final = LAT(min(row_lon):max(row_lon),:);
SALT_final = SALT(min(row_lon):max(row_lon),:,:);
TEMP_final = TEMP(min(row_lon):max(row_lon),:,:);
clear LON LAT SALT TEMP col_lat row_lon

%%% Formatting LLC Data for Interpolation
LON_final = LON_final .* ones(size(LON_final,1), size(LON_final,2), size(LLC.depth,1));
LAT_final = LAT_final .* ones(size(LAT_final,1), size(LAT_final,2), size(LLC.depth,1));
DEPTH(1,1,:) = -1 .* LLC.depth;
DEPTH_final = DEPTH .* ones(size(LAT_final,1), size(LAT_final,2), size(LLC.depth,1));
clear DEPTH

%%% Interpolating Data
ts_time = datenum(qc_ts(tag_no).time);
[lat,lon,depth] = meshgrid(qc_ts(tag_no).lat, qc_ts(tag_no).lon, depth_grid);
%salt = interp3(LAT_final, LON_final, DEPTH_final, SALT_final, lat, lon, depth);
temp = interp3(LAT_final, LON_final, DEPTH_final, TEMP_final, lat, lon, depth);
%clear lat lon depth DEPTH_final SALT_final

salt_interp = NaN(size(qc_ts(tag_no).salt));
temp_interp = NaN(size(qc_ts(tag_no).temp));
for i = 1:length(qc_ts(tag_no).lat)
    %salt_interp(:,i) = salt(i,i,:);
    temp_interp(:,i) = temp(i,i,:);
end
%clear salt temp i

parfor i = 1:length(qc_ts(tag_no).lat)
    salt_interp(:,i) = interp3(LAT_final, LON_final, DEPTH_final, SALT_final,qc_ts(tag_no).lat(i), qc_ts(tag_no).lon(i),depth_grid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Visualizing LLC Data Along Seal Track %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
delz = 5;
figure('Renderer', 'painters', 'Position', [0 0 1000 850])

%%% MEOP Subplot
ax1 = subplot(4,1,1);
hold on
pp = pcolor(ax1, ts_time, depth_grid(1:delz:end), qc_ts(tag_no).salt(1:delz:end,:));
%[IA,IB] = unique(ts_time);
%[C,h] = contour(ax1, unique(ts_time), depth_grid, ts_density(:,IB), round(min(min(ts_density)):0.1:max(max(ts_density)), 1), 'k');
hold off
%clabel(C,h,'LabelSpacing',500);
set(pp, 'EdgeColor', 'none');
set(gca, 'YDir','reverse');
set(gca, 'Layer','top');
cmap = cmocean('haline');
colormap(ax1, cmap);
colorbar;
caxis([min(min(qc_ts(tag_no).salt)) max(max(qc_ts(tag_no).salt))])
xticks(linspace(ts_time(1), ts_time(end), (ts_time(end) - ts_time(1)) / 5));
datetick('x', 'mm/dd', 'keepticks');
ylabel('Pressure (dbar)');
title('MEOP Salinity');

%%% WOA Subplot
ax2 = subplot(4,1,2);
hold on
pp = pcolor(ax2, ts_time, depth_grid, salt_interp);
%[IA,IB] = unique(ts_time);
%[C,h] = contour(ax2, unique(ts_time), depth_grid, woa_density(:,IB), round(min(min(woa_density)):0.1:max(max(woa_density)), 1), 'k');
hold off
%clabel(C,h,'LabelSpacing',500);
set(pp, 'EdgeColor', 'none');
set(gca, 'YDir','reverse');
set(gca, 'Layer','top');
cmap = cmocean('haline');
colormap(ax2, cmap);
colorbar;
caxis([min(min(salt_interp)) max(max(salt_interp))])
xticks(linspace(ts_time(1), ts_time(end), (ts_time(end) - ts_time(1)) / 5));
datetick('x', 'mm/dd', 'keepticks');
ylabel('Pressure (dbar)');
title('LLC Salinity');

%%% MEOP Subplot
ax3 = subplot(4,1,3);
hold on
pp = pcolor(ax3, ts_time, depth_grid(1:delz:end), qc_ts(tag_no).temp(1:delz:end,:));
%[IA,IB] = unique(ts_time);
%[C,h] = contour(ax3, unique(ts_time), depth_grid, ts_density(:,IB), round(min(min(ts_density)):0.1:max(max(ts_density)), 1), 'k');
hold off
%clabel(C,h,'LabelSpacing',500);
set(pp, 'EdgeColor', 'none');
set(gca, 'YDir','reverse');
set(gca, 'Layer','top');
cmap = cmocean('thermal');
colormap(ax3, cmap);
colorbar;
caxis([min(min(qc_ts(tag_no).temp)) max(max(qc_ts(tag_no).temp))])
xticks(linspace(ts_time(1), ts_time(end), (ts_time(end) - ts_time(1)) / 5));
datetick('x', 'mm/dd', 'keepticks');
ylabel('Pressure (dbar)');
title('MEOP Temperature');

%%% WOA Subplot
ax4 = subplot(4,1,4);
hold on
pp = pcolor(ax4, ts_time, depth_grid, temp_interp);
%[IA,IB] = unique(ts_time);
%[C,h] = contour(ax4, unique(ts_time), depth_grid, woa_density(:,IB), round(min(min(woa_density)):0.1:max(max(woa_density)), 1), 'k');
hold off
%clabel(C,h,'LabelSpacing',500);
set(pp, 'EdgeColor', 'none');
set(gca, 'YDir','reverse');
set(gca, 'Layer','top');
cmap = cmocean('thermal');
colormap(ax4, cmap);
colorbar;
caxis([min(min(qc_ts(tag_no).temp)) max(max(qc_ts(tag_no).temp))])
xticks(linspace(ts_time(1), ts_time(end), (ts_time(end) - ts_time(1)) / 5));
datetick('x', 'mm/dd', 'keepticks');
ylabel('Pressure (dbar)');
title('LLC Temperature');

