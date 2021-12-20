%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Formatting LLC Data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% LLC Data
LLC_data_1 = matfile('SO_snapshots_NSF_LLC4320_k1-86_face1_01-Dec-2011.mat');
LLC_data_2 = matfile('SO_snapshots_NSF_LLC4320_k1-86_face2_01-Dec-2011.mat');
LLC_data_4 = matfile('SO_snapshots_NSF_LLC4320_k1-86_face4_01-Dec-2011.mat');
LLC_data_5 = matfile('SO_snapshots_NSF_LLC4320_k1-86_face5_01-Dec-2011.mat');

%%% Longitude Data
xc_1 = double(LLC_data_1.xc);
xc_2 = double(LLC_data_2.xc);
xc_4 = double(LLC_data_4.xc);
xc_5 = double(LLC_data_5.xc);

%%% Latitude Data
yc_1 = double(LLC_data_1.yc);
yc_2 = double(LLC_data_2.yc);
yc_4 = double(LLC_data_4.yc);
yc_5 = double(LLC_data_5.yc);

%%% Depth Data
depth_4 = LLC_data_4.rc(1:40,:);

%%% Salinity Data
salt_1 = LLC_data_1.s(:,:,1:40);
salt_2 = LLC_data_2.s(:,:,1:40);
salt_4 = LLC_data_4.s(:,:,1:40);
salt_5 = LLC_data_5.s(:,:,1:40);

%%% Temperature Data
temp_1 = LLC_data_1.t(:,:,1:40);
temp_2 = LLC_data_2.t(:,:,1:40);
temp_4 = LLC_data_4.t(:,:,1:40);
temp_5 = LLC_data_5.t(:,:,1:40);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Visualizing Quad Locations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure();
hold on
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -57],'MLineLocation', 30, 'PLineLocation', 10, 'FontSize',10);
axis off;
framem on;
gridm on;
mlabel on;
plabel on;
plotm(yc_1, xc_1, 'b');
plotm(yc_2, xc_2, 'r');
plotm(yc_4, xc_4, 'k');
plotm(yc_5, xc_5, 'y');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LLC Data Near Sear Track %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Cutting LLC Data by Latitude
[~,col_lat] = find(yc_4 < max(qc_ts(tag_no).lat) & yc_4 > min(qc_ts(tag_no).lat));
LON = xc_4(:,min(col_lat):max(col_lat));
LAT = yc_4(:,min(col_lat):max(col_lat));
SALT = salt_4(:,min(col_lat):max(col_lat),:);
TEMP = temp_4(:,min(col_lat):max(col_lat),:);

%%% Cutting LLC Data by Longitude
[row_lon,~] = find(LON < max(qc_ts(tag_no).lon) & LON > min(qc_ts(tag_no).lon));
LON_final = LON(min(row_lon):max(row_lon),:);
LAT_final = LAT(min(row_lon):max(row_lon),:);
SALT_final = SALT(min(row_lon):max(row_lon),:,:);
TEMP_final = TEMP(min(row_lon):max(row_lon),:,:);
clear LON LAT SALT TEMP col_lat row_lon

%%% Formatting LLC Data for Interpolation
LON_final = LON_final .* ones(size(LON_final,1), size(LON_final,2), size(depth_4,1));
LAT_final = LAT_final .* ones(size(LAT_final,1), size(LAT_final,2), size(depth_4,1));
DEPTH(1,1,:) = -1 .*depth_4;
DEPTH_final = DEPTH .* ones(size(LAT_final,1), size(LAT_final,2), size(depth_4,1));
clear DEPTH

%%% Interpolating Data
tag_no = 131;
ts_time = datenum(qc_ts(tag_no).time);
[lat,lon,depth] = meshgrid(qc_ts(tag_no).lat, qc_ts(tag_no).lon, depth_grid);
salt = interp3(LAT_final, LON_final, DEPTH_final, SALT_final, lat, lon, depth);
temp = interp3(LAT_final, LON_final, DEPTH_final, TEMP_final, lat, lon, depth);
clear lat lon depth LAT_final LON_final DEPTH_final SALT_final

salt_interp = NaN(size(qc_ts(tag_no).salt));
temp_interp = NaN(size(qc_ts(tag_no).temp));
for i = 1:length(qc_ts(tag_no).lat)
    salt_interp(:,i) = salt(i,i,:);
    temp_interp(:,i) = temp(i,i,:);
end
clear salt temp i

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

