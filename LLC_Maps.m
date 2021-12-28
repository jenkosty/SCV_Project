level = 25;

% LLC_1.salt(LLC_1.salt == 0) = NaN;
% LLC_2.salt(LLC_2.salt == 0) = NaN;
% LLC_4.salt(LLC_4.salt == 0) = NaN;
% LLC_5.salt(LLC_5.salt == 0) = NaN;

figure('Renderer', 'painters', 'Position', [0 0 1300 700]);
sgtitle('Depth: ' + string(LLC_1.depth(level)));

ax1 = subplot(1,2,1);
hold on
load coastlines
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -60],'MLineLocation', 30, 'PLineLocation', 10, 'FontSize',10);
axis off;
framem on;
gridm on;
mlabel on;
plabel on;
cmap = cmocean('haline');
colormap(ax1, cmap);
hold on
%pcolorm(LLC_1.lat(1:5:end, 1:5:end), LLC_1.lon(1:5:end,1:5:end), LLC_1.salt(1:5:end,1:5:end,level));
%pcolorm(LLC_2.lat(1:5:end,1:5:end), LLC_2.lon(1:5:end,1:5:end), LLC_2.salt(1:5:end,1:5:end,level));
%pcolorm(LLC_4.lat(1:5:end,1:5:end), LLC_4.lon(1:5:end,1:5:end), LLC_4.salt(1:5:end,1:5:end,level));
pcolorm(LLC_5.lat(1:5:end,1:5:end), LLC_5.lon(1:5:end,1:5:end), LLC_5.salt(1:5:end,1:5:end,level));
geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor','white')
hold off
colorbar;

ax2 = subplot(1,2,2);
hold on
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -60],'MLineLocation', 30, 'PLineLocation', 10, 'FontSize',10);
axis off;
framem on;
gridm on;
mlabel on;
plabel on;
cmap = cmocean('thermal');
colormap(ax2, cmap);
hold on
%pcolorm(LLC_1.lat(1:5:end,1:5:end), LLC_1.lon(1:5:end,1:5:end), LLC_1.temp(1:5:end,1:5:end,level));
%pcolorm(LLC_2.lat(1:5:end,1:5:end), LLC_2.lon(1:5:end,1:5:end), LLC_2.temp(1:5:end,1:5:end,level));
%pcolorm(LLC_4.lat(1:5:end,1:5:end), LLC_4.lon(1:5:end,1:5:end), LLC_4.temp(1:5:end,1:5:end,level));
pcolorm(LLC_5.lat(1:5:end,1:5:end), LLC_5.lon(1:5:end,1:5:end), LLC_5.temp(1:5:end,1:5:end,level));
geoshow(coastlat, coastlon, 'DisplayType','polygon','FaceColor','white');
hold off
colorbar;


