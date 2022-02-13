%%% Load data
% RTOPO_lat = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'lat'));
% RTOPO_lon = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'lon'))';
% RTOPO_bedrock_topography = double(ncread('/Users/jenkosty/Downloads/Research/detectSCV-main/RTOPO2.nc', 'bedrock_topography'))';

for tag_no = 1
    
    depth_grid = 1:500;

    %%% Converting time data into datenum
    ts_time = datetime(qc_ts(tag_no).time);
    ts_months = month(ts_time);

    %%% Getting the first and last months
    first_month = sprintf('%02d', min(ts_months));
    last_month = sprintf('%02d', max(ts_months));

    %%% Sea Ice Info
    first_sea_ice = shaperead('/Users/jenkosty/Downloads/Research/Sea_Ice_Extent/median_extent_S_' + string(first_month) + '_1981-2010_polyline_v3.0/median_extent_S_' + string(first_month) + '_1981-2010_polyline_v3.0.shp');
    first_proj = shapeinfo('/Users/jenkosty/Downloads/Research/Sea_Ice_Extent/median_extent_S_' + string(first_month) + '_1981-2010_polyline_v3.0/median_extent_S_' + string(first_month) + '_1981-2010_polyline_v3.0.shp').CoordinateReferenceSystem;
    [first_sea_ice_lat, first_sea_ice_lon] = projinv(first_proj, [first_sea_ice.X], [first_sea_ice.Y]);

    last_sea_ice = shaperead('/Users/jenkosty/Downloads/Research/Sea_Ice_Extent/median_extent_S_' + string(last_month) + '_1981-2010_polyline_v3.0/median_extent_S_' + string(last_month) + '_1981-2010_polyline_v3.0.shp');
    last_proj = shapeinfo('/Users/jenkosty/Downloads/Research/Sea_Ice_Extent/median_extent_S_' + string(last_month) + '_1981-2010_polyline_v3.0/median_extent_S_' + string(last_month) + '_1981-2010_polyline_v3.0.shp').CoordinateReferenceSystem;
    [last_sea_ice_lat, last_sea_ice_lon] = projinv(last_proj, [last_sea_ice.X], [last_sea_ice.Y]);

    %%% Pcolor vertical resolution
    delz = 5;
    
    fig = figure('Renderer', 'painters', 'Position', [0 0 1000 950]);

    %%% Full Map Subplot
    subplot(4,2,1)
    hold on
    load coastlines
    axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -57],'MLineLocation', 30, 'PLineLocation', 10, 'FontSize',10);
    axis off;
    framem on;
    gridm on;
    mlabel on;
    plabel on;
    plotm(coastlat, coastlon, 'k');
    a = plotm(first_sea_ice_lat, first_sea_ice_lon, 'g', 'DisplayName', 'Median ' + string(datestr(datetime(1,min(ts_months),1),'mmmm')) + ' Sea Ice');
    b = plotm(last_sea_ice_lat, last_sea_ice_lon, 'r', 'DisplayName', 'Median ' + string(datestr(datetime(1,max(ts_months),1),'mmmm')) + ' Sea Ice');
    c = plotm(qc_ts(tag_no).lat, qc_ts(tag_no).lon, '.', 'DisplayName','Seal Track');
    plotm(qc_ts(tag_no).lat(1), qc_ts(tag_no).lon(1), 'linestyle', 'none', 'Marker','o', 'MarkerFaceColor','[0 1 0]', 'MarkerSize', 7, 'DisplayName', 'Track Beginning');
    plotm(qc_ts(tag_no).lat(end), qc_ts(tag_no).lon(end), 'linestyle', 'none', 'Marker','o', 'MarkerFaceColor','[1 0 0]', 'MarkerSize', 7, 'DisplayName', 'Track Beginning');
    hold off
    legend([a,b,c], 'Location', 'best')

    %%% Zoomed Map Subplot
    subplot(4,2,2)
    hold on
    axesm('lambertstd','MapParallels',[-75 -15],'MapLatLimit',[min(qc_ts(tag_no).lat)-1 max(qc_ts(tag_no).lat)+1],'MapLonLimit',[min(qc_ts(tag_no).lon)-1 max(qc_ts(tag_no).lon)+1], 'MLineLocation', 10, 'PLineLocation', 2, 'FontSize',10);
    axis off;
    framem on;
    gridm on;
    mlabel on;
    plabel on;
    plotm(coastlat, coastlon, 'k');
    a = plotm(first_sea_ice_lat, first_sea_ice_lon, 'g', 'DisplayName', 'Median ' + string(datestr(datetime(1,min(ts_months),1),'mmmm')) + ' Sea Ice');
    b = plotm(last_sea_ice_lat, last_sea_ice_lon, 'r', 'DisplayName', 'Median ' + string(datestr(datetime(1,max(ts_months),1),'mmmm')) + ' Sea Ice');
    c = scatterm(qc_ts(tag_no).lat, qc_ts(tag_no).lon,10,datenum(qc_ts(tag_no).time),'filled', 'DisplayName','Seal Track');
    plotm(qc_ts(tag_no).lat(1), qc_ts(tag_no).lon(1), 'linestyle', 'none', 'Marker','o', 'MarkerFaceColor','[0 1 0]', 'MarkerSize', 7, 'DisplayName', 'Track Beginning');
    plotm(qc_ts(tag_no).lat(end), qc_ts(tag_no).lon(end), 'linestyle', 'none', 'Marker','o', 'MarkerFaceColor','[1 0 0]', 'MarkerSize', 7, 'DisplayName', 'Track Beginning');
    hold off
    a = colorbar;
    cbdate;

    clear ts_months first_month last_month first_sea_ice first_proj first_sea_ice_lat first_sea_ice_lon last_sea_ice last_proj last_sea_ice_lat last_sea_ice_lon a b c ts_time

    ts_time = datenum(qc_ts(tag_no).time);

    %%% Bathymetry Calculation
    ts_bathymetry = interp2(RTOPO_lon, RTOPO_lat', RTOPO_bedrock_topography, qc_ts(tag_no).lon, qc_ts(tag_no).lat);

    %%% Density Calculation
    ts_density = gsw_rho(qc_ts(tag_no).salt, qc_ts(tag_no).temp, zeros(size(qc_ts(tag_no).salt,1), size(qc_ts(tag_no).salt,2)));

    %%% Bathymetry Subplot
    ax1 = subplot(4,2,3:4);
    plot(ts_time, ts_bathymetry)
    xticks(linspace(ts_time(1), ts_time(end), (ts_time(end) - ts_time(1)) / 5))
    datetick('x', 'mm/dd/yy', 'keepticks');
    ylabel('Depth (m)','FontSize',12);
    title('Bedrock Topography','FontSize', 12);
    
    %%% MEOP Subplot
    ax2 = subplot(4,2,5:6);
    hold on
    pp = pcolor(ax2, ts_time, depth_grid(1:delz:end), qc_ts(tag_no).salt(1:delz:end,:));
    [~,IB] = unique(ts_time);
    [C,h] = contour(ax2, unique(ts_time), depth_grid, ts_density(:,IB), round(min(min(ts_density)):0.04:max(max(ts_density)), 2), 'k');
    hold off
    clabel(C,h,'LabelSpacing',500);
    set(pp, 'EdgeColor', 'none');
    set(gca, 'YDir','reverse');
    set(gca, 'Layer','top');
    cmap = cmocean('haline');
    colormap(ax2, cmap);
    a = colorbar;
    caxis([min(min(qc_ts(tag_no).salt)) max(max(qc_ts(tag_no).salt))])
    %ylabel(a,'Salinity (g/kg)','FontSize',16,'Rotation',270);
    xticks(linspace(ts_time(1), ts_time(end), (ts_time(end) - ts_time(1)) / 5));
    datetick('x', 'mm/dd/yy', 'keepticks');
    ylabel('Pressure (dbar)', 'FontSize', 12);
    title('Salinity', 'FontSize', 12);
    
    %%% MEOP Subplot
    ax3 = subplot(4,2,7:8);
    hold on
    pp = pcolor(ax3, ts_time, depth_grid(1:delz:end), qc_ts(tag_no).temp(1:delz:end,:));
    [~,IB] = unique(ts_time);
    [C,h] = contour(ax3, unique(ts_time), depth_grid, ts_density(:,IB), round(min(min(ts_density)):0.04:max(max(ts_density)), 2), 'k');
    hold off
    clabel(C,h,'LabelSpacing',500);
    set(pp, 'EdgeColor', 'none');
    set(gca, 'YDir','reverse');
    set(gca, 'Layer','top');
    cmap = cmocean('thermal');
    colormap(ax3, cmap);
    a = colorbar;
    caxis([min(min(qc_ts(tag_no).temp)) max(max(qc_ts(tag_no).temp))])
    %ylabel(a,'Temperature (' + string(char(176)) + 'C)','FontSize',16,'Rotation',270);
    xticks(linspace(ts_time(1), ts_time(end), (ts_time(end) - ts_time(1)) / 5));
    datetick('x', 'mm/dd/yy', 'keepticks');
    ylabel('Pressure (dbar)', 'FontSize', 12);
    title('Temperature', 'FontSize', 12);
     
    p1 = get(ax1, 'Position');
    p2 = get(ax2, 'Position');
    p3 = get(ax3, 'Position');

    p2(3) = p1(3);
    p3(3) = p1(3);

    set(ax2, 'Position', p2);
    set(ax3, 'Position', p3);

    saveas(fig, '/Users/jenkosty/Downloads/Research/detectSCV-main/MEOP Time Series/Seal_' + string(qc_ts(tag_no).tag), 'fig')
    
    clear ax1 ax2 ax3 IB C h pp p1 p2 p3 delz ts_density ts_time cmap ts_bathymetry a
end
%%
% figure();
% hold on
% load coastlines
% axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -57],'MLineLocation', 30, 'PLineLocation', 10, 'FontSize',10);
% axis off;
% framem on;
% gridm on;
% mlabel on;
% plabel on;
% plotm(coastlat, coastlon, 'k');
% plotm(qc_ts(tag_no).lat, qc_ts(tag_no).lon, '.', 'DisplayName','Seal Track');
% plotm(qc_ts(tag_no).lat(1), qc_ts(tag_no).lon(1), 'linestyle', 'none', 'Marker','o', 'MarkerFaceColor','[0 1 0]', 'MarkerSize', 7, 'DisplayName', 'Track Beginning');
% plotm(qc_ts(tag_no).lat(end), qc_ts(tag_no).lon(end), 'linestyle', 'none', 'Marker','o', 'MarkerFaceColor','[1 0 0]', 'MarkerSize', 7, 'DisplayName', 'Track Beginning');
% hold off

    