for tag_no = 200
    
    %%% Converting time data into datenum
    ts_time = datenum(qc_ts(tag_no).time);
    ts_months = month(datetime(qc_ts(tag_no).time));
    
    ts_density = gsw_rho(qc_ts(tag_no).salt, qc_ts(tag_no).temp, zeros(size(qc_ts(tag_no).salt,1), size(qc_ts(tag_no).salt,2)));
    
    [LAT,LON,DEPTH] = meshgrid(woa_lats(1:120), woa_lons, woa_depths);
    
    salt = [];
    temp = [];
    for i = unique(ts_months)'
        meop_section = find(ts_months == i);
        [lat,lon,depth] = meshgrid(qc_ts(tag_no).lat(meop_section), qc_ts(tag_no).lon(meop_section), depth_grid);
        salt(meop_section, meop_section',:) = interp3(LAT, LON, DEPTH, woa_salt(:,1:120,:,i), lat, lon, depth);
        temp(meop_section, meop_section',:) = interp3(LAT, LON, DEPTH, woa_temp(:,1:120,:,i), lat, lon, depth);
    end

%     [lat,lon,depth] = meshgrid(qc_ts(tag_no).lat, qc_ts(tag_no).lon, depth_grid);
%     salt = interp3(LAT, LON, DEPTH, woa_salt(:,1:120,:,6), lat, lon, depth);
%     temp = interp3(LAT, LON, DEPTH, woa_temp(:,1:120,:,6), lat, lon, depth);    
    
    salt_interp = [];
    temp_interp = [];
    for i = 1:length(qc_ts(tag_no).lat)
        salt_interp(:,i) = salt(i,i,:);
        temp_interp(:,i) = temp(i,i,:);
    end
    
    woa_density = gsw_rho(salt_interp, temp_interp, zeros(size(salt_interp,1), size(salt_interp,2)));
    
    delz = 5;
    
    figure('Renderer', 'painters', 'Position', [0 0 1000 850])
    
    %%% MEOP Subplot
    ax1 = subplot(4,1,1);
    hold on
    pp = pcolor(ax1, ts_time, depth_grid(1:delz:end), qc_ts(tag_no).salt(1:delz:end,:));
    [IA,IB] = unique(ts_time);
    [C,h] = contour(ax1, unique(ts_time), depth_grid, ts_density(:,IB), round(min(min(ts_density)):0.1:max(max(ts_density)), 1), 'k');
    hold off
    clabel(C,h,'LabelSpacing',500);
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
    [IA,IB] = unique(ts_time);
    [C,h] = contour(ax2, unique(ts_time), depth_grid, woa_density(:,IB), round(min(min(woa_density)):0.1:max(max(woa_density)), 1), 'k');
    hold off
    clabel(C,h,'LabelSpacing',500);
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
    title('WOA Salinity');
    
    %%% MEOP Subplot
    ax3 = subplot(4,1,3);
    hold on
    pp = pcolor(ax3, ts_time, depth_grid(1:delz:end), qc_ts(tag_no).temp(1:delz:end,:));
    [IA,IB] = unique(ts_time);
    [C,h] = contour(ax3, unique(ts_time), depth_grid, ts_density(:,IB), round(min(min(ts_density)):0.1:max(max(ts_density)), 1), 'k');
    hold off
    clabel(C,h,'LabelSpacing',500);
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
    [IA,IB] = unique(ts_time);
    [C,h] = contour(ax4, unique(ts_time), depth_grid, woa_density(:,IB), round(min(min(woa_density)):0.1:max(max(woa_density)), 1), 'k');
    hold off
    clabel(C,h,'LabelSpacing',500);
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
    title('WOA Temperature');
end