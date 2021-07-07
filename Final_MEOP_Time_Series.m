%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Final MEOP Time Series %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tag_no = 70;
    
    day_average = 1; % half-range for two-day running mean
    ID = raw_tag_final(tag_no); % MEOP Tag ID
    
    %%% creating time matrix for seal profiles (year, month, day, hour, minute, and second separated in different columns)
    ts_time = [];
    for i = 1:length(interp_meop_data(tag_elements_final{1, tag_no}))
        ts_time(i,:) = interp_meop_data(tag_elements_final{1, tag_no}(i)).time;
    end
    
    %%% converting seal profile times into serial date numbers (used as x-axis of
    %%% time series)
    ts_dates = [];
    for i = 1:length(tag_elements_final{1, tag_no})
        ts_dates(i) = datenum([ts_time(i,1),ts_time(i,2),ts_time(i,3),ts_time(i,4),ts_time(i,5),ts_time(i,6)]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Pressure Space Calculations %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% creating matrices to hold the time series (colorbar) salinity, temperature, and density data
    %%% Note: dates must be unique for the time series
    u = 1;
    ts_salt = [];
    ts_temp = [];
    ts_density = [];
    for i = unique(ts_dates)
        l = find(ts_dates == i, 1);
        ts_salt(:,u) = interp_meop_data(tag_elements_final{1, tag_no}(l)).salt;
        ts_temp(:,u) = interp_meop_data(tag_elements_final{1, tag_no}(l)).temp;
        ts_density(:,u) = gsw_rho(ts_salt(:,u), ts_temp(:,u), zeros(1, length(ts_salt(:,u)))) - 1000;
        u = u + 1;
    end
    
    %%% Computing averaging windows for two-day running mean
    %%% Note: number of profiles collected every two days is not consistent
    j = 1;
    ts_dates_tda = {};
    for i = unique(ts_dates)
        ts_dates_tda(1,j) = {find(ts_dates > i & ts_dates < i+day_average)};
        ts_dates_tda(2,j) = {find(ts_dates < i & ts_dates > i-day_average)};
        j = j + 1;
    end
    
    %%% Computing temperature, salinity, and density data for two-day running mean
    ts_salt_runmean = [];
    ts_salt_runmean_initial = [];
    ts_temp_runmean = [];
    ts_temp_runmean_initial = [];
    ts_density_runmean = [];
    ts_density_runmean_initial = [];
    
    for i = 1:length(woa_depths(1:41))
        for j = 1:length(ts_dates_tda)
            ts_salt_runmean_initial(i,:) = movmean(ts_salt(i,:), [length(ts_dates_tda{1,j}) length(ts_dates_tda{2,j})], 'omitnan');
            ts_salt_runmean(i,j) = ts_salt_runmean_initial(i,j);
            ts_temp_runmean_initial(i,:) = movmean(ts_temp(i,:), [length(ts_dates_tda{1,j}) length(ts_dates_tda{2,j})], 'omitnan');
            ts_temp_runmean(i,j) = ts_temp_runmean_initial(i,j);
            ts_density_runmean_initial(i,:) = movmean(ts_density(i,:), [length(ts_dates_tda{1,j}) length(ts_dates_tda{2,j})], 'omitnan');
            ts_density_runmean(i,j) = ts_density_runmean_initial(i,j);
        end
    end
    
    %%% Computing distance between isopycnals for two-day running mean
    isopycnal_dist_runmean = [];
    for i = 1:length(unique(ts_dates))
        isopycnal_dist_runmean(1,i) = length(find(ts_density_runmean(:,i) >= 27.5 & ts_density_runmean(:,i) <= 27.6));
        isopycnal_dist_runmean(2,i) = length(find(ts_density_runmean(:,i) >= 27.4 & ts_density_runmean(:,i) <= 27.6));
        isopycnal_dist_runmean(3,i) = length(find(ts_density_runmean(:,i) >= 27.3 & ts_density_runmean(:,i) <= 27.6));
    end
        
    %%% creating figure and general title
    fig = figure('Renderer', 'painters', 'Position', [0 0 700 850]);
    sgtitle('MEOP Tag ' + string(ID), 'FontSize', 18, 'FontWeight', 'bold')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Creating Isopycnal Separation Part of Figure %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1 = subplot(5,1,1);
    
    hold on
    plot(unique(ts_dates), isopycnal_dist_runmean(1,:) * 5, 'DisplayName', '27.5 - 27.6')
    plot(unique(ts_dates), isopycnal_dist_runmean(2,:) * 5, 'DisplayName', '27.4 - 27.6')
    plot(unique(ts_dates), isopycnal_dist_runmean(3,:) * 5, 'DisplayName', '27.3 - 27.6')
    hold off
    
    title('Running Mean Isopycnal Separation')
    ylabel('Isopycnal Separation (m)')
    set(gca, 'Layer','top')
    xlim([min(unique(ts_dates)), max(unique(ts_dates))]);
    xticks(linspace(min(unique(ts_dates)), max(unique(ts_dates)), (max(unique(ts_dates)) - min(unique(ts_dates))) / 3))
    datetick('x', 'mm/dd/yy', 'keepticks');
    
    ax = gca;
    labels_1 = string(ax.XAxis.TickLabels);
    labels = string(ax.XAxis.TickLabels); % extract
    labels(1:5:end) = nan; % remove every other one
    
    for i = 1:length(labels)
        if strcmp(labels_1(i), labels(i)) == 0
            labels(i) = labels_1(i);
        else
            labels(i) = nan;
        end
    end
    
    ax.XAxis.TickLabels = labels; % set
    
    legend();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Creating Pressure Space Part of Figure %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% salinity time series
    ax2 = subplot(5,1,2);
    time_series_gen(unique(ts_dates), woa_depths(1:41),'pressure', ts_salt_runmean, 'salinity', ax2)
    contours(unique(ts_dates), woa_depths(1:41), 'density', ts_salt_runmean, ts_temp_runmean)
    title('Running Mean Salinity Time Series - Pressure Space')
    
    %%% temperature time series
    ax4 = subplot(5,1,4);
    time_series_gen(unique(ts_dates), woa_depths(1:41),'pressure', ts_temp_runmean, 'temperature', ax4)
    contours(unique(ts_dates), woa_depths(1:41), 'density', ts_salt_runmean, ts_temp_runmean)
    title('Running Mean Temperature Time Series - Pressure Space')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Density Space Calculations %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% creating matrices to hold the time series (colorbar) salinity, temperature, and density data
    %%% Note: dates must be unique for the time series
    u = 1;
    ts_salt = [];
    ts_temp = [];
    ts_y_range = [];
    for i = unique(ts_dates)
        l = find(ts_dates == i, 1);
        ts_salt(:,u) = ds_interp_meop_data(tag_elements_final{1, tag_no}(l)).salt;
        ts_temp(:,u) = ds_interp_meop_data(tag_elements_final{1, tag_no}(l)).temp;
        ts_y_range(:,u) = isnan(ds_interp_meop_data(tag_elements_final{1, tag_no}(l)).salt);
        u = u + 1;
    end
    
    %%% Calculating y-axis range for time series
    density_elements = [];
    density_elements = find(all(ts_y_range, 2) == 0);
    
    density_range = [];
    density_range = [density_grid(min(density_elements)), density_grid(max(density_elements))];
    
    %%% Computing salinty and temperature data for running mean
    ts_salt_runmean = [];
    ts_salt_runmean_initial = [];
    ts_temp_runmean = [];
    ts_temp_runmean_initial = [];
    for i = 1:length(density_grid)
        for j = 1:length(ts_dates_tda)
            ts_salt_runmean_initial(i,:) = movmean(ts_salt(i,:), [length(ts_dates_tda{1,j}) length(ts_dates_tda{2,j})], 'omitnan');
            ts_salt_runmean(i,j) = ts_salt_runmean_initial(i,j);
            ts_temp_runmean_initial(i,:) = movmean(ts_temp(i,:), [length(ts_dates_tda{1,j}) length(ts_dates_tda{2,j})], 'omitnan');
            ts_temp_runmean(i,j) = ts_temp_runmean_initial(i,j);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Creating Density Space Part of Figure %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% salinity time series
    ax3 = subplot(5,1,3);
    time_series_gen(unique(ts_dates), density_grid, 'density', ts_salt_runmean, 'salinity', ax3)
    contours(unique(ts_dates), density_grid, 'spice', ts_salt_runmean, ts_temp_runmean)
    title('Running Mean Salinity Time Series - Density Space')
    
    %%% temperature time series
    ax5 = subplot(5,1,5);
    time_series_gen(unique(ts_dates), density_grid, 'density', ts_temp_runmean, 'temperature', ax5)
    contours(unique(ts_dates), density_grid, 'spice', ts_salt_runmean, ts_temp_runmean)
    title('Running Mean Temperature Time Series - Density Space')
    
    pos1 = get(ax1, 'Position');
    pos2 = get(ax2, 'Position');
    pos3 = get(ax3, 'Position');
    pos4 = get(ax4, 'Position');
    pos5 = get(ax5, 'Position');
    pos2(3) = pos1(3);
    pos3(3) = pos1(3);
    pos4(3) = pos1(3);
    pos5(3) = pos1(3);
    set(ax2, 'Position', pos2);
    set(ax3, 'Position', pos3);
    set(ax4, 'Position', pos4);
    set(ax5, 'Position', pos5);
    
    %%% assigning variables for map function
    lat = double([interp_meop_data(tag_elements_final{1, tag_no}).lat]);
    lon = double([interp_meop_data(tag_elements_final{1, tag_no}).lon]);

    subplot(6,2,11);
    map(lat, lon, ID, 'meop');
    
    subplot(6,2,12)
    map_zoom(lat, lon, ID, 'meop');

    %saveas(fig, '/Users/jenkosty/Downloads/detectSCV-main/MEOP Time Series/MEOP TS ' + string(tag_no), 'png')
end
