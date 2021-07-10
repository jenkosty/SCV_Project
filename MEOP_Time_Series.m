%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pressure-Space MEOP Time Series %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Script creates 2 pressure-space time series figures (raw data, two-day
%%% running mean)

tag_no = 70; % seal tag index
day_average = 1; % half-range for two-day running mean

%%% assigning variables for map function
lat = double([interp_meop_data(tag_elements_final{1, tag_no}).lat]);
lon = double([interp_meop_data(tag_elements_final{1, tag_no}).lon]);
ID = raw_tag_final(tag_no);

%%% map function displays tag movement over time
figure();
map(lat, lon, ID, 'meop');

figure();
map_zoom(lat, lon, ID, 'meop');

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

%%% Stratification Calculation (partial derivative of density with respect
%%% to depth)
M_2 = [];
N_2 = [];
[M_2, N_2] = gradient(ts_density, 5);

%%% Identifying indicies above different isopycnals and calculating the
%%% mean N^2 above each isopycnal
low_strat = {};
isopycnals = {};
isopycnal_dist = [];
mean_N_2 = [];

isopycnal_dist_y_axis = [];

isopycnal_dist_y_axis = linspace(min(min(ts_density)), max(max(ts_density)), 20);

for i = 1:length(unique(ts_dates))
    low_strat(i) = {find(N_2(:,i) <= 0.001)}; %%% removing stratification less than 0.001
    N_2(low_strat{1,i}, i) = NaN;
    isopycnals(1,i) = {find(ts_density(:,i) < 27.6)};
    isopycnals(2,i) = {find(ts_density(:,i) < 27.4)};
    isopycnals(3,i) = {find(ts_density(:,i) < 27.2)};
    isopycnals(4,i) = {find(ts_density(:,i) < 27.0)};
   
    mean_N_2(1,i) = mean(N_2(isopycnals{1,i}, i), 'omitnan');
    mean_N_2(2,i) = mean(N_2(isopycnals{2,i}, i), 'omitnan');
    mean_N_2(3,i) = mean(N_2(isopycnals{3,i}, i), 'omitnan');
    mean_N_2(4,i) = mean(N_2(isopycnals{4,i}, i), 'omitnan');
    
    
    isopycnal_dist(1,i) = length(find(ts_density(:,i) >= 27.4 & ts_density(:,i) <= 27.6));
    isopycnal_dist(2,i) = length(find(ts_density(:,i) >= 27.3 & ts_density(:,i) <= 27.5));
    isopycnal_dist(3,i) = length(find(ts_density(:,i) >= 27.2 & ts_density(:,i) <= 27.4));
    isopycnal_dist(4,i) = length(find(ts_density(:,i) >= 27.1 & ts_density(:,i) <= 27.3));
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
N_2_runmean_initial = [];
N_2_runmean = [];

for i = 1:length(woa_depths(1:41))
    for j = 1:length(ts_dates_tda)
        ts_salt_runmean_initial(i,:) = movmean(ts_salt(i,:), [length(ts_dates_tda{1,j}) length(ts_dates_tda{2,j})], 'omitnan');
        ts_salt_runmean(i,j) = ts_salt_runmean_initial(i,j);
        ts_temp_runmean_initial(i,:) = movmean(ts_temp(i,:), [length(ts_dates_tda{1,j}) length(ts_dates_tda{2,j})], 'omitnan');
        ts_temp_runmean(i,j) = ts_temp_runmean_initial(i,j);
        ts_density_runmean_initial(i,:) = movmean(ts_density(i,:), [length(ts_dates_tda{1,j}) length(ts_dates_tda{2,j})], 'omitnan');
        ts_density_runmean(i,j) = ts_density_runmean_initial(i,j);
        N_2_runmean_initial(i,:) = movmean(N_2(i,:), [length(ts_dates_tda{1,j}) length(ts_dates_tda{2,j})], 'omitnan');
        N_2_runmean(i,j) = N_2_runmean_initial(i,j);
    end
end

isopycnals_runmean = {};
isopycnal_dist_runmean = [];
mean_N_2_runmean = [];
for i = 1:length(unique(ts_dates))
    isopycnals_runmean(1,i) = {find(ts_density_runmean(:,i) < 27.6)};
    isopycnals_runmean(2,i) = {find(ts_density_runmean(:,i) < 27.4)};
    isopycnals_runmean(3,i) = {find(ts_density_runmean(:,i) < 27.2)};
    isopycnals_runmean(4,i) = {find(ts_density_runmean(:,i) < 27.0)};
    
    mean_N_2_runmean(1,i) = mean(N_2(isopycnals_runmean{1,i}, i), 'omitnan');
    mean_N_2_runmean(2,i) = mean(N_2(isopycnals_runmean{2,i}, i), 'omitnan');
    mean_N_2_runmean(3,i) = mean(N_2(isopycnals_runmean{3,i}, i), 'omitnan');
    mean_N_2_runmean(4,i) = mean(N_2(isopycnals_runmean{4,i}, i), 'omitnan');
    
    isopycnal_dist_runmean(1,i) = length(find(ts_density_runmean(:,i) >= 27.4 & ts_density_runmean(:,i) <= 27.6));
    isopycnal_dist_runmean(2,i) = length(find(ts_density_runmean(:,i) >= 27.3 & ts_density_runmean(:,i) <= 27.5));
    isopycnal_dist_runmean(3,i) = length(find(ts_density_runmean(:,i) >= 27.2 & ts_density_runmean(:,i) <= 27.4));
    isopycnal_dist_runmean(4,i) = length(find(ts_density_runmean(:,i) >= 27.1 & ts_density_runmean(:,i) <= 27.3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% creating raw time series figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);
sgtitle('MEOP Tag ' + string(ID), 'FontSize', 18, 'FontWeight', 'bold')

%%% Ploting mean Stratification above different isopycnals
ax1 = subplot(5,1,1);

hold on
plot(unique(ts_dates), mean_N_2(1,:), 'LineWidth', 2, 'DisplayName', '27.6')
plot(unique(ts_dates), mean_N_2(2,:), 'LineWidth', 2, 'DisplayName', '27.4')
plot(unique(ts_dates), mean_N_2(3,:), 'LineWidth', 2, 'DisplayName', '27.2')
plot(unique(ts_dates), mean_N_2(4,:), 'LineWidth', 2, 'DisplayName', '27.0')
hold off

datetick('x', 'mm/dd/yyyy');
xlim([min(unique(ts_dates)), max(unique(ts_dates))]);
title('Mean N^2 Above Different Isopycnals')
ylabel('Mean N^2')
legend();

%%% Plotting distance between two isopycnals
ax2 = subplot(5,1,2);

hold on
plot(unique(ts_dates), isopycnal_dist(1,:) * 5, 'DisplayName', '27.4 - 27.6')
plot(unique(ts_dates), isopycnal_dist(2,:) * 5, 'DisplayName', '27.3 - 27.5')
plot(unique(ts_dates), isopycnal_dist(3,:) * 5, 'DisplayName', '27.2 - 27.4')
plot(unique(ts_dates), isopycnal_dist(4,:) * 5, 'DisplayName', '27.1 - 27.3')
hold off

datetick('x', 'mm/dd/yyyy');
xlim([min(unique(ts_dates)), max(unique(ts_dates))]);
title('Distance Between Isopycnals')
ylabel('Isopycnal Separation (m)')
legend();

%%% salinity time series
ax3 = subplot(5,1,3);
time_series_gen(unique(ts_dates), woa_depths(1:41),'pressure', ts_salt, 'salinity', ax3)
contours(unique(ts_dates), woa_depths(1:41), 'density', ts_salt, ts_temp)
title('Salinity Time Series')

%%% temperature time series
ax4 = subplot(5,1,4);
time_series_gen(unique(ts_dates), woa_depths(1:41),'pressure', ts_temp, 'temperature', ax4)
contours(unique(ts_dates), woa_depths(1:41), 'density', ts_salt, ts_temp)
title('Temperature Time Series')

%%% density time series
ax5 = subplot(5,1,5);
time_series_gen(unique(ts_dates), woa_depths(1:41),'pressure', ts_density, 'density', ax5)
contours(unique(ts_dates), woa_depths(1:41), 'density', ts_salt, ts_temp)
title('Density Time Series')

%%% Realigning subplot 0 to match color-bar axes
pos1 = get(ax1, 'Position');
pos3 = get(ax3, 'Position'); 
pos4 = get(ax4, 'Position');
pos5 = get(ax5, 'Position');
pos3(3) = pos1(3);
pos4(3) = pos1(3);
pos5(3) = pos1(3);
set(ax3, 'Position', pos3);
set(ax4, 'Position', pos4);
set(ax5, 'Position', pos5);

%saveas(gcf,'/Users/jenkosty/Downloads/detectSCV-main/Plots/WOA, Argo, and MEOP Plots/6-25-2021/' + string(ID) + ' (ps)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% creating running mean time series figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);
sgtitle('MEOP Tag ' + string(ID), 'FontSize', 18, 'FontWeight', 'bold')

%%% Ploting mean Stratification above different isopycnals
ax1 = subplot(5,1,1);

hold on
plot(unique(ts_dates), mean_N_2_runmean(1,:), 'LineWidth', 2, 'DisplayName', '27.6')
plot(unique(ts_dates), mean_N_2_runmean(2,:), 'LineWidth', 2, 'DisplayName', '27.4')
plot(unique(ts_dates), mean_N_2_runmean(3,:), 'LineWidth', 2, 'DisplayName', '27.2')
plot(unique(ts_dates), mean_N_2_runmean(4,:), 'LineWidth', 2, 'DisplayName', '27.0')
hold off;

datetick('x', 'mm/dd/yyyy');
xlim([min(ts_dates), max(ts_dates)]);
title('Running Mean Average N^2 Above Different Isopycnals')
ylabel('Mean N^2')
legend();

%%% Plotting distance between two isopycnals
ax2 = subplot(5,1,2);

hold on
plot(unique(ts_dates), isopycnal_dist_runmean(1,:)* 5, 'DisplayName', '27.4 - 27.6')
plot(unique(ts_dates), isopycnal_dist_runmean(2,:)* 5, 'DisplayName', '27.3 - 27.5')
plot(unique(ts_dates), isopycnal_dist_runmean(3,:)* 5, 'DisplayName', '27.2 - 27.4')
plot(unique(ts_dates), isopycnal_dist_runmean(4,:)* 5, 'DisplayName', '27.1 - 27.3')
hold off

datetick('x', 'mm/dd/yyyy');
xlim([min(unique(ts_dates)), max(unique(ts_dates))]);
title('Running Mean Distance Between Isopycnals')
ylabel('Isopycnal Separation (m)')
legend();

%%% salinity time series
ax3 = subplot(5,1,3);
time_series_gen(unique(ts_dates), woa_depths(1:41),'pressure', ts_salt_runmean, 'salinity', ax3)
contours(unique(ts_dates), woa_depths(1:41), 'density', ts_salt_runmean, ts_temp_runmean)
title('Running Mean Salinity Time Series')

%%% temperature time series
ax4 = subplot(5,1,4);
time_series_gen(unique(ts_dates), woa_depths(1:41),'pressure', ts_temp_runmean, 'temperature', ax4)
contours(unique(ts_dates), woa_depths(1:41), 'density', ts_salt_runmean, ts_temp_runmean)
title('Running Mean Temperature Time Series')

%%% density time series
ax5 = subplot(5,1,5);
time_series_gen(unique(ts_dates), woa_depths(1:41),'pressure', ts_density_runmean, 'density', ax5)
contours(unique(ts_dates), woa_depths(1:41), 'spice', ts_salt_runmean, ts_temp_runmean)
title('Running Mean Density Time Series')

%%% Realigning subplot 0 to match color-bar axes
pos1 = get(ax1, 'Position');
pos3 = get(ax3, 'Position'); 
pos4 = get(ax4, 'Position');
pos5 = get(ax5, 'Position');
pos3(3) = pos1(3);
pos4(3) = pos1(3);
pos5(3) = pos1(3);
set(ax3, 'Position', pos3);
set(ax4, 'Position', pos4);
set(ax5, 'Position', pos5);

%saveas(gcf,'/Users/jenkosty/Downloads/detectSCV-main/Plots/WOA, Argo, and MEOP Plots/6-25-2021/' + string(ID) + ' (ps, rm)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Stratification Figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1 = figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);

time_series_gen(unique(ts_dates), woa_depths(1:41), 'pressure', N_2_runmean, 'stratification', ax1)
contours(unique(ts_dates), woa_depths(1:41), 'density', ts_salt_runmean, ts_temp_runmean)

%%%%%%%%%%%%%%%%%%%%%
%%%%%% Set Up %%%%%%%
%%%%%%%%%%%%%%%%%%%%%

% %% creating variable with the MEOP tag IDs
% tag = [];
% for i = 1:length(interp_meop_data)
%     tag(i) = interp_meop_data(i).tag;
% end
% 
% %%% identifying the unique MEOP tag IDs and finding their associated
% %%% profiles
% raw_tag = unique(tag);
% u = 1;
% for i = raw_tag
%     tag_elements(u) = {find(tag == i)};
%     u = u + 1;
% end
% 
% 
% %%% identifying tags with > 40 profiles
% u = 1;
% for i = 1:length(tag_elements)
%     if length(tag_elements{1, i}) > 40
%         raw_tag_final(u) = raw_tag(i);
%         tag_elements_final(u) = {tag_elements{1, i}};
%         u = u + 1;
%     end
% end