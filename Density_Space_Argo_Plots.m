%%% input month (1 - 12)
month = input('Enter a month: ');

%%% input latitude (-90 - 90)
lat_input = input('Enter a latitude: ');

%%% input longitude (-180 - 180)
lon_input = input('Enter a longitude: ');

earthRad = 6371; % Earth radius in kilometers
no_profiles = 10; % Number of profiles 

%%% initializing matrices to hold the distance from each Argo profile to
%%% the input corrdinate and the month of each Argo profile
dist_argo = zeros(1,length(ds_interp_argo_data));
month_argo_gen = zeros(1,length(ds_interp_argo_data));

%%% calculating the distance from each Argo dataset to the input
%%% latitude and longitude; identifying the month in which each Argo
%%% dataset was collected
for i = 1:length(ds_interp_argo_data)
    dist_argo(i) = distance(lat_input, lon_input, ds_interp_argo_data(i).lat, ds_interp_argo_data(i).lon, earthRad);
    month_argo_gen(i) = ds_interp_argo_data(i).month;
end

locations_argo = [];

%%% identifying "no_profiles" number of Argo profiles closest to the
%%% input latitude and longitude and that is within the correct month
[locations_argo(:,1)] = find(month_argo_gen == month);
[sorted_dist, sorted_dist_index] = sort(dist_argo(locations_argo));
coordinates_argo = sorted_dist_index(1:no_profiles);

%%% identifying all the profiles within the correct month
[monthly_argo_data] = ds_interp_argo_data(locations_argo);

%%% creating matrices to hold Argo salinity and temperature data
salt_argo = [];
temp_argo = [];
argo_years = [];

for i = 1:length(coordinates_argo)
    %%% creating matrices of Argo salinity and temperature with depth
    %%% for the input month
    salt_argo(:,i) = monthly_argo_data(coordinates_argo(i)).salt;
    temp_argo(:,i) = monthly_argo_data(coordinates_argo(i)).temp;
    argo_years(i) = monthly_argo_data(coordinates_argo(i)).year;
end

%%% creating vectors to hold min/max/mean/iqr data for each depth of the
%%% Argo salinity and temperature data
min_argo_salt = zeros(length(density_grid),1);
max_argo_salt = zeros(length(density_grid),1);
min_argo_temp = zeros(length(density_grid),1);
max_argo_temp = zeros(length(density_grid),1);
mean_argo_salt = zeros(length(density_grid),1);
mean_argo_temp = zeros(length(density_grid),1);
iqr_argo_salt = zeros(length(density_grid),1);
iqr_argo_temp = zeros(length(density_grid),1);

%%% calculating the min/max/mean/iqr for each depth of the Argo salinity
%%% and temperature data
for i = 1:length(density_grid)
    min_argo_salt(i) = min(salt_argo(i,:));
    max_argo_salt(i) = max(salt_argo(i,:));
    min_argo_temp(i) = min(temp_argo(i,:));
    max_argo_temp(i) = max(temp_argo(i,:));
    mean_argo_salt(i) = mean(salt_argo(i,:), 'omitnan');
    mean_argo_temp(i) = mean(temp_argo(i,:), 'omitnan'); 
    iqr_argo_salt(i) = iqr(salt_argo(i,:), 'all');
    iqr_argo_temp(i) = iqr(temp_argo(i,:), 'all');
end

%%% creating anomaly plot
figure('Renderer', 'painters', 'Position', [10 10 1500 800]);
sgtitle({string(no_profiles) + ' Argo Profiles Closest to Lat = ' + string(lat_input) + ' and Lon = ' + string(lon_input) ' for Month ' + string(month) + '; Years: ' + string(min(argo_years)) + '-' + string(max(argo_years))});

%%% salinity vs depths subplot
subplot(1,4,1);
xlabel('Salinity')
ylabel('Density')
grid on, grid minor
hold on
for i = 1:length(coordinates_argo)
    plot(salt_argo(:,i), density_grid)
end
plot(mean_argo_salt,density_grid, 'r', 'LineWidth', 4);
plot(mean_argo_salt + 0.5 * iqr_argo_salt, density_grid, 'b', 'LineWidth', 4);
plot(mean_argo_salt - 0.5 * iqr_argo_salt, density_grid, 'b', 'LineWidth', 4);
set(gca, 'YDir','reverse')
ylim([27.5,27.75])
hold off

%%% salinity anomaly vs depths subplot
subplot(1,4,2);
xlabel('Salinity Anomaly');
ylabel('Density');
grid on, grid minor
hold on
%%% calculating and plotting anomalies and iqr
anom_argo_salt = [];
iqr_argo_salt_anom = [];

for i = 1:length(coordinates_argo)
    anom_argo_salt(:,i) = salt_argo(:,i) - mean_argo_salt;
    plot(anom_argo_salt(:,i), density_grid)
end

for i = 1:length(density_grid)
    iqr_argo_salt_anom(i) = iqr(anom_argo_salt(i,:));
end
plot(0.5 * iqr_argo_salt_anom, density_grid,'b', 'LineWidth',3)
plot(-0.5 * iqr_argo_salt_anom, density_grid,'b', 'LineWidth',3)
set(gca, 'YDir','reverse')
ylim([27.5,27.75])
hold off

%%% temperature vs depths subplot
subplot(1,4,3);
xlabel('Temperature')
ylabel('Density')
grid on, grid minor
hold on
for i = 1:length(coordinates_argo)
    plot(temp_argo(:,i), density_grid)
end
plot(mean_argo_temp, density_grid, 'r', 'LineWidth', 3)
plot(mean_argo_temp + 0.5 * iqr_argo_temp, density_grid, 'b', 'LineWidth', 3)
plot(mean_argo_temp - 0.5 * iqr_argo_temp, density_grid, 'b', 'LineWidth', 3)
set(gca, 'YDir','reverse')
ylim([27.5,27.75])
hold off

%%% temperature anomaly vs depths subplot
subplot(1,4,4);
xlabel('Temperature Anomaly');
ylabel('Density');
grid on, grid minor
hold on

anom_argo_temp = [];
iqr_argo_temp_anom = [];

for i = 1:length(coordinates_argo)
    anom_argo_temp(:,i) = temp_argo(:,i) - mean_argo_temp;
    plot(anom_argo_temp(:,i), density_grid)
end

for i = 1:length(density_grid)
    iqr_argo_temp_anom(i) = iqr(anom_argo_temp(i,:));
end
plot(0.5 * iqr_argo_temp_anom, density_grid, 'b', 'LineWidth',3)
plot(-0.5 * iqr_argo_temp_anom, density_grid, 'b', 'LineWidth',3)
set(gca, 'YDir','reverse')
ylim([27.5,27.75])
hold off
drawnow
saveas(gcf,'/Users/jenkosty/Downloads/detectSCV-main/Plots/WOA, Argo, and MEOP Plots/4-2-2021/' + string(month) + ',' + string(lat_input) + ',' + string(lon_input) + '.fig');

%%% creating T/S plot with isopycnals
figure()
sgtitle({string(no_profiles) + ' Argo Profiles Closest to Lat = ' + string(lat_input) + ' and Lon = ' + string(lon_input) ' for Month ' + string(month) + '; Years: ' + string(min(argo_years)) + '-' + string(max(argo_years))});
set(gcf,'Position',[400 400 700 700]);
xlabel('Salinity')
ylabel('Temperature')
hold on
for i = 1:length(coordinates_argo)
    plot(salt_argo(:,i), temp_argo(:,i))
end

%%% plotting isopycnals
pt_max = max(max(temp_argo)) + 0.25;
pt_min = min(min(temp_argo)) - 0.25;
ss_max = max(max(salt_argo)) + 0.1;
ss_min = min(min(salt_argo)) - 0.1;

%%% Grid for contouring density
pt_step = (pt_max-pt_min)/399;
pt_grid = pt_min:pt_step:pt_max;
ss_step = (ss_max-ss_min)/399;
ss_grid = ss_min:ss_step:ss_max;
[PT_grid,SS_grid] = meshgrid(pt_grid,ss_grid);

%%% calculating density
pd = gsw_rho(SS_grid,PT_grid,density_grid) - 1000;

%%% plotting curves
[C,h] = contour(SS_grid, PT_grid, pd, 15,'EdgeColor','k');
clabel(C,h);
hold off
drawnow
saveas(gcf,'/Users/jenkosty/Downloads/detectSCV-main/Plots/WOA, Argo, and MEOP Plots/4-2-2021/' + string(month) + ',' + string(lat_input) + ',' + string(lon_input) + ' (TS).fig');

%%% creating spice plot
figure()
sgtitle({string(no_profiles) + ' Argo Profiles Closest to Lat = ' + string(lat_input) + ' and Lon = ' + string(lon_input) ' for Month ' + string(month) + '; Years: ' + string(min(argo_years)) + '-' + string(max(argo_years))});
set(gcf,'Position',[400 400 700 700]);
subplot(1,2,1)
spice_argo = [];
hold on
for i = 1:length(coordinates_argo)
    spice_argo(:,i) = gsw_spiciness0(salt_argo(:,i), temp_argo(:,i));
    plot(spice_argo(:,i), density_grid)
end

%%% creating matrices to hold the min/max/mean/iqr at each depth for the
%%% spice data
min_argo_spice = zeros(length(density_grid),1);
max_argo_spice = zeros(length(density_grid),1);
mean_argo_spice = zeros(length(density_grid),1);
iqr_argo_spice = zeros(length(density_grid),1);

%%% calculating the min/max/mean/iqr at each depth for the spice data
for i = 1:length(density_grid)
    min_argo_spice(i) = min(spice_argo(i,:));
    max_argo_spice(i) = max(spice_argo(i,:));
    mean_argo_spice(i) = mean(spice_argo(i,:), 'omitnan');
    iqr_argo_spice(i) = iqr(spice_argo(i,:));
end

plot(mean_argo_spice, density_grid, 'r', 'LineWidth', 3)
plot(mean_argo_spice + 0.5 * iqr_argo_spice, density_grid,'b', 'LineWidth',3)
plot(mean_argo_spice - 0.5 * iqr_argo_spice, density_grid,'b', 'LineWidth',3)
hold off
xlabel('Spice')
ylabel('Density')
set(gca, 'YDir','reverse')
ylim([27.5,27.75])
grid on, grid minor

%%% creating spice anomaly subplot
subplot(1,2,2)
anom_spice_argo = [];
hold on

%%% calculating and plotting the spice anomaly for each profile
for i = 1:length(coordinates_argo)
    anom_spice_argo(:,i) = spice_argo(:,i) - mean_argo_spice;
    plot(anom_spice_argo(:,i), density_grid)
end

iqr_argo_spice_anom = zeros(length(density_grid),1);

%%% calculating the spice anomaly iqr at each depth
for i = 1:length(density_grid)
    iqr_argo_spice_anom(i) = iqr(anom_spice_argo(i,:));
end

plot(0.5 * iqr_argo_spice_anom, density_grid,'b', 'LineWidth',3);
plot(-0.5 * iqr_argo_spice_anom, density_grid,'b', 'LineWidth',3);
hold off
xlabel('Spice Anomaly')
ylabel('Density')
set(gca, 'YDir','reverse')
ylim([27.5,27.75])
grid on, grid minor
drawnow
saveas(gcf,'/Users/jenkosty/Downloads/detectSCV-main/Plots/WOA, Argo, and MEOP Plots/4-2-2021/' + string(month) + ',' + string(lat_input) + ',' + string(lon_input) + ' (Spice).fig');

%%% creating map to display the profiles
figure();
load coastlines
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -60],'FontSize',10)
set(gcf,'Position',[400 400 700 700]);
axis off;
framem on;
gridm on;
mlabel on;
plabel on;
plotm(double([monthly_argo_data.lat]),double([monthly_argo_data.lon]),'m.')
plotm(double([monthly_argo_data(coordinates_argo).lat]),double([monthly_argo_data(coordinates_argo).lon]),'co')
plotm(lat_input, lon_input,'g*')
plotm(coastlat, coastlon, 'k')