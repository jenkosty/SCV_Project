%%% creating temperature, salinity plots for Southern Ocean MEOP, Argo, and WOA data

%%% input month (1 - 12)
month = input('Enter a month: ');

%%% input latitude (-90 - 90)
lat_input = input('Enter a latitude: ');

%%% input longitude (-180 - 180)
lon_input = input('Enter a longitude: ');

%%% finding Southern Ocean (under 60S) WOA Data
woa_lats_SO = woa_lats(find(woa_lats < -60));

earthRad = 6371; % Earth radius in kilometers
radius_dist = 50; % Search radius in kilometers

%%% initializing WOA distance matrix
dist_woa = zeros(length(woa_lats_SO),length(woa_lons));

%%% calculating the distance from each WOA entry to the input latitude and
%%% longitude
for i = 1:length(woa_lats_SO)
    for f = 1:length(woa_lons)
        dist_woa(i,f) = distance(lat_input,lon_input,woa_lats_SO(i),woa_lons(f),earthRad);
    end
end

%%% initializing month and distance matrices for Argo data
dist_argo = zeros(1,length(SO_argo));
month_argo_gen = zeros(1,length(SO_argo));

%%% calculating the distance from each Argo dataset to the input
%%% latitude and longitude ; indentifying the month in which each Argo
%%% daraset was connected
for i = 1:length(SO_argo)
    dist_argo(i) = distance(lat_input, lon_input, SO_argo(i).lat, SO_argo(i).lon, earthRad);
    month_argo_gen(i) = SO_argo(i).time(2);
end

%%% initializing month and distance matrices for MEOP data
dist_meop = zeros(1,length(SO_sealdata_qc));
month_meop_gen = zeros(1,length(SO_sealdata_qc));

%%% calculating the distance from each MEOP dataset to the input
%%% latitude and longitude; identifying the month in which each MEOP
%%% dataset was collected
for i = 1:length(SO_sealdata_qc)
    dist_meop(i) = distance(lat_input,lon_input,SO_sealdata_qc(i).LAT,SO_sealdata_qc(i).LON,earthRad);
    month_meop_gen(i) = SO_sealdata_qc(i).TIME(2);
end

%%% initializing coordinate matrices for the WOA, MEOP, and Argo data
coordinates_woa = [];
coordinates_meop = [];
locations_meop = [];
coordinates_argo = [];
locations_argo = [];

%%% identifying the WOA, MEOP, and Argo data that is within "radius_dist" kilometers of the
%%% input latitude and longitude (and is within the correct month for the
%%% MEOP/Argo data)
[coordinates_woa(:,1), coordinates_woa(:,2)] = find(dist_woa < radius_dist);
[coordinates_meop(:,1)] = find(dist_meop < radius_dist & month_meop_gen == month);
[locations_meop(:,1)] = find(month_meop_gen == month);
[coordinates_argo(:,1)] = find(dist_argo < radius_dist & month_argo_gen == month);
[locations_argo(:,1)] = find(month_argo_gen == month);

%%% creating figure to display 3 subplots
figure('Renderer', 'painters', 'Position', [10 10 1300 700]);

%%% creating subplot 1: temperature vs salinity
subplot(1,13,[1,6]);
xlabel('Salinity');
ylabel('Temperature')
hold on

%%% initializing matrices for the WOA salinity and temperature matrices
salt_woa = zeros(length(coordinates_woa),length(woa_depths));
temp_woa = zeros(length(coordinates_woa),length(woa_depths));

%%% creating WOA salinity and temperature matrices and plotting WOA
%%% salinity vs temperature
for i = 1:length(coordinates_woa)
    %%% creating matrices of WOA salinity and temperature with depth
    %%% for the input month
    salt_woa(i,:) = squeeze(woa_salt(coordinates_woa(i,2), coordinates_woa(i,1),:, month));
    temp_woa(i,:) = squeeze(woa_temp(coordinates_woa(i,2), coordinates_woa(i,1),:, month));
    
    %%% plotting WOA temperature vs salinity in yellow
    plot(salt_woa(i,1:41),temp_woa(i,1:41),'y');
end

%%% plotting MEOP temperature and salinity in blue
for i = 1:length(coordinates_meop)
    idz = find(SO_sealdata_qc(coordinates_meop(i,1)).PRES <= 700);
    plot(SO_sealdata_qc(coordinates_meop(i,1)).SALT(idz), SO_sealdata_qc(coordinates_meop(i,1)).TEMP(idz),'b'); 
end

%%% plotting Argo temperature and salinity in magenta
for i = 1:length(coordinates_argo)
    idx = find(SO_argo(coordinates_argo(i,1)).pres <= 700);
    plot(SO_argo(coordinates_argo(i,1)).salt(idx), SO_argo(coordinates_argo(i,1)).temp(idx),'m');
end

%%% plotting isopycnals
pt_max = max(max(temp_woa)) + 1;
pt_min = min(min(temp_woa)) - 1;
ss_max = max(max(salt_woa)) + 0.1;
ss_min = min(min(salt_woa)) - 0.1;

%%% Grid for contouring density
pt_step = (pt_max-pt_min)/56;
pt_grid = pt_min:pt_step:pt_max;
ss_step = (ss_max-ss_min)/56;
ss_grid = ss_min:ss_step:ss_max;
[PT_grid,SS_grid] = meshgrid(pt_grid,ss_grid);

%%% calculating density
pd = gsw_rho(SS_grid,PT_grid,woa_depths) - 1000;

%%% plotting curves
[C,h] = contour(SS_grid, PT_grid, pd, 15,'EdgeColor','k');
clabel(C,h);
hold off

%%% creating subplot 2: salinity vs depths
subplot(1,14,[8,10]);
xlabel('Salinity');
ylabel('Depths');
hold on
%%% plotting WOA salinity vs depths in yellow
for i = 1:length(coordinates_woa)
    plot(salt_woa, woa_depths,'y');
end
    
%%% plotting MEOP salinity vs depths in blue
for i = 1:length(coordinates_meop)
    plot(SO_sealdata_qc(coordinates_meop(i,1)).SALT, SO_sealdata_qc(coordinates_meop(i,1)).PRES, 'b');
end

%%% plotting Argo salinity vs depths in magenta
for i = 1:length(coordinates_argo)
    plot(SO_argo(coordinates_argo(i,1)).salt, SO_argo(coordinates_argo(i,1)).pres,'m');
end
set(gca, 'YDir','reverse')
ylim([0,700])
grid on, grid minor
hold off

%%% creating subplot 3: temperature vs depths
subplot(1,14,[12,14]);
xlabel('Temperature');
ylabel('Depths');
hold on
%%% plotting WOA temperature vs depths in yellow
for i = 1:length(coordinates_woa)
    plot(temp_woa, woa_depths, 'y');
end

%%% plotting MEOP temperature vs depths in blue
for i = 1:length(coordinates_meop)
    plot(SO_sealdata_qc(coordinates_meop(i,1)).TEMP, SO_sealdata_qc(coordinates_meop(i,1)).PRES, 'b');
end

%%% plotting Argo temperature vs depths in magenta
for i = 1:length(coordinates_argo)
    plot(SO_argo(coordinates_argo(i,1)).temp, SO_argo(coordinates_argo(i,1)).pres,'m');
end
set(gca, 'YDir','reverse')
ylim([0,700])
grid on, grid minor
hold off

%%% calculating the years for the Argo and Meop profiles
meop_years = [];
argo_years = [];
for i = 1:length(coordinates_meop)
    meop_years(i) = SO_sealdata_qc(coordinates_meop(i)).TIME(1);
end

if isempty(coordinates_meop)
    meop_years = [0,0,0];
end

for i = 1:length(coordinates_argo)
    argo_years(i) = SO_argo(coordinates_argo(i)).time(1);
end

if isempty(coordinates_argo)
    argo_years = [0,0,0];
end

%%% adding plot title
sgtitle({'WOA (yellow), MEOP (blue), and Argo (magenta) Profiles Within ' + string(radius_dist) + ' km of Latitude = ' + string(lat_input) + ' and Longitude = ' + string(lon_input) + ' for Month #' + string(month); 'Number of MEOP Profiles = ' + string(length(coordinates_meop)) + '; Number of Argo Profiles = ' + string(length(coordinates_argo)); 'MEOP: ' + string(min(meop_years)) + '-' + string(max(meop_years)) + '; Argo: ' + string(min(argo_years)) + '-' + string(max(argo_years))});
drawnow

%%% creating map of MEOP and Argo profiles for the input month
figure();
load coastlines
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -60],'FontSize',10)
set(gcf,'Position',[400 400 700 700]);
axis off;
framem on;
gridm on;
mlabel on;
plabel on;
plotm([SO_sealdata_qc(locations_meop(:,1)).LAT],[SO_sealdata_qc(locations_meop(:,1)).LON],'b.')
plotm(double([SO_argo(locations_argo(:)).lat]),double([SO_argo(locations_argo(:)).lon]),'m.')
plotm([SO_sealdata_qc(coordinates_meop(:,1)).LAT],[SO_sealdata_qc(coordinates_meop(:,1)).LON],'ro')
plotm(double([SO_argo(coordinates_argo(:,1)).lat]),double([SO_argo(coordinates_argo(:,1)).lon]),'co')
plotm(lat_input, lon_input,'g*')
plotm(coastlat, coastlon, 'k')
title({'MEOP (blue) and Argo (magenta) Profiles for Month ' + string(month); 'Number of MEOP Profiles = ' + string(length(locations_meop)) + '; Number of Argo Profiles = ' + string(length(locations_argo))}, 'FontSize',14);

%%% creating figure to display anomalies
figure('Renderer', 'painters', 'Position', [10 10 900 600]);
sgtitle({' MEOP (blue), and Argo (magenta) Anomalies Within ' + string(radius_dist) + ' km of Latitude = ' + string(lat_input) + ' and Longitude = ' + string(lon_input) + ' for Month #' + string(month); 'Number of MEOP Profiles = ' + string(length(coordinates_meop)) + '; Number of Argo Profiles = ' + string(length(coordinates_argo)); 'MEOP: ' + string(min(meop_years)) + '-' + string(max(meop_years)) + '; Argo: ' + string(min(argo_years)) + '-' + string(max(argo_years))});

min_woa_salt = zeros(length(woa_depths));
max_woa_salt = zeros(length(woa_depths));
min_woa_temp = zeros(length(woa_depths));
max_woa_temp = zeros(length(woa_depths));
mean_woa_salt = zeros(length(woa_depths));
mean_woa_temp = zeros(length(woa_depths));

for i = 1:length(woa_depths)
    min_woa_salt(i) = min(salt_woa(:,i));
    max_woa_salt(i) = max(salt_woa(:,i));
    min_woa_temp(i) = min(temp_woa(:,i));
    max_woa_temp(i) = max(temp_woa(:,i));
    mean_woa_salt(i) = mean(salt_woa(:,i), 'omitnan');
    mean_woa_temp(i) = mean(temp_woa(:,i), 'omitnan');  
end

interpolated_salt_meop = [];
interpolated_salt_argo = [];
anom_salt_meop = [];
anom_salt_argo = [];

%%% creating salinity anomaly subplot
subplot(1,7,[1,3]);
xlabel('Salinity Anomaly');
ylabel('Depths');

hold on
%%% calculating and plotting MEOP salinity anomalies (from mean WOA
%%% salinity)
for i = 1:length(coordinates_meop)
    interpolated_salt_meop{i} = interp1(woa_depths, mean_woa_salt(:,1), SO_sealdata_qc(coordinates_meop(i,:)).PRES);
    anom_salt_meop{i} = SO_sealdata_qc(coordinates_meop(i)).SALT - interpolated_salt_meop{i};
    %plot(interpolated_salt_meop{i},SO_sealdata_qc(coordinates_meop(i,:)).PRES, 'm');
    plot(anom_salt_meop{i},SO_sealdata_qc(coordinates_meop(i,:)).PRES, 'b')
    %plot(salt_woa(i,:),woa_depths, 'r')
    %plot(mean_woa_salt(:,1), woa_depths, 'g')
    %plot(SO_sealdata_qc(coordinates_meop(i,:)).SALT, SO_sealdata_qc(coordinates_meop(i,:)).PRES, 'b')
end

%%% calculating and plotting Argo salinity anomalies (from mean WOA
%%% salinity)
for i = 1:length(coordinates_argo)
    interpolated_salt_argo{i} = interp1(woa_depths, mean_woa_salt(:,1),SO_argo(coordinates_argo(i,:)).pres);
    anom_salt_argo{i} = SO_argo(coordinates_argo(i,:)).salt - interpolated_salt_argo{1,i};
    plot(anom_salt_argo{i}, SO_argo(coordinates_argo(i,:)).pres, 'm')
    %plot(SO_argo(coordinates_argo(i,:)).salt, SO_argo(coordinates_argo(i,:)).pres, 'y')
end
set(gca, 'YDir','reverse')
ylim([0,700])
grid on, grid minor
hold off

interpolated_temp_meop = [];
interpolated_temp_argo = [];
anom_temp_meop = [];
anom_temp_argo = [];

%%% creating temperature anomaly subplot
subplot(1,7,[5,7]);
xlabel('Temperature Anomaly');
ylabel('Depths');

hold on
%%% calculating and plotting MEOP temperature anomalies (from mean WOA
%%% temperature)
for i = 1:length(coordinates_meop)
    interpolated_temp_meop{i} = interp1(woa_depths, mean_woa_temp(:,1),SO_sealdata_qc(coordinates_meop(i,:)).PRES);
    anom_temp_meop{i} = SO_sealdata_qc(coordinates_meop(i,:)).TEMP - interpolated_temp_meop{1,i};
    plot(anom_temp_meop{i},SO_sealdata_qc(coordinates_meop(i,:)).PRES, 'b')
    %plot(mean_woa_temp(:,1),woa_depths, 'r')
    %plot(SO_sealdata_qc(coordinates_meop(i,:)).TEMP, SO_sealdata_qc(coordinates_meop(i,:)).PRES, 'g')
end

%%% calculating and plotting Argo temperature anomalies (from mean WOA
%%% temperature)
for i = 1:length(coordinates_argo)
    interpolated_temp_argo{i} = interp1(woa_depths, mean_woa_temp(:,1),SO_argo(coordinates_argo(i,:)).pres);
    anom_temp_argo{i} = SO_argo(coordinates_argo(i,:)).temp - interpolated_temp_argo{1,i};
    plot(anom_temp_argo{i}, SO_argo(coordinates_argo(i,:)).pres, 'm')
    %plot(SO_argo(coordinates_argo(i,:)).temp, SO_argo(coordinates_argo(i,:)).pres, 'y')
end
set(gca, 'YDir','reverse')
grid on, grid minor
hold off