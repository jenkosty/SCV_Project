%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Density-Space MEOP Time Series %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tag_no = 103;
day_average = 1;
runmean_window_length = 5;

lat = double([ds_interp_meop_data(tag_elements_final{1, tag_no}).lat]);
lon = double([ds_interp_meop_data(tag_elements_final{1, tag_no}).lon]);
ID = raw_tag_final(tag_no);

%%% creating a map showing the movement of the tag over time
map(lat, lon, ID);
%saveas(gcf,'/Users/jenkosty/Downloads/detectSCV-main/Plots/WOA, Argo, and MEOP Plots/6-4-2021/' + string(raw_tag_final(tag_no)) + ' (Map)')

ts_time = [];
for i = 1:length(ds_interp_meop_data(tag_elements_final{1, tag_no}))
    ts_time(i,:) = ds_interp_meop_data(tag_elements_final{1, tag_no}(i)).time;
end

ts_dates = [];
for i = 1:length(tag_elements_final{1, tag_no})
    ts_dates(i) = datenum([ts_time(i,1),ts_time(i,2),ts_time(i,3),ts_time(i,4),ts_time(i,5),ts_time(i,6)]);
end

%%% creating matrices to hold the time series temperature and salinity data
u = 1;
ts_salt = [];
ts_temp = [];
ts_years = [];
ts_months = [];
ts_density = [];
for i = unique(ts_dates)
    l = find(ts_dates == i, 1);
    ts_salt(:,u) = ds_interp_meop_data(tag_elements_final{1, tag_no}(l)).salt;
    ts_temp(:,u) = ds_interp_meop_data(tag_elements_final{1, tag_no}(l)).temp;
    ts_years(:,u) = ds_interp_meop_data(tag_elements_final{1, tag_no}(l)).year;
    ts_months(:,u) = ds_interp_meop_data(tag_elements_final{1, tag_no}(l)).month;
    ts_density(:,u) = isnan(ds_interp_meop_data(tag_elements_final{1, tag_no}(l)).salt);
    u = u + 1;
end  

ts_salt_nu = [];
ts_temp_nu = [];
ts_years_nu = [];
ts_months_nu = [];

for i = 1:length(ts_dates)
    ts_salt_nu(:,i) = ds_interp_meop_data(tag_elements_final{1, tag_no}(i)).salt;
    ts_temp_nu(:,i) = ds_interp_meop_data(tag_elements_final{1, tag_no}(i)).temp;
    ts_years_nu(:,i) = ds_interp_meop_data(tag_elements_final{1, tag_no}(i)).year;
    ts_months_nu(:,i) = ds_interp_meop_data(tag_elements_final{1, tag_no}(i)).month;
end  

%%% creating range for density
density_elements = [];
density_elements = find(all(ts_density, 2) == 0);

density_range = [];
density_range = [density_grid(min(density_elements)), density_grid(max(density_elements))]; 

%%% Computing ranges to 2-day average
i = ts_dates(1);
j = 1;
ts_dates_tda_axis = [];
ts_dates_tda(:) = [];
while i < ts_dates(end)
ts_dates_tda(j) = {find(ts_dates >= i & ts_dates < i+day_average)};
ts_dates_tda_axis(j) = i;
i = i + day_average;
j = j + 1;
end

%%% Computing temperature, salinity, and density data for 2-day average
ts_salt_tda = [];
ts_temp_tda = [];

for j = 1:length(ts_dates_tda)
    for i = 1:length(density_grid)
        ts_salt_tda(i,j) = mean(ts_salt_nu(i,ts_dates_tda{j}), 'omitnan');
        ts_temp_tda(i,j) = mean(ts_temp_nu(i,ts_dates_tda{j}), 'omitnan');
    end
end

%%% Computing temperature, salinity, and density data for running mean
ts_salt_runmean = [];
ts_temp_runmean = [];
ts_density_runmean = [];
for i = 1:length(density_grid)
    ts_salt_runmean(i,:) = movmean(ts_salt(i,:), runmean_window_length, 'omitnan');
    ts_temp_runmean(i,:) = movmean(ts_temp(i,:), runmean_window_length, 'omitnan');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% creating salinity time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);
ax1 = subplot(2,1,1);
hold on
cmap = cmocean('haline');
p = pcolor(unique(ts_dates), density_grid, ts_salt);
colormap(ax1, cmap);
colorbar;
caxis([min(min(ts_salt)),max(max(ts_salt))]);

pd = [];

%%% calculating spice
for i = 1:length(density_grid)
    pd(i,:) = gsw_spiciness0(ts_salt(i,:), ts_temp(i,:));
end

%%% plotting curves
[C,h] = contour(unique(ts_dates),density_grid, pd, [0,-.1,-.2,-.3,-.4,-.5,-.6,-.7,-.8,-.9,-1,-1.1,-1.2,-1.3,-1.4],'EdgeColor','k');
clabel(C,h,'LabelSpacing',250);
hold off

set(p, 'EdgeColor', 'none')
set(gca, 'YDir','reverse');
datetick('x', 'mm/dd/yyyy');
ylim(density_range);
xlim([min(ts_dates), max(ts_dates)]);
title('Salinity Time Series; Tag ' + string(raw_tag_final(tag_no)))
ylabel('Density')
xlabel('Date')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% creating temperature time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax2 = subplot(2,1,2);
cmap = cmocean('thermal');
p = pcolor(ax2, unique(ts_dates), density_grid, ts_temp);
colormap(ax2, cmap);
colorbar;
caxis([min(min(ts_temp)),max(max(ts_temp))]);

pd = [];

%%% calculating spice
for i = 1:length(density_grid)
    pd(i,:) = gsw_spiciness0(ts_salt(i,:), ts_temp(i,:));
end

%%% plotting curves
hold on
[C,h] = contour(unique(ts_dates), density_grid, pd,[0,-.1,-.2,-.3,-.4,-.5,-.6,-.7,-.8,-.9,-1,-1.1,-1.2,-1.3,-1.4],'EdgeColor','k');
clabel(C,h, 'LabelSpacing', 250);
hold off

set(p, 'EdgeColor', 'none')
set(gca, 'YDir','reverse');
datetick('x', 'mm/dd/yyyy');
ylim(density_range);
xlim([min(ts_dates), max(ts_dates)]);
title('Temperature Time Series; Tag ' + string(raw_tag_final(tag_no)));
ylabel('Density')
xlabel('Date')
%saveas(gcf,'/Users/jenkosty/Downloads/detectSCV-main/Plots/WOA, Argo, and MEOP Plots/6-4-2021/' + string(raw_tag_final(tag_no)) + ' (ds)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% creating running mean salinity time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);
ax1 = subplot(2,1,1);
hold on
cmap = cmocean('haline');
p = pcolor(unique(ts_dates), density_grid, ts_salt_runmean);
colormap(ax1, cmap);
colorbar;
caxis([min(min(ts_salt_runmean)),max(max(ts_salt_runmean))]);

pd = [];

%%% calculating spice
for i = 1:length(density_grid)
    pd(i,:) = gsw_spiciness0(ts_salt_runmean(i,:), ts_temp_runmean(i,:));
end

%%% plotting curves
[C,h] = contour(unique(ts_dates), density_grid, pd, [0,-.1,-.2,-.3,-.4,-.5,-.6,-.7,-.8,-.9,-1,-1.1,-1.2,-1.3,-1.4],'EdgeColor','k');
clabel(C,h,'LabelSpacing',500);
hold off

set(p, 'EdgeColor', 'none')
set(gca, 'YDir','reverse');
datetick('x', 'mm/dd/yyyy');
ylim(density_range);
xlim([min(ts_dates), max(ts_dates)]);
title('Running Mean Salinity Time Series; Tag ' + string(raw_tag_final(tag_no)))
ylabel('Density')
xlabel('Time')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% creating running mean temperature time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax2 = subplot(2,1,2);
cmap = cmocean('thermal');
p = pcolor(unique(ts_dates), density_grid, ts_temp_runmean);
colormap(ax2, cmap);
colorbar;
caxis([min(min(ts_temp_runmean)),max(max(ts_temp_runmean))]);

pd = [];

%%% calculating spice
for i = 1:length(density_grid)
    pd(i,:) = gsw_spiciness0(ts_salt_runmean(i,:), ts_temp_runmean(i,:));
end

hold on
%%% plotting curves
[C,h] = contour(unique(ts_dates), density_grid, pd,  [0,-.1,-.2,-.3,-.4,-.5,-.6,-.7,-.8,-.9,-1,-1.1,-1.2,-1.3,-1.4],'EdgeColor','k');
clabel(C,h,'LabelSpacing',500);
hold off

set(p, 'EdgeColor', 'none')
set(gca, 'YDir','reverse');
datetick('x', 'mm/dd/yyyy');
ylim(density_range);
xlim([min(ts_dates), max(ts_dates)]);
title('Running Mean Temperature Time Series; Tag ' + string(raw_tag_final(tag_no)));
ylabel('Density')
xlabel('Time')
%saveas(gcf,'/Users/jenkosty/Downloads/detectSCV-main/Plots/WOA, Argo, and MEOP Plots/6-4-2021/' + string(raw_tag_final(tag_no)) + ' (ds, rm)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% creating 2-day average salinity time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);
ax1 = subplot(2,1,1);
hold on
cmap = cmocean('haline');
p = pcolor(unique(ts_dates_tda_axis), density_grid, ts_salt_tda);
colormap(ax1, cmap);
colorbar;
caxis([min(min(ts_salt_tda)),max(max(ts_salt_tda))]);

pd = [];

%%% calculating spice
for i = 1:length(density_grid)
    pd(i,:) = gsw_spiciness0(ts_salt_tda(i,:), ts_temp_tda(i,:));
end

%%% plotting curves
[C,h] = contour(unique(ts_dates_tda_axis),density_grid, pd, [0,-.1,-.2,-.3,-.4,-.5,-.6,-.7,-.8,-.9,-1,-1.1,-1.2,-1.3,-1.4],'EdgeColor','k');
clabel(C,h,'LabelSpacing',250);
hold off

set(p, 'EdgeColor', 'none')
set(gca, 'YDir','reverse');
datetick('x', 'mm/dd/yyyy');
ylim(density_range);
xlim([min(ts_dates_tda_axis), max(ts_dates_tda_axis)]);
title(string(day_average) + '-Day Average Salinity Time Series; Tag ' + string(raw_tag_final(tag_no)))
ylabel('Density')
xlabel('Date')

%%% creating 2-day average temperature time series
ax2 = subplot(2,1,2);
cmap = cmocean('thermal');
p = pcolor(unique(ts_dates_tda_axis), density_grid, ts_temp_tda);
colormap(ax2, cmap);
colorbar;
caxis([min(min(ts_temp_tda)),max(max(ts_temp_tda))]);

pd = [];

%%% calculating spice
for i = 1:length(density_grid)
    pd(i,:) = gsw_spiciness0(ts_salt_tda(i,:), ts_temp_tda(i,:));
end

%%% plotting curves
hold on
[C,h] = contour(unique(ts_dates_tda_axis), density_grid, pd,[0,-.1,-.2,-.3,-.4,-.5,-.6,-.7,-.8,-.9,-1,-1.1,-1.2,-1.3,-1.4],'EdgeColor','k');
clabel(C,h, 'LabelSpacing', 250);
hold off

set(p, 'EdgeColor', 'none')
set(gca, 'YDir','reverse');
datetick('x', 'mm/dd/yyyy');
ylim(density_range);
xlim([min(ts_dates_tda_axis), max(ts_dates_tda_axis)]);
title(string(day_average) + '-Day Average Temperature Time Series; Tag ' + string(raw_tag_final(tag_no)));
ylabel('Density')
xlabel('Date')
%saveas(gcf,'/Users/jenkosty/Downloads/detectSCV-main/Plots/WOA, Argo, and MEOP Plots/6-4-2021/' + string(raw_tag_final(tag_no)) + ' (ds, ' + string(day_average) + '-day average)')

%%%%%%%%%%%%%%%%%%%%%
%%%%%% Set Up %%%%%%%
%%%%%%%%%%%%%%%%%%%%%

% %%% creating variable with the MEOP tag IDs
% tag = [];
% for i = 1:length(ds_interp_meop_data)
%     tag(i) = ds_interp_meop_data(i).tag;
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
% 
% function datatip_z2cdata(h)
%     % h is a graphics object with a default X/Y/Z datatip and a 
%     % CData property (e.g. pcolor).  The Z portion of the datatip
%     % will be relaced with CData
%     % Generate an invisible datatip to ensure that DataTipTemplate is generated
%     dt = datatip(h,h.XData(1),h.YData(1),'Visible','off'); 
%     % Replace datatip row labeled Z with CData
%     idx = strcmpi('Z',{h.DataTipTemplate.DataTipRows.Label});
%     newrow = dataTipTextRow('C',h.CData);
%     h.DataTipTemplate.DataTipRows(idx) = newrow;
%     % Remove invisible datatip
%     delete(dt)
% end