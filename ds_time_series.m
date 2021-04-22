%%% creating density space time-series plots and maps for ARGO floats

float_no = 218;

%%% creating a map showing the movement of the float over time
figure();
hold on
load coastlines
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -60],'FontSize',10)
set(gcf,'Position',[400 400 700 700]);
axis off;
framem on;
gridm on;
mlabel on;
plabel on;
plotm(double([ds_interp_argo_data(ds_float_elements_final{1, float_no}).lat]),double([ds_interp_argo_data(ds_float_elements_final{1, float_no}).lon]), '-s');
plotm(double([ds_interp_argo_data(ds_float_elements_final{1, float_no}(1)).lat]),double([ds_interp_argo_data(ds_float_elements_final{1, float_no}(1)).lon]), 'g*');
plotm(double([ds_interp_argo_data(ds_float_elements_final{1, float_no}(length(ds_float_elements_final{1, float_no}))).lat]),double([ds_interp_argo_data(ds_float_elements_final{1, float_no}(length(ds_float_elements_final{1, float_no}))).lon]), 'r*');
plotm(coastlat, coastlon, 'k')
title('Float ' + string(ds_raw_float_final(float_no)))
hold off
%saveas(gcf,'/Users/jenkosty/Downloads/detectSCV-main/Plots/WOA, Argo, and MEOP Plots/4-16-2021/' + string(raw_float(float_no)) + ' (Map)')

%%% creating matrices to hold the time series temperature and salinity data
ts_salt = [];
ts_temp = [];
ts_years = [];
ts_months = [];
ts_density = [];
for i = 1:length(ds_interp_argo_data(ds_float_elements_final{1, float_no}))
    ts_salt(:,i) = ds_interp_argo_data(ds_float_elements_final{1, float_no}(i)).salt;
    ts_temp(:,i) = ds_interp_argo_data(ds_float_elements_final{1, float_no}(i)).temp;
    ts_years(:,i) = ds_interp_argo_data(ds_float_elements_final{1, float_no}(i)).year;
    ts_months(:,i) = ds_interp_argo_data(ds_float_elements_final{1, float_no}(i)).month;
    ts_density(:,i) = isnan(ds_interp_argo_data(ds_float_elements_final{1, float_no}(i)).salt);
end

%%% creating range for density
density_elements = [];
density_elements = find(all(ts_density, 2) == 0);

density_range = [];
density_range = [density_grid(min(density_elements)), density_grid(max(density_elements))];

%%% finding xtick locations (first occurance of each year)
u = 1;
ts_xticks = [];
for i = unique(ts_years)
    ts_xticks(u) = find(ts_years == i, 1);
    u = u + 1;
end    

%%% creating salinity time series
figure();
hold on
p = pcolor(double([ds_interp_argo_data(ds_float_elements_final{1, float_no}).cycle]), density_grid, ts_salt);
colorbar;
set(p, 'EdgeColor', 'none')
set(gca, 'YDir','reverse');
ylim([27.4,27.7]);
xticks(double([ds_interp_argo_data(ds_float_elements_final{1, float_no}(ts_xticks)).cycle]))
xticklabels(string(ts_years(ts_xticks)));
title('Salinity Time Series; Float ' + string(ds_raw_float_final(float_no)))
ylabel('Density')
xlabel('Year')
%saveas(gcf,'/Users/jenkosty/Downloads/detectSCV-main/Plots/WOA, Argo, and MEOP Plots/4-16-2021/' + string(raw_float(float_no)) + ' (Salt)');

%%% creating temperature time series
figure();
p = pcolor(double([ds_interp_argo_data(ds_float_elements_final{1, float_no}).cycle]), density_grid, ts_temp);
colorbar;
set(p, 'EdgeColor', 'none')
set(gca, 'YDir','reverse');
ylim(density_range)
xticks(double([ds_interp_argo_data(ds_float_elements_final{1, float_no}(ts_xticks)).cycle]));
xticklabels(string(ts_years(ts_xticks)));
title('Temperature Time Series; Float ' + string(ds_raw_float_final(float_no)));
ylabel('Density')
xlabel('Years')
%saveas(gcf,'/Users/jenkosty/Downloads/detectSCV-main/Plots/WOA, Argo, and MEOP Plots/4-16-2021/' + string(raw_float(float_no)) + ' (Temp)')

%%%%%%%%%%%%%%%%%%%%%
%%%%%% Set Up %%%%%%%
%%%%%%%%%%%%%%%%%%%%%
% 
% %%% creating variable with the ARGO float IDs
% ds_float = [];
% for i = 1:length(ds_interp_argo_data)
%     ds_float(i) = str2double(ds_interp_argo_data(i).float);
% end
% 
% %%% identifying the unique ARGO float IDs and finding their associated
% %%% profiles
% ds_raw_float = unique(ds_float);
% u = 1;
% for i = ds_raw_float
%     ds_float_elements(u) = {find(ds_float == i)};
%     u = u + 1;
% end
% 
% %%% identifying floats with > 19 profiles
% u = 1;
% for i = 1:length(ds_float_elements)
%     if length(ds_float_elements{1, i}) > 19
%         ds_raw_float_final(u) = ds_raw_float(i);
%         ds_float_elements_final(u) = {ds_float_elements{1, i}};
%         u = u + 1;
%     end
% end

