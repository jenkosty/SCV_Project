%%% creating time-series plots and maps for ARGO floats

%%% creating variable with the ARGO float IDs
float = [];
for i = 1:length(interp_argo_data)
    float(i) = str2double(interp_argo_data(i).float);
end

%%% identifying the unique ARGO float IDs and finding their associated
%%% profiles
raw_float = unique(float);
u = 1;
for i = raw_float
    float_elements(u) = {find(float == i)};
    u = u + 1;
end

for i = 154
float_no = i;

%%% creating a map showing the movement of the float over time
figure();
hold on
for i = float_no
    load coastlines
    axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -60],'FontSize',10)
    set(gcf,'Position',[400 400 700 700]);
    axis off;
    framem on;
    gridm on;
    mlabel on;
    plabel on;
    plotm(double([interp_argo_data(float_elements{1, i}).lat]),double([interp_argo_data(float_elements{1, i}).lon]), '-s');
    plotm(coastlat, coastlon, 'k')
    title('Float ' + string(raw_float(float_no)))
    %saveas(gcf,'/Users/jenkosty/Downloads/detectSCV-main/Plots/WOA, Argo, and MEOP Plots/4-16-2021/' + string(raw_float(float_no)) + ' (Map)')
end
hold off
end

%%% creating matrices to hold the time series temperature and salinity data
ts_salt = [];
ts_temp = [];
ts_years = [];
ts_months = [];
for i = 1:length(interp_argo_data(float_elements{1, float_no}))
    ts_salt(:,i) = interp_argo_data(float_elements{1, float_no}(i)).salt;
    ts_temp(:,i) = interp_argo_data(float_elements{1, float_no}(i)).temp;
    ts_years(:,i) = interp_argo_data(float_elements{1, float_no}(i)).year;
    ts_months(:,i) = interp_argo_data(float_elements{1, float_no}(i)).month;
end

%%% finding xtick locations (first occurance of each year)
u = 1;
ts_xticks = [];
for i = unique(ts_years)
    ts_xticks(u) = find(ts_years == i, 1);
    u = u + 1;
end

    

%%% creating salinity time series
figure();
p = pcolor(double([interp_argo_data(float_elements{1, float_no}).cycle]), woa_depths(1:41), ts_salt);
colorbar;
set(p, 'EdgeColor', 'none')
set(gca, 'YDir','reverse');
ylim([0,500]);
xticks(double([interp_argo_data(float_elements{1, float_no}(ts_xticks)).cycle]))
xticklabels(string(ts_years(ts_xticks)));
title('Salinity Time Series; Float ' + string(raw_float(float_no)))
ylabel('Depths')
xlabel('Year')
%saveas(gcf,'/Users/jenkosty/Downloads/detectSCV-main/Plots/WOA, Argo, and MEOP Plots/4-16-2021/' + string(raw_float(float_no)) + ' (Salt)');

%%% creating temperature time series
figure();
p = pcolor(double([interp_argo_data(float_elements{1, float_no}).cycle]), woa_depths(1:41), ts_temp);
colorbar;
set(p, 'EdgeColor', 'none')
set(gca, 'YDir','reverse');
ylim([0,500]);
xticks(double([interp_argo_data(float_elements{1, float_no}(ts_xticks)).cycle]));
xticklabels(string(ts_years(ts_xticks)));
title('Temperature Time Series; Float ' + string(raw_float(float_no)));
ylabel('Depths')
xlabel('Years')
%saveas(gcf,'/Users/jenkosty/Downloads/detectSCV-main/Plots/WOA, Argo, and MEOP Plots/4-16-2021/' + string(raw_float(float_no)) + ' (Temp)')