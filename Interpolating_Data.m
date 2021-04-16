% for i = 1:length(SO_sealdata_qc)
%     interp_meop_salt(i,:) = interp1(SO_sealdata_qc(i).PRES, SO_sealdata_qc(i).SALT, woa_depths(1:41));
% end

% for i = 1:length(SO_sealdata_qc)
%     interp_meop_temp(i,:) = interp1(SO_sealdata_qc(i).PRES, SO_sealdata_qc(i).TEMP, woa_depths(1:41));
% end

% u = 1;
% for i = 1:length(SO_argo)
%     [a,b] = unique(SO_argo(i).pres);
%     if length(b) > 2
%         interp_argo_temp(u,:) = interp1(SO_argo(i).pres(b), SO_argo(i).temp(b), woa_depths(1:41));
%         u = u + 1;
%     end
% end

% u = 1;
% for i = 1:length(SO_argo)
%     [a,b] = unique(SO_argo(i).pres);
%     if length(b) > 2
%         interp_argo_temp(u,:) = interp1(SO_argo(i).pres(b), SO_argo(i).temp(b), woa_depths(1:41));
%         interp_argo_salt(u,:) = interp1(SO_argo(i).pres(b), SO_argo(i).salt(b), woa_depths(1:41));
%         argo_lats_lons(u,:) = [SO_argo(i).lat, SO_argo(i).lon];
%         argo_months(u) = SO_argo(i).time(2);
%         argo_years(u) = SO_argo(i).time(1);
%         argo_float(u) = SO_argo(i).float;
%         argo_cycle(u) = SO_argo(i).cycle;
%         interp_argo_data(u) = struct('float', SO_argo(i).float, 'cycle', SO_argo(i).cycle, 'lat',SO_argo(i).lat,'lon',SO_argo(i).lon,'month',SO_argo(i).time(2),'year',SO_argo(i).time(1), 'salt', interp1(SO_argo(i).pres(b), SO_argo(i).salt(b), woa_depths(1:41)), 'temp', interp1(SO_argo(i).pres(b), SO_argo(i).temp(b), woa_depths(1:41)));
%         u = u + 1;
%     end
% end
% 
% for i = 1:length(argo_lats_lons)
%     interp_argo_data(i) = struct('float', argo_float, 'cycle', argo_cycle, 'lat',argo_lats_lons(i,1),'lon',argo_lats_lons(i,2),'month',argo_months(i),'year',argo_years(i), 'salt', interp_argo_salt(i,:), 'temp', interp_argo_temp(i,:));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolating into Density Space %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
salt_argo = [];
temp_argo = [];
max_argo_density = [];
min_argo_density = [];
u = 1;
for i = 1:length(SO_argo)
    if isempty(SO_argo(i).salt) == 0
        argo_data_w_den(u) = struct('float', SO_argo(i).float,'cycle', SO_argo(i).cycle,'lat', SO_argo(i).lat, 'lon', SO_argo(i).lon,'den',gsw_rho(SO_argo(i).salt,SO_argo(i).temp,zeros(length(SO_argo(i)),1)) - 1000, 'salt', SO_argo(i).salt, 'temp', SO_argo(i).temp, 'month', SO_argo(i).time(2), 'year', SO_argo(i).time(1), 'time',SO_argo(i).time);
        max_argo_density(u) = max(argo_data_w_den(u).den);
        min_argo_density(u) = min(argo_data_w_den(u).den);
        u = u + 1;
    end
end

density_grid = [];
density_grid(:,1) = linspace(26, 28, 400);
u = 1;
for i = 1:length(argo_data_w_den)
    [a,b] = unique(argo_data_w_den(i).den);
    if length(b) > 2
        ds_interp_argo_data(u) = struct('float', argo_data_w_den(i).float,'cycle', argo_data_w_den(i).cycle, 'lat', argo_data_w_den(i).lat, 'lon', argo_data_w_den(i).lon,'salt',interp1(argo_data_w_den(i).den(b), argo_data_w_den(i).salt(b),density_grid),'temp', interp1(argo_data_w_den(i).den(b), argo_data_w_den(i).temp(b),density_grid), 'month', argo_data_w_den(i).month, 'year', argo_data_w_den(i).year, 'time',argo_data_w_den(i).time);
        u = u + 1;
    end
end