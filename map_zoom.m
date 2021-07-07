%%% function to create zoomed in location map (accompanies "map" function)
%%% Note: lat, lon must be vectors

function map_zoom(lat, lon, identifier, dataset)

hold on
load coastlines

%%% creating map figure settings
axesm('mercator','MapLonLimit', [min(lon) max(lon)],'MapLatLimit', [min(lat) max(lat)],'MLineLocation', 5, 'PLineLocation', 2, 'FontSize',10);
axis off;
framem on;
gridm on;
mlabel on;
plabel on;

%%% plotting profile location data
plotm(lat,lon, '-s');                                       % displaying all profiles
s = plotm(lat(1),lon(1), 'g*', 'DisplayName', 'Start');     % highlighting first profile
e = plotm(lat(end),lon(end), 'r*', 'DisplayName', 'End');   % highlighting last profile
plotm(coastlat, coastlon, 'k')                              % displaying coastlines

% %%% creating appropriate title for the dataset (meop tags vs argo floats)
% if strcmp(dataset, 'meop')
%     title('MEOP Tag ' + string(identifier))
% elseif strcmp(dataset, 'argo')
%     title('ARGO Float ' + string(identifier))
% end

hold off
end