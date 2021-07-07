%%% function to create location map
%%% Note: lat, lon must be vectors

function map(lat, lon, identifier, dataset)
hold on
load coastlines

%%% creating map figure settings
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -60],'MLineLocation', 30, 'PLineLocation', 10, 'FontSize',10)
%set(gcf,'Position',[400 400 700 700]);
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
%legend([s e]);                                              % creating a legend

% %%% creating appropriate title for the dataset (meop tags vs argo floats)
% if strcmp(dataset, 'meop')
%     title('MEOP Tag ' + string(identifier))
% elseif strcmp(dataset, 'argo')
%     title('ARGO Float ' + string(identifier))
% end

hold off
end

