%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loads LLC map and seal track data from specified snapshot %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snapshot_dates = {'01-Oct-2011','01-Nov-2011', '01-Dec-2011', '01-Jan-2012', '01-Feb-2012', '01-Mar-2012', '01-Apr-2012', '01-May-2012', ...
    '01-Jun-2012', '01-Jul-2012', '01-Aug-2012', '01-Sep-2012'};

%%% Snapshot date
date = '01-May-2012';

%%% Loading LLC Okubo Weiss Data
LLC = cell(4,1);
sectors = {'LLC_1', 'LLC_2', 'LLC_4', 'LLC_5'};
for i = 1:4
    load('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/' + string(sectors{i}) + '/lats.mat');
    LLC{i}.lats = lats;
    load('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/' + string(sectors{i}) + '/lons.mat');
    LLC{i}.lons = lons;
    load('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/' + string(sectors{i}) + '/OW.mat');
    LLC{i}.OW = OW;
    load('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/' + string(sectors{i}) + '/vort.mat');
    LLC{i}.vort = vort;
    load('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/' + string(sectors{i}) + '/depth.mat')
    LLC{i}.depth = depth;
    LLC{i}.date = date;
    disp(i)
    clear lats lons OW vort depth
end
LLC_1 = LLC{1};
LLC_2 = LLC{2};
LLC_4 = LLC{3};
LLC_5 = LLC{4};
clear LLC i sectors

%%% Loading LLC seal track data
load('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/LLCsealdata.mat')

%%% Extracting SCVs from seal tracks
LLCscvs = [];
u = length(LLCscvs)+1;
for tag_no = 1:length(LLCsealdata)
    flagged_scvs = find(LLCsealdata(tag_no).scv);
    for i = flagged_scvs
        LLCscvs(u).tag = LLCsealdata(tag_no).tag;
        LLCscvs(u).tag_no = tag_no;
        LLCscvs(u).cast = i;
        LLCscvs(u).lat = LLCsealdata(tag_no).lat(i);
        LLCscvs(u).lon = LLCsealdata(tag_no).lon(i);
        LLCscvs(u).date = LLCsealdata(tag_no).date;
        if mean(LLCsealdata(tag_no).vort(:,i), 1, 'omitnan') > 0
            LLCscvs(u).type = "Cyclonic";
        else
            LLCscvs(u).type = "Anticyclonic";
        end
        u = u + 1;
    end
end

%%% Saving SCV time series
tags = unique([LLCscvs.tag_no]);
dstr = erase(date, '-');
assignin('base', 'SCVs_' + string(dstr), [LLCsealdata(tags)])
% save('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/scvs', 'SCVs_' + string(dstr));

clear u tag_no snapshot snapshot_no LLCsnapshots i flagged_scvs u dstr tags

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DETECTION ALGORITHM START
%%

test_prof = 1:length(LLCsealdata);

for tag_no = test_prof

    %%% Calculating iqr, lower-threshold, and upper threshold
    LLCsealdata(tag_no).ds.iqrs = calc_iqr(LLCsealdata, tag_no);
    
end

%%% Loading algorithm settings
run('LLCseals_algorithm_settings.m')

isa_amplitude = 7;
dha_amplitude = 0; %0:0.02:0.1;

uu = 1;
for isa = 1:length(isa_amplitude)
    %prms.isa_gauss.amplitude = isa_amplitude(isa);

    for dha = 1:length(dha_amplitude)
        %prms.dha.amplitude = dha_amplitude(dha);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Creating Arrays for Detection Results %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp('-------------------------')
        disp('Creating Detection Arrays')
        disp('-------------------------')

        for tag_no = test_prof

            %%% Creating array to hold anticyclonic detections
            anticyclones.isa_iqr = [];
            anticyclones.isa_gaussian = [];
            anticyclones.dha = [];
            anticyclones.isopycnal_stability = [];
            anticyclones.bathymetric_stability = [];
            anticyclones.final = [];
            LLCsealdata(tag_no).anticyclones = anticyclones;

            %%% Creating array to hold cyclonic detections
            cyclones.spice_iqr = [];
            cyclones.spice_gaussian = [];
            cyclones.dha = [];
            cyclones.isopycnal_stability = [];
            cyclones.bathymetric_stability = [];
            cyclones.final = [];
            LLCsealdata(tag_no).cyclones = cyclones;

            %%% Creating arrays to hold profile rejections and reasons
            rejected.anticyclones = zeros(size(LLCsealdata(tag_no).cast)); 
            reason.anticyclones = strings(size(LLCsealdata(tag_no).cast));
            rejected.cyclones = zeros(size(LLCsealdata(tag_no).cast));
            reason.cyclones = strings(size(LLCsealdata(tag_no).cast));
            LLCsealdata(tag_no).reason = reason;
            LLCsealdata(tag_no).rejected = rejected;

        end

        clear anticyclones cyclones

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Anticyclonic Detections %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        disp('---------------------------------')
        disp('Beginning Anticyclonic Detections')
        disp('---------------------------------')

        LLCsealdata = anticyclonic_checks(LLCsealdata, test_prof);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Cyclonic Detection Checks %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp('-----------------------------')
        disp('Beginning Cyclonic Detections')
        disp('-----------------------------')

        %LLCsealdata = cyclonic_checks(LLCsealdata, test_prof);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Isopycnal Stability Check (anticyclones AND cyclones) %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp('-------------------------')
        disp('Isopycnal Stability Check')
        disp('-------------------------')

        for tag_no = test_prof

            %%% Creating array to hold isopycnal separation variance data
            LLCsealdata(tag_no).isopycnal_stability = NaN(size(LLCsealdata(tag_no).cast));

            %%% Isopycnal stability check (anticyclones AND cyclones)
            [LLCsealdata(tag_no).anticyclones, LLCsealdata(tag_no).cyclones, ...
                LLCsealdata(tag_no).isopycnal_stability, LLCsealdata(tag_no).rejected, ...
                LLCsealdata(tag_no).reason] = isopycnal_stability_check(LLCsealdata, tag_no, prms);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Bathymetric Stability Check (anticyclones AND cyclones) %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp('---------------------------');
        disp('Bathymetric Stability Check');
        disp('---------------------------');

        for tag_no = test_prof

            %%% Creating array to hold bathymetric variance data
            LLCsealdata(tag_no).bathymetric_var = NaN(size(LLCsealdata(tag_no).cast));

            %%% Bathymetric stability check (anticyclones AND cyclones)
            [LLCsealdata(tag_no).anticyclones, LLCsealdata(tag_no).cyclones,...
                LLCsealdata(tag_no).bathymetric_var, LLCsealdata(tag_no).rejected,...
                LLCsealdata(tag_no).reason] = bathymetry_check(LLCsealdata, tag_no, prms);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Removing Edge Cases %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp('-------------------');
        disp('Removing Edge Cases');
        disp('-------------------');

        for tag_no = test_prof

            %%% Removing edge cases (anticyclones AND cyclones)
            [LLCsealdata(tag_no).anticyclones, LLCsealdata(tag_no).cyclones,...
                LLCsealdata(tag_no).rejected, LLCsealdata(tag_no).reason]...
                = removing_edge_cases(LLCsealdata, tag_no);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Summarizing Detected SCVs %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp('----------------------')
        disp('Summarizing Detections')
        disp('----------------------')

        detected_scvs = [];
        u = 1;
        for tag_no = test_prof
            for i = [LLCsealdata(tag_no).anticyclones.final LLCsealdata(tag_no).cyclones.final]

                detected_scvs(u).tag = LLCsealdata(tag_no).tag;
                detected_scvs(u).tag_no = tag_no;
                detected_scvs(u).cast = i;
                detected_scvs(u).lat = LLCsealdata(tag_no).lat(i);
                detected_scvs(u).lon = LLCsealdata(tag_no).lon(i);
                if ismember(i, LLCsealdata(tag_no).anticyclones.final)
                    detected_scvs(u).type = "Anticyclonic";
                else
                    detected_scvs(u).type = "Cyclonic";
                end
                detected_scvs(u).isopycnal_stability = LLCsealdata(tag_no).isopycnal_stability(i);
                detected_scvs(u).bathymetric_variance = LLCsealdata(tag_no).bathymetric_var(i);
                detected_scvs(u).OW = LLCsealdata(tag_no).OW(i);

                %%% Checking to see if profile is flagged as an eddy by OW
                if LLCsealdata(tag_no).scv(i) == 1
                    detected_scvs(u).True_SCV = "Yes";
                else
                    detected_scvs(u).True_SCV = "No";
                end

                u = u + 1;

            end
        end

        %%% Noting region of detection
        for i = 1:length(detected_scvs)
            if detected_scvs(i).lon > -60 && detected_scvs(i).lon < 0
                detected_scvs(i).region = "Weddell";
            elseif detected_scvs(i).lon <= -60 && detected_scvs(i).lon > -120
                detected_scvs(i).region = "WAP";
            elseif detected_scvs(i).lon <= -120 || detected_scvs(i).lon > 170
                detected_scvs(i).region = "Ross";
            else
                detected_scvs(i).region = "East";
            end
        end
       
        %%% Calculating number of true/false positives
        if ~isempty(detected_scvs)
            true_positives = size(detected_scvs(strcmp([detected_scvs.True_SCV], "Yes")), 2);
            false_positives = size(detected_scvs(strcmp([detected_scvs.True_SCV], "No")), 2);
        else
            true_positives = 0;
            false_positives = 0;
        end
        disp('True Positives: ' + string(true_positives))
        disp('False Positives: ' + string(false_positives))

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Summarizing Missed SCVs %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        missed_scvs = [];
        u = 1;
        for tag_no = test_prof

            detected = [LLCsealdata(tag_no).anticyclones.final LLCsealdata(tag_no).cyclones.final];
            OW_flagged = find(LLCsealdata(tag_no).scv);
            ind = setdiff(OW_flagged, detected);

            for i = ind

                missed_scvs(u).tag = LLCsealdata(tag_no).tag;
                missed_scvs(u).tag_no = tag_no;
                missed_scvs(u).cast = i;
                missed_scvs(u).lat = LLCsealdata(tag_no).lat(i);
                missed_scvs(u).lon = LLCsealdata(tag_no).lon(i);

                %%% Noting type of SCV and reason for exclusion
                if mean(LLCsealdata(tag_no).vort(:,i), 1, 'omitnan') > 0
                    missed_scvs(u).type = "Cyclonic";
                    missed_scvs(u).reason = LLCsealdata(tag_no).reason.cyclones(i);
                else
                    missed_scvs(u).type = "Anticyclonic";
                    missed_scvs(u).reason = LLCsealdata(tag_no).reason.anticyclones(i);
                end

                missed_scvs(u).OW = LLCsealdata(tag_no).OW(i);
                missed_scvs(u).bathymetric_variance = LLCsealdata(tag_no).bathymetric_var(i);
                missed_scvs(u).isopycnal_stability = LLCsealdata(tag_no).isopycnal_stability(i);

                u = u + 1;

            end
        end

        %%% Noting region of occurance
        for i = 1:length(missed_scvs)
            if missed_scvs(i).lon > -60 && missed_scvs(i).lon < 0
                missed_scvs(i).region = "Weddell";
            elseif missed_scvs(i).lon <= -60 && missed_scvs(i).lon > -120
                missed_scvs(i).region = "WAP";
            elseif missed_scvs(i).lon <= -120 || missed_scvs(i).lon > 170
                missed_scvs(i).region = "Ross";
            else
                missed_scvs(i).region = "East";
            end
        end

        clear u OW_flagged detected ind i

        %%% Calculating number of true/false negatives
        false_negatives = length(missed_scvs);
        disp('False Negatives: ' + string(false_negatives))

        true_negatives = length([LLCsealdata(test_prof).cast]) - false_negatives - false_positives - true_positives; 
        disp('True Negatives: ' + string(true_negatives))

        %%% Summarizing results
        results(uu).isa_amplitude = prms.isa_gauss.amplitude;
        results(uu).dha_amplitude = prms.dha.amplitude;
        results(uu).false_positives = false_positives;
        results(uu).true_positives = true_positives;
        results(uu).false_negatives = false_negatives;
        results(uu).true_negatives = true_negatives;
        results(uu).sensitivity = true_positives / (true_positives + false_negatives);
        results(uu).specificity = true_positives / (true_negatives + false_positives);
        results(uu).precision = true_positives / (true_positives + false_positives);
        results(uu).missed_scvs = missed_scvs;
        results(uu).detected_scvs = detected_scvs;
        uu = uu + 1;

        clear false_positives false_negatives true_positives true_negatives
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure of Optimization Results %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position', [100 100 800 800])
sgtitle('Anticyclone Optimization')
m = 6;

subplot(m,1,1)
hold on
for i = 1:length(detected_scvs)
    if strcmp(detected_scvs(i).True_SCV, 'Yes')
        clr = 'g';
        z = 1;
    else
        clr = 'r';
        z = 0;
    end
    plot(detected_scvs(i).isopycnal_stability,z, 'LineStyle','none', 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', clr)
end
for i = 1:length(missed_scvs)
    if strcmp(missed_scvs(i).type, "Anticyclonic")
        plot(missed_scvs(i).isopycnal_stability, -1, 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c')
    end
end
xlabel('Isopycnal Stability')
ylim([-1 1])

subplot(m,1,2)
hold on
for i = 1:length(detected_scvs)
    if strcmp(detected_scvs(i).True_SCV, 'Yes')
        clr = 'g';
        z = 1;
    else
        clr = 'r';
        z = 0;
    end
    plot(detected_scvs(i).bathymetric_variance,z, 'LineStyle','none', 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', clr)
end
for i = 1:length(missed_scvs)
    if strcmp(missed_scvs(i).type, "Anticyclonic")
        plot(missed_scvs(i).bathymetric_variance, -1, 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c')
    end
end
xlabel('Bathymetric Stability')
ylim([-1 1])

subplot(m,1,3)
hold on
for j = 1:length(detected_scvs)
    if strcmp(detected_scvs(j).True_SCV, 'Yes')
        clr = 'g';
        z = 1;
    else
        clr = 'r';
        z = 0;
    end
    tag_no = detected_scvs(j).tag_no;
    i = detected_scvs(j).cast;
    plot(LLCsealdata(tag_no).gauss(i).A, z, 'LineStyle','none', 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', clr)
end
xlabel('Max Isopycnal Separation Anomaly (Normalized)')
ylim([-1 1])

subplot(m,1,4)
hold on
for j = 1:length(detected_scvs)
    if strcmp(detected_scvs(j).True_SCV, 'Yes')
        clr = 'g';
        z = 1;
    else
        clr = 'r';
        z = 0;
    end
    tag_no = detected_scvs(j).tag_no;
    i = detected_scvs(j).cast;
    plot(LLCsealdata(tag_no).gauss(i).H, z, 'LineStyle','none', 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', clr)
end
xlabel('Isopycnal Separation Anomaly Width')
ylim([-1 1])

subplot(m,1,5)
hold on
for j = 1:length(detected_scvs)
    if strcmp(detected_scvs(j).True_SCV, 'Yes')
        clr = 'g';
        z = 1;
    else
        clr = 'r';
        z = 0;
    end
    tag_no = detected_scvs(j).tag_no;
    i = detected_scvs(j).cast;
    plot(LLCsealdata(tag_no).gauss(i).minlse, z, 'LineStyle','none', 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', clr)
end
xlabel('Min LSE')
ylim([-1 1])

subplot(m,1,6)
hold on
for j = 1:length(detected_scvs)
    if strcmp(detected_scvs(j).True_SCV, 'Yes')
        clr = 'g';
        z = 1;
    else
        clr = 'r';
        z = 0;
    end
    tag_no = detected_scvs(j).tag_no;
    i = detected_scvs(j).cast;
    plot(max(LLCsealdata(tag_no).ps.anoms.dyn_height_anom(:,i)), z, 'LineStyle','none', 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', clr)
end
xlabel('Max Dynamic Height Anomaly')
ylim([-1 1])

clear clr z
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Summarizing Figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

isopycnals = 0.01;
lw = 2; % linewidth for plotting

j = 0;
for tag_no = vertcat(detected_scvs.tag_no)'
        j = j + 1;
        for i = detected_scvs(j).cast

        casts = [LLCsealdata(tag_no).ref_ind{1,i} LLCsealdata(tag_no).ref_ind{2,i} i];
        ind = min(casts):max(casts);

        fig = figure('Position', [0 0 1300 950]);
        sgtitle('LLC Seal Track, Tag = ' + string(LLCsealdata(tag_no).tag) + ', Tag # = ' + string(tag_no) + ', Cast = ' + string(i) + ', True SCV: ' + string(detected_scvs(j).True_SCV) + ', Type: ' + string(detected_scvs(j).type), 'FontSize', 18, 'FontWeight', 'bold')
        %sgtitle('LLC Seal Track, Tag = ' + string(LLCsealdata(tag_no).tag) + ', Tag # = ' + string(tag_no) + ', Cast = ' + string(i) + ', Date = ' + string(LLC_1.date), 'FontSize', 18, 'FontWeight', 'bold')

        %%% Bathymetry Subplot
        ax0 = subplot(5,3, [1 2]);
        hold on
        a = plot(LLCsealdata(tag_no).cast(ind), LLCsealdata(tag_no).bathymetry(ind), 'k-o', 'MarkerFaceColor', 'k', 'LineWidth', 2, 'DisplayName', 'RTOPO2');
        b = plot(LLCsealdata(tag_no).cast(ind), LLCsealdata(tag_no).BedMachine(ind), 'r-o', 'MarkerFaceColor', 'r', 'LineWidth', 2, 'DisplayName', 'BedMachine');
        xline(LLCsealdata(tag_no).cast(i), 'r', 'LineWidth', 1.5)
        hold off
        xlim([LLCsealdata(tag_no).cast(ind(1)) LLCsealdata(tag_no).cast(ind(end))])
        ylabel('Bedrock Topography (m)', 'FontSize', 12)
        title('Bathymetry', 'FontSize', 12)
        legend([a b])

        %%% Temperature Subplot
        ax1 = subplot(5,3,[4 5]);
        hold on
        pp = pcolor(LLCsealdata(tag_no).cast(ind), LLCsealdata(tag_no).ps.pres(:,1), LLCsealdata(tag_no).temp(:,ind));
        set(pp, 'EdgeColor', 'none');
        [C,h] = contour(ax1, LLCsealdata(tag_no).cast(ind), depth_grid, LLCsealdata(tag_no).ps.sigma0(:,ind), round(min(min(LLCsealdata(tag_no).ps.sigma0)):isopycnals:max(max(LLCsealdata(tag_no).ps.sigma0)), 2), 'k');
        clabel(C,h,'LabelSpacing',500);
        xline(LLCsealdata(tag_no).cast(i), 'r', 'LineWidth', 1.5)
        hold off
        cmap = cmocean('thermal'); colormap(ax1, cmap); colorbar;
        clim([min(min(LLCsealdata(tag_no).temp(:,ind))) max(max(LLCsealdata(tag_no).temp(:,ind)))])
        set(gca, 'YDir','reverse'); set(gca, 'Layer','top');
        ylabel('Pressure (dbar)', 'FontSize', 12); ylim([0 500])
        title('Temperature', 'FontSize', 12);

        %%% Salinity Subplot
        ax2 = subplot(5,3,[7 8]);
        hold on
        pp = pcolor(LLCsealdata(tag_no).cast(ind), LLCsealdata(tag_no).ps.pres(:,1), LLCsealdata(tag_no).salt(:,ind));
        set(pp, 'EdgeColor', 'none');
        [C,h] = contour(ax2, LLCsealdata(tag_no).cast(ind), depth_grid, LLCsealdata(tag_no).ps.sigma0(:,ind), round(min(min(LLCsealdata(tag_no).ps.sigma0)):isopycnals:max(max(LLCsealdata(tag_no).ps.sigma0)), 2), 'k');
        clabel(C,h,'LabelSpacing',500);
        xline(LLCsealdata(tag_no).cast(i), 'r', 'LineWidth', 1.5)
        hold off
        cmap = cmocean('haline'); colormap(ax2, cmap); colorbar;
        clim([min(min(LLCsealdata(tag_no).salt(:,ind))) max(max(LLCsealdata(tag_no).salt(:,ind)))])
        set(gca, 'YDir','reverse'); set(gca, 'Layer','top');
        ylabel('Pressure (dbar)', 'FontSize', 12); ylim([0 500]);
        title('Salinity', 'FontSize', 12);

        %%% Realigning Plots
        p0 = get(ax0, 'Position');
        p1 = get(ax1, 'Position');
        p2 = get(ax2, 'Position');
        p1(3) = p0(3);
        p2(3) = p0(3);
        set(ax1, 'Position', p1);
        set(ax2, 'Position', p2);
        clear p0 p1 p2

        %%% Map Subplot
        subplot(5,3,3)
        load coastlines
        axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -57]);
        axis off; framem on; 
        geoshow(coastlat, coastlon, 'Color', 'k')
        scatterm(LLCsealdata(tag_no).lat(i), LLCsealdata(tag_no).lon(i), 20, 'ko', 'MarkerFaceColor', 'g')

        %%% Zoomed Map Subplot
        ax3 = subplot(5,3,[6 9]);
        hold on

            %%% Plotting basemap of LLC time series
            depth_ind = 34;
            if ismember(4, unique(LLCsealdata(tag_no).sector(ind)))
                pp = pcolorps(LLC_4.lats(:,:), LLC_4.lons(:,:), LLC_4.vort(:,:,depth_ind));
                set(pp, 'EdgeColor', 'none');
            end
            if ismember(1, unique(LLCsealdata(tag_no).sector(ind)))
                pp = pcolorps(LLC_1.lats(:,:), LLC_1.lons(:,:), LLC_1.vort(:,:,depth_ind));
                set(pp, 'EdgeColor', 'none');
            end
            if ismember(2, unique(LLCsealdata(tag_no).sector(ind)))
                pp = pcolorps(LLC_2.lats(:,:), LLC_2.lons(:,:), LLC_2.vort(:,:,depth_ind));
                set(pp, 'EdgeColor', 'none');
            end
            if ismember(5, unique(LLCsealdata(tag_no).sector(ind)))
                pp = pcolorps(LLC_5.lats(:,:), LLC_5.lons(:,:), LLC_5.vort(:,:,depth_ind));
                set(pp, 'EdgeColor', 'none');
            end

            %%% Plotting data
            colormap(ax3, cmocean('balance')); colorbar; clim([-0.4 0.4]);
            plotps(LLCsealdata(tag_no).lat(ind), LLCsealdata(tag_no).lon(ind), '-s', 'Color', [0.5 0.5 0.5],...
                'MarkerSize', 4, 'MarkerFaceColor', [0.5 0.5 0.5], 'LineWidth', 2)
            scatterps(LLCsealdata(tag_no).lat(i), LLCsealdata(tag_no).lon(i), 6, 'gs', 'MarkerFaceColor', 'g')
            mapzoomps(LLCsealdata(tag_no).lat(i), LLCsealdata(tag_no).lon(i), 'size', [75 75])
            title('Depth = ' + string(round(LLC_1.depth(depth_ind))) + 'm')

        %%% Spice Profile
        subplot(5,4,13);
        hold on
        plot(LLCsealdata(tag_no).ds.spice(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', lw)
        plot(LLCsealdata(tag_no).ds.ref_spice(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'k', 'DisplayName', 'Reference', 'LineWidth', lw)
        hold off
        xlabel('Spice');
        set(gca, 'YDir', 'reverse');   
        ylim([0 500]);
        legend('Location', 'best')

        %%% Isopycnal Separation Profile
        subplot(5,4,14);
        hold on
        plot(LLCsealdata(tag_no).ds.isopycnal_separation(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', lw)
        plot(LLCsealdata(tag_no).ds.ref_isopycnal_separation(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'k','DisplayName', 'Reference', 'LineWidth', lw)
        hold off
        xlabel('Isopycnal Separation');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);
        legend('Location', 'best')
        
        %%% N2 Profile
        subplot(5,4,15);
        hold on
        plot(LLCsealdata(tag_no).ds.N2(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', lw)
        plot(LLCsealdata(tag_no).ds.ref_N2(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'k','DisplayName', 'Reference', 'LineWidth', lw)
        hold off
        xlabel('N^2');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);
        legend('Location', 'best')

        %%% DHA Profile
        subplot(5,4,16)
        hold on
        plot(LLCsealdata(tag_no).ps.dyn_height_anom(:,i), LLCsealdata(tag_no).ps.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', lw)
        plot(LLCsealdata(tag_no).ps.ref_dyn_height_anom(:,i), LLCsealdata(tag_no).ps.pres(:,i), 'k','DisplayName', 'Reference', 'LineWidth', lw)
        hold off
        xlabel('Dynamic Height Anomaly');
        set(gca, 'YDir', 'reverse');
        legend('Location', 'best')

        %%% Spice Anomaly Profile
        subplot(5,4,17);
        hold on
        plot(LLCsealdata(tag_no).ds.anoms.spice(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'b', 'DisplayName', 'Profile', 'LineWidth', lw)
        xline(0, '--k', 'LineWidth', 1);
        hold off
        xlabel('Spice Anomaly');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);

        %%% Isopycnal Separation Anomaly Profile
        subplot(5,4,18);
        hold on
        plot(LLCsealdata(tag_no).ds.anoms.isopycnal_separation(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'b', 'DisplayName', 'Profile', 'LineWidth', lw);
        a = plot(LLCsealdata(tag_no).ds.anoms.isopycnal_separation_normalized(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'm', 'DisplayName', 'Normalized', 'LineWidth', lw);
        xline(0, '--k', 'LineWidth', 1);
        hold off
        xlabel('Isopycnal Separation Anomaly');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);
        legend(a, 'Location', 'best')

        clear a 

        %%% N2 Anomaly Profile
        subplot(5,4,19);
        hold on
        plot(LLCsealdata(tag_no).ds.anoms.N2(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'b', 'DisplayName', 'Profile', 'LineWidth', lw)
        xline(0, '--k', 'LineWidth', 1);
        hold off
        xlabel('N^2 Anomaly');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);

        %%% Dynamic Height Anomaly Profile
        subplot(5,4,20);
        hold on
        xline(0, '--k', 'LineWidth', 1);
        plot(LLCsealdata(tag_no).ps.anoms.dyn_height_anom(:,i), LLCsealdata(tag_no).ps.pres(:,i),'b', 'DisplayName', 'PS','LineWidth',lw);
        hold off
        xlabel('Dynamic Height Anomaly Anomaly');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);
        
        %saveas(fig, '/Volumes/Elements/LLCsealtrack Figures/8Mar2023/' + string(LLCsealdata(tag_no).tag) + '_' + string(i), 'png')

        end
end

clear ax0 ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 C h cmap fig h IB isopycnals p1 p2 p3 p4 p5 pp casts lon_max...
    lon_min lat_max lat_min coastlat coastlon i j a b vort_var OW_var excluded_lats excluded_lons excluded lw

