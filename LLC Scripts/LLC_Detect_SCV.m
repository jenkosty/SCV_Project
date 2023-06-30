%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading non-LLC data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading MEOP seal data
load("qc_ts.mat")

%%% Loading RTOPO data
RTOPO.lat = double(ncread('/Volumes/Elements/Raw Data/RTOPO2.nc', 'lat'));
RTOPO.lon = double(ncread('/Volumes/Elements/Raw Data/RTOPO2.nc', 'lon'))';
RTOPO.bedrock_topography = double(ncread('/Volumes/Elements/Raw Data/RTOPO2.nc', 'bedrock_topography'))';

%%% Loading BedMachineAntarctica data
BedMachineAntarctica.x = double(ncread('/Volumes/Elements/Raw Data/BedMachineAntarctica-v3.nc', 'x'));
BedMachineAntarctica.y = double(ncread('/Volumes/Elements/Raw Data/BedMachineAntarctica-v3.nc', 'y'));
%[BedMachineAntarctica.lat, BedMachineAntarctica.lon] = xy2ll(BedMachineAntarctica.x, BedMachineAntarctica.y, -1);
BedMachineAntarctica.bedrock_topography = double(ncread('/Volumes/Elements/Raw Data/BedMachineAntarctica-v3.nc', 'bed'));

%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading LLC data %%%
%%%%%%%%%%%%%%%%%%%%%%%%

%%% Clearing structures
clear LLC_1 LLC_2 LLC_4 LLC_5

%%% Date of LLC snapshot
date = '01-Jan-2012';

%%% Loading mat files
LLC_1.mat = matfile('/Volumes/Elements/LLC_Snapshots/SO_snapshots_NSF_LLC4320_k1-86_face1_' + string(date) + '.mat');
LLC_2.mat = matfile('/Volumes/Elements/LLC_Snapshots/SO_snapshots_NSF_LLC4320_k1-86_face2_' + string(date) + '.mat');
LLC_4.mat = matfile('/Volumes/Elements/LLC_Snapshots/SO_snapshots_NSF_LLC4320_k1-86_face4_' + string(date) + '.mat');
LLC_5.mat = matfile('/Volumes/Elements/LLC_Snapshots/SO_snapshots_NSF_LLC4320_k1-86_face5_' + string(date) + '.mat');

%%% Face 1
LLC_1.edge_lats = [double(LLC_1.mat.yc(1,:)), double(LLC_1.mat.yc(:,end))', flip(double(LLC_1.mat.yc(end,:))), flip(double(LLC_1.mat.yc(:,1)))'];
LLC_1.edge_lons = [double(LLC_1.mat.xc(1,:)), double(LLC_1.mat.xc(:,end))', flip(double(LLC_1.mat.xc(end,:))), flip(double(LLC_1.mat.xc(:,1)))'];
LLC_1.polygon = geopolyshape(LLC_1.edge_lats, LLC_1.edge_lons);
LLC_1.lats = double(LLC_1.mat.yc);
LLC_1.lons = double(LLC_1.mat.xc);
LLC_1.salt = double(LLC_1.mat.s);
LLC_1.temp = double(LLC_1.mat.t);
LLC_1.vort = double(LLC_1.mat.Ro);
LLC_1.depth = LLC_1.mat.rc;
LLC_1.grid = load('/Volumes/Elements/LLC_Snapshots/SO_grids_NSF_LLC4320_k1-86_face1.mat');
LLC_1.date = date;

%%% Face 2
LLC_2.edge_lats = [double(LLC_2.mat.yc(1,:)), double(LLC_2.mat.yc(:,end))', flip(double(LLC_2.mat.yc(end,:))), flip(double(LLC_2.mat.yc(:,1)))'];
LLC_2.edge_lons = [double(LLC_2.mat.xc(1,:)), double(LLC_2.mat.xc(:,end))', flip(double(LLC_2.mat.xc(end,:))), flip(double(LLC_2.mat.xc(:,1)))'];
LLC_2.polygon = geopolyshape(LLC_2.edge_lats, LLC_2.edge_lons);
LLC_2.lats = double(LLC_2.mat.yc);
LLC_2.lons = double(LLC_2.mat.xc);
LLC_2.salt = double(LLC_2.mat.s);
LLC_2.temp = double(LLC_2.mat.t);
LLC_2.vort = double(LLC_2.mat.Ro);
LLC_2.depth = LLC_2.mat.rc;
LLC_2.grid = load('/Volumes/Elements/LLC_Snapshots/SO_grids_NSF_LLC4320_k1-86_face2.mat');
LLC_2.date = date;

%%% Face 4
LLC_4.edge_lats = [double(LLC_4.mat.yc(1,:)), double(LLC_4.mat.yc(:,end))', flip(double(LLC_4.mat.yc(end,:))), flip(double(LLC_4.mat.yc(:,1)))'];
LLC_4.edge_lons = [double(LLC_4.mat.xc(1,:)), double(LLC_4.mat.xc(:,end))', flip(double(LLC_4.mat.xc(end,:))), flip(double(LLC_4.mat.xc(:,1)))'];
LLC_4.polygon = geopolyshape(LLC_4.edge_lats, LLC_4.edge_lons);
LLC_4.lats = double(LLC_4.mat.yc);
LLC_4.lons = double(LLC_4.mat.xc);
LLC_4.salt = double(LLC_4.mat.s);
LLC_4.temp = double(LLC_4.mat.t);
LLC_4.vort = double(LLC_4.mat.Ro);
LLC_4.depth = LLC_4.mat.rc;
LLC_4.grid = load('/Volumes/Elements/LLC_Snapshots/SO_grids_NSF_LLC4320_k1-86_face4.mat');
LLC_4.date = date;

%%% Face 5
LLC_5.edge_lats = [double(LLC_5.mat.yc(1,:)), double(LLC_5.mat.yc(:,end))', flip(double(LLC_5.mat.yc(end,:))), flip(double(LLC_5.mat.yc(:,1)))'];
LLC_5.edge_lons = [double(LLC_5.mat.xc(1,:)), double(LLC_5.mat.xc(:,end))', flip(double(LLC_5.mat.xc(end,:))), flip(double(LLC_5.mat.xc(:,1)))'];
LLC_5.polygon = geopolyshape(LLC_5.edge_lats, LLC_5.edge_lons);
LLC_5.lats = double(LLC_5.mat.yc);
LLC_5.lons = double(LLC_5.mat.xc);
LLC_5.salt = double(LLC_5.mat.s);
LLC_5.temp = double(LLC_5.mat.t);
LLC_5.vort = double(LLC_5.mat.Ro);
LLC_5.depth = LLC_5.mat.rc;
LLC_5.grid = load('/Volumes/Elements/LLC_Snapshots/SO_grids_NSF_LLC4320_k1-86_face5.mat');
LLC_5.date = date;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating the Okubo Weiss parameter %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Depth Averaged Okubo Weiss calculation
LLC_1.OW = [];
LLC_2.OW = [];
LLC_4.OW = [];
LLC_5.OW = [];
LLC_1.OW(:,:) = Okubo_Weiss(LLC_1.grid.dxg(:,:), LLC_1.grid.dyg(:,:), LLC_1.grid.dxc(:,:), LLC_1.grid.dyc(:,:), LLC_1.mat.u(:,:,:), LLC_1.mat.v(:,:,:), LLC_1.depth);
LLC_2.OW(:,:) = Okubo_Weiss(LLC_2.grid.dxg(:,:), LLC_2.grid.dyg(:,:), LLC_2.grid.dxc(:,:), LLC_2.grid.dyc(:,:), LLC_2.mat.u(:,:,:), LLC_2.mat.v(:,:,:), LLC_2.depth);
LLC_4.OW(:,:) = Okubo_Weiss(LLC_4.grid.dxg(:,:), LLC_4.grid.dyg(:,:), LLC_4.grid.dxc(:,:), LLC_4.grid.dyc(:,:), LLC_4.mat.u(:,:,:), LLC_4.mat.v(:,:,:), LLC_4.depth);
LLC_5.OW(:,:) = Okubo_Weiss(LLC_5.grid.dxg(:,:), LLC_5.grid.dyg(:,:), LLC_5.grid.dxc(:,:), LLC_5.grid.dyc(:,:), LLC_5.mat.u(:,:,:), LLC_5.mat.v(:,:,:), LLC_5.depth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Saving LLC Output Data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LLC = {LLC_1, LLC_2, LLC_4, LLC_5};
sectors = {'LLC_1', 'LLC_2', 'LLC_4', 'LLC_5'};
for i = 1:4
   lats = LLC{i}.lats;
   save('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/' + string(sectors{i}) + '/lats', 'lats');
   lons = LLC{i}.lons;
   save('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/' + string(sectors{i}) + '/lons', 'lons');
   grid = LLC{i}.grid;
   save('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/' + string(sectors{i}) + '/grid', 'grid');
   % temp = LLC{i}.temp;
   % save('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/' + string(sectors{i}) + '/temp', 'temp');
   % salt = LLC{i}.salt;
   % save('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/' + string(sectors{i}) + '/salt', 'salt');
   vort = LLC{i}.vort;
   save('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/' + string(sectors{i}) + '/vort', 'vort');
   depth = LLC{i}.depth;
   save('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/' + string(sectors{i}) + '/depth', 'depth');
   OW = LLC{i}.OW;
   save('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/' + string(sectors{i}) + '/OW', 'OW');
   polygon = LLC{i}.polygon;
   save('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/' + string(sectors{i}) + '/polygon', 'polygon');
   save('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/' + string(sectors{i}) + '/date', 'date');
   disp(i);

   clear lats lons grid temp salt vort depth OW polygon
end
clear LLC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolating LLC data to MEOP seal tracks %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Seal tags used for interpolation
test_prof = 1:467;

%%% Pressure grid for interpolation
depth_grid = 0:800;
depth_grid = depth_grid';

for tag_no = test_prof

    %%% Grabbing data on MEOP seal time series
    LLCsealdata(tag_no).tag = qc_ts(tag_no).tag;
    LLCsealdata(tag_no).cast = qc_ts(tag_no).cast;
    LLCsealdata(tag_no).lat = qc_ts(tag_no).lat;
    LLCsealdata(tag_no).lon = qc_ts(tag_no).lon;
    LLCsealdata(tag_no).time = qc_ts(tag_no).time;
    LLCsealdata(tag_no).date = LLC_1.date;

    sector = NaN(1,length(LLCsealdata(tag_no).cast));
    %%% Identifying the LLC sector of each MEOP profile
    for i = 1:length(LLCsealdata(tag_no).cast)
        if isinterior(LLC_1.polygon, geopointshape(LLCsealdata(tag_no).lat(i), LLCsealdata(tag_no).lon(i)))
            sector(i) = 1;
        elseif isinterior(LLC_2.polygon, geopointshape(LLCsealdata(tag_no).lat(i), LLCsealdata(tag_no).lon(i)))
            sector(i) = 2;
        elseif isinterior(LLC_5.polygon, geopointshape(LLCsealdata(tag_no).lat(i), LLCsealdata(tag_no).lon(i)))
            sector(i) = 5;
        else
            sector(i) = 4;
        end 
    end

    LLCsealdata(tag_no).sector = sector;

    %%% Initializing matrices for interpolated data
    LLCseal_salt = NaN(length(depth_grid), length(LLCsealdata(tag_no).cast));
    LLCseal_temp = NaN(length(depth_grid), length(LLCsealdata(tag_no).cast));
    LLCseal_vort = NaN(length(depth_grid), length(LLCsealdata(tag_no).cast));
    LLCseal_OW = NaN(size(LLC_1.OW, 3), length(LLCsealdata(tag_no).cast));

    for i = 1:length(LLCsealdata(tag_no).cast)

        %%% Grabbing LLC data for profile sector
        if sector(i) == 1
            LLC = LLC_1;
        elseif sector(i) == 2
            LLC = LLC_2;
        elseif sector(i) == 4
            LLC = LLC_4;
        elseif sector(i) == 5
            LLC = LLC_5;
        end

        %%% Finding LLC points close to MEOP profile
        LLC_lats = LLC.lats;
        LLC_lons = LLC.lons;
        nan_ind = LLC_lats > LLCsealdata(tag_no).lat(i) + 0.025 | LLC_lats < LLCsealdata(tag_no).lat(i) - 0.025 | LLC_lons > LLCsealdata(tag_no).lon(i) + 0.04 | LLC_lons < LLCsealdata(tag_no).lon(i) - 0.04;
        LLC_lats(nan_ind) = NaN;
        close_ind = find(~isnan(LLC_lats));
        [rows, cols] = ind2sub(size(LLC_lats), close_ind);

        lats_unfmt = NaN(1,length(rows));
        lons_unfmt = NaN(1,length(rows));
        salt = NaN(86, length(rows));
        temp = NaN(86, length(rows));
        vort = NaN(86, length(rows));
        okubo_weiss = NaN(size(LLC.OW, 3), length(rows));

        %%% Extracting LLC data close to MEOP profile
        for j = 1:length(rows)
            lats_unfmt(j) = double(LLC.lats(rows(j), cols(j)));
            lons_unfmt(j) = double(LLC.lons(rows(j), cols(j)));
            salt(:,j) = squeeze(LLC.salt(rows(j), cols(j), :));
            temp(:,j) = squeeze(LLC.temp(rows(j), cols(j), :));
            vort(:,j) = squeeze(LLC.vort(rows(j), cols(j), :));
            okubo_weiss(:,j) = squeeze(LLC.OW(rows(j), cols(j),:));
        end
        temp(salt == 0) = NaN;
        vort(salt == 0) = NaN;
        salt(salt == 0) = NaN;

        lats = double(lats_unfmt .* ones(size(salt)));
        lons = double(lons_unfmt .* ones(size(salt)));
        depths = LLC.depth(1:size(salt, 1), 1 ) .* ones(size(salt));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Interpolating LLC data %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% Salinity
        S = scatteredInterpolant(lats(:), lons(:), depths(:), salt(:),'linear','none');
        S_prof = S(LLCsealdata(tag_no).lat(i).*ones(size(depth_grid)), LLCsealdata(tag_no).lon(i).*ones(size(depth_grid)), -depth_grid);
        if isempty(S_prof)
            LLCseal_salt(:,i) = NaN(size(depth_grid));
        else
            LLCseal_salt(:,i) = S_prof;
        end

        %%% Temperature
        T = scatteredInterpolant(lats(:), lons(:), depths(:), temp(:),'linear','none');
        T_prof = T(LLCsealdata(tag_no).lat(i).*ones(size(depth_grid)), LLCsealdata(tag_no).lon(i).*ones(size(depth_grid)), -depth_grid);
        if isempty(T_prof)
            LLCseal_temp(:,i) = NaN(size(depth_grid));
        else
            LLCseal_temp(:,i) = T_prof;
        end

        %%% Vorticity
        V = scatteredInterpolant(lats(:), lons(:), depths(:), vort(:),'linear','none');
        V_prof = V(LLCsealdata(tag_no).lat(i).*ones(size(depth_grid)), LLCsealdata(tag_no).lon(i).*ones(size(depth_grid)), -depth_grid);
        if isempty(V_prof)
            LLCseal_vort(:,i) = NaN(size(depth_grid));
        else
            LLCseal_vort(:,i) = V_prof;
        end

        %%% Okubo Weiss
        for j = 1:size(LLC.OW, 3)
            OW = scatteredInterpolant(lats_unfmt', lons_unfmt', okubo_weiss(j,:)');
            OW_prof = OW(LLCsealdata(tag_no).lat(i), LLCsealdata(tag_no).lon(i));
            if isempty(OW_prof)
                LLCseal_OW(j,i) = NaN;
            else
                LLCseal_OW(j,i) = OW_prof;
            end
        end

        clear lats_unfmt lons_unfmt salt vort temp okubo_weiss cols rows nan_ind lats lons depths S T V OW j close_ind S_prof T_prof V_prof OW_prof
    end

    LLCsealdata(tag_no).salt = LLCseal_salt;
    LLCsealdata(tag_no).temp = LLCseal_temp;
    LLCsealdata(tag_no).vort = LLCseal_vort;
    LLCsealdata(tag_no).OW = LLCseal_OW;

    clear sector i LLCseal LLCseal_salt LLCseal_temp LLCseal_vort LLCseal_OW LLC_lats LLC_lons
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Flagging Profiles as SCVs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

for tag_no = test_prof
    disp(tag_no)
    for i = find(LLCsealdata(tag_no).scv == 1)
        [scv, note] = trueSCVdetections(LLCsealdata, tag_no, i, LLC_1, LLC_2, LLC_4, LLC_5);
    end
end

%%
for tag_no = test_prof

    disp(tag_no)

    %%% Creating array to hold notes about SCV classification
    %LLCsealdata(tag_no).scv_reason = strings(size(LLCsealdata(tag_no).cast));

    for i = find(LLCsealdata(tag_no).scv == 1)%1:length(LLCsealdata(tag_no).cast)

        %%% Extracting profile sector
        sector = LLCsealdata(tag_no).sector(i);
        if sector == 1
            LLC = LLC_1;
        elseif sector == 2
            LLC = LLC_2;
        elseif sector == 4
            LLC = LLC_4;
        elseif sector == 5
            LLC = LLC_5;
        end

        %%% Extracting grid cell closest to profile
        scv_center_lat = LLCsealdata(tag_no).lat(i);
        scv_center_lon = LLCsealdata(tag_no).lon(i);
        dist = distance(scv_center_lat, scv_center_lon, LLC.lats, LLC.lons);
        [j,k] = find(dist == min(dist(:)));
        %clear dist scv_center_lat scv_center_lon

        %%% Extracting a square grid surrounding the profile
        no_cells = 30;
        j_ind = j-no_cells:j+no_cells;
        if max(j_ind) > 4320 || min(j_ind) < 1
            LLCsealdata(tag_no).scv(i) = 0;
            continue
        end
        k_ind = k-no_cells:k+no_cells;
        if max(k_ind) > 3295 || min(k_ind) < 1
            LLCsealdata(tag_no).scv(i) = 0;
            continue
        end
        %clear no_cells
    
        %%% Extracting lat/lon/OW data in square grid
        lat = LLC.lats(j_ind, k_ind);
        lon = LLC.lons(j_ind, k_ind);
        [x,y] = ll2xy(lat, lon, -1);
        bathymetry = interpBedmachineAntarctica(x, y, 'bed');
        OW = LLC.OW(j_ind, k_ind);
        OWnans = OW;
        OWnans(OW == 0) = NaN;

        % figure()
        % subplot(121)
        % pcolor(x, y, OW)
        % shading interp
        % colormap(cmocean('balance')); colorbar; clim([-2e-9 2e-9]);
        % subplot(122)
        % pcolor(x, y, bathymetry)
        % shading interp
        % colorbar;

        %%% Checking bathymetry
        shallow = sum(double(abs(bathymetry) <= 1000), 'all');
        deep = sum(double(abs(bathymetry) > 1000), 'all');
        if (shallow ~= 0) && (deep ~= 0)
            bathy_ratio = min([shallow deep]) / length(bathymetry(:));
            if bathy_ratio > 0.25
                LLCsealdata(tag_no).scv(i) = 0;
                LLCsealdata(tag_no).scv_reason(i) = 'Proximity to Shelf Break';
                continue
            end
        end

        %%% Checking number of profiles below 350
        if (shallow ~= 0)
            too_shallow = sum(double(abs(bathymetry(:))) < 350);
            shallow_ratio = too_shallow / length(bathymetry(:));
            if shallow_ratio > 0.25
                LLCsealdata(tag_no).scv(i) = 0;
                LLCsealdata(tag_no).scv_reason(i) = 'Too Shallow';
                continue
            end
        end

        %%% Checking depth of profile of interest
        if abs(LLCsealdata(tag_no).bathymetry(i)) < 350
            LLCsealdata(tag_no).scv(i) = 0;
            LLCsealdata(tag_no).scv_reason(i) = 'Profile Too Shallow';
            continue
        end

        %%% Checking standard deviation of OW field
        %disp(std(OWnans(:), 'omitnan'))
        if std(OWnans(:), 'omitnan') > 2e-9
            LLCsealdata(tag_no).scv(i) = 0;
            LLCsealdata(tag_no).scv_reason(i) = 'STD of OW field';
            continue
        end

        %%% Calculating closed curves of low Okubo Weiss
        contour_level = -2*std(OWnans(:), 'omitnan');
        if isnan(contour_level)
            LLCsealdata(tag_no).scv(i) = 0;
            LLCsealdata(tag_no).scv_reason(i) = 'NaNs in OW field';
            continue
        end
        [xc, yc] = closedcurves(lat, lon, OW, -2*std(OWnans(:), 'omitnan'));

        %%% Skipping if no closed curves found
        if isempty(xc) || isempty(yc)
            LLCsealdata(tag_no).scv(i) = 0;
            LLCsealdata(tag_no).scv_reason(i) = 'No closed contours';
            continue
        end

        %%% Checking if SCV lies in one of the closed curves
        for ii = 1:length(xc)
            contour_shape = geopolyshape(xc{ii}, yc{ii});
            point = geopointshape(LLCsealdata(tag_no).lat(i), LLCsealdata(tag_no).lon(i));
            [xo,yo,~,~,~,a,b,theta]=curvemoments(xc{ii},yc{ii});
            A = 1.25 * a;
            B = 1.25 * b;
            z = LLCsealdata(tag_no).lat(i) + (sqrt(-1) * LLCsealdata(tag_no).lon(i));
            wgs84 = wgs84Ellipsoid("km");
            area = areaint(xc{ii},yc{ii},wgs84);
            if isinterior(contour_shape, point) && area > 15 %%% Normal Contour
                ind = ii;
                LLCsealdata(tag_no).scv(i) = 1;
                LLCsealdata(tag_no).scv_reason(i) = 'Normal';
                break
            elseif inellipse(z, sqrt((a^2+b^2)/2),(a^2-b^2)/(a^2+b^2),theta,xo+(sqrt(-1)*yo)) && area > 15 %%% Fitted Ellipse
                ind = ii;
                LLCsealdata(tag_no).scv(i) = 1;
                LLCsealdata(tag_no).scv_reason(i) = 'Fitted';
                break
            elseif inellipse(z, sqrt((A^2+B^2)/2),(A^2-B^2)/(A^2+B^2),theta,xo+(sqrt(-1)*yo)) && area > 15 %%% Fitted and Expanded Ellipse
                ind = ii;
                LLCsealdata(tag_no).scv(i) = 1;
                LLCsealdata(tag_no).scv_reason(i) = 'Expanded + Fitted';
                break
            else %%% Not in Contour
                LLCsealdata(tag_no).scv(i) = 0;
                LLCsealdata(tag_no).scv_reason(i) = 'Not in Contour';
            end
        end
        clear A B a b theta xo yo z ii point contour_shape

        %%% Skipping if point is not enclosed in any contour
        if LLCsealdata(tag_no).scv(i) == 0
            continue
        end

        % %%% Checking circularity of contour
        % [~,~,~,~,~,a,b,~]=curvemoments(xc{ind},yc{ind});
        % ellipse_area = pi * a * b;
        % perimeter = 2*pi*sqrt((a^2 + b^2) / 2);
        % circularity = ((4 * pi * ellipse_area) / perimeter^2);
        % if circularity < 0.4
        %     LLCsealdata(tag_no).scv(i) = 0;
        %     LLCsealdata(tag_no).scv_reason(i) = 'Circularity';
        %     continue
        % end
        %clear ellipse_area perimeter a b

        %%% Checking overlap betweem contour and fitted ellipse
        [xo,yo,~,~,~,a,b,theta]=curvemoments(xc{ind},yc{ind});
        [k,l]=ab2kl(a,b);
        [ellxc,ellyc] = ellcurves(k,l,theta,xo + sqrt(-1)*yo,'npoints',64);
        ellipse_area = areaint(ellxc,ellyc,wgs84);
        ratio = area / ellipse_area;
        if ratio < 0.5
            LLCsealdata(tag_no).scv(i) = 0;
            LLCsealdata(tag_no).scv_reason(i) = 'Ratio';
            continue
        end
        ellpt = b / a;
        clear k l a b theta ellipse_area

        %%% Checking magnitude of contour interior
        points = geopointshape(lat, lon);
        incontour = isinterior(geopolyshape(xc{ind}, yc{ind}), points);
        OW_in_contour = OW(incontour);
        if mean(OW_in_contour) > -0.5e-9
            LLCsealdata(tag_no).scv(i) = 0;
            LLCsealdata(tag_no).scv_reason(i) = 'Low OW';
            continue
        end
        clear points incontour
        
        %%% Saving contour data
        contourdata.xc = xc{ind}; % Raw contour latitude
        contourdata.yc = yc{ind}; % Raw contour longitude
        contourdata.ellxc = ellxc; % Ellipse latitude
        contourdata.ellyc = ellyc; % Ellipse longitude
        contourdata.centerlat = xo; % Center latitude
        contourdata.centerlon = yo; % Center longitude
        contourdata.ellpt = ellpt;
        contourdata.yc_constantlat = unique([linspace(contourdata.centerlon - 0.75, contourdata.centerlon, 13), linspace(contourdata.centerlon + 0.75, contourdata.centerlon, 13)]);
        contourdata.xc_constantlat = contourdata.centerlat * ones(size(contourdata.yc_constantlat));
        contourdata.xc_constantlon = unique([linspace(contourdata.centerlat - 0.25, contourdata.centerlat, 13), linspace(contourdata.centerlat + 0.25, contourdata.centerlat, 13)]);
        contourdata.yc_constantlon = contourdata.centerlon * ones(size(contourdata.xc_constantlon));
        contourdata.lat = lat;
        contourdata.lon = lon;
        contourdata.OW = OW;
        LLCsealdata(tag_no).contourdata(i) = contourdata;

        %%%%%%%%%%%%%%
        %%% Figure %%%
        %%%%%%%%%%%%%%

        % fig = figure('Position', [100 100 1500 500]);
        % sgtitle('Tag #: ' + string(tag_no) + ', Cast: ' + string(i) + ', Contour: ' + string(LLCsealdata(tag_no).scv_reason(i)) + ', Mean OW: ' + string(mean(OW_in_contour) * 1e9) + ', Area: ' + string(area) + ', Ratio: ' + string(ratio) + ' STD: ' + string(std(OWnans(:), 'omitnan')));
        % 
        % %%% Okubo-Weiss subplot
        % subplot(131)
        % h = pcolor(lat, lon, OW); %%% Colorplot
        % set(h, 'EdgeColor', 'none');
        % hold on
        % plot(LLCsealdata(tag_no).lat(i), LLCsealdata(tag_no).lon(i), 'Marker', 'o', 'MarkerSize', 8, 'Color', 'g', 'MarkerFaceColor', 'g')
        % colormap(cmocean('balance')); colorbar; clim([-2e-9 2e-9]);
        % %title('Depth-Averaged Okubo-Weiss')
        % cellplot(xc,yc,'2b') %%% Plotting contours
        % plot(contourdata.centerlat, contourdata.centerlon, 'Marker', 'o', 'MarkerSize', 8, 'Color', 'm', 'MarkerFaceColor', 'm')
        % plot(contourdata.xc_constantlat, contourdata.yc_constantlat, '-o')
        % [xo,yo,~,~,~,a,b,theta]=curvemoments(xc,yc);
        % [k,l]=ab2kl(a,b);
        % ellipseplot(k,l,theta,xo+sqrt(-1)*yo,'2r') %%% Plotting best-fit ellipses
        % [k,l]=ab2kl(1.25*a,1.25*b);
        % ellipseplot(k,l,theta,xo+sqrt(-1)*yo,'2--r') %%% Plotting expanded best-fit ellipses
        % daspect('auto')
        % 
        % %%% Bathymetry subplot
        % ax2 = subplot(132);
        % h = pcolor(lat, lon, bathymetry); %%% Colorplot
        % set(h, 'EdgeColor', 'none');
        % hold on
        % contour(lat, lon, bathymetry, [-1000 -1000], 'k')
        % plot(LLCsealdata(tag_no).lat(i), LLCsealdata(tag_no).lon(i), 'Marker', 'o', 'MarkerSize', 8, 'Color', 'g', 'MarkerFaceColor', 'g')
        % colormap(ax2, flipud(cmocean('deep'))); colorbar;
        % 
        % %%% Map subplot
        % subplot(133)
        % load coastlines
        % axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -57]);
        % axis off; framem on;
        % geoshow(coastlat, coastlon, 'Color', 'k')
        % scatterm(LLCsealdata(tag_no).lat(i), LLCsealdata(tag_no).lon(i), 20, 'ko', 'MarkerFaceColor', 'g')

        %saveas(fig, '/Users/jenkosty/Documents/Research/SCV_Project/Figures/22June2023/Tag' + string(tag_no) + 'Cast' + string(i), 'png')

    end

    clear sector scv_center_lat scv_center_lon dist j k j_ind k_ind lat lon vort Zr R vort_all ...
        LLC area xc yc k l a b xo yo kappa x y z wgs84 theta ratio fig circularity OWnans OW ...
        ind bathymetry coastlat coastlon h ax2 OW_in_contour ellxc ellyc bathy_ratio contourdata ...
        contour_level no_vells shallow shallow_ratio too_shallow 

end

%save('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/LLCsealdata', 'LLCsealdata', 'depth_grid');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Start of Detection Algorithm %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_prof = 1:467;

%%% Loading algorithm settings
load("LLCseals_algorithm_settings.mat")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------------')
disp('Calculating Variables')
disp('---------------------')

for tag_no = test_prof

    %%% Creating arrays to hold rejection and justification data
    rejected.anticyclones = zeros(size(LLCsealdata(tag_no).cast)); % Anticyclones
    reason.anticyclones = strings(size(LLCsealdata(tag_no).cast));
    rejected.cyclones = zeros(size(LLCsealdata(tag_no).cast));     % Cyclones
    reason.cyclones = strings(size(LLCsealdata(tag_no).cast));
    LLCsealdata(tag_no).rejected = rejected;                       
    LLCsealdata(tag_no).reason = reason;                           
    
    %%% Calculating RTOPO2 bathymetry
    LLCsealdata(tag_no).bathymetry = interp2(RTOPO.lon, RTOPO.lat', RTOPO.bedrock_topography, LLCsealdata(tag_no).lon, LLCsealdata(tag_no).lat);

    %%% Calculating BedMachineAntarctica bathymetry
    [x,y] = ll2xy(LLCsealdata(tag_no).lat, LLCsealdata(tag_no).lon, -1);
    LLCsealdata(tag_no).BedMachine = interpBedmachineAntarctica(x, y, 'bed');

    %%% Calculating pressure space variables
    LLCsealdata(tag_no).ps = calc_pres_space_vars(LLCsealdata, tag_no, depth_grid, 1);

end

clear rejected reason x y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Checking for Density Inversions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-------------------------------')
disp('Checking for Density Inversions')
disp('-------------------------------')

for tag_no = test_prof

    %%% Checking for density inversions. Small inversions are corrected.
    %%% Profiles with large inversions will be excluded.
    [LLCsealdata(tag_no).ps.sigma0, LLCsealdata(tag_no).ps.max_sigma0_inversion,...
        LLCsealdata(tag_no).rejected, LLCsealdata(tag_no).reason]...
        = checking_for_density_inversions(LLCsealdata, tag_no);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolating to Density Grid %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('------------------------------')
disp('Interpolating to Density Space')
disp('------------------------------')

%%% Creating density grid
density_grid = (26.0:0.001:28.5)';

for tag_no = test_prof
    
    %%% Interpolating pressure space variables to density grid
    LLCsealdata(tag_no).ds = interp_to_sigma0_space(LLCsealdata, tag_no, density_grid);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Isopycnal Separation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--------------------------------')
disp('Calculating Isopycnal Separation')
disp('--------------------------------')

for tag_no = test_prof
    
    %%% Calculating isopycnal separation
    LLCsealdata(tag_no).ds.isopycnal_separation = calc_isopycnal_separation(LLCsealdata, tag_no);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Finding Indices to Build Reference Profiles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--------------------------------------')
disp('Getting Indices for Reference Profiles')
disp('--------------------------------------')
 
for tag_no = test_prof
        
    %%% Finding indices to build reference profiles
    LLCsealdata(tag_no).ref_ind = indices_for_ref_profiles(LLCsealdata,tag_no, refprof);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Building Reference Profiles %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------------------')
disp('Building Reference Profiles')
disp('---------------------------')

for tag_no = test_prof

    %%% Building density-space reference profiles
    [LLCsealdata(tag_no).ds.ref_salt, LLCsealdata(tag_no).ds.ref_temp,...
        LLCsealdata(tag_no).ds.ref_N2, LLCsealdata(tag_no).ds.ref_spice,...
        LLCsealdata(tag_no).ds.ref_isopycnal_separation] = build_density_space_ref_profiles(LLCsealdata, tag_no);

    %%% Building pressure-space reference profiles
    [LLCsealdata(tag_no).ps.ref_dyn_height_anom, LLCsealdata(tag_no).ps.ref_N2]...
        = build_pres_space_ref_profiles(LLCsealdata, tag_no);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating Anomalies %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------------')
disp('Calculating Anomalies')
disp('---------------------')

for tag_no = test_prof

    %%% Calculating density-space anomalies
    LLCsealdata(tag_no).ds.anoms = calc_density_space_anomalies(LLCsealdata, tag_no);

    %%% Calculating pressure-space anomalies
    LLCsealdata(tag_no).ps.anoms = calc_pres_space_anomalies(LLCsealdata, tag_no);

end

%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculating IQR %%%
%%%%%%%%%%%%%%%%%%%%%%%

disp('---------------')
disp('Calculating IQR')
disp('---------------')

for tag_no = test_prof

    %%% Calculating iqr, lower-threshold, and upper threshold
    LLCsealdata(tag_no).ds.iqrs = calc_iqr(LLCsealdata, tag_no);
    
end

% save('/Volumes/Elements/LLCsealdata/Snapshot_' + string(date) + '/LLCsealdata', 'LLCsealdata', 'depth_grid');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    anticyclones.MLD = [];
    anticyclones.final = [];
    LLCsealdata(tag_no).anticyclones = anticyclones;

    %%% Creating array to hold cyclonic detections
    cyclones.spice_iqr = [];
    cyclones.spice_gaussian = [];
    cyclones.dha = [];
    cyclones.isopycnal_stability = [];
    cyclones.bathymetric_stability = [];
    cyclones.MLD = [];
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

LLCsealdata = anticyclonic_checks(LLCsealdata, test_prof, prms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cyclonic Detection Checks %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-----------------------------')
disp('Beginning Cyclonic Detections')
disp('-----------------------------')

LLCsealdata = cyclonic_checks(LLCsealdata, test_prof, prms);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mixed Layer Depth Check %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('-----------------------')
disp('Mixed Layer Depth Check')
disp('-----------------------')

for tag_no = test_prof
    
    LLCsealdata(tag_no).MLD = NaN(size(LLCsealdata(tag_no).cast));

    [LLCsealdata(tag_no).anticyclones, LLCsealdata(tag_no).cyclones,...
        LLCsealdata(tag_no).MLD, LLCsealdata(tag_no).rejected,...
        LLCsealdata(tag_no).reason] = MLD_check(LLCsealdata, tag_no, prms);

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

scvs = [];
u = 1;
for tag_no = test_prof
    for i = [LLCsealdata(tag_no).anticyclones.final LLCsealdata(tag_no).cyclones.final]

        scvs(u).tag = LLCsealdata(tag_no).tag;
        scvs(u).tag_no = tag_no;
        scvs(u).cast = i;
        scvs(u).time = LLCsealdata(tag_no).time(i,:);
        scvs(u).lat = LLCsealdata(tag_no).lat(i);
        scvs(u).lon = LLCsealdata(tag_no).lon(i);
        if ismember(i, LLCsealdata(tag_no).anticyclones.final)
            scvs(u).type = "Anticyclonic";
        else
            scvs(u).type = "Cyclonic";
        end
        scvs(u).isopycnal_variance = LLCsealdata(tag_no).isopycnal_var(i);
        scvs(u).bathymetric_variance = LLCsealdata(tag_no).bathymetric_var(i);
        scvs(u).OW = LLCsealdata(tag_no).OW(i);

        %%% Checking to see if profile is flagged as an eddy by OW
        if LLCsealdata(tag_no).scv(i) == 1
            scvs(u).OW_check = "Passed";
        else
            scvs(u).OW_check = "Failed";
        end

        u = u + 1;

    end
end

%%% Noting region of detection
for i = 1:length(scvs)
    if scvs(i).lon > -60 && scvs(i).lon < 0
        scvs(i).region = "Weddell";
    elseif scvs(i).lon <= -60 && scvs(i).lon > -120
        scvs(i).region = "WAP";
    elseif scvs(i).lon <= -120 || scvs(i).lon > 170
        scvs(i).region = "Ross";
    else
        scvs(i).region = "East";
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Summarizing Missed SCVs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

missed_scvs = [];
u = 1;
for tag_no = test_prof

    detected = [LLCsealdata(tag_no).anticyclones.final];
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
       
        missed_scvs(u).min_density = min(LLCsealdata(tag_no).ps.sigma0(:), [],'omitnan');
        missed_scvs(u).max_density = max(LLCsealdata(tag_no).ps.sigma0(:), [],'omitnan');

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

% %%% Adjusted DHA calculation for missed SCVs
% if ~isempty(missed_scvs)
%     j = 0;
%     for tag_no = vertcat(missed_scvs.tag_no)'
% 
%         j = j+1;
%         i = missed_scvs(j).cast;
% 
%         %%% DHA check
%         LLCsealdata(tag_no).dha_check{i} = dha_check_modes(LLCsealdata, tag_no, i);
% 
%     end
% end

missed_scvs(strcmp([missed_scvs.type], "Cyclonic")) = [];
failed_scvs = scvs;
failed_scvs(strcmp([failed_scvs.OW_check], "Passed")) = [];

isa = [2, 3, 4];
missed_isa = [34, 34, 35];
false_detections_isa = [115, 101, 88];
real_detections_isa = [6, 6, 5];

clear u OW_flagged detected ind i

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Summarizing Figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

isopycnals = 0.01;

j = 0;
for tag_no = 205
         j = j + 1;
        for i = 184 %find(LLCsealdata(tag_no).scv == 1)

        casts = [LLCsealdata(tag_no).ref_ind{1,i} LLCsealdata(tag_no).ref_ind{2,i} i];
        ind = min(casts):max(casts);

        fig = figure('Position', [0 0 1300 950]);
        sgtitle('LLC Seal Track, Tag = ' + string(LLCsealdata(tag_no).tag) + ', Tag # = ' + string(tag_no) + ', Cast = ' + string(i) + ', Date = ' + string(LLC_1.date), 'FontSize', 18, 'FontWeight', 'bold')

        %%% Bathymetry Subplot
        ax0 = subplot(5,3, [1 2]);
        hold on
        plot(LLCsealdata(tag_no).cast(ind), LLCsealdata(tag_no).bathymetry(ind), 'k-o', 'MarkerFaceColor', 'k', 'LineWidth', 2)
        xline(LLCsealdata(tag_no).cast(i), 'r', 'LineWidth', 1.5)
        hold off
        xlim([LLCsealdata(tag_no).cast(ind(1)) LLCsealdata(tag_no).cast(ind(end))])
        ylabel('Bedrock Topography (m)', 'FontSize', 12)
        title('Bathymetry', 'FontSize', 12)

        %%% Temperature Subplot
        ax1 = subplot(5,3, [4 5]);
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
            depth_ind = 30;
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
        plot(LLCsealdata(tag_no).ds.spice(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        plot(LLCsealdata(tag_no).ds.ref_spice(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'k', 'DisplayName', 'Reference', 'LineWidth', 1.5)
        hold off
        xlabel('Spice');
        set(gca, 'YDir', 'reverse');   
        ylim([0 500]);
        legend('Location', 'best')

        %%% Isopycnal Separation Profile
        subplot(5,4,14);
        hold on
        plot(LLCsealdata(tag_no).ds.isopycnal_separation(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        plot(LLCsealdata(tag_no).ds.ref_isopycnal_separation(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'k','DisplayName', 'Reference', 'LineWidth', 1.5)
        hold off
        xlabel('Isopycnal Separation');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);
        legend('Location', 'best')
        
        %%% N2 Profile
        subplot(5,4,15);
        hold on
        plot(LLCsealdata(tag_no).ds.N2(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        plot(LLCsealdata(tag_no).ds.ref_N2(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'k','DisplayName', 'Reference', 'LineWidth', 1.5)
        hold off
        xlabel('N^2');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);
        legend('Location', 'best')

        %%% DHA Profile
        subplot(5,4,16)
        hold on
        plot(LLCsealdata(tag_no).ps.dyn_height_anom(:,i), LLCsealdata(tag_no).ps.pres(:,i), 'r', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        plot(LLCsealdata(tag_no).ps.ref_dyn_height_anom(:,i), LLCsealdata(tag_no).ps.pres(:,i), 'k','DisplayName', 'Reference', 'LineWidth', 1.5)
        hold off
        xlabel('Dynamic Height Anomaly');
        set(gca, 'YDir', 'reverse');
        legend('Location', 'best')

        %%% Spice Anomaly Profile
        subplot(5,4,17);
        hold on
        plot(LLCsealdata(tag_no).ds.anoms.spice(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'b', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        xline(0, '--k', 'LineWidth', 1);
        hold off
        xlabel('Spice Anomaly');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);

        %%% Isopycnal Separation Anomaly Profile
        subplot(5,4,18);
        hold on
        plot(LLCsealdata(tag_no).ds.anoms.isopycnal_separation_normalized(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'b', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        xline(0, '--k', 'LineWidth', 1);
        hold off
        xlabel('Isopycnal Separation Anomaly');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);

        %%% N2 Anomaly Profile
        subplot(5,4,19);
        hold on
        plot(LLCsealdata(tag_no).ds.anoms.N2(:,i), LLCsealdata(tag_no).ds.pres(:,i), 'b', 'DisplayName', 'Profile', 'LineWidth', 1.5)
        xline(0, '--k', 'LineWidth', 1);
        hold off
        xlabel('N^2 Anomaly');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);

        %%% Dynamic Height Anomaly Profile
        subplot(5,4,20);
        hold on
        xline(0, '--k', 'LineWidth', 1);
        a = plot(LLCsealdata(tag_no).ps.anoms.dyn_height_anom(:,i), LLCsealdata(tag_no).ps.pres(:,i),'b', 'DisplayName', 'PS','LineWidth',1.5);
%         if LLCsealdata(tag_no).dha_check{i}.rejected == 0
%             b = plot(LLCsealdata(tag_no).dha_check{i}.dyn_height_anom_BC1,LLCsealdata(tag_no).dha_check{i}.dyn_height_pres_BC1,'m','LineWidth',1.5, 'DisplayName', 'Adjusted');
%             legend([a b], 'Location', 'best')
%         end
        hold off
        xlabel('Dynamic Height Anomaly Anomaly');
        set(gca, 'YDir', 'reverse');
        ylim([0 500]);
        
        %saveas(fig, '/Volumes/Elements/LLCsealtrack Figures/8Mar2023/' + string(LLCsealdata(tag_no).tag) + '_' + string(i), 'png')

        end
end

clear ax0 ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8 ax9 C h cmap fig h IB isopycnals p1 p2 p3 p4 p5 pp casts lon_max...
    lon_min lat_max lat_min coastlat coastlon i j a b vort_var OW_var excluded_lats excluded_lons excluded

%%
%%%%%%%%%%%%%%%%%%
%%% Map Figure %%%
%%%%%%%%%%%%%%%%%%

test_prof = 1:467;

%%% Extracting SCVs from the LLC seal data
LLCscvs = [];
u = 1;
for tag_no = test_prof
    flagged_scvs = find(LLCsealdata(tag_no).scv);
    for i = flagged_scvs
        LLCscvs(u).tag = LLCsealdata(tag_no).tag;
        LLCscvs(u).tag_no = tag_no;
        LLCscvs(u).cast = i;
        LLCscvs(u).lat = LLCsealdata(tag_no).lat(i);
        LLCscvs(u).lon = LLCsealdata(tag_no).lon(i);
        LLCscvs(u).date = LLCsealdata(tag_no).date;
        LLCscvs(u).vort = mean(LLCsealdata(tag_no).vort(:,i), 'omitnan');
        LLCscvs(u).OW = LLCsealdata(tag_no).OW(i);
        LLCscvs(u).label = string(LLCscvs(u).tag_no) + ', ' + string(LLCscvs(u).cast);
        u = u + 1;
    end
end

%%% Figure settings
load coastlines
figure('Position', [500 100 1000 850])
axesm('stereo','Origin',[-90 0],'MapLatLimit',[-90 -57],'MLineLocation', 30, 'PLineLocation', 10, 'FontSize',10);
axis off; framem on; gridm on; mlabel on; plabel on;
hold on

%%% Plotting LLC data (Choose either Okubo Weiss or Vorticity)
% colormap(cmocean('balance')); colorbar; clim([-1e-8 1e-8]);
% pcolorm(LLC_1.lats(:,:), LLC_1.lons(:,:), LLC_1.OW(:,:));
% pcolorm(LLC_2.lats(:,:), LLC_2.lons(:,:), LLC_2.OW(:,:));
% pcolorm(LLC_4.lats(:,:), LLC_4.lons(:,:), LLC_4.OW(:,:));
% pcolorm(LLC_5.lats(:,:), LLC_5.lons(:,:), LLC_5.OW(:,:));
% title('Depth-Averaged Okubo Weiss for ' + string(LLC_1.date), 'FontSize', 20)
% geoshow(coastlat, coastlon, 'Color', 'k')

%%% Plotting LLC data (Choose either Okubo Weiss or Vorticity)
colormap(cmocean('balance')); colorbar; clim([-0.4 0.4]);
depth_ind = 28:37;
pcolorm(LLC_1.lats(:,:), LLC_1.lons(:,:), mean(LLC_1.vort(:,:,depth_ind), 3, 'includenan'));
pcolorm(LLC_2.lats(:,:), LLC_2.lons(:,:), mean(LLC_2.vort(:,:,depth_ind), 3, 'includenan'));
pcolorm(LLC_4.lats(:,:), LLC_4.lons(:,:), mean(LLC_4.vort(:,:,depth_ind), 3, 'includenan'));
pcolorm(LLC_5.lats(:,:), LLC_5.lons(:,:), mean(LLC_5.vort(:,:,depth_ind), 3, 'includenan'));
title('Depth-Averaged Vortex Rossby Number for ' + string(LLC_1.date), 'FontSize', 20)
geoshow(coastlat, coastlon, 'Color', 'k')

%%% Plotting and labeling location of SCVs in the LLC seal data
scatterm([LLCscvs.lat], [LLCscvs.lon], 20, 'ks', 'markerfacecolor','g')
textm([LLCscvs.lat], [LLCscvs.lon], [LLCscvs.label], 'color', 'white')

clear coastlat coastlon u flagged_scvs
