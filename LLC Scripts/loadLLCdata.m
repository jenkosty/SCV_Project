%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function to extract LLC data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% INPUTS:
    %%% sector - (1,2,4,5)
        %%% 1 = Eastern Weddell, Northern East
        %%% 2 = Eastern East
        %%% 4 = Southern East, Ross
        %%% 5 = WAP, Western Weddell
    %%% depth_levels

%%% OUTPUTS:
    %%% LLC - struct with lat, lon, temp, salt, vorticity data 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LLC = loadLLCdata(sector, depth_levels)

    %%% Loading matfile
    LLC_data = matfile('/Volumes/Elements/LLC_Snapshots/SO_snapshots_NSF_LLC4320_k1-86_face' + string(sector) + '_01-Dec-2011.mat');
    
    %%% Longitude Data
    xc = double(LLC_data.xc);
    
    %%% Latitude Data
    yc = double(LLC_data.yc);
    
    %%% Depth Data
    depth = LLC_data.rc(depth_levels,:);
    
    %%% Salinity Data
    salt = LLC_data.s(:,:,depth_levels);
    
    %%% Temperature Data
    temp = LLC_data.t(:,:,depth_levels);
    
    %%% Vorticity Data
    vort = LLC_data.Ro(:,:,depth_levels);

    %%% Adding NaNs where there is no data
    temp(salt == 0) = NaN;
    vort(salt == 0) = NaN;
    salt(salt == 0) = NaN;
    
    %%% Combining data into one structure
    LLC = struct('lat', yc, 'lon', xc, 'depth', depth, 'salt', salt, 'temp', temp,'vort', vort);

end