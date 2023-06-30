%%% 
%%% combineSealData.m
%%%
%%% Loads all seal data as a single .mat file
%%% 

countries = { ...
   % 'AUSTRALIA', ...
  'BRAZIL', ...
  % 'CANADA', ...
  % 'CHINA', ...
  % 'DENMARK', ...
  % 'FRANCE', ...
  % 'GERMANY', ...
  % 'NORWAY', ...
  % 'SOUTH-AFRICA', ...
  % 'SWEDEN', ...
  % 'UK', ...
  % 'USA' ...
};

%%% First pass - count number of casts
ncasts = 0;
for m=1:length(countries)
  
  %%% Read all files in the argo data directory for this country and loop
  %%% through them
  dirpath = ['/Volumes/Elements/MEOP-CTD_2021-11-26/',countries{m},'/DATA'];
  ncfiles = dir(dirpath);
  for n=1:length(ncfiles)
    
    %%% Only look at NetCDF files
    fname = ncfiles(n).name;    
    if ((length(fname)>=3) && endsWith(fname,'.nc'))
      
      %%% Check that the latitude lies in the Southern Hemisphere
      %%% N.B. Here we insist that all casts taken by a given seal lie in
      %%% the Southern Hemisphere
      LAT = ncread(fullfile(dirpath,fname),'LATITUDE');
      if (isempty(find(LAT>=0)))
        ncasts = ncasts + length(LAT);
      end
      
    end
    
  end
  
end

%%% Allocate data structure for seal data
sealdata_all = struct(...
  'TAG', cell(1,ncasts), ...
  'LAT', cell(1,ncasts), ...
  'LON', cell(1,ncasts), ...
  'TIME', cell(1,ncasts), ...
  'TIME_QC', cell(1,ncasts), ...
  'PRES', cell(1,ncasts), ...
  'PRES_QC', cell(1,ncasts), ...
  'PRES_RAW', cell(1,ncasts), ...
  'PRES_RAW_QC', cell(1,ncasts), ...
  'PRES_ERROR', cell(1,ncasts), ...
  'TEMP', cell(1,ncasts), ...
  'TEMP_QC', cell(1,ncasts), ...
  'TEMP_RAW', cell(1,ncasts), ...
  'TEMP_RAW_QC', cell(1,ncasts), ...
  'TEMP_ERROR', cell(1,ncasts), ...
  'SALT', cell(1,ncasts), ...
  'SALT_QC', cell(1,ncasts), ...
  'SALT_RAW', cell(1,ncasts), ...
  'SALT_RAW_QC', cell(1,ncasts), ...
  'SALT_ERROR', cell(1,ncasts)...
);

%%
%%% Second pass - extract cast data
castcntr = 0;
for m= 1:length(countries)
  
  %%% Read all files in the argo data directory for this country and loop
  %%% through them
  dirpath = ['/Volumes/Elements/MEOP-CTD_2021-11-26/',countries{m},'/DATA'];
  ncfiles = dir(dirpath);
  for n = 35 %1:length(ncfiles)
    
    %%% Only look at NetCDF files
    fname = ncfiles(n).name;    
    if ((length(fname)>=3) && endsWith(fname,'all_prof.nc'))
      
      %%% Check that the latitude lies in the Southern Hemisphere
      LAT = ncread(fullfile(dirpath,fname),'LATITUDE');
      if (~isempty(find(LAT>=0,1)))
        continue;
      end
      
      %%% Load the rest of the required data. Use adjusted data if
      %%% available.
      LON = ncread(fullfile(dirpath,fname),'LONGITUDE');
      TIME = ncread(fullfile(dirpath,fname),'JULD');
      TIME_QC = ncread(fullfile(dirpath,fname),'JULD_QC');

      %%% Pressure
      PRES = ncread(fullfile(dirpath,fname),'PRES_ADJUSTED');     
      PRES_QC = ncread(fullfile(dirpath,fname),'PRES_ADJUSTED_QC');     
      if (isempty(PRES) || isempty(PRES_QC))
        PRES = ncread(fullfile(dirpath,fname),'PRES');     
        PRES_QC = ncread(fullfile(dirpath,fname),'PRES_QC');   
      end
      PRES_RAW = ncread(fullfile(dirpath,fname),'PRES');
      PRES_RAW_QC = ncread(fullfile(dirpath,fname),'PRES_QC'); 
      PRES_ERROR = ncread(fullfile(dirpath,fname),'PRES_ADJUSTED_ERROR');   

      %%% Temperature
      TEMP = ncread(fullfile(dirpath,fname),'TEMP_ADJUSTED');
      TEMP_QC = ncread(fullfile(dirpath,fname),'TEMP_ADJUSTED_QC');
      if (isempty(TEMP) || isempty(TEMP_QC))
        TEMP = ncread(fullfile(dirpath,fname),'TEMP');
        TEMP_QC = ncread(fullfile(dirpath,fname),'TEMP_QC');
      end
      TEMP_RAW = ncread(fullfile(dirpath,fname),'TEMP');
      TEMP_RAW_QC = ncread(fullfile(dirpath,fname),'TEMP_QC'); 
      TEMP_ERROR = ncread(fullfile(dirpath,fname),'TEMP_ADJUSTED_ERROR');

      %%% Salinity
      SALT = ncread(fullfile(dirpath,fname),'PSAL_ADJUSTED');
      SALT_QC = ncread(fullfile(dirpath,fname),'PSAL_ADJUSTED_QC');
      if (isempty(SALT) || isempty(SALT_QC))
        SALT = ncread(fullfile(dirpath,fname),'PSAL');
        SALT_QC = ncread(fullfile(dirpath,fname),'PSAL_QC');
      end
      SALT_RAW = ncread(fullfile(dirpath,fname),'PSAL');
      SALT_RAW_QC = ncread(fullfile(dirpath,fname),'PSAL_QC'); 
      SALT_ERROR = ncread(fullfile(dirpath,fname),'PSAL_ADJUSTED_ERROR'); 

      %%% Tag number
      TAG = ncread(fullfile(dirpath,fname),'PLATFORM_NUMBER');
      
      %%% Saving all profiles into main structure
      for i=1:length(LAT)
        castcntr = castcntr + 1;
        %disp(castcntr)
        sealdata_all(castcntr).TAG = string(TAG(:,i)');
        sealdata_all(castcntr).LAT = LAT(i);
        sealdata_all(castcntr).LON = LON(i);
        sealdata_all(castcntr).TIME = datevec(datenum('1950-01-01 00:00:00') + TIME(i));
        sealdata_all(castcntr).TIME_QC = TIME_QC(i);

        %%% Pressure
        sealdata_all(castcntr).PRES = PRES(:,i);
        sealdata_all(castcntr).PRES_QC = PRES_QC(:,i);
        sealdata_all(castcntr).PRES_RAW = PRES_RAW(:,i);
        sealdata_all(castcntr).PRES_RAW_QC = PRES_RAW_QC(:,i);
        sealdata_all(castcntr).PRES_ERROR = PRES_ERROR(:,i);

        %%% Temperature
        sealdata_all(castcntr).TEMP = TEMP(:,i);
        sealdata_all(castcntr).TEMP_QC = TEMP_QC(:,i);
        sealdata_all(castcntr).TEMP_RAW = TEMP_RAW(:,i);
        sealdata_all(castcntr).TEMP_RAW_QC = TEMP_RAW_QC(:,i);
        sealdata_all(castcntr).TEMP_ERROR = TEMP_ERROR(:,i);

        %%% Salinity
        sealdata_all(castcntr).SALT = SALT(:,i);
        sealdata_all(castcntr).SALT_QC = SALT_QC(:,i);
        sealdata_all(castcntr).SALT_RAW = SALT_RAW(:,i);
        sealdata_all(castcntr).SALT_RAW_QC = SALT_RAW_QC(:,i);
        sealdata_all(castcntr).SALT_ERROR = SALT_ERROR(:,i);
      end
      
    end
    
  end
  
end

%%% Write to a .mat file
% save('/Volumes/Elements/Raw Data/sealdata_all.mat', 'sealdata_all')

