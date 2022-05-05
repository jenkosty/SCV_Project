%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Filtering Out "Bad" Data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Identifying MEOP Profiles in the Southern Ocean (Below 60S)
idy = find([sealdata_all.LAT] < -60);
SO_sealdata = sealdata_all(idy);
clear idy

%%% Quality Controlling Southern Ocean MEOP Profiles (Identifying profiles
%%% with consistent temperature and salinity quality ratings of 1 - i.e. the best rating)

%%% Extracting rating data
time_qc = NaN(length(SO_sealdata),1);
pres_qc = cell(length(SO_sealdata),1);
salt_qc = cell(length(SO_sealdata),1);
temp_qc = cell(length(SO_sealdata),1);

for i = 1:length(SO_sealdata)
    time_qc(i) = str2double(SO_sealdata(i).TIME_QC);
    pres_qc{i,:} = str2num(SO_sealdata(i).PRES_QC);
    salt_qc{i,:} = str2num(SO_sealdata(i).SALT_QC);
    temp_qc{i,:} = str2num(SO_sealdata(i).TEMP_QC);
end

%%% Isolating profiles with the best data
for i = 1:length(SO_sealdata)
    if time_qc(i) == 1
        good_data{i,:} = find((pres_qc{i,1} == 1) & (salt_qc{i,1} == 1) & (temp_qc{i,1} == 1));
    else 
        good_data{i,:} = [];
    end
end

%%% Creating new structure with QC profiles
clear SO_sealdata_qc
u = 1;
for i = 1:length(SO_sealdata)
    if isempty(good_data{i,1}) == 0
        SO_sealdata_qc(u).lat = SO_sealdata(i).LAT;
        SO_sealdata_qc(u).lon = SO_sealdata(i).LON;
        SO_sealdata_qc(u).time = SO_sealdata(i).TIME;
        SO_sealdata_qc(u).pres = SO_sealdata(i).PRES(good_data{i,1});
        SO_sealdata_qc(u).salt = SO_sealdata(i).SALT(good_data{i,1});
        SO_sealdata_qc(u).temp = SO_sealdata(i).TEMP(good_data{i,1});
        SO_sealdata_qc(u).tag = SO_sealdata(i).TAG;
        u = u + 1;
    end
end

clear time_qc pres_qc salt_qc temp_qc u i SO_sealdata good_data