function infotable_melt = combineinfo_PTSD(dir_gardel, stringToBeFound, clinical_info)

% to combine: anat labels in gardel folder (per chan), clinical info from
% Lisa's table (PTSD score, anxiety, depression, etiology, bilateral, age,
% sex)

% find anatomical labels in bids and put them in one structure that will be used afterwards for combining chan info in each function

subdirs = strcat(dir_gardel, '\**\ses-postimp\anat\*.xlsx');
files = dir(subdirs);
infotable = {};
idx=0;

for i=1:length(files)
    if contains(files(i).name, stringToBeFound)
        idx=idx+1;
        infotable{idx,1} = extractBefore(files(i).name,"_");
        cd(files(i).folder);
        temp = readtable(files(i).name);
        temp.channels = string(temp.channels);        
        infotable{idx,2} = temp(:,1:2);
        clearvars temp
    end
end

colnames = clinical_info.Properties.VariableNames;
all = table2cell(clinical_info);

for sub=1:length(infotable(:,1))
    code = string(infotable(sub,1));
    idxsub = strcmp(all(:,1), code);
    all(idxsub, 17) = infotable(sub,2);
    clearvars idxsub code
end


%% create a big table with all the channels and subject code for each channel
row = 1;
infotable_melt = {};
for i = 1:length(all)
    l = size(all{i,17}, 1);
    infotable_melt(row:row+l-1,1) = repelem(all(i,1), l)';
    infotable_melt(row:row+l-1,2:3) = table2cell(all{i,17});
    infotable_melt(row:row+l-1,4:18) = repmat(all(i,2:16), l,1);
    row = row+l;
end

infotable_melt = cell2table(infotable_melt, 'VariableNames',{'subject','channel', 'brain_area', colnames{2:end}});


end