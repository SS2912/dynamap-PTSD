function [results_all, delphos_melttable] = delphos_PTSD(input_folder, mtgs_folder, sections, varargin)

%% Analyze DELPHOS .mat results and create a table with all results for all subjects
% 
% Syntax:  
%    [results_all, delphos_melttable] = delphos_PTSD(input_folder, mtgs, sections, varargin)
%    
% Input format: 
%    1 folder per section, 1 file per subject with subj code, all subjects in the same folder (ex: sub-0b12810f4f6d_postSTIM_B.mat)
%
% Inputs:
%   input_folder    - input adress of folder containing all different sections' folders (ex:  "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\Delphos\mat_threshold80";)
%   mtgs_folder     - adress of raw montages folders (ex: "\\dynaserv\Galvani_ps2\montages\delphos\*.mtg");
%   sections        - cell array of all input sections (==folders). Ex: {"baseline", "sham", "postA", "postB", "postall"}
% 
% Outputs:
%   mean_window - [Nch x Nw] Mean rate at each window and for each channel
%   spikes      - All the detections (timepoints) for each channel
%
% Required functions: 
%   - dyn_read_mtg.m,
%   - mean_rate_delphos.m, 
%
% Authors: Sara Simula (original: taken from delphos_PS2, last version: Sept 2023)


%% 1a. prepare to read the .mat files in the folders 
% EXAMPLE input directory (with delphos results): 

% clear

% input_folder  = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\Delphos\mat_threshold80";
% mtgs_folder = "\\dynaserv\Galvani_ps2\montages\delphos\*.mtg";
% sections = {"baseline", "sham", "postA", "postB", "postall"};

%% 1b. read the montages of each patient, in which you will insert the calculated SR 
% (step which is necessary here cause delphos only keeps the non-zero detections channels in its results tables)

simEFdir = 0;
tabledir = "\\dynaserv\meg\nicolas\PTSD\PTSD_meltinfo.xlsx";

for ii = 1:2:nargin-3
        if strcmp('tabledir', varargin{ii})
            tabledir = varargin{ii+1}; 
        elseif strcmp('simEFdir', varargin{ii})
            simEFdir = varargin{ii+1};
        end
end

mtgs = dir(mtgs_folder);
SR_diff_table = table();
results_all = struct();
delphos_melttable = table();
info_tot = {};
events_all_str = ["Alpha","Beta","Fast Ripple","Gamma","Ripple","Spike"];
events_all = cellstr(["Alpha","Beta","Fast Ripple","Gamma","Ripple","Spike"]);

detections = struct();

for i = 1:size(mtgs,1) % for each subject
    delphos_subtab = table();
    sub = string(extractBetween(mtgs(i).name, '-', '_'));
    montage_struct = dyn_read_mtg(strcat(mtgs(i).folder, '\', mtgs(i).name));
    results_all(i).subj = sub(1);

    % create table with bipolar channels and sub code where you'll insert SR info
    for j = 1:length(montage_struct)
        delphos_subtab.subj(j)= sub(1);  % sub(1) because in old file names there were different string delimited by - and _ -> I only took the first one (subj code)
        delphos_subtab.chan(j) = strcat(string(montage_struct(j).name), '-', string(montage_struct(j).reference));
    end
    

    %% 2. find the same subject in delphos input folder and read + write the results into SR_diff_table
    % 2a. for each section(folder), find the current subject and load their delphos results
    idxev = 0;

    for ii = 1:length(sections)
        curr_folder = strcat(input_folder, "\", sections{ii});
        cd(curr_folder)
        myfiles  = dir(strcat(curr_folder, "\*.mat"));
        match    = contains(string({myfiles.name}), sub(1), 'IgnoreCase', true);
       
        if sum(match)
           load(myfiles(match).name);
           windowsize    = 60; % now working only with windowsize of 60 sec (1 min), need to change mean_rate_delphos.m to make it work with different window sizes
           events_labels = string({results.markers.label});
           events_list   = cellstr(unique(events_labels));
           if size(results.labels,1)<size(results.labels,2)
              chan_list   = string(results.labels)';
           else 
              chan_list = string(results.labels);
           end
           
           % create big table where to insert matching channels from delphos results to raw channel list
           subtab_event         = table();           
           subtab_event.chan    = chan_list;
           subtab_event.section = repelem(sections{ii}, length(chan_list), 1);
           colstart             = size(subtab_event,2);

    % 2b. for each subject, and for each section, calculate spike/osc rate via the "mean_rate_delphos.m" function
        for ev = 1:length(events_all)

           event = events_all{ev};
           detections(ev +idxev).subj                 = sub(1);
           detections(ev +idxev).section              = sections{ii};
           detections(ev +idxev).event                = event; 
           detections(ev +idxev).channels             = chan_list;  
           detections(ev +idxev).windowsize           = windowsize;
           [detections(ev +idxev).meanrate_window, ~] = mean_rate_delphos(results, windowsize, event); % mean rate per window (mean_window) and timepoints of all detections for each channel

           subtab_event(:,ev+colstart) = table(mean(detections(ev+idxev).meanrate_window, 2)*(60/windowsize));
         
        end  
       
        idxev = 1+idxev+length(sections);
        results_all(i).detections = detections;
        
        startcol = 3;

        % now match the channels in delphos with the ones in raw mtg
        subtab_event.Properties.VariableNames = ["chan", "section", events_all{:}];
        for chan=1:length(delphos_subtab.chan)
            chan_match = strcmp(string(delphos_subtab.chan(chan)), string(subtab_event.chan));
            
            if ~isempty(subtab_event(chan_match, :))
                delphos_subtab(chan, startcol:size(subtab_event,2)+startcol-1) = subtab_event(chan_match, :);
            else
                delphos_subtab(chan, startcol:size(subtab_event,2)+startcol-1) = {NaN};
            end
        end

            info = readinfo_PTSD(sub, delphos_subtab.chan, tabledir);
            varnames_info = info.Properties.VariableNames;
            info_tot = [info_tot; table2cell(info)];

        delphos_melttable = [delphos_melttable; delphos_subtab];
        end
        
    end

    clearvars delphos_subtab detections colstart sub results info
end

%delphos_melttable.Properties.VariableNames = ["subject", "chanRaw", "chanDelphos", "section", events_list{:}];
delphos_melttable.Properties.VariableNames = ["subject_raw", "chanRaw", "chanDelphos", "section", events_all_str];

info_tot = cell2table(info_tot, 'VariableNames', varnames_info);
delphos_melttable = [delphos_melttable info_tot];
end

% In delphos github:
% mrk_chn_bln = strcmp({detection_markers_chn(:).label}, handles.markers2display{mrk_idx});
% handles.rate_markers(i,mrk_idx) = sum(mrk_chn_bln)/(handles.duration/60);
            