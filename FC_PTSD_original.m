function [FCtable, h2_amyhpc] = FC_PTSD(dir_data, dir_info, varargin)
%
% Calculates node strength for each channel and subj in dir_data
% Syntax:  
%    FCtable = FC_PTSD(dir_data, dir_info, varargin)
%
% Inputs:
%   dir_data    - input adress of folder containing h2 results from bids pipeline with sub-code at beginning 1 file per subject with subj code, all subjects in the same folder 
%   dir_info    - dir of info table with info on each channel
% 
% Varargin:
%   thr_h2, thr_lag - optional thershold for h2 (thr_h2=0 for strength, >0 for degrees) and lag. Default: thr_h2 = 0; thr_lag = 0
%   norm            - 1 if you want to nromalise each subj by number of channels (Default= 1)
%   graph           - 1 (default) to show connectivity matrices of zvalues compared to baseline chosen in base, 0 to not show any graph
%   roi             - "all" (default, h2 between all channels) or "EZ", "EZPZ", "NI" to calculate node strength only in subset of EZ channels (or non-inv channels)
%
% Output:
%   FCtable     - 1 row per channel (and per subject), mean node stength OUT and TOT
%
% Required functions: 
%   - ins_countlinks.m,
% Authors: Sara Simula (original: Aug 2023. Last version: )

% Example/debug:
% clear
% dir_data = "\\dynaserv\meg\nicolas\PTSD\analysis\FC\broad";
% dir_info = "\\dynaserv\meg\nicolas\PTSD\PTSD_meltinfo.xlsx";

% 1. Optional variables: default values
thr_h2 = 0;  % thr_h2=0 for strength, >0 for degrees. Default: strength
thr_lag = 0; % value in ms, optional input. Default: 0 % TO CHECK???
norm =1; %% normalise strength by number of channels
graph = 0;
roi = "all";

for ii = 1:2:nargin-3
        if strcmp('thr_h2', varargin{ii})
            thr_h2 = varargin{ii+1}; 
        elseif strcmp('thr_lag', varargin{ii})
            thr_lag = varargin{ii+1};
        elseif strcmp('norm', varargin{ii})
            norm = varargin{ii+1};
        elseif strcmp('graph', varargin{ii})
            graph = varargin{ii+1};
        elseif strcmp('roi', varargin{ii})
            roi = varargin{ii+1};
        end
end

% 2. set directory for h2 files to analyse
cd(dir_data)
myfiles = dir('*.mat');

% 3. Calculate the strength/degrees for each channel and each patient and each channel
subj_info= [];
widedata=[];
FCtable = [];
h2_amyhpc = struct();

chan_infoAll = readtable(dir_info);
varnames = chan_infoAll.Properties.VariableNames;

for index=1:length(myfiles)
    rest = load(myfiles(index).name);
    subj = string(extractBetween(myfiles(index).name, 'sub-', '_ses'));
    channels = string(rest.electrode_names);
   
    % read and add info on subject and channels
    
    chan_infoSub = chan_infoAll(chan_infoAll.subject == subj, :);

    b = table2cell(chan_infoSub);
    subj_info = string(b);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3.2. Reduce matrix of h2 and lag in case "roi" is different from "all"
    if ~strcmp("all", roi)
    switch roi
        case "EZ"
            good_labels = "EZ";
        case "EZPZ"
            good_labels = '[EP.]Z'; %any character in [] followed by Z (matches "EZ" or "PZ")
        case "NI"
            good_labels = "NI";
        otherwise
            disp("incorrect roi input, possible values are: EZ, EZPZ or NI")
    end
            idx = ~cellfun(@isempty,regexp(subj_info(:,TOCHANGE), good_labels));
            rest.aw_h2 = rest.aw_h2(idx,idx,:);
            rest.aw_lag = rest.aw_lag(idx,idx,:);
            rest.electrode_names = rest.electrode_names(idx);
            
%         channels = string(rest.electrode_names); %update the selected channels
        subj_info = subj_info(idx,:);
        clearvars idx
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 3.3 calculate the in, out, tot strengths/degrees using the func countlinks in graphcompare
    if thr_h2==0 && norm==1
        
        [linksrest,~]=ins_countlinks(rest,thr_h2);        % linkspre contains 1 row per channel, 1 column per window of h2
        OUTrest = mean(linksrest.outstrength_norm,2);  % median across windows of calculation h2
        TOTrest = mean(linksrest.totstrength_norm,2);  % median across windows of calculation h2
  
    else %degrees YET TO BE NORMALIZED!!!!! NOT TO USE FOR NOW
        [linksrest,~]=ins_countlinks(rest,thr_h2);
        OUTrest = mean(linksrest.outdegree,2);
        TOTrest = mean(linksrest.totdegree,2);

    end

    tmp = [string(rest.electrode_names)', OUTrest, TOTrest];
    for chan = 1:length(tmp(:,1))
        match = sum(contains(subj_info(:,2), tmp(chan,1)));
        if match % does not always contain all elecs in bipolar selection montage cause some of them were WM or excluded from labels table
            subj_info(contains(subj_info(:,2), tmp(chan,1)), 18:19) = tmp(chan,2:3);
        else 
            subj_info(contains(subj_info(:,2), tmp(chan,1)), 18:19) = NaN;
        end
    end
   
    clearvars tmp match
    subj_info = subj_info(~ismissing(subj_info(:,18)),:);
        
    % now we have one strength value per node, and we need to also make a
    % selection of the areas we are interested in to average values inside

    % regions to select: amygdala, hpc, orbito-frontal, insula (ant),
    % cingulate gyrus --> not all subjects have at least one electrodes in
    % each of these regions, so we need to average them all into the
    % category "involved" vs "non-involed"

    % a. look for involved regions labels
    inv_labels = ["Amygdala", "Hippocampus", "Orbito-frontal-cortex", "Insula", "Rhinal-cortex"];
    ni_labels  = ["ITS-anterior", "Thalamus", "F2", "STS-anterior", "T1-lateral-posterior", "T3", "Precuneus", "T2-posterior"]; % ok

    chan_inv    = chan_infoSub.channel(contains(chan_infoSub.brain_area, inv_labels));
    chan_ni     = chan_infoSub.channel(contains(chan_infoSub.brain_area, ni_labels));

    idx_inv     = contains(rest.electrode_names, chan_inv);
    idx_ni      = contains(rest.electrode_names, chan_ni);
    idx_invni   = logical(idx_ni + idx_inv); % involved and non-involved selected regions

    %% calculate node strength only between involved nodes
    rest_inv.aw_h2 = rest.aw_h2(idx_inv,idx_inv,:);
    rest_inv.aw_lag = rest.aw_lag(idx_inv,idx_inv,:);
    rest_inv.electrode_names = rest.electrode_names(idx_inv);

    if thr_h2==0 && norm==1
        [linksrestINV,~]=ins_countlinks(rest_inv,thr_h2);        % linkspre contains 1 row per channel, 1 column per window of h2
        OUTrest_inv = mean(linksrestINV.outstrength_norm,2);
        TOTrest_inv = mean(linksrestINV.totstrength_norm,2);
   
        tmp = [string(rest_inv.electrode_names)', OUTrest_inv, TOTrest_inv];
        for chan = 1:length(tmp(:,1))
            match = sum(contains(subj_info(:,2), tmp(chan,1)));
            if match % does not always contain all elecs in bipolar selection montage cause some of them were WM or excluded from labels table
                subj_info(contains(subj_info(:,2), tmp(chan,1)), 20:21) = tmp(chan,2:3);
            else 
                subj_info(contains(subj_info(:,2), tmp(chan,1)), 20:21) = NaN;
            end
        end

    else
    error("edit function to calculate degrees")
    end

    clearvars tmp match

    %% calculate node strength only between non-involved nodes
    rest_ni.aw_h2 = rest.aw_h2(idx_ni,idx_ni,:);
    rest_ni.aw_lag = rest.aw_lag(idx_ni,idx_ni,:);
    rest_ni.electrode_names = rest.electrode_names(idx_ni);
    
    if thr_h2==0 && norm==1
        [linksrestNI,~]=ins_countlinks(rest_ni,thr_h2);        % linkspre contains 1 row per channel, 1 column per window of h2
        OUTrest_ni = mean(linksrestNI.outstrength_norm,2);
        TOTrest_ni = mean(linksrestNI.totstrength_norm,2);

        tmp = [string(rest_ni.electrode_names)', OUTrest_ni, TOTrest_ni];
        for chan = 1:length(tmp(:,1))
            match = sum(contains(subj_info(:,2), tmp(chan,1)));
            if match % does not always contain all elecs in bipolar selection montage cause some of them were WM or excluded from labels table
                subj_info(contains(subj_info(:,2), tmp(chan,1)), 22:23) = tmp(chan,2:3);
            else 
                subj_info(contains(subj_info(:,2), tmp(chan,1)), 22:23) = NaN;
            end
        end

    else
    error("edit function to calculate degrees")
    end

    clearvars tmp match

%      %% calculate node strength between amygdala and hpc   
%     
%     left_amy  = chan_infoSub.channel(contains(chan_infoSub.brain_area, "Left-Amygdala"));
%     left_hpc  = chan_infoSub.channel(contains(chan_infoSub.brain_area, ["Left-Hippocampus-anterior", "Left-Hippocampus-posterior"]));
%     right_amy = chan_infoSub.channel(contains(chan_infoSub.brain_area, "Right-Amygdala"));
%     right_hpc = chan_infoSub.channel(contains(chan_infoSub.brain_area, ["Right-Hippocampus-anterior", "Right-Hippocampus-posterior"]));
% 
%     sel = {left_amy, left_hpc, right_amy, right_hpc};
%     for x = 1:2:size(sel,2)
%         h2_amyhpc.subject(index)  = subj;
%         if ~isempty(sel{x}) && ~isempty(sel{x+1})
%             idx_amy  = contains(rest.electrode_names, sel{x});
%             idx_hpc  = contains(rest.electrode_names, sel{x+1});
%             [h2_array,~]=ins_findmaxh2(rest); % whole matrix
% 
%             linksrest_amyhpc = struct();
%         
%             % sum on columns (in-strength)
%             linksrest_amy_hpc.instrength  = squeeze(sum(h2_array(idx_hpc, idx_amy,:), 1));
%             linksrest_hpc_amy.instrength  = squeeze(sum(h2_array(idx_amy, idx_hpc,:), 1));
% 
% %             % sum on rows (out-strength)
% %             linksrest_amy_hpc.outstrength = squeeze(sum(h2_array(idx_amy, idx_hpc,:), 2));
% %             linksrest_hpc_amy.outstrength = squeeze(sum(h2_array(idx_hpc, idx_amy,:), 2));
% %             
% %             numchan_amy = sum(idx_hpc);
% %             linksrest_amy_hpc.instrength_norm  = linksrest_amy_hpc.instrength/(numchan_amy-1);
% %             linksrest_amy_hpc.outstrength_norm = linksrest_amy_hpc.outstrength/(numchan_amy-1);
% %             linksrest_amy_hpc.totstrength_norm = linksrest_amy_hpc.instrength_norm + linksrest_amy_hpc.outstrength_norm;
% %             OUTrest_inv_ni = mean(linksrest_amy_hpc.outstrength_norm,2);  % median across windows of calculation h2
% %             TOTrest_inv_ni = mean(linksrest_amy_hpc.totstrength_norm,2);
% 
% % symmetrize matrix
%             for ii = 1 : size(h2_array,1)
%                     for jj = 1 : size(h2_array,2)
%                         if abs(h2_array(jj,ii))>=abs(h2_array(ii,jj))
%                             h2_array(ii,jj) = h2_array(jj,ii);
%                         end
%                     end
%             end
%             h2_amyhpc.h2matrix(index) = {mean(h2_array(idx_amy, idx_hpc, :),3)}; 
%             linksrest_amy_hpc.ychannels = rest.electrode_names(idx_amy);
%             linksrest_amy_hpc.xchannels = rest.electrode_names(idx_hpc);
%             %%% TO CHANGE !!! match anat label swith channels to show in
%             %%% matrix h2 graphs
% %             linksrest_amy_hpc.channels = rest.electrode_names(idx_amy);
% %             tmp = string(linksrest_amy_hpc.channels)';
% %             for chan = 1:length(tmp(:,1))
% %                 match = sum(contains(subj_info(:,2), tmp(chan,1)));
% %                 if match % does not always contain all elecs in bipolar selection montage cause some of them were WM or excluded from labels table
% %                     h2_amyhpc.anat = subj_info(contains(subj_info(:,2), tmp(chan,1)), 3);
% %                 else 
% %                     subj_info(contains(subj_info(:,2), tmp(chan,1)), 24:25) = NaN;
% %                 end
% %             end
% 
% %             h2_amyhpc.brain_areas     = 
% 
%             clearvars tmp match h2_array
%        end
%     end
    %% calculate node strength between inv and non-inv channels (not inside)   
    if thr_h2==0 && norm==1
        [h2_array,~]=ins_findmaxh2(rest); 
        linksrest_inv_ni = struct();
        linksrest_ni_inv = struct();
    
        % sum on columns (in-strength)
        linksrest_ni_inv.instrength  = squeeze(sum(h2_array(idx_inv,idx_ni,:), 1));
        linksrest_inv_ni.instrength  = squeeze(sum(h2_array(idx_ni,idx_inv,:), 1));
    
        % sum on rows (out-strength)
        linksrest_inv_ni.outstrength  = squeeze(sum(h2_array(idx_inv,idx_ni,:), 2));
        linksrest_ni_inv.outstrength  = squeeze(sum(h2_array(idx_ni,idx_inv,:), 2));

%% inv --> ni 
        numchan_ni = sum(idx_ni);
        linksrest_inv_ni.instrength_norm  = linksrest_inv_ni.instrength/(numchan_ni-1);
        linksrest_inv_ni.outstrength_norm = linksrest_inv_ni.outstrength/(numchan_ni-1);
        linksrest_inv_ni.totstrength_norm = linksrest_inv_ni.instrength_norm + linksrest_inv_ni.outstrength_norm;
        OUTrest_inv_ni = mean(linksrest_inv_ni.outstrength_norm,2);  % median across windows of calculation h2
        TOTrest_inv_ni = mean(linksrest_inv_ni.totstrength_norm,2);
        linksrest_inv_ni.channels = rest.electrode_names(idx_inv);

        tmp = [string(linksrest_inv_ni.channels)', OUTrest_inv_ni, TOTrest_inv_ni];
        for chan = 1:length(tmp(:,1))
            match = sum(contains(subj_info(:,2), tmp(chan,1)));
            if match % does not always contain all elecs in bipolar selection montage cause some of them were WM or excluded from labels table
                subj_info(contains(subj_info(:,2), tmp(chan,1)), 24:25) = tmp(chan,2:3);
            else 
                subj_info(contains(subj_info(:,2), tmp(chan,1)), 24:25) = NaN;
            end
        end

        clearvars tmp match

%% ni --> inv    
        numchan_inv = sum(idx_inv);
        linksrest_ni_inv.instrength_norm  = linksrest_ni_inv.instrength/(numchan_inv-1);
        linksrest_ni_inv.outstrength_norm = linksrest_ni_inv.outstrength/(numchan_inv-1);
        linksrest_ni_inv.totstrength_norm = linksrest_ni_inv.instrength_norm + linksrest_ni_inv.outstrength_norm;
        OUTrest_ni_inv = mean(linksrest_ni_inv.outstrength_norm,2);  % median across windows of calculation h2
        TOTrest_ni_inv = mean(linksrest_ni_inv.totstrength_norm,2);
        linksrest_ni_inv.channels = rest.electrode_names(idx_ni);

        tmp = [string(linksrest_ni_inv.channels)', OUTrest_ni_inv, TOTrest_ni_inv];
        for chan = 1:length(tmp(:,1))
            match = sum(contains(subj_info(:,2), tmp(chan,1)));
            if match % does not always contain all elecs in bipolar selection montage cause some of them were WM or excluded from labels table
                subj_info(contains(subj_info(:,2), tmp(chan,1)), 26:27) = tmp(chan,2:3);
            else 
                subj_info(contains(subj_info(:,2), tmp(chan,1)), 26:27) = NaN;
            end
        end

    else
        error("edit function to calculate degrees")
    end
  
    clearvars tmp match h2_array
    
    FCtable = [FCtable; subj_info];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% graphs
   if graph
       lim_max = max(max(max(mean(rest_inv.aw_h2,3))), max(max(mean(rest_ni.aw_h2,3))));
       
       color = viridis;
       figure('Name', subj)
    
       subplot(2,1,1)
       imagesc(mean(rest_inv.aw_h2,3))
       colormap(color)
       title("Involved regions")
       colorbar
       xticks(1:length(string(rest_inv.electrode_names)));
       xticklabels(rest_inv.electrode_names);
       xtickangle(60)
       yticks(1:length(string(rest_inv.electrode_names)));
       yticklabels(rest_inv.electrode_names);
       ylabel = "channel";
       xlabel = "channel";
       clim([0, lim_max])
    
       subplot(2,1,2)
       imagesc(mean(rest_ni.aw_h2,3))
       colormap(color)
       title("Non-involved regions")
       colorbar
       xticks(1:length(string(rest_ni.electrode_names)));
       xticklabels(string(rest_ni.electrode_names));
       xtickangle(60)
       yticks(1:length(string(rest_ni.electrode_names)));
       yticklabels(string(rest_ni.electrode_names));
       ylabel = "channel";
       xlabel = "channel";
       clim([0, lim_max])
    
       print(strcat('\\dynaserv\meg\nicolas\PTSD\analysis\FC\figures\', subj, 'inv_ni_SAMESCALE'), '-dpdf', '-fillpage')
    
   end   
    clearvars  subj channels check linksrest OUTrest TOTrest subj_info min max rest h2_array
end

%% save the table

    varnames = {varnames{:}, 'OUTrest_all', 'TOTrest_all', 'OUTrest_inv_inv','TOTrest_inv_inv','OUTrest_ni_ni', 'TOTrest_ni_ni', 'OUTrest_inv_ni', 'TOTrest_inv_ni', 'OUTrest_ni_inv', 'TOTrest_ni_inv'};
    FCtable = array2table(FCtable, 'VariableNames', varnames);

    
    FCtable.involved(~ismissing(FCtable.TOTrest_inv_inv)) = "inv";
    FCtable.involved(~ismissing(FCtable.TOTrest_ni_ni)) = "non-inv";
   

end   
