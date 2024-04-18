function [FCtable, h2_amyhpc] = FC_PTSD(dir_data, dir_info, category, varargin)
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
%   meth            - how to obtain signle value of h2 for each region ("median", "mean". default = "median")
%
% Output:
%   FCtable     - 1 row per channel (and per subject), mean node stength OUT and TOT
%
% Required functions: 
%   - ins_countlinks.m,
%   - FCmatrix_no_rep.m (TOOLBOX_Sara, Onedrive Phd projects)
%   - optional: reduce_h2matrix.m (if intraconnectivity only in "roi" regions)
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
        elseif strcmp('meth', varargin{ii})
            meth = varargin{ii+1};
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

    % reduce subj_info to only the channels in the montage
    anat_column = contains(varnames,"brain_area");
    channel_col = contains(varnames,"channel");

    channels = upper(channels);
    chan_subj = upper(subj_info(:, channel_col));
    idx = zeros(length(chan_subj), 1);    
    for chan=1:length(channels)
        idx = idx + contains(chan_subj, channels(chan));
    end
    idx = logical(idx);
    subj_info_new = subj_info(idx,:);
    subj_info = [];
    subj_info = subj_info_new;

    %% 3.1 Only take 1 value per each structure (varragin for max, median or mean

    subj_anat = subj_info(:, anat_column);
    subj_chan = subj_info(:, channel_col);
   
    if exist("meth", "var")
        [h2_mean, h2_median] = FCmatrix_no_rep(rest.aw_h2, rest.aw_lag, subj_anat, subj_chan, channels);
       
        % substitute new matrices and chan names instead of old rest.aw_h2 etc
        rest.aw_h2 = [];
        rest.aw_lag = [];
        rest.electrode_names = {};
    
        switch meth
            case "mean"
                rest.aw_h2  = h2_mean.h2;
                rest.aw_lag = h2_mean.lag;
            case "median"
                rest.aw_h2  = h2_median.h2;
                rest.aw_lag = h2_median.lag;
        end
    
        rest.electrode_names = h2_mean.chan; %mean or median channels are the same
        channels = upper(string(rest.electrode_names));

    end
    %% 3.2. Reduce matrix of h2 and lag in case "roi" is different from "all"
%     roi_col = ????
    if ~strcmp("all", roi)
        [h2_roi, subj_info_new] = reduce_h2matrix(rest, roi, subj_info, roi_col);
        rest = [];
        rest = h2_roi;
        clearvars subj_info
        subj_info = subj_info_new;
    end
    
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

    tmp = [upper(string(rest.electrode_names))', OUTrest, TOTrest];
    subj_channels = upper(subj_info(:,channel_col));
    for chan = 1:length(tmp(:,1))
        match = sum(contains(subj_channels, tmp(chan,1)));
        if match % does not always contain all elecs in bipolar selection montage cause some of them were WM or excluded from labels table
            subj_info(contains(subj_channels, tmp(chan,1)), 18:19) = tmp(chan,2:3);
        else 
            subj_info(contains(subj_channels, tmp(chan,1)), 18:19) = NaN;
        end
    end
   
    clearvars tmp match
    subj_info = subj_info(~ismissing(subj_info(:,18)),:);
    subj_channels = upper(subj_info(:,channel_col));

        
    % now we have one strength value per node, and we need to also make a
    % selection of the areas we are interested in to average values inside

    % regions to select: amygdala, hpc, orbito-frontal, insula (ant),
    % cingulate gyrus --> not all subjects have at least one electrodes in
    % each of these regions, so we need to average them all into the
    % category "involved" vs "non-involed"

    % a. look for involved regions labels or ipsi/contra
    side = string(unique(chan_infoSub.EpilepsySide));
    if category == "ipsi_contra"
        
        if side == "R"
            ipsi_labels = "Right";  % ipsi, called inv only to adapt to the script and choose
            contra_labels = "Left";   % contra side
        elseif side == "L"
            ipsi_labels = "Left";
            contra_labels = "Right";
        else                            % bilateral epi
            ipsi_labels = ["Left", "Right"];
            contra_labels = "nothing";  % need to have no match with contra when epi is bilateral 
        end

    elseif category == "inv_ni"
            ipsi_labels = ["Amygdala", "Hippocampus"]; % calculate only between amygdala and hippocampus
            %ipsi_labels = ["Amygdala", "Hippocampus", "Orbito-frontal-cortex", "Insula", "Rhinal-cortex", "Anterior-cingulate-cortex"];
            contra_labels  = ["ITS-anterior", "Thalamus", "F2", "STS-anterior", "T1-lateral-posterior", "T3", "Precuneus", "T2-posterior"]; % ok
    elseif category == "DMN_CEN"
            ipsi_labels = ["Parahippocampal", "Temporal", "Posterior-cingulate-cortex", "Precuneus", "prefrontal", "Angular"]; % DMN
            contra_labels  = ["F2", "SFS"]; % CEN

    else 
        error("choose between ""ipsi/contra"" and ""inv/ni""");
    end

    chan_inv    = chan_infoSub.channel(contains(chan_infoSub.brain_area, ipsi_labels));
    chan_ni     = chan_infoSub.channel(contains(chan_infoSub.brain_area, contra_labels));

    idx_inv     = contains(rest.electrode_names, upper(string(chan_inv)));
    idx_ni      = contains(rest.electrode_names, upper(string(chan_ni)));
    idx_invni   = logical(idx_ni + idx_inv); % involved and non-involved selected regions

    %% calculate node strength only between involved nodes
    rest_inv.aw_h2 = rest.aw_h2(idx_inv,idx_inv,:);
    rest_inv.aw_lag = rest.aw_lag(idx_inv,idx_inv,:);
    rest_inv.electrode_names = rest.electrode_names(idx_inv);

    if thr_h2==0 && norm==1
        [linksrestINV,~]=ins_countlinks(rest_inv,thr_h2);        % linkspre contains 1 row per channel, 1 column per window of h2
        OUTrest_inv = mean(linksrestINV.outstrength_norm,2);
        TOTrest_inv = mean(linksrestINV.totstrength_norm,2);
   
        if length(string(rest_inv.electrode_names))~=1
           tmp = [upper(string(rest_inv.electrode_names))', OUTrest_inv, TOTrest_inv];
        else
            tmp = [upper(string(rest_inv.electrode_names))', NaN(1), NaN(1)];
        end
        for chan = 1:length(tmp(:,1))
            match = sum(contains(subj_channels, tmp(chan,1)));
            if match % does not always contain all elecs in bipolar selection montage cause some of them were WM or excluded from labels table
                subj_info(contains(subj_channels, tmp(chan,1)), 20:21) = tmp(chan,2:3);
            else 
                subj_info(contains(subj_channels, tmp(chan,1)), 20:21) = NaN;
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
        
        if length(string(rest_ni.electrode_names))~=1
            tmp = [upper(string(rest_ni.electrode_names))', OUTrest_ni, TOTrest_ni];
        else
            tmp = [upper(string(rest_ni.electrode_names))', NaN(1), NaN(1)];
        end

        for chan = 1:length(tmp(:,1))
            match = sum(contains(subj_channels, tmp(chan,1)));
            if match % does not always contain all elecs in bipolar selection montage cause some of them were WM or excluded from labels table
                subj_info(contains(subj_channels, tmp(chan,1)), 22:23) = tmp(chan,2:3);
            else 
                subj_info(contains(subj_channels, tmp(chan,1)), 22:23) = NaN;
            end
        end

    else
    error("edit function to calculate degrees")
    end

    clearvars tmp match


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

         if length(string(linksrest_inv_ni.channels))~=1
            tmp = [upper(string(linksrest_inv_ni.channels))', OUTrest_inv_ni, TOTrest_inv_ni];       
         else
            tmp = [upper(string(linksrest_inv_ni.channels))', NaN(1), NaN(1)];
         end
         
        if isempty(tmp)
            subj_info(:, 24:25) = NaN;
        else
        for chan = 1:length(tmp(:,1))
            match = sum(contains(subj_channels, tmp(chan,1)));
            if match % does not always contain all elecs in bipolar selection montage cause some of them were WM or excluded from labels table
                subj_info(contains(subj_channels, tmp(chan,1)), 24:25) = tmp(chan,2:3);
            else 
                subj_info(contains(subj_channels, tmp(chan,1)), 24:25) = NaN;
            end
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

        if length(string(linksrest_ni_inv.channels))~=1
            tmp = [upper(string(linksrest_ni_inv.channels))', OUTrest_ni_inv, TOTrest_ni_inv];       
         else
            tmp = [upper(string(linksrest_ni_inv.channels))', NaN(1), NaN(1)];
         end        
        if isempty(tmp)
            subj_info(:, 26:27) = NaN;
        else
            for chan = 1:length(tmp(:,1))
                match = sum(contains(subj_channels, tmp(chan,1)));
                if match 
                    subj_info(contains(subj_channels, tmp(chan,1)), 26:27) = tmp(chan,2:3);
                else
                    subj_info(contains(subj_channels, tmp(chan,1)), 26:27) = NaN;
                end
            end

        end
  
    clearvars tmp match h2_array
    end
    FCtable = [FCtable; subj_info];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% graphs
   channels_inv = rest_inv.electrode_names;
   channels_ni  = rest_ni.electrode_names;

    for j=1:length(channels_inv)
       match = strcmp(upper(subj_info(:,channel_col)), channels_inv(j));
       if sum(match)
           anat_inv(j) = subj_info(match, strcmp(varnames,"brain_area"));
       else
           anat_inv(j) = "missing in excel";
       end
    end

    clearvars match
    for j=1:length(channels_ni)
        match = strcmp(upper(subj_info(:,channel_col)), channels_ni(j));
        if sum(match)
           anat_ni(j) = subj_info(match, strcmp(varnames,"brain_area"));
        else
           anat_ni(j) = "missing in excel";
        end
    end

   if graph
       if ~isempty(rest_inv.aw_h2) && ~isempty(rest_ni.aw_h2)
        lim_max = max(max(max(mean(rest_inv.aw_h2,3))), max(max(mean(rest_ni.aw_h2,3))));
       elseif ~isempty(rest_inv.aw_h2) && isempty(rest_ni.aw_h2)
           lim_max = max(max(mean(rest_inv.aw_h2,3)));
       else
           lim_max = max(max(mean(rest_ni.aw_h2,3)));
       end
       
       color = viridis;
       figure('Name', strcat(subj, "- side: ", side))
    
       subplot(2,1,1)
       imagesc(mean(rest_inv.aw_h2,3))
       colormap(color)
       switch category
           case "ipsi/contra"
               title("Ipsi")
           case "inv/ni"
               title("Involved regions")
       end
       colorbar
       xticks(1:length(string(rest_inv.electrode_names)));
       xticklabels(rest_inv.electrode_names);
       xtickangle(60)
       yticks(1:length(string(rest_inv.electrode_names)));
       yticklabels(anat_inv);
       ylabel = "channel";
       xlabel = "channel";
       clim([0, lim_max])
    
       subplot(2,1,2)
       imagesc(mean(rest_ni.aw_h2,3))
       colormap(color)
       switch category
           case "ipsi/contra"
               title("Contra")
           case "inv/ni"
                title("Non-involved regions")
       end
       colorbar
       xticks(1:length(string(rest_ni.electrode_names)));
       xticklabels(string(rest_ni.electrode_names));
       xtickangle(60)
       yticks(1:length(string(rest_ni.electrode_names)));
       yticklabels(anat_ni);
       ylabel = "channel";
       xlabel = "channel";
       clim([0, lim_max])
       
       switch category
           case "ipsi/contra"
               print(strcat('\\dynaserv\meg\nicolas\PTSD\analysis\FC\figures\', subj, 'ipsi_contra_SAMESCALE'), '-dpng');
%                print(strcat('\\dynaserv\meg\nicolas\PTSD\analysis\FC\figures\', subj, 'ipsi_contra_SAMESCALE'), '-dpdf', '-fillpage')
           case "inv/ni"
               print(strcat('\\dynaserv\meg\nicolas\PTSD\analysis\FC\figures\', subj, 'inv_ni_SAMESCALE'), '-dpng');
%                print(strcat('\\dynaserv\meg\nicolas\PTSD\analysis\FC\figures\', subj, 'inv_ni_SAMESCALE'), '-dpdf', '-fillpage')
       end
    
   end   
    clearvars  subj channels check linksrest OUTrest TOTrest subj_info min max rest h2_array
    
end

%% save the table
% switch category
%    case "ipsi/contra"
%        varnames = {varnames{:}, 'OUTrest_all', 'TOTrest_all', 'OUTrest_ipsi_ipsi','TOTrest_ipsi_ipsi','OUTrest_contra_contra', 'TOTrest_contra_contra', 'OUTrest_ipsi_contra', 'TOTrest_ipsi_contra', 'OUTrest_contra_ipsi', 'TOTrest_contra_ipsi'};
%    case "inv/ni"
%        varnames = {varnames{:}, 'OUTrest_all', 'TOTrest_all', 'OUTrest_inv_inv','TOTrest_inv_inv','OUTrest_ni_ni', 'TOTrest_ni_ni', 'OUTrest_inv_ni', 'TOTrest_inv_ni', 'OUTrest_ni_inv', 'TOTrest_ni_inv'};
% end
    
varnames = {varnames{:}, 'OUTrest_all', 'TOTrest_all', 'OUTrest_inv_inv','TOTrest_inv_inv','OUTrest_ni_ni', 'TOTrest_ni_ni', 'OUTrest_inv_ni', 'TOTrest_inv_ni', 'OUTrest_ni_inv', 'TOTrest_ni_inv'};
FCtable = array2table(FCtable, 'VariableNames', varnames);

switch category
    case "ipsi_contra"
        FCtable.involved(~ismissing(FCtable.TOTrest_inv_inv)) = "ipsi";
        FCtable.involved(~ismissing(FCtable.TOTrest_ni_ni)) = "contra";           
    case "inv_ni"
        FCtable.involved(~ismissing(FCtable.TOTrest_inv_inv)) = "inv";
        FCtable.involved(~ismissing(FCtable.TOTrest_ni_ni)) = "non-inv";
    case "DMN_CEN"
        FCtable.involved(~ismissing(FCtable.TOTrest_inv_inv)) = "DMN";
        FCtable.involved(~ismissing(FCtable.TOTrest_ni_ni)) = "CEN";
end


 % if also inv-ni or ipsi-contra
%        varnames = {varnames{:}, 'OUTrest_all', 'TOTrest_all', 'OUTrest_ipsi_ipsi','TOTrest_ipsi_ipsi','OUTrest_contra_contra', 'TOTrest_contra_contra', 'OUTrest_ipsi_contra', 'TOTrest_ipsi_contra', 'OUTrest_contra_ipsi', 'TOTrest_contra_ipsi'};
%        varnames = {varnames{:}, 'OUTrest_all', 'TOTrest_all', 'OUTrest_inv_inv','TOTrest_inv_inv','OUTrest_ni_ni', 'TOTrest_ni_ni', 'OUTrest_inv_ni', 'TOTrest_inv_ni', 'OUTrest_ni_inv', 'TOTrest_ni_inv'};


end   
