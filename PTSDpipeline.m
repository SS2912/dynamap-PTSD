%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%     PTSD pipeline     %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sara simula, August 2023
% pipeline of analysis for PTSD project with Lisa S. (PSD, FC, IEDs)

%% initialization and settings
clear
clc
date = datestr(clock,'YYYY-mm-dd_HH.MM');

%% 1. create table with anat info for each channel and clincal info

dir_gardel = "\\dynaserv\meg\nicolas\PTSD\dataset_PTSD\derivatives\Gardel";
stringToBeFound = 'SEEG_anatomical_labels_bipolar';
clinical_info = readtable("\\dynaserv\meg\nicolas\PTSD\Selection_info.xlsx"); 

infotable_melt = combineinfo_PTSD(dir_gardel, stringToBeFound, clinical_info);
name = strcat('\\dynaserv\meg\nicolas\PTSD\', 'PTSD_meltinfo', "-", date, '.xlsx'); % when saving to csv, it doesnt save variable names
writetable(infotable_melt, name)  

%% %%%%%%%%%%%%%  2. FC analysis  %%%%%%%%%%%%%%%%%%%
clear
close all

% Varargin:
%   thr_h2, thr_lag - optional thershold for h2 (thr_h2=0 for strength, >0 for degrees) and lag. Default: thr_h2 = 0; thr_lag = 0
%   norm            - 1 if you want to nromalise each subj by number of channels (Default= 1)
%   graph           - 1 (default) to show connectivity matrices of zvalues compared to baseline chosen in base, 0 to not show any graph
%   roi             - "all" (default, h2 between all channels) or "EZ", "EZPZ", "NI" to calculate node strength only in subset of EZ channels (or non-inv channels)

dir_info = "\\dynaserv\meg\nicolas\PTSD\PTSD_meltinfo.xlsx";
category = "inv_ni";  % ipsi_contra or inv_ni or DMN_CEN
method = "mean";

%% Broad %%
dir_data = "\\dynaserv\meg\nicolas\PTSD\analysis\FC\h2matrices\broad";
[FCtable, h2_amyhpc] = FC_PTSD(dir_data, dir_info, category, 'graph', 0, 'meth', method); 

name = strcat('\\dynaserv\meg\nicolas\PTSD\analysis\FC\', 'FCtable-', category, "-broad-AMYHPC", method, "-", date, '.xlsx'); % when saving to csv, it doesnt save variable names
writetable(FCtable, name)  

%% Delta %%
dir_data = "\\dynaserv\meg\nicolas\PTSD\analysis\FC\h2matrices\delta";
[FCtable, h2_amyhpc] = FC_PTSD(dir_data, dir_info, category, 'graph', 0);

name = strcat('\\dynaserv\meg\nicolas\PTSD\analysis\FC\', 'FCtable-', category, "-delta-AMYHPC", method, "-", date, '.xlsx'); % when saving to csv, it doesnt save variable names
writetable(FCtable, name)  

%% Theta %%
dir_data = "\\dynaserv\meg\nicolas\PTSD\analysis\FC\h2matrices\theta";
[FCtable, h2_amyhpc] = FC_PTSD(dir_data, dir_info, category, 'graph', 0);

name = strcat('\\dynaserv\meg\nicolas\PTSD\analysis\FC\', 'FCtable-', category, "-theta-AMYHPC", method, "-", date, '.xlsx'); % when saving to csv, it doesnt save variable names
writetable(FCtable, name)  

%% Alpha %%
dir_data = "\\dynaserv\meg\nicolas\PTSD\analysis\FC\h2matrices\alpha";
[FCtable, h2_amyhpc] = FC_PTSD(dir_data, dir_info, category, 'graph', 0);

name = strcat('\\dynaserv\meg\nicolas\PTSD\analysis\FC\', 'FCtable-', category, "-alpha-AMYHPC", method, "-", date, '.xlsx'); % when saving to csv, it doesnt save variable names
writetable(FCtable, name)  

%% Beta %%
dir_data = "\\dynaserv\meg\nicolas\PTSD\analysis\FC\h2matrices\beta";
[FCtable, h2_amyhpc] = FC_PTSD(dir_data, dir_info, category, 'graph', 0);

name = strcat('\\dynaserv\meg\nicolas\PTSD\analysis\FC\', 'FCtable-', category, "-beta-AMYHPC", method, "-", date, '.xlsx'); % when saving to csv, it doesnt save variable names
writetable(FCtable, name)  

%% Lowgamma %%
dir_data = "\\dynaserv\meg\nicolas\PTSD\analysis\FC\h2matrices\lowgamma";
[FCtable, h2_amyhpc] = FC_PTSD(dir_data, dir_info, category, 'graph', 0);

name = strcat('\\dynaserv\meg\nicolas\PTSD\analysis\FC\', 'FCtable-', category, "-lowgamma-AMYHPC", method, "-", date, '.xlsx'); % when saving to csv, it doesnt save variable names
writetable(FCtable, name)  

%% 3. Delphos analysis

% needed to change script from PS2 cause patietns have different events
% list!! changes can be applied to PS2 too but not needed for now

% % change names to montages to add subject code to each mtg file
% folder_raw = "\\dynaserv\meg\nicolas\PTSD\montages"; % one folder per subject with 3 different mtgs files
% folder_output = "\\dynaserv\meg\nicolas\PTSD\analysis\delphos\matresults\mtgs";
% myfiles  = dir(strcat(folder_raw, "\**\*raw.mtg")); % find the delphos one (bipolar_raw)
% 
% for i=1:length(myfiles)
%     subj = extractAfter(myfiles(i).folder, 'sub-');
%     newname = strcat("sub-", subj, "_bipolar-raw.mtg");
%     movefile(strcat(myfiles(i).folder, '\', myfiles(i).name), strcat(folder_output, '\', newname))
% end
clear
clc

input_folder  = "\\dynaserv\meg\nicolas\PTSD\analysis\delphos\matresults";
mtgs_folder = "\\dynaserv\meg\nicolas\PTSD\analysis\delphos\mtgs\*.mtg";
dir_info = "\\dynaserv\meg\nicolas\PTSD\PTSD_meltinfo.xlsx";
sections = "rest";

[results_all, delphos_melttable] = delphos_PTSD(input_folder, mtgs_folder, sections, 'tabledir', dir_info);

name = strcat('\\dynaserv\meg\nicolas\PTSD\analysis\delphos\', 'Delphos_table', date, '.xlsx'); % when saving to csv, it doesnt save variable names
writetable(delphos_melttable, name)  

%% 4. PSD analysis
clear; clc

dir_data = "\\dynaserv\meg\nicolas\PTSD\analysis\psd\bids_results";
dir_info = "\\dynaserv\meg\nicolas\PTSD\PTSD_meltinfo.xlsx";
sections = "rest";

bands    = [1 4; 4 8; 8 15; 15 30; 30 45]; % this is the default in spectral_PTSD.m
band_name = ["delta", "theta", "alpha", "beta", "gamma"]; % default in function

PSDtable = spectral_PTSD(dir_data, dir_info, 'window', 4, 'bands', bands);

PSD_firsthalf = PSDtable(1:size(PSDtable,1)/2, :);
name = strcat('\\dynaserv\meg\nicolas\PTSD\analysis\psd\Tables\', 'PSD_fooof_table_1_', date, '.xlsx'); % when saving to csv, it doesnt save variable names
writetable(PSD_firsthalf, name)  

PSD_sechalf = PSDtable((size(PSDtable,1)/2)+1:end, :);
name = strcat('\\dynaserv\meg\nicolas\PTSD\analysis\psd\Tables\', 'PSD_fooof_table_2_', date, '.xlsx'); % when saving to csv, it doesnt save variable names
writetable(PSD_sechalf, name) 
