function info = readinfo_PTSD(subj_code, chan_names, tabledir, varargin)
%
% combine anatomical labels to subject's clincal info in one melt table (1 row/channel) where you will fill other measures in
% inputs
%   subj_code: bids code of subject to match for adding info
%   chan_names: list of electrodes of that subject
%   tabledir: path of table with info on channels and subjects (for now: roi, anat, group)


% % debug and test: 
% tabledir = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\ChanInfo_PS2.xlsx";
% simEFdir = "C:\Users\simula\OneDrive - Aix-Marseille Université\PhD\MyProjects\Galvani\analysis\EFtable2023-04-28_12.00.csv";
% subj_code = "2305a0fd66e8";
% load('\\dynaserv\Galvani_ps2\analysis\delphos\all_mat-thr80\baseline\sub-2305a0fd66e8_baseline.mat');
% chan_names = results.labels;

%% simEF was made for PS2 study but no need for PTSD study
%     for ii = 1:2:nargin-3
%             if strcmp('simEFdir', varargin{ii})
%                 simEFdir = varargin{ii+1}; 
%             end
%     end

%% Read table and select info
    chan_info = readtable(tabledir);
    
    %Associate each channel to its ROI (or other info) from the info_table
    info = {};
    subchan = chan_info(chan_info.subject == string(subj_code), :); %extract subtable of subject's data from the excel sheet
    varnames = subchan.Properties.VariableNames;

    % for each channel in 'channels', add info present on the excel table
    for i=1:length(chan_names)
       check = string(chan_names(i)); % extracts info from monopolar channel list since montages can be diff in chaninfo table and in analysed file (mono vs bipolar montage)
       
       if ~isempty(string(subchan.brain_area(subchan.channel == check)))
%            anat = string(subchan.brain_area(subchan.channel == check));
           toadd = table2cell(subchan(subchan.channel == check, :));
       else
           toadd = array2table(nan(1,size(subchan,2)));
           toadd = cellfun(@(x) nan(1,size(x,2)),table2cell(toadd),'UniformOutput',false);
       end
           info = [info; toadd];
       
           
    end
     
    info = cell2table(info);
    allVars = 1:width(info);
    newNames = {varnames{:}};
    info = renamevars(info,allVars,newNames);

    clearvars check subchan 


end