close all
clear all
clc
%%
% This code adds the ECoG channels to EEG:
%
% See: README.txt and [1] for more info.

% [1] Alireza Chamanzar, Jonathan Elmer, Lori Shutter, Jed Hartings, Pulkit Grover,
%  "Noninvasive and reliable automated detection of spreading depolarization in severe traumatic brain injury using scalp EEG",
%   Submitted to Nature Comms Med, 2022.

% Author: Alireza 	Date: 2022/01/01 09:00:00 	Revision: 0.1
% Copyright: This project is licensed - see the LICENSE.md file for details

%% Preprocessing steps for EEG data:
eegpath = 'C:\Users\Alireza\Downloads\eeglab_current\eeglab2020_0';
% Configure local paths so that eeglab can be loaded and run
current_path = pwd;
cd(eegpath)
eeglab % Run eeglab for supporting functions 

% Enter directory containing Subject ID specific data
cd(current_path)

% Load all Sessions of chosen Patient ID
Session_names = dir('Patient*');
Sessions_size = size(Session_names,1);

% Loop through all sessions of Patient ID
for ss = 1:Sessions_size
    % Enter Session-specific directory 
    cd(Session_names(ss).name)
    current_path = pwd;

    % Load 'wo_gamma' file which contains Session-wide filtered EEG data
    Part_names = dir('*_woGamma_preica_withDC_vII.set');
    Part_size = size(Part_names,1);
    
    EEG = pop_loadset('filename',Part_names(1).name,'filepath',current_path);
    % Filter out non-standard EEG channels and Physio channels
    EEG = pop_select( EEG,'nochannel',{'ECG' 'PLETH' 'RESP' 'EOGv' 'EOGh' 'X1' 'X2' 'Y1' 'Y2'});
    % Load Session-wide Impedance file 
    load([Session_names(ss).name,'_woGamma_preica_Imp.mat'])
    
    % Load Session-wide ECoG file
    Part_names = dir('*_woGamma_ECoG_withDC_vII.set');
    Part_size = size(Part_names,1);
    ECoG = pop_loadset('filename',Part_names(1).name,'filepath',current_path);
    
    srate = EEG.srate; % Read sampling rate from EEG struct
    
    % In the next block of code we number our CSD events
    CSD_Labels = []; % Initialize a vector for CSD Labels
    CSD_ind = 1; % Start CSD counting at 1
    for i=1:size(EEG.event,2) % Loop through all events
        % find 'Spreading Depolarization' events
        if(strcmp(EEG.event(i).type,'CSD') || strcmp(EEG.event(i).type,'CSD/ISD') ||...
                strcmp(EEG.event(i).type,'ISD') || strcmp(EEG.event(i).type,'scCSD'))
            % At each 'SD' event, ennumerate beginging of event label
            EEG.event(i).type = [sprintf('%d_',CSD_ind) EEG.event(i).type];
            CSD_Labels = cat(1,CSD_Labels,EEG.event(i));
            CSD_ind = CSD_ind+1;
        end
    end
    
    EEG.data(measurement_data_Imp_temp_copy<=0) = 0;
    
    figure;plot(EEG.data(1,:))
   
    
     %%
     % This code block is a 'dead-end' block, thus it is commented out 
     
    % CSD_Labels = [];
    % Other_Labels = [];
    % WL = 5*60*srate;
    % CSD_GT_total = zeros(1,size(EEG.data,2));
    % Event_total = zeros(1,size(EEG.data,2));
    % CSD_Labels_total = [];
    % Event_Labels_total = [];
    % for i=1:size(EEG.event,2)
    %     if(contains(EEG.event(i).type,'CSD') || contains(EEG.event(i).type,'CSD/ISD') ||...
    %             contains(EEG.event(i).type,'ISD') || contains(EEG.event(i).type,'scCSD'))
    %         EEG.event(i).latency = EEG.event(i).latency;
    %         CSD_Labels = cat(1,CSD_Labels,EEG.event(i));
    %         CSD_GT_total(round(max(EEG.event(i).latency,1)))=1;
    %         CSD_Labels_total = cat(1,CSD_Labels_total,{EEG.event(i).type});
    %     else%if(~(strcmp(EEG.event(i).type,'boundary')))% || strcmp(EEG.event(i).type,'Data Quality Normal')))
    %         EEG.event(i).latency = EEG.event(i).latency;
    %         Other_Labels = cat(1,Other_Labels,EEG.event(i));
    %         Event_total(round(max(EEG.event(i).latency,1)))=1;
    %         Event_Labels_total = cat(1,Event_Labels_total,{EEG.event(i).type});
    %     end
    % end
    
    % sensor_locs = [EEG.chanlocs.X;EEG.chanlocs.Y;EEG.chanlocs.Z]';
    
    % Add the ECoG channels to EEG:
    for ECoG_ind=1:6
        EEG.data(end+1,:) = ECoG.data(ECoG_ind,:);
    end
    EEG.nbchan = size(EEG.data,1);
    if ~isempty(EEG.chanlocs)
        EEG.chanlocs(end+1).labels = 'ECoG_1';
        EEG.chanlocs(end+1).labels = 'ECoG_2';
        EEG.chanlocs(end+1).labels = 'ECoG_3';
        EEG.chanlocs(end+1).labels = 'ECoG_4';
        EEG.chanlocs(end+1).labels = 'ECoG_5';
        EEG.chanlocs(end+1).labels = 'ECoG_6';
    end
    
    eeglab redraw
    pop_eegplot( EEG, 1, 1, 1);
    
    EEG = pop_saveset( EEG, 'filename',[Session_names(ss).name,'_Visualization_withDC_vII.set'],'filepath',current_path);

    cd ..
end



