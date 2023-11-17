close all
clear all
clc
%%
% This code downsamples the dataset, concatenates small parts of each EEG session, prunes EEG based on impedance:
%
% See: README.txt and [1] for more info.

% [1] Alireza Chamanzar, Jonathan Elmer, Lori Shutter, Jed Hartings, Pulkit Grover,
%  "Noninvasive and reliable automated detection of spreading depolarization in severe traumatic brain injury using scalp EEG",
%   Submitted to Nature Comms Med, 2022.

% Author: Alireza 	Date: 2021/12/19 12:00:48 	Revision: 0.1
% Copyright: This project is licensed - see the LICENSE.md file for details

%% Preprocessing steps for EEG data:
% First enter eeglab directory (change with local eeglab path)
eegpath = 'C:\Users\Alireza\Downloads\eeglab_current\eeglab2020_0';
cd(eegpath) 
eeglab % Run eeglab for supporting functions 
% Jump back to Patient ID directory containing Session Folders
cd(current_path) 
% Find all Session Folders in given Patient ID
Session_names = dir('Patient*');
Sessions_size = size(Session_names,1); 
% Load Population-Wide median Impedance (modify path to local path)
load D:\SD-II_POST_ICA\SD-II\median_Imp.mat
median_Imp_tot = median_Imp;
% Loop through all recording sessions of given Patient ID
for ss = 1:Sessions_size
    % Enter Session specific folder for each Session loop
    cd(Session_names(ss).name)
    % Load Session-specific Impdence Data 
    EEG_Imp = dir('EEG,Composite,Impedance,*.mat')
    % Begin by loading first EEG Impedence File
    load(EEG_Imp(1).name)
    % Assign loaded variables as 'Impedance' Variables 
    comp_elements_Imp = comp_elements;
    measurement_data_Imp = measurement_data;
    time_vector_Imp = time_vector;
    Imp_srate = 1/(time_vector_Imp(2) - time_vector_Imp(1));
    % If more than one Impedance file is present, append all information
    for Imp_ind = 2:size(EEG_Imp,1)
        load(EEG_Imp(Imp_ind).name)
        measurement_data_Imp = cat(2,measurement_data_Imp,measurement_data);
        time_vector_Imp = cat(2,time_vector_Imp,time_vector);
    end 
    current_path = pwd; % Ensure pathing is still correct
    % Find all ECoG files for given Session
    ECoG_Part_names = dir('*_ECoG_filtered_withDC.set');
    
    current_path = pwd;
    % Find all 'preICA w DC' files which represent our EEG data
    Part_names = dir('*_preica_withDC.set');
    Part_name_char = string({Part_names.name});
    % Remove possible 'extra' files in the directory to prevent nested merging process:
    part_woGamma_ind = find(contains(Part_name_char,'_woGamma_preica_withDC.set'));
    Part_name_char(part_woGamma_ind) = [];
    Part_names(part_woGamma_ind) = [];
    % Measure number of Parts for given Session
    Part_size = size(Part_names,1);
    % Find indices of specific parts of file name
    part_ind = strfind(Part_name_char,'part');
    of_ind = strfind(Part_name_char,'of');
    preica_ind = strfind(Part_name_char,'_preica');
    Part_num = [];  % Initialize list of Part Numbers
    for i=1:Part_size
        % Using string locations from above, create list of part numbers
        Part_num = cat(1,Part_num,str2double(Part_name_char{i}(part_ind{i}+4:of_ind{i}-1)));
    end
    AllEEG = []; % Initilize 'ALLEEG' variable for merging parts
    [~,Part_ind] = sort(Part_num); % Ensure Part numbers are in correct order
    measurement_data_Imp_temp = []; % Initialize list for imp_data
    % Load a sample EEG data file for channel order information
    EEG_test = pop_loadset('filename',Part_names(Part_ind(1)).name,'filepath',current_path);
    EEG_test = pop_select( EEG_test,'nochannel',{'Depth1' 'Depth2' 'Depth3' 'Depth4' 'Depth5' 'Depth6'});
    Imp_ind = [];
    % Find channels which overlap b/w EEG data and Impedence Data
    for el = 1:19
        Imp_ind = cat(1,Imp_ind,find(contains(comp_elements_Imp,EEG_test.chanlocs(el).labels)));
    end
    % Only keep Impedence Channels which overlap with EEG Channls
    measurement_data_Imp = measurement_data_Imp(Imp_ind,:);
    % Check for any overlapping time index differences
    [time_vector_Imp,ind] = unique(time_vector_Imp);
    % Remove any overlapping time points from Impedance data
    measurement_data_Imp = measurement_data_Imp(:,ind);
    % Loop through all 'parts' of Session
    for pp = 1: Part_size
        Part_names(Part_ind(pp)).name
        % Load specific EEG part file 'pp'
        EEG = pop_loadset('filename',Part_names(Part_ind(pp)).name,'filepath',current_path);
        % Remove some unnecessary channels from EEG data
        EEG = pop_select( EEG,'nochannel',{'Depth1' 'Depth2' 'Depth3' 'Depth4' 'Depth5' 'Depth6'});
        % Load specific ECoG part file 'pp'
        ECoG = pop_loadset('filename',ECoG_Part_names(Part_ind(pp)).name,'filepath',current_path);
        % Perfrom lowpass Filtering on EEG Data below 30Hz
        EEG = pop_eegfiltnew(EEG, [],30,1000,0,[],0);
        % Resample EEG data from 256 Hz down to 64 Hz
        [EEG_resample] = pop_resample(EEG, 64);
        % Save resampled EEG data back to main EEG variable
        EEG = EEG_resample;
        % Perform lowpass Filtering on ECoG Data below 30Hz 
        ECoG = pop_eegfiltnew(ECoG, [],30,1000,0,[],0);
        % Resample ECoG Data from 256 to 64 Hz
        [ECoG_resample] = pop_resample(ECoG, 64);
        % Save resampled ECoG data back to main ECoG variable
        ECoG = ECoG_resample;

        %% Remove poor data equality events:
        % Commented out code below no longer used for output generation
        % Remove_ind = find(strcmp({EEG.event(:).type},'Electrode Impedance Uncertain') | ...
        %     strcmp({EEG.event(:).type},'Electrode Impedance High') | ...
        %     strcmp({EEG.event(:).type},'Electrode Impedance Abnormal Low') | ...
        %     strcmp({EEG.event(:).type},'Electrode Impedance Imbalance'));
        % Norm_ind = find(strcmp({EEG.event(:).type},'Data Quality Normal'));
        % 
        % Remove_latency = [EEG.event(Remove_ind).latency];
        % Norm_latency= [EEG.event(Norm_ind).latency];
        % 
        % eerej_vec = [-1,-1];%dummy value
        % for ee=1:size(Remove_latency,2)
        %     if(Remove_latency(ee)>eerej_vec(end,2))
        %         ind_temp = find(Norm_latency>Remove_latency(ee));
        %         if(~isempty(ind_temp))
        %             eerej_vec = cat(1,eerej_vec,[Remove_latency(ee),Norm_latency(ind_temp(1))]);
        %         else
        %             eerej_vec = cat(1,eerej_vec,[Remove_latency(ee),EEG.pnts]);
        %         end
        %     end
        % end
        % eerej_vec(1,:) = [];%dummy value 
        % Find where Imp vector begins and ends relative to specific part
        Imp_strt = find(time_vector_Imp>=EEG.xmin-10);
        Imp_strt = Imp_strt(1);
        Imp_end = find(time_vector_Imp<=EEG.xmax+10);
        Imp_end = Imp_end(end);
        
        % This next line does several things at once, it first uses
        % imp_start and imp_end to trim Impedance data to appropriate time
        % range. It then interpolates the impedance data so that it is of
        % equal length to the EEG data vector, finally it appends the 
        % part specific data to a Session-specific 'storage' variable 
        measurement_data_Imp_temp = cat(2,measurement_data_Imp_temp,interp1(time_vector_Imp(Imp_strt:Imp_end)*1000,measurement_data_Imp(:,Imp_strt:Imp_end)',EEG.times,'pchip')');
        
        EEG.urevent = {};
        if(pp==1)
            AllEEG = EEG;
        else
            AllEEG = pop_mergeset(AllEEG,EEG);
        end
        
        if(pp==1)
            AllECoG = ECoG;
        else
            AllECoG = pop_mergeset(AllECoG,ECoG);
        end
    end
    
    
    measurement_data_Imp_temp = cat(2,measurement_data_Imp_temp(:,320:end),100*ones(19,319));
    data_temp = AllEEG.data(1:19,:);
    measurement_data_Imp_temp_copy = measurement_data_Imp_temp;
    for rep_i=1:1
    median_Imp = median(measurement_data_Imp_temp_copy,2);
    quantile_Imp = quantile(measurement_data_Imp_temp_copy,0.9,2);
    thr_Imp = 2;
    for ch_i=1:19
        eerej_vec = measurement_data_Imp_temp_copy(ch_i,:)>thr_Imp*min(median_Imp(ch_i),median_Imp_tot);
        eerej_vec_comp = bwconncomp(eerej_vec);
        for comp_ind=1:eerej_vec_comp.NumObjects
            data_temp(ch_i,max(1,eerej_vec_comp.PixelIdxList{1,comp_ind}(1))...
                :min(size(data_temp,2),eerej_vec_comp.PixelIdxList{1,comp_ind}(end))) = 0;
            measurement_data_Imp_temp_copy(ch_i,max(1,eerej_vec_comp.PixelIdxList{1,comp_ind}(1))...
                :min(size(measurement_data_Imp_temp_copy,2),eerej_vec_comp.PixelIdxList{1,comp_ind}(end))) = 0;
        end
    end
    data_temp(measurement_data_Imp_temp_copy<=0) = 0;
    end
    AllEEG.data(1:19,:) = data_temp;

%% 
    pop_eegplot(AllEEG, 1, 1, 1)

    save([Session_names(ss).name,'_woGamma_preica_Imp.mat'],'measurement_data_Imp_temp_copy','-v7.3')
    EEG = pop_saveset( AllEEG, 'filename',[Session_names(ss).name,'_woGamma_preica_withDC_vII.set'],'filepath',current_path);
    ECoG = pop_saveset( AllECoG, 'filename',[Session_names(ss).name,'_woGamma_ECoG_withDC_vII.set'],'filepath',current_path);

    cd ..
end





