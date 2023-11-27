close all
clear all
clc
%%
% This code does the EEG outlier detection and rejection step in WAVEFRONT
%
% See: README.txt and [1] for more info.

% [1] Alireza Chamanzar, Jonathan Elmer, Lori Shutter, Jed Hartings, Pulkit Grover,
%  "Noninvasive and reliable automated detection of spreading depolarization in severe traumatic brain injury using scalp EEG",
%   Submitted to Nature Comms Med, 2022.

% Author: Alireza 	Date: 2022/01/03 12:00:00 	Revision: 0.1
% Copyright: This project is licensed - see the LICENSE.md file for details

%% Preprocessing steps for EEG data:
eegpath = 'C:\Users\Alireza\Downloads\eeglab_current\eeglab2020_0';
current_path = pwd;
cd(eegpath)
eeglab
cd(current_path)
% Similar to previous steps, ensure eeglab is running and initialize path

% Find all Sessions of chosen Patient ID
Session_names = dir('Patient*');
Sessions_size = size(Session_names,1);

% Loop through all sessions of Patient ID
for ss = 1: Sessions_size
  
    % Enter Session-specific directory 
    cd(Session_names(ss).name)
    current_path = pwd;

    % Load 'Visualization' file which contains Session-wide EEG & ECoG data
    EEG_vis = pop_loadset('filename',[Session_names(ss).name,'_Visualization_Delta_ECoGcleaned_vII.set'],'filepath',current_path);
    EEG_vis = pop_select( EEG_vis,'nochannel',{'Depth1' 'Depth2' 'Depth3' 'Depth4' 'Depth5' 'Depth6'});
    
    load([Session_names(ss).name,'_woGamma_preica_Imp.mat'])
    
    % Collection of commented out code blocks
    % X = [EEG_vis.chanlocs.X];
    % Y = [EEG_vis.chanlocs.Y];
    % Z = [EEG_vis.chanlocs.Z];
    % EEG = EEG_vis;
    % Add the Impedance channels to EEG:
    % for Imp_ind=1:19
        % EEG.data(end+1,:) = measurement_data_Imp_temp_copy(Imp_ind,:);
        % EEG.chanlocs(end+1).labels = sprintf('Imp_%d',Imp_ind);
    % end
    % EEG.nbchan = size(EEG.data,1);    
    % Remove poor data equality events:
    % Remove_ind = find(strcmp({EEG.event(:).type},'Electrode Impedance High'));
    % Norm_ind = find(strcmp({EEG.event(:).type},'Data Quality Normal'));
    % 
    % Remoxve_latency = [EEG.event(Remove_ind).latency];
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
    % measurement_data_Imp_temp_copy = EEG.data(end-18:end,:);
    % EEG = pop_select( EEG,'nochannel',size(EEG.data,1)-18:size(EEG.data,1));    
    % EEG_vis = EEG;
    %
    % WL = 10*60*srate; 
    
    srate = EEG_vis.srate;
    
    %% Begin 'Part' processing

    % Set Part size of 4 hours
    Block_size = 4*60*60*srate+1;

    % Create an array of the EEG data from EEG vis data field 
    EEG = EEG_vis.data(1:19,:);
    
    % Check for edge cases of full data less than 4 hours (or Part length)
    if(size(EEG,2)<Block_size)
        Block_size = size(EEG,2);
    end
    
    % M_global_tot = zeros(19,1);
    % M_global_Count = zeros(19,1);

    % Create starting index array of 3.5 hr intervals
    % This creates Parts of 4 hrs with 30 min overlap of adjacent Parts 
    for CSD_ind=1:3.5*60*60*srate:size(EEG,2)-Block_size+1
        
        % Set Part start index as 'CSD_ind' from for loop
        T_start = CSD_ind;
        % Set Part end as CSD + 4 hrs, or as end of data in edge-cases
        T_end = min(CSD_ind+Block_size-1,size(EEG,2));
        
        % Check CSD_ind is less than 3.5 hrs away from end of data
        % if small amount of data, extend 'final' Part beyond Part size
        if(abs(size(EEG,2)-Block_size+1-CSD_ind)<3.5*60*60*srate)
            T_end = size(EEG,2);
        end
        
        %Create logic array for EEG mask
        EEG_mask = (measurement_data_Imp_temp_copy>0);

        %Trim data & logic array from full Session to specific Part
        EEG_sub = EEG(:,T_start:T_end);
        EEG_mask = EEG_mask(:,T_start:T_end);
        
        % Anywhere EEG_mask is zero, mask from EEG data 
        EEG_sub(EEG_mask==0) = 0;
        
        % Create copies of EEG data for subsequent processing 
        EEG_data = EEG_sub;
        EEG_temp = EEG_sub;
        % Initialize a struct for bad indicies 
        Delete_ind_tot = struct;

        % Loop through all 19 EEG channels
        parfor j=1:19 % This part can either be run in parallel or series
            j;
            % Create data array of all non-zero data points at channel j  
            temp = EEG_sub(j,EEG_mask(j,:)>0); 
            % Find indices of EEG sub which are not zero
            TF_temp = find(EEG_mask(j,:)>0);

            % std_max = (std(temp,[],2));
            % if(std_max<=0 || isempty(temp))
                % std_max=30;
            % end

            % Find outliers which 3 IQRs above the upper quartile of data
            [~,TF] = rmoutliers(temp,'quartiles','ThresholdFactor',3);
            T_ind = TF_temp(TF);
            Delete_ind = [];
            Delete_ind = cat(2,Delete_ind,T_ind);
            % Store channel specific 'bad' indices in Part wide struct
            Delete_ind_tot(j).Delete = unique(Delete_ind);
        end
        
        % Set window length of 1 sec (in samps) for outlier removal padding
        w_length = 1*64;
        for j=1:19 % Loop through each channel 
             % Remove window length before and after each outlier index
            for i=-w_length:w_length 
                Delete_ind_delayed = Delete_ind_tot(j).Delete+i;
                % Ensure no delayed indices are negative or above upper lim
                Delete_ind_delayed(Delete_ind_delayed<=0) = [];
                Delete_ind_delayed(Delete_ind_delayed>size(EEG_data,2)) = [];
                % Remove outlier delayed indices from data & mask
                EEG_data(j,Delete_ind_delayed) = 0;
                EEG_mask(j,Delete_ind_delayed) = 0;
            end
        end
        
        % Set EEG sub as newly cleaned data
        EEG_sub = EEG_data;
       
        for j=1:19 % For each of the 19 channels
            % Extract data portions which are not equal to zero
            EEG_temp = EEG_data(j,EEG_mask(j,:)>0);
            % If temp is not zero, normalize channel by channel's std
            if(~isempty(EEG_temp))
                EEG_sub(j,:) = EEG_sub(j,:)./std(EEG_temp,0,2);
            end
        end
        
        % Ensure no data points are NaN
        EEG_sub(isnan(EEG_sub)) = 0;
        % Final masking of EEG data 
        EEG_sub(EEG_mask==0) = 0;
        
        % 'save' normalized and cleaned EEG data back to EEG_vis struct for
        % appropriate time region [T_start:T_end] 
        EEG_vis.data(1:19,T_start:T_end) = EEG_sub;
        
        % Save EEG mask as Impedance variable for proper time region
        measurement_data_Imp_temp_copy(1:19,T_start:T_end) = EEG_mask;
        
        clear CSD_type Events_type
    end
    
    EEG_vis = pop_saveset(EEG_vis, 'filename',[Session_names(ss).name,'_Visualization_Delta.set'],'filepath',current_path);
    save([Session_names(ss).name,'_woGamma_preica_Imp_Delta.mat'],'measurement_data_Imp_temp_copy','-v7.3')
    
    cd ..
end



