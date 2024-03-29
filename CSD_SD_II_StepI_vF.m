clear all
close all
clc
%%
% This code does the ECoG outlier detection and rejection step in WAVEFRONT
%
% See: README.txt and [1] for more info.

% [1] Alireza Chamanzar, Jonathan Elmer, Lori Shutter, Jed Hartings, Pulkit Grover,
%  "Noninvasive and reliable automated detection of spreading depolarization in severe traumatic brain injury using scalp EEG",
%   Submitted to Nature Comms Med, 2022.

% Author: Alireza 	Date: 2022/01/02 12:00:00 	Revision: 0.1
% Copyright: This project is licensed - see the LICENSE.md file for details

%% Initialization:
% cd 'H:\SD-II_POST_ICA\SD-II\04-1203\04-1203'

eegpath = 'C:\Users\Alireza\Downloads\eeglab_current\eeglab2020_0';
current_path = pwd;
cd(eegpath)
eeglab  
% Similar to previous steps, ensure eeglab is running and initialize path

% Find all Sessions of chosen Patient ID
Session_names = dir('Patient*');
Sessions_size = size(Session_names,1);

% Loop through all sessions of Patient ID
for ss = 1:Sessions_size
    
    % Enter Session-specific directory 
    cd(Session_names(ss).name)
    current_path = pwd;
    % Load 'Visualization' file which contains Session-wide EEG & ECoG data
    EEG_vis = pop_loadset('filename',[Session_names(ss).name,'_Visualization_withDC_vII.set'],'filepath',current_path);
    
    % Find any instances of 'Not a Number' [NaN] within the data
    % [ind_r,ind_c] = find(isnan(EEG_vis.data));
    % Zero anywhere the data is 'NaN'
    % EEG_vis.data(ind_r,ind_c) = 0;
    
    EEG_data = EEG_vis.data;
    EEG_data(isnan(EEG_data)) = 0;
    EEG_vis.data = EEG_data;

    % Filter EEG data to extract the delta frequency band 
    EEG_vis = pop_eegfiltnew(EEG_vis, 0.5,4,10000,0,[],1);
    
    % Extract ECoG channels from EEG struct and save as "EEG_sub"
    EEG_sub = EEG_vis.data(end-5:end,:);
    Limit_vec = ones(size(EEG_sub,1),2);%[1e5,1.5e5];
    Limit_vec(:,2) = Limit_vec(:,2)*size(EEG_sub,2);
    
    % Save copies of ECoG data
    EEG_data = EEG_sub; % Copy for masking 
    EEG_temp = EEG_sub; % Copy for plotting pre outlier removed data
    Delete_ind_tot = struct; % Initialize a struct for bad indicies

    % Loop through each of the 6 ECoG channels
    for j=1:6
        j;
        % Create data array of all non-zero data points at channel j 
        temp = EEG_sub(j,EEG_sub(j,:)~=0); 
        % Find indices of EEG sub which are not zero
        TF_temp = find(EEG_sub(j,:)~=0);
        % Find outliers 10 IQRs above the upper quartile of the data
        [~,TF] = rmoutliers(temp,'quartiles','ThresholdFactor',10);
        % Map outlier indices back to non-zeroed space
        T_ind = TF_temp(TF);
        Delete_ind = [];
        Delete_ind = cat(2,Delete_ind,T_ind);
        % Store channel specific mask indices in Session wide struct
        Delete_ind_tot(j).Delete = unique(Delete_ind);
    end
    
    % Set window length for outlier removal padding  
    w_length = 100;
    for j=1:6 % Loop through each channel 
        % Remove window lengths before and after each outlier index
        for i=-w_length:w_length
            Delete_ind_delayed = Delete_ind_tot(j).Delete+i;
            % Ensure no delayed indices are negative or above upper limit
            Delete_ind_delayed(Delete_ind_delayed<=0) = [];
            Delete_ind_delayed(Delete_ind_delayed>size(EEG_data,2)) = [];
            % Remove outlier delayed indices
            EEG_data(j,Delete_ind_delayed) = 0;
            % EEG_mask(j,Delete_ind_delayed) = 0;
        end
    end
    
    EEG_sub = EEG_data(:,:); % Assign masked EEG_data back to EEG_sub
    
    close all % Call for closing any previous open MATLAB figures
    % For each channel plot premasked and masked data for visualization
    for i=1:6
        figure;plot(EEG_temp(i,:)) % EEG temp represents 'precleaned' data
        hold on
        plot(EEG_sub(i,:)) % EEG sub represents 'cleaned' data
        hold off
    end
    
    % Normalize each ECoG channel by channels' standard deviations 
    EEG_sub = EEG_sub./std(EEG_sub,0,2);
    % Ensure no indices are NaN
    EEG_sub(isnan(EEG_sub)) = 0;
    % Add ECoG channels back to EEG_vis variable 
    EEG_vis.data(20:25,:) = EEG_sub;
    % Save EEG_Vis with cleaned and normalized ECoG channels as .set file
    EEG_vis = pop_saveset(EEG_vis, 'filename',[Session_names(ss).name,'_Visualization_Delta_ECoGcleaned_vII.set'],'filepath',current_path);
    
    % For visualization purposes only, zero EEG chans and plot only ECoG
    EEG_vis.data(1:19,:) = 0;
    pop_eegplot(EEG_vis, 1, 1, 1);

    % Exit Session Specific directory
    cd ..
end


