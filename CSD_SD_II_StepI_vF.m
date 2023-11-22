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
    
    EEG_sub = EEG_vis.data(end-5:end,:);
    Limit_vec = ones(size(EEG_sub,1),2);%[1e5,1.5e5];
    Limit_vec(:,2) = Limit_vec(:,2)*size(EEG_sub,2);
    
    EEG_data = EEG_sub;
    EEG_temp = EEG_sub;
    Delete_ind_tot = struct;
    for j=1:6
        j
        temp = EEG_sub(j,EEG_sub(j,:)~=0);
        TF_temp = find(EEG_sub(j,:)~=0);
        [~,TF] = rmoutliers(temp,'quartiles','ThresholdFactor',10);
        T_ind = TF_temp(TF);
        Delete_ind = [];
        Delete_ind = cat(2,Delete_ind,T_ind);
        Delete_ind_tot(j).Delete = unique(Delete_ind);
    end
    
    
    w_length = 100;
    for j=1:6
        for i=-w_length:w_length
            Delete_ind_delayed = Delete_ind_tot(j).Delete+i;
            Delete_ind_delayed(Delete_ind_delayed<=0) = [];
            Delete_ind_delayed(Delete_ind_delayed>size(EEG_data,2)) = [];
            EEG_data(j,Delete_ind_delayed) = 0;
            EEG_mask(j,Delete_ind_delayed) = 0;
        end
    end
    
    EEG_sub = EEG_data(:,:);
    
    close all
    for i=1:6
        figure;plot(EEG_temp(i,:))
        hold on
        plot(EEG_sub(i,:))
        hold off
    end
    
    
    EEG_sub = EEG_sub./std(EEG_sub,0,2);
    EEG_sub(isnan(EEG_sub)) = 0;
    
    EEG_vis.data(20:25,:) = EEG_sub;
    
    EEG_vis = pop_saveset(EEG_vis, 'filename',[Session_names(ss).name,'_Visualization_Delta_ECoGcleaned_vII.set'],'filepath',current_path);
    EEG_vis.data(1:19,:) = 0;
    pop_eegplot(EEG_vis, 1, 1, 1);
    cd ..
end


