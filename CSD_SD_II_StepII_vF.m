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

Session_names = dir('Patient*');
Sessions_size = size(Session_names,1);


for ss = 1: Sessions_size
    cd(Session_names(ss).name)
    current_path = pwd;
    
    
    EEG_vis = pop_loadset('filename',[Session_names(ss).name,'_Visualization_Delta_ECoGcleaned_vII.set'],'filepath',current_path);
    EEG_vis = pop_select( EEG_vis,'nochannel',{'Depth1' 'Depth2' 'Depth3' 'Depth4' 'Depth5' 'Depth6'});
    
       
    
    X = [EEG_vis.chanlocs.X];
    Y = [EEG_vis.chanlocs.Y];
    Z = [EEG_vis.chanlocs.Z];
    
    
    EEG = EEG_vis;
    load([Session_names(ss).name,'_woGamma_preica_Imp.mat'])
    
    % Add the ECoG channels to EEG:
    for Imp_ind=1:19
        EEG.data(end+1,:) = measurement_data_Imp_temp_copy(Imp_ind,:);
        EEG.chanlocs(end+1).labels = sprintf('Imp_%d',Imp_ind);
    end
    EEG.nbchan = size(EEG.data,1);
    
    %% Remove poor data equality events:
    
    Remove_ind = find(strcmp({EEG.event(:).type},'Electrode Impedance High'));
    
    
    Norm_ind = find(strcmp({EEG.event(:).type},'Data Quality Normal'));
    
    Remove_latency = [EEG.event(Remove_ind).latency];
    Norm_latency= [EEG.event(Norm_ind).latency];
    
    eerej_vec = [-1,-1];%dummy value
    for ee=1:size(Remove_latency,2)
        if(Remove_latency(ee)>eerej_vec(end,2))
            ind_temp = find(Norm_latency>Remove_latency(ee));
            if(~isempty(ind_temp))
                eerej_vec = cat(1,eerej_vec,[Remove_latency(ee),Norm_latency(ind_temp(1))]);
            else
                eerej_vec = cat(1,eerej_vec,[Remove_latency(ee),EEG.pnts]);
            end
        end
    end
    

    measurement_data_Imp_temp_copy = EEG.data(end-18:end,:);
    EEG = pop_select( EEG,'nochannel',size(EEG.data,1)-18:size(EEG.data,1));
    
    srate = EEG.srate;
    
    EEG_vis = EEG;
    
    %%
    
    WL = 10*60*srate;
    
    Block_size = 4*60*60*srate+1;
    EEG = EEG.data(1:19,:);
    
    if(size(EEG,2)<Block_size)
        Block_size = size(EEG,2);
    end
    
    M_global_tot = zeros(19,1);
    M_global_Count = zeros(19,1);
    for CSD_ind=1:3.5*60*60*srate:size(EEG,2)-Block_size+1
        
        T_start = CSD_ind;
        T_end = min(CSD_ind+Block_size-1,size(EEG,2));
        
        if(abs(size(EEG,2)-Block_size+1-CSD_ind)<3.5*60*60*srate)
            T_end = size(EEG,2);
        end
        
        EEG_mask = (measurement_data_Imp_temp_copy>0);
        EEG_sub = EEG(:,T_start:T_end);
        EEG_mask = EEG_mask(:,T_start:T_end);
        
        
        EEG_sub(EEG_mask==0) = 0;
        
        EEG_data = EEG_sub;
        EEG_temp = EEG_sub;
        Delete_ind_tot = struct;
        parfor j=1:19
            j
            temp = EEG_sub(j,EEG_mask(j,:)>0);
            TF_temp = find(EEG_mask(j,:)>0);
            std_max = (std(temp,[],2));
            if(std_max<=0 || isempty(temp))
                std_max=30;
            end
            [~,TF] = rmoutliers(temp,'quartiles','ThresholdFactor',3);
            T_ind = TF_temp(TF);
            Delete_ind = [];
            Delete_ind = cat(2,Delete_ind,T_ind);
            Delete_ind_tot(j).Delete = unique(Delete_ind);
        end
        
        w_length = 1*64;
        for j=1:19
            for i=-w_length:w_length
                Delete_ind_delayed = Delete_ind_tot(j).Delete+i;
                Delete_ind_delayed(Delete_ind_delayed<=0) = [];
                Delete_ind_delayed(Delete_ind_delayed>size(EEG_data,2)) = [];
                EEG_data(j,Delete_ind_delayed) = 0;
                EEG_mask(j,Delete_ind_delayed) = 0;
            end
        end
        
        EEG_sub = EEG_data;
        for j=1:19
            EEG_temp = EEG_data(j,EEG_mask(j,:)>0);
            if(~isempty(EEG_temp))
                EEG_sub(j,:) = EEG_sub(j,:)./std(EEG_temp,0,2);
            end
        end
        EEG_sub(isnan(EEG_sub)) = 0;
        EEG_sub(EEG_mask==0) = 0;
        
        
        EEG_vis.data(1:19,T_start:T_end) = EEG_sub;
        measurement_data_Imp_temp_copy(1:19,T_start:T_end) = EEG_mask;
        
        clear CSD_type Events_type
    end
    
    EEG_vis = pop_saveset(EEG_vis, 'filename',[Session_names(ss).name,'_Visualization_Delta.set'],'filepath',current_path);
    save([Session_names(ss).name,'_woGamma_preica_Imp_Delta.mat'],'measurement_data_Imp_temp_copy','-v7.3')
    
    cd ..
end



