close all
clear all
clc
%%
% This code downsamples the dataset, concatenates small parts of each EEG session, prunes EEG based on impedance:
%
% See: README.txt and [1] for more info.

% [1] Alireza Chamanzar, Jonathan Elmer, Lori Shutter, Jed Hartings, Pulkit Grover,
%  "Noninvasive, automated and reliable detection of spreading depolarizations in 
%                                   severe traumatic brain injury using scalp EEG",
%   Submitted to Nature Comms Med, 2022.

% Author: Alireza 	Date: 2021/12/19 12:00:48 	Revision: 0.1
% Copyright: This project is licensed - see the LICENSE.md file for details

%% Preprocessing steps for EEG data:
eegpath = 'C:\Users\Alireza\Downloads\eeglab_current\eeglab2020_0';
% cd(eegpath)
% eeglab
% cd(current_path)

Session_names = dir('Patient*');
Sessions_size = size(Session_names,1);

load D:\SD-II_POST_ICA\SD-II\median_Imp.mat
median_Imp_tot = median_Imp;


for ss = 1:Sessions_size
    cd(Session_names(ss).name)
    EEG_Imp = dir('EEG,Composite,Impedance,*.mat')
    load(EEG_Imp(1).name)
    
    comp_elements_Imp = comp_elements;
    measurement_data_Imp = measurement_data;
    time_vector_Imp = time_vector;
    Imp_srate = 1/(time_vector_Imp(2) - time_vector_Imp(1));
    
    for Imp_ind = 2:size(EEG_Imp,1)
    load(EEG_Imp(Imp_ind).name)
    measurement_data_Imp = cat(2,measurement_data_Imp,measurement_data);
    time_vector_Imp = cat(2,time_vector_Imp,time_vector);
    end 
    current_path = pwd;
    ECoG_Part_names = dir('*_ECoG_filtered_withDC.set');
    
    current_path = pwd;
    Part_names = dir('*_preica_withDC.set');
    Part_name_char = string({Part_names.name});
    %Remove possible files in the directory to prevent nested merging process:
    part_woGamma_ind = find(contains(Part_name_char,'_woGamma_preica_withDC.set'));
    Part_name_char(part_woGamma_ind) = [];
    Part_names(part_woGamma_ind) = [];
    
    Part_size = size(Part_names,1);
    part_ind = strfind(Part_name_char,'part');
    of_ind = strfind(Part_name_char,'of');
    preica_ind = strfind(Part_name_char,'_preica');
    Part_num = [];
    for i=1:Part_size
        Part_num = cat(1,Part_num,str2double(Part_name_char{i}(part_ind{i}+4:of_ind{i}-1)));
    end
    AllEEG = [];
    [~,Part_ind] = sort(Part_num);
    measurement_data_Imp_temp = [];
    
    EEG_test = pop_loadset('filename',Part_names(Part_ind(1)).name,'filepath',current_path);
    EEG_test = pop_select( EEG_test,'nochannel',{'Depth1' 'Depth2' 'Depth3' 'Depth4' 'Depth5' 'Depth6'});
    Imp_ind = [];
    for el = 1:19
        Imp_ind = cat(1,Imp_ind,find(contains(comp_elements_Imp,EEG_test.chanlocs(el).labels)));
    end
    measurement_data_Imp = measurement_data_Imp(Imp_ind,:);
    [time_vector_Imp,ind] = unique(time_vector_Imp);
    measurement_data_Imp = measurement_data_Imp(:,ind);
    
    for pp = 1: Part_size
        Part_names(Part_ind(pp)).name
        
        EEG = pop_loadset('filename',Part_names(Part_ind(pp)).name,'filepath',current_path);
        EEG = pop_select( EEG,'nochannel',{'Depth1' 'Depth2' 'Depth3' 'Depth4' 'Depth5' 'Depth6'});
        ECoG = pop_loadset('filename',ECoG_Part_names(Part_ind(pp)).name,'filepath',current_path);
     
        EEG = pop_eegfiltnew(EEG, [],30,1000,0,[],0);
        [EEG_resample] = pop_resample(EEG, 64);
        EEG = EEG_resample;
        
        ECoG = pop_eegfiltnew(ECoG, [],30,1000,0,[],0);
        [ECoG_resample] = pop_resample(ECoG, 64);
        ECoG = ECoG_resample;

        %% Remove poor data equality events:
        Remove_ind = find(strcmp({EEG.event(:).type},'Electrode Impedance Uncertain') | ...
            strcmp({EEG.event(:).type},'Electrode Impedance High') | ...
            strcmp({EEG.event(:).type},'Electrode Impedance Abnormal Low') | ...
            strcmp({EEG.event(:).type},'Electrode Impedance Imbalance'));
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
        eerej_vec(1,:) = [];%dummy value        
        Imp_strt = find(time_vector_Imp>=EEG.xmin-10);
        Imp_strt = Imp_strt(1);
        Imp_end = find(time_vector_Imp<=EEG.xmax+10);
        Imp_end = Imp_end(end);
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




