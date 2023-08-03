close all
clear all
clc
%%
% This code does the impedace masking, interpolation, power envelope extraction, 
%                   and depression detection steps for the EEG data in WAVEFRONT
%
% See: README.txt and [1] for more info.

% [1] Alireza Chamanzar, Jonathan Elmer, Lori Shutter, Jed Hartings, Pulkit Grover,
%  "Noninvasive and reliable automated detection of spreading depolarization in severe traumatic brain injury using scalp EEG",
%   Submitted to Nature Comms Med, 2022.

% Author: Alireza 	Date: 2022/01/04 12:00:00 	Revision: 0.1
% Copyright: This project is licensed - see the LICENSE.md file for details

%% Preprocessing steps for EEG data:
eegpath = 'C:\Users\Alireza\Downloads\eeglab_current\eeglab2020_0';
current_path = pwd;
cd(eegpath)
eeglab
cd(current_path)

Session_names = dir('Patient*');
Sessions_size = size(Session_names,1);


for ss = 1:Sessions_size
    cd(Session_names(ss).name)
    current_path = pwd;

    EEG_vis = pop_loadset('filename',[Session_names(ss).name,'_Visualization_Delta.set'],'filepath',current_path);
    
    event_delete_ind = [];
    for event_ind=1:size(EEG_vis.event,2)
        test_ind = strfind(string(EEG_vis.event(event_ind).type),'Data Quality');
        test_ind = cat(1,test_ind,strfind(string(EEG_vis.event(event_ind).type),'Detection'));
        test_ind = cat(1,test_ind,strfind(string(EEG_vis.event(event_ind).type),'boundary'));
        test_ind = cat(1,test_ind,strfind(string(EEG_vis.event(event_ind).type),'Electrode'));
        
        target_detected = 0;
        for ss_ind = 1:size(test_ind,1)
            if(size(test_ind,1)>1)
                if(~isempty(test_ind{ss_ind,1}))
                    target_detected = 1;
                end
            else
                if(~isempty(test_ind))
                    target_detected = 1;
                end
            end
        end
        
        if(target_detected==1)
            event_delete_ind = cat(1,event_delete_ind,event_ind);
            %             EEG_vis.event(event_ind) = [];
            continue;
        end
        
    end
    
%     pop_eegplot(EEG_vis, 1, 1, 1)
    
    EEG_vis = pop_editeventvals(EEG_vis,'delete',event_delete_ind);
    
    
    
    X = [EEG_vis.chanlocs.X];
    Y = [EEG_vis.chanlocs.Y];
    Z = [EEG_vis.chanlocs.Z];
    
    
    EEG = EEG_vis;
    EEG_vis_PW = EEG_vis;
    %loading the corresponding impedance data:
    load([Session_names(ss).name,'_woGamma_preica_Imp_Delta.mat'])
    
    srate = EEG.srate;
    
    %%
    WL = 5*60*srate;
    CSD_GT_total = zeros(1,size(EEG.data,2));
    Event_total = zeros(1,size(EEG.data,2));
    CSD_Labels_total = [];
    Event_Labels_total = [];
    for i=1:size(EEG.event,2)
        if(contains(EEG.event(i).type,'CSD') || contains(EEG.event(i).type,'CSD/ISD') ||...
                contains(EEG.event(i).type,'ISD') || contains(EEG.event(i).type,'scCSD'))
            EEG.event(i).latency = EEG.event(i).latency;% + WL-1;
            CSD_GT_total(round(max(EEG.event(i).latency,1)))=1;
            CSD_Labels_total = cat(1,CSD_Labels_total,{EEG.event(i).type});
        else
            EEG.event(i).latency = EEG.event(i).latency;% + WL-1;
            Event_total(round(max(EEG.event(i).latency,1)))=1;
            Event_Labels_total = cat(1,Event_Labels_total,{EEG.event(i).type});
        end
    end
    
    Block_size = 240*60*srate+1;
    Step_size = 180*60*srate;
    EEG = EEG.data(1:19,:);
    
    if(size(EEG,2)<Block_size)
        Block_size = size(EEG,2);
    end
    
    M_global_tot = zeros(19,1);
    M_global_Count = zeros(19,1);
    valid_ind_Prev = zeros(19,1);
    EEG_PW_Prev = zeros(19,1);
    for CSD_ind=1:Step_size:size(EEG,2)-Block_size+1
        
        T_start = CSD_ind;
        T_end = min(CSD_ind+Block_size-1,size(EEG,2));
        
        if(abs(size(EEG,2)-Block_size+1-CSD_ind)<180*60*srate)
            T_end = size(EEG,2);
        end
        
        EEG_mask = (measurement_data_Imp_temp_copy>0);
        CSD_GT = CSD_GT_total;
        Events = Event_total;

        
        EEG_sub = EEG(:,T_start:T_end);
        EEG_mask = EEG_mask(:,T_start:T_end);
        CSD_GT = CSD_GT(1,T_start:T_end);
        Events = Events(1,T_start:T_end);
        

        EEG_sub(EEG_mask==0) = 0;
        
        
        CSD_type_counter = find(CSD_GT_total);
        CSD_type_ind_strt = min(find(CSD_type_counter>=T_start));
        CSD_type_ind_Endd = max(find(CSD_type_counter<=T_end));
        CSD_type = [];
        if(~isempty(CSD_type_ind_strt) && ~isempty(CSD_type_ind_Endd))
            if(CSD_type_ind_strt<=CSD_type_ind_Endd)
                CSD_type = CSD_Labels_total(CSD_type_ind_strt:CSD_type_ind_Endd);
            end
        end
        
        
        Event_type_counter = find(Event_total);
        Event_type_ind_strt = min(find(Event_type_counter>=T_start));
        Event_type_ind_Endd = max(find(Event_type_counter<=T_end));
        Events_type = [];
        if(~isempty(Event_type_ind_strt) && ~isempty(Event_type_ind_Endd))
            if(Event_type_ind_strt<=Event_type_ind_Endd)
                Events_type = Event_Labels_total(Event_type_ind_strt:Event_type_ind_Endd);
            end
        end
        
        EEG_data = EEG_sub;
        EEG_temp = EEG_sub;
        EEG_sub = EEG_data(1:19,:);
        
        
        %Extract large connected zero intervals in the data with normal
        %impedance:
        
        CC_WL = 10*srate;
        for Ch = 1:19
            EEG_data_ConComp = bwconncomp(EEG_sub(Ch,:)==0);
            Strt = [];
            Endd = [];
            
            for i=1:size(EEG_data_ConComp.PixelIdxList,2)
                if(size(EEG_data_ConComp.PixelIdxList{1,i},1)>CC_WL)
                     EEG_mask(Ch,EEG_data_ConComp.PixelIdxList{1,i}) = 0;
                end
            end
        end
        
        %Remove small and isolated connected valid portions:
        EEG_mask_gaps = EEG_mask; 
        CC_WL = 20*60*srate;
        for Ch = 1:19
            EEG_mask_ConComp = bwconncomp(EEG_mask(Ch,:));
            Strt = [];
            Endd = [];
            
            for i=2:size(EEG_mask_ConComp.PixelIdxList,2)-1
                if(size(EEG_mask_ConComp.PixelIdxList{1,i},1)<CC_WL)
                    if(EEG_mask_ConComp.PixelIdxList{1,i+1}(1)-EEG_mask_ConComp.PixelIdxList{1,i}(end)<0.05*CC_WL)
                        EEG_mask_gaps(Ch,EEG_mask_ConComp.PixelIdxList{1,i}(end):EEG_mask_ConComp.PixelIdxList{1,i+1}(1)) = 1;
                    end
                    if(EEG_mask_ConComp.PixelIdxList{1,i}(1)-EEG_mask_ConComp.PixelIdxList{1,i-1}(end)<0.05*CC_WL)
                        EEG_mask_gaps(Ch,EEG_mask_ConComp.PixelIdxList{1,i-1}(end):EEG_mask_ConComp.PixelIdxList{1,i}(1)) = 1;
                    end
                end
            end
            
            if(size(EEG_mask_ConComp.PixelIdxList,2)<2)
                continue;
            end
            
            if(size(EEG_mask_ConComp.PixelIdxList{1,1},1)<CC_WL)
                if(EEG_mask_ConComp.PixelIdxList{1,2}(1)-EEG_mask_ConComp.PixelIdxList{1,1}(end)<0.05*CC_WL)
                    EEG_mask_gaps(Ch,EEG_mask_ConComp.PixelIdxList{1,1}(end):EEG_mask_ConComp.PixelIdxList{1,2}(1)) = 1;
                end
            end
            if(size(EEG_mask_ConComp.PixelIdxList{1,end},1)<CC_WL)
                if(EEG_mask_ConComp.PixelIdxList{1,end}(1)-EEG_mask_ConComp.PixelIdxList{1,end-1}(end)<0.05*CC_WL)
                    EEG_mask_gaps(Ch,EEG_mask_ConComp.PixelIdxList{1,end-1}(end):EEG_mask_ConComp.PixelIdxList{1,end}(1)) = 1;
                end
            end
        end
        
        
        %Remove small and isolated connected valid portions:
        CC_WL = 20*60*srate;
        for Ch = 1:19
            EEG_mask_ConComp = bwconncomp(EEG_mask_gaps(Ch,:));
            Strt = [];
            Endd = [];
            
            for i=1:size(EEG_mask_ConComp.PixelIdxList,2)
                if(size(EEG_mask_ConComp.PixelIdxList{1,i},1)<CC_WL)
                        EEG_mask(Ch,EEG_mask_ConComp.PixelIdxList{1,i}) = 0;
                        EEG_mask_gaps(Ch,EEG_mask_ConComp.PixelIdxList{1,i}) = 0;
                end
            end
        end
        
        
        EEG_sub(EEG_mask==0) = 0;
        
        %figure;stackedplot(1:size(EEG_sub,2),EEG_sub')
        
        
        EEG_PW_tot = zeros(size(EEG_sub));
        for Ch=1:19
            
            EEG_mask_ConComp = bwconncomp(EEG_mask_gaps(Ch,:));
            
            for i=1:size(EEG_mask_ConComp.PixelIdxList,2)
                
                EEG_sub_sub = EEG_sub(Ch,EEG_mask_ConComp.PixelIdxList{1,i});
                EEG_mask_sub_sub = EEG_mask(Ch,EEG_mask_ConComp.PixelIdxList{1,i});
                EEG_PW_tot_temp = EEG_PW_tot(Ch,EEG_mask_ConComp.PixelIdxList{1,i});
                
                
                EEG_PW = EEG_sub_sub(EEG_mask_sub_sub==1).^2;
                EEG_mask_sub = EEG_mask_sub_sub(EEG_mask_sub_sub==1);
                if(size(EEG_PW,2)<3*WL)
                    EEG_PW_tot(Ch,EEG_mask_ConComp.PixelIdxList{1,i}) = 0;
                    EEG_mask(Ch,EEG_mask_ConComp.PixelIdxList{1,i}) = 0;
                    continue;
                end
                
                EEG_PW_temp = EEG_PW;
                [u,l] = envelope(EEG_PW_temp',WL,'rms');
                EEG_PW = u';
                EEG_PW(:,size(EEG_PW,2)-WL+1:end) = 0;
                EEG_PW(:,1:WL-1) = 0;
                EEG_mask_sub(:,size(EEG_mask_sub,2)-WL+1:end) = 0;
                EEG_mask_sub(:,1:WL-1) = 0;
                EEG_PW_temp = EEG_PW(:,WL:size(EEG_PW,2)-WL);
                
                EEG_PW_tot_temp(EEG_mask_sub_sub==1) = EEG_PW;
                EEG_mask_sub_sub(EEG_mask_sub_sub==1) = EEG_mask_sub;
                
                EEG_PW_tot(Ch,EEG_mask_ConComp.PixelIdxList{1,i}) = EEG_PW_tot_temp;
                EEG_mask(Ch,EEG_mask_ConComp.PixelIdxList{1,i}) = EEG_mask_sub_sub;
            end
        end
        
%         figure;stackedplot(1:size(EEG_sub,2),EEG_sub')
%         figure;stackedplot(1:size(EEG_sub,2),EEG_PW_tot')

        
        long_invalid_WL = 1*60*srate;
   
        EEG_mask_long_invalid = double(EEG_mask);
        for Ch = 1:19
            EEG_mask_ConComp = bwconncomp(EEG_mask_long_invalid(Ch,:)==0);
            Strt = [];
            Endd = [];
            for i=1:size(EEG_mask_ConComp.PixelIdxList,2)
                if(size(EEG_mask_ConComp.PixelIdxList{1,i},1)>long_invalid_WL)
                    EEG_mask_long_invalid(Ch,EEG_mask_ConComp.PixelIdxList{1,i}) = -1;
                end
            end
        end
%             
        
        for Ch = 1:19
            valid_ind = find(EEG_mask(Ch,:)==1);
            invalid_ind = find(EEG_mask_long_invalid(Ch,:)==0);
            %Mean adjust due to introduced offset in RMS envelope:
            Overlap_CK = find(ismember(valid_ind_Prev(Ch)-T_start+1,valid_ind));
            if(CSD_ind>1 && ~isempty(Overlap_CK))
                Adjust_mean = EEG_PW_Prev(Ch) - EEG_PW_tot(Ch,valid_ind_Prev(Ch)-T_start+1);
                EEG_PW_tot(Ch,:) = EEG_PW_tot(Ch,:) + Adjust_mean;
            end
            if(~isempty(valid_ind))
                valid_ind_Prev(Ch) = valid_ind(end)+T_start-1;
                EEG_PW_Prev(Ch) = EEG_PW_tot(Ch,valid_ind(end));
                EEG_PW_tot(Ch,invalid_ind) = ...
                    interp1(valid_ind,EEG_PW_tot(Ch,valid_ind),invalid_ind,'nearest');
                EEG_mask(Ch,:) = 1;
                EEG_mask(Ch,isnan(EEG_PW_tot(Ch,:))) = 0;
                EEG_mask(Ch,EEG_mask_long_invalid(Ch,:)==-1) = 0;
                EEG_PW_tot(Ch,isnan(EEG_PW_tot(Ch,:))) = 0;
                EEG_PW_tot(Ch,EEG_mask_long_invalid(Ch,:)==-1) = 0;
            else
                EEG_mask(Ch,:) = 0;
                valid_ind_Prev(Ch) = 0;
                EEG_PW_Prev(Ch) = 0;
            end
        end
        
        
      
        %figure;stackedplot(1:size(EEG_sub,2),EEG_PW_tot')


%       Cross Correlation:
        WL_xcrr = 5*60*srate;
        Ptrn = [(1/WL_xcrr)*ones(1,WL_xcrr),zeros(1,WL),-(1/WL_xcrr)*ones(1,WL_xcrr)];%*ones(1,WL_xcrr)];
        for i=1:size(EEG_PW_tot,1)
            
            EEG_PW_tot_ConComp = bwconncomp(EEG_mask(i,:)==0);
            Strt = [];
            Endd = [];
            for j=1:size(EEG_PW_tot_ConComp.PixelIdxList,2)
                if(size(EEG_PW_tot_ConComp.PixelIdxList{1,j},1)>0.1*WL_xcrr)
                    Strt = cat(1,Strt,EEG_PW_tot_ConComp.PixelIdxList{1,j}(1,1));
                    Endd = cat(1,Endd,EEG_PW_tot_ConComp.PixelIdxList{1,j}(end,1));
                end
            end

            [c,lags] = xcorr(EEG_PW_tot(i,:),Ptrn);
            EEG_PW_tot(i,:) = c(find(lags==0):end);
            EEG_PW_tot(i,:) = [zeros(1,2*WL_xcrr+WL),EEG_PW_tot(i,1:size(EEG_PW_tot,2)-(2*WL_xcrr+WL))];
            
            for j=1:size(Strt,1)
%                 EEG_PW_tot(i,max(Strt(j)-0.5*size(Ptrn,2),1):min(Strt(j)+0.5*size(Ptrn,2)-1,size(EEG_PW_tot,2))) = 0;
                EEG_PW_tot(i,max(Strt(j)-(2*WL_xcrr+WL),1):min(Endd(j)+(2*WL_xcrr+WL),size(EEG_PW_tot,2))) = 0;
                EEG_mask(i,max(Strt(j)-(2*WL_xcrr+WL),1):min(Endd(j)+(2*WL_xcrr+WL),size(EEG_PW_tot,2))) = 0;
            end 
        end
        
        
        %Remove small and isolated connected valid portions:
        for Ch = 1:19
            EEG_mask_ConComp = bwconncomp(EEG_mask(Ch,:));
            Strt = [];
            Endd = [];
            
            for i=2:size(EEG_mask_ConComp.PixelIdxList,2)-1
                if(size(EEG_mask_ConComp.PixelIdxList{1,i},1)<CC_WL)
                    if(EEG_mask_ConComp.PixelIdxList{1,i+1}(1)-EEG_mask_ConComp.PixelIdxList{1,i}(end)>0.05*CC_WL...
                            && EEG_mask_ConComp.PixelIdxList{1,i}(1)-EEG_mask_ConComp.PixelIdxList{1,i-1}(end)>0.05*CC_WL )
                        EEG_mask(Ch,EEG_mask_ConComp.PixelIdxList{1,i}) = 0;
                    end
                end
            end
            
            if(size(EEG_mask_ConComp.PixelIdxList,2)<2)
                continue;
            end
            
            if(size(EEG_mask_ConComp.PixelIdxList{1,1},1)<CC_WL)
                if(EEG_mask_ConComp.PixelIdxList{1,2}(1)-EEG_mask_ConComp.PixelIdxList{1,1}(end)>0.05*CC_WL)
                    EEG_mask(Ch,EEG_mask_ConComp.PixelIdxList{1,1}) = 0;
                end
            end
            if(size(EEG_mask_ConComp.PixelIdxList{1,end},1)<CC_WL)
                if(EEG_mask_ConComp.PixelIdxList{1,end}(1)-EEG_mask_ConComp.PixelIdxList{1,end-1}(end)>0.05*CC_WL)
                    EEG_mask(Ch,EEG_mask_ConComp.PixelIdxList{1,end}) = 0;
                end
            end
        end
        
        EEG_PW_tot(EEG_mask==0) = 0;
        
%         figure;stackedplot(1:size(EEG_sub,2),EEG_PW_tot')


        EEG_vis_temp = EEG_PW_tot;
        for i=1:19
            Ch_inter_dist = (X(i)-X).^2+(Y(i)-Y).^2+(Z(i)-Z).^2;
            [~,dist_ind] = sort(Ch_inter_dist);
            EEG_PW_temp =  EEG_PW_tot(i,EEG_mask(i,:)==1);
            M_global_tot(i) = M_global_tot(i) + sum(EEG_PW_temp,2);
            M_global_Count(i) = M_global_Count(i) + size(EEG_PW_temp,2);
        end
        
        
        
        EEG_PW_tot(EEG_mask==0) = 0;
% %         
%         figure;stackedplot(1:size(EEG_sub,2),EEG_sub')
%         figure;stackedplot(1:size(EEG_sub,2),EEG_PW_tot') 

        EEG_PW = EEG_PW_tot;
           
        EEG_elec = 19;
%         close all
%         
%         figure;stackedplot(1:size(EEG_sub,2),EEG_sub') 
%         figure;stackedplot(1:size(EEG_sub,2),EEG_PW')
        
        save([Session_names(ss).name ,'_',sprintf('CSD_ind_%d',CSD_ind), '_Delta_SDII_vConfusion_240minOv180min_CCWL20min_Xcrr_Pruned_vVII.mat'], 'EEG_PW', 'EEG_sub', 'EEG_mask', 'CSD_type', 'Events_type', 'Events', 'CSD_GT', 'srate', 'CSD_GT_total','CSD_Labels_total','T_start','T_end', '-v7.3');
        close all
        
        clear CSD_type Events_type
    end
    Global_M = M_global_tot./M_global_Count;
    save([Session_names(ss).name ,'GlobalMean_Delta_SDII_vConfusion_240minOv180min_CCWL20min_Xcrr_Improved_Pruned_vVII.mat'],'Global_M')
    cd ..
end



