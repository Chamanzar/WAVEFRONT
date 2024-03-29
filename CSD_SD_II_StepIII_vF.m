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

% Similar to previous steps, ensure eeglab is running and initialize path
eegpath = 'C:\Users\Alireza\Downloads\eeglab_current\eeglab2020_0';
current_path = pwd;
cd(eegpath)
eeglab
cd(current_path)

% Find all Sessions of chosen Patient ID
Session_names = dir('Patient*');
Sessions_size = size(Session_names,1);

% Loop through all sessions of Patient ID
for ss = 1:Sessions_size

    % Enter Session-specific directory 
    cd(Session_names(ss).name)
    current_path = pwd;
    
    % Load 'Visualization' file which contains Session-wide EEG & ECoG data
    EEG_vis = pop_loadset('filename',[Session_names(ss).name,'_Visualization_Delta.set'],'filepath',current_path);
    
    % Code block for processing events
    event_delete_ind = []; % Initialize array for storing indices to delete
    for event_ind=1:size(EEG_vis.event,2) % Loop through all events
        % Check each string for keywords indicating events to be filtered
        test_ind = strfind(string(EEG_vis.event(event_ind).type),'Data Quality');
        test_ind = cat(1,test_ind,strfind(string(EEG_vis.event(event_ind).type),'Detection'));
        test_ind = cat(1,test_ind,strfind(string(EEG_vis.event(event_ind).type),'boundary'));
        test_ind = cat(1,test_ind,strfind(string(EEG_vis.event(event_ind).type),'Electrode'));
        
        target_detected = 0; % Set target variable to zero
        for ss_ind = 1:size(test_ind,1)
            % If any target keywords are found, set target detected to 1
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
        
        % If keyword found & target set to 1, store event ind to be deleted
        if(target_detected==1)
            event_delete_ind = cat(1,event_delete_ind,event_ind);
            %             EEG_vis.event(event_ind) = [];
            continue;
        end
        % If no keyword found, 'keep' event by not deleting
    end
    
    %pop_eegplot(EEG_vis, 1, 1, 1) % optional call for plotting data
    
    % Delete events at every index found above
    EEG_vis = pop_editeventvals(EEG_vis,'delete',event_delete_ind);
     
    % X = [EEG_vis.chanlocs.X];
    % Y = [EEG_vis.chanlocs.Y];
    % Z = [EEG_vis.chanlocs.Z];
    
    % Create copies of EEG_vis 
    EEG = EEG_vis;
    EEG_vis_PW = EEG_vis;

    % Loading the corresponding impedance data
    load([Session_names(ss).name,'_woGamma_preica_Imp_Delta.mat'])
   
    srate = EEG.srate; % Set srate as defined from EEG (from EEG_vis)
    WL = 5*60*srate; % Set window length of 5 minutes 
    
    %%
    % This code seperates CSD events from other non-removed events  

    % Initalize vectors for storing event times
    CSD_GT_total = zeros(1,size(EEG.data,2)); 
    Event_total = zeros(1,size(EEG.data,2));

    % Initalize vectors for storing event strings 
    CSD_Labels_total = [];
    Event_Labels_total = [];

    % Loop through all events which were not filtered out in previous step
    for i=1:size(EEG.event,2)
        % If event contains any of keywords, store info to CSD variables
        if(contains(EEG.event(i).type,'CSD') || contains(EEG.event(i).type,'CSD/ISD') ||...
                contains(EEG.event(i).type,'ISD') || contains(EEG.event(i).type,'scCSD'))
            % EEG.event(i).latency = EEG.event(i).latency;% + WL-1;
            CSD_GT_total(round(max(EEG.event(i).latency,1)))=1;
            CSD_Labels_total = cat(1,CSD_Labels_total,{EEG.event(i).type});
        else
            % EEG.event(i).latency = EEG.event(i).latency;% + WL-1;
            Event_total(round(max(EEG.event(i).latency,1)))=1;
            Event_Labels_total = cat(1,Event_Labels_total,{EEG.event(i).type});
        end
    end
    
    % Beging 'Part-specific' data processing block (similar to StepII)

    Block_size = 240*60*srate+1; % Set block length to 4 hours
    Step_size = 180*60*srate; % Set step length to 3 hours
    EEG = EEG.data(1:19,:); % Set EEG to data array of only EEG data
    
    % Check for edge cases of full data less than 4 hours (or Block size)
    if(size(EEG,2)<Block_size)
        Block_size = size(EEG,2);
    end
    
    % Initalize Session wide variables for channel dependent information
    M_global_tot = zeros(19,1);
    M_global_Count = zeros(19,1);
    valid_ind_Prev = zeros(19,1);
    EEG_PW_Prev = zeros(19,1);
    
    % Create starting indices of 3 hour intervals 
    for CSD_ind=1:Step_size:size(EEG,2)-Block_size+1
        
        % Set Part Start as CSD_ind 
        T_start = CSD_ind;
        % Set Part End as CSD + 4 hrs, or, as 'end of data' in edge-cases
        T_end = min(CSD_ind+Block_size-1,size(EEG,2));
        
        % Check if CSD_ind is less than 3 hrs away from end of data
        % if this is true, extend 'final' Part beyond "Part_size"
        if(abs(size(EEG,2)-Block_size+1-CSD_ind)<180*60*srate)
            T_end = size(EEG,2);
        end
        
        % Set EEG mask as nonzero Imp data
        EEG_mask = (measurement_data_Imp_temp_copy>0);

        % Copy 'total' event data to trim in next step
        CSD_GT = CSD_GT_total;
        Events = Event_total;

        % Trim data, mask, and events to 'Part' time interval
        EEG_sub = EEG(:,T_start:T_end);
        EEG_mask = EEG_mask(:,T_start:T_end);
        CSD_GT = CSD_GT(1,T_start:T_end);
        Events = Events(1,T_start:T_end);
        
        % Mask data anywhere masking variable is equal to zero
        EEG_sub(EEG_mask==0) = 0;
        
        % Find the time indices of the CSD events
        CSD_type_counter = find(CSD_GT_total);

        % Determine if any of the CSD events are within the Part range
        CSD_type_ind_strt = min(find(CSD_type_counter>=T_start));
        CSD_type_ind_Endd = max(find(CSD_type_counter<=T_end));
        
        % If the CSD Event is in range, set CSD_type from CSD_Labels_tot
        CSD_type = [];
        if(~isempty(CSD_type_ind_strt) && ~isempty(CSD_type_ind_Endd))
            if(CSD_type_ind_strt<=CSD_type_ind_Endd)
                CSD_type = CSD_Labels_total(CSD_type_ind_strt:CSD_type_ind_Endd);
            end
        end
        
        % Analogous to above CSD event block, find non-CSD events in range
        Event_type_counter = find(Event_total);
        Event_type_ind_strt = min(find(Event_type_counter>=T_start));
        Event_type_ind_Endd = max(find(Event_type_counter<=T_end));
        Events_type = [];
        if(~isempty(Event_type_ind_strt) && ~isempty(Event_type_ind_Endd))
            if(Event_type_ind_strt<=Event_type_ind_Endd)
                Events_type = Event_Labels_total(Event_type_ind_strt:Event_type_ind_Endd);
            end
        end
        
        % Redeclare copies of EEG data 
        EEG_data = EEG_sub;
        EEG_temp = EEG_sub;
        EEG_sub = EEG_data(1:19,:);
        
        % Begining of channel dependent Connected Component Processing
        % The folowing steps are similar but each perform unique task
        
        %Extract large connected 'zero-intervals' in the data with normal
        %impedance:
        
        % Set 'Connect Component Window Length' as 10 seconds
        CC_WL = 10*srate;
        for Ch = 1:19 
            % At each channel, find intervals of comtinuous zeros
            EEG_data_ConComp = bwconncomp(EEG_sub(Ch,:)==0);
            Strt = [];
            Endd = [];
            % Loop through all regions and if region is longer than CC_WL
            % add indices to EEG_mask for eventual filtering 
            for i=1:size(EEG_data_ConComp.PixelIdxList,2)
                if(size(EEG_data_ConComp.PixelIdxList{1,i},1)>CC_WL)
                     EEG_mask(Ch,EEG_data_ConComp.PixelIdxList{1,i}) = 0;
                end
            end
        end
        
        %Remove small and isolated valid portions:
        % Save a copy of EEG_mask
        EEG_mask_gaps = EEG_mask; 
        % Reset CC_WL as 20 minutes for 'small valid portion' removal 
        CC_WL = 20*60*srate;
        for Ch = 1:19
            % Find regions of 'connectivity' in EEG_mask
            EEG_mask_ConComp = bwconncomp(EEG_mask(Ch,:));
            Strt = [];
            Endd = [];

            % Check if connected components regions are atleast 20 mins 
            for i=2:size(EEG_mask_ConComp.PixelIdxList,2)-1
                if(size(EEG_mask_ConComp.PixelIdxList{1,i},1)<CC_WL)
                    % If interval is less than 20 mins & not w/n a 1 min
                    % interval from neighboring regions, add region to mask
                    if(EEG_mask_ConComp.PixelIdxList{1,i+1}(1)-EEG_mask_ConComp.PixelIdxList{1,i}(end)<0.05*CC_WL)
                        EEG_mask_gaps(Ch,EEG_mask_ConComp.PixelIdxList{1,i}(end):EEG_mask_ConComp.PixelIdxList{1,i+1}(1)) = 1;
                    end
                    if(EEG_mask_ConComp.PixelIdxList{1,i}(1)-EEG_mask_ConComp.PixelIdxList{1,i-1}(end)<0.05*CC_WL)
                        EEG_mask_gaps(Ch,EEG_mask_ConComp.PixelIdxList{1,i-1}(end):EEG_mask_ConComp.PixelIdxList{1,i}(1)) = 1;
                    end
                end
            end
            
            % Check if more than 2 ROIs are found for Channel Ch 
            if(size(EEG_mask_ConComp.PixelIdxList,2)<2)
                continue;
            end
            
            % If >2 ROIs, perform edge case processing as above  
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
        
        CC_WL = 20*60*srate;
        for Ch = 1:19
            % Find connected regions of new mask from above
            EEG_mask_ConComp = bwconncomp(EEG_mask_gaps(Ch,:));
            Strt = [];
            Endd = [];
            
            % If ROI is less than 20 mins in length set region to be masked
            for i=1:size(EEG_mask_ConComp.PixelIdxList,2)
                if(size(EEG_mask_ConComp.PixelIdxList{1,i},1)<CC_WL)
                        EEG_mask(Ch,EEG_mask_ConComp.PixelIdxList{1,i}) = 0;
                        EEG_mask_gaps(Ch,EEG_mask_ConComp.PixelIdxList{1,i}) = 0;
                end
            end
        end
        
        % Perform masking using mask generated above
        EEG_sub(EEG_mask==0) = 0;
        
        %figure;stackedplot(1:size(EEG_sub,2),EEG_sub')
        
        % Initalize matrix for creating Power Enveloped data
        EEG_PW_tot = zeros(size(EEG_sub));
        for Ch=1:19
            
            % Find valid portions of data
            EEG_mask_ConComp = bwconncomp(EEG_mask_gaps(Ch,:));
            
            % Find the power of the signal at each valid region
            for i=1:size(EEG_mask_ConComp.PixelIdxList,2)
                
                EEG_sub_sub = EEG_sub(Ch,EEG_mask_ConComp.PixelIdxList{1,i});
                EEG_mask_sub_sub = EEG_mask(Ch,EEG_mask_ConComp.PixelIdxList{1,i});
                EEG_PW_tot_temp = EEG_PW_tot(Ch,EEG_mask_ConComp.PixelIdxList{1,i});
                
                
                EEG_PW = EEG_sub_sub(EEG_mask_sub_sub==1).^2;
                EEG_mask_sub = EEG_mask_sub_sub(EEG_mask_sub_sub==1);

                % If region is less than 15 mins, mask region
                if(size(EEG_PW,2)<3*WL)
                    EEG_PW_tot(Ch,EEG_mask_ConComp.PixelIdxList{1,i}) = 0;
                    EEG_mask(Ch,EEG_mask_ConComp.PixelIdxList{1,i}) = 0;
                    continue;
                end
                
                % Calculate the envelope from power of singal 
                EEG_PW_temp = EEG_PW;
                [u,l] = envelope(EEG_PW_temp',WL,'rms');
                EEG_PW = u';

                % Trim edges of Power Envelope
                EEG_PW(:,size(EEG_PW,2)-WL+1:end) = 0;
                EEG_PW(:,1:WL-1) = 0;
                EEG_mask_sub(:,size(EEG_mask_sub,2)-WL+1:end) = 0;
                EEG_mask_sub(:,1:WL-1) = 0;
                EEG_PW_temp = EEG_PW(:,WL:size(EEG_PW,2)-WL);
                
                % Assign mask and pwr env back into full Part array
                EEG_PW_tot_temp(EEG_mask_sub_sub==1) = EEG_PW;
                EEG_mask_sub_sub(EEG_mask_sub_sub==1) = EEG_mask_sub;
                
                % Set ch specific info back to Part specific data
                EEG_PW_tot(Ch,EEG_mask_ConComp.PixelIdxList{1,i}) = EEG_PW_tot_temp;
                EEG_mask(Ch,EEG_mask_ConComp.PixelIdxList{1,i}) = EEG_mask_sub_sub;
            end
        end
        
%         figure;stackedplot(1:size(EEG_sub,2),EEG_sub')
%         figure;stackedplot(1:size(EEG_sub,2),EEG_PW_tot')

        % Set param long invalid window length as 1 min
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


%%      Cross Correlation:
        % Cross correlate the Power Envs w first-derivative kernal 'Ptrn'
        WL_xcrr = 5*60*srate;
        Ptrn = [(1/WL_xcrr)*ones(1,WL_xcrr),zeros(1,WL),-(1/WL_xcrr)*ones(1,WL_xcrr)];%*ones(1,WL_xcrr)];
        for i=1:size(EEG_PW_tot,1)
            
            EEG_PW_tot_ConComp = bwconncomp(EEG_mask(i,:)==0);
            Strt = [];
            Endd = [];

            % If mask is greater than 0.5 mins save as region to mask
            for j=1:size(EEG_PW_tot_ConComp.PixelIdxList,2)
                if(size(EEG_PW_tot_ConComp.PixelIdxList{1,j},1)>0.1*WL_xcrr)
                    Strt = cat(1,Strt,EEG_PW_tot_ConComp.PixelIdxList{1,j}(1,1));
                    Endd = cat(1,Endd,EEG_PW_tot_ConComp.PixelIdxList{1,j}(end,1));
                end
            end
            
            % Perform cross-correlation
            [c,lags] = xcorr(EEG_PW_tot(i,:),Ptrn);
            EEG_PW_tot(i,:) = c(find(lags==0):end);
            EEG_PW_tot(i,:) = [zeros(1,2*WL_xcrr+WL),EEG_PW_tot(i,1:size(EEG_PW_tot,2)-(2*WL_xcrr+WL))];
            
            % Mask regions found from previous masking loop from PW Env
            for j=1:size(Strt,1)
%                 EEG_PW_tot(i,max(Strt(j)-0.5*size(Ptrn,2),1):min(Strt(j)+0.5*size(Ptrn,2)-1,size(EEG_PW_tot,2))) = 0;
                EEG_PW_tot(i,max(Strt(j)-(2*WL_xcrr+WL),1):min(Endd(j)+(2*WL_xcrr+WL),size(EEG_PW_tot,2))) = 0;
                EEG_mask(i,max(Strt(j)-(2*WL_xcrr+WL),1):min(Endd(j)+(2*WL_xcrr+WL),size(EEG_PW_tot,2))) = 0;
            end 
        end
        
        
        % Analgous to prev code block remove isolated portions from EEG PW  
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
        
        % Perform masking on mask calculated above
        EEG_PW_tot(EEG_mask==0) = 0;
        
        %figure;stackedplot(1:size(EEG_sub,2),EEG_PW_tot')
        % EEG_vis_temp = EEG_PW_tot;

        % Find channel specific Session-wide mean of power envelope
        for i=1:19
            % Ch_inter_dist = (X(i)-X).^2+(Y(i)-Y).^2+(Z(i)-Z).^2;
            % [~,dist_ind] = sort(Ch_inter_dist);
            EEG_PW_temp =  EEG_PW_tot(i,EEG_mask(i,:)==1);
            M_global_tot(i) = M_global_tot(i) + sum(EEG_PW_temp,2);
            M_global_Count(i) = M_global_Count(i) + size(EEG_PW_temp,2);
        end
        
        % Final masking of data 
        EEG_PW_tot(EEG_mask==0) = 0;

%         figure;stackedplot(1:size(EEG_sub,2),EEG_sub')
%         figure;stackedplot(1:size(EEG_sub,2),EEG_PW_tot') 
        
        EEG_PW = EEG_PW_tot;
           
        % EEG_elec = 19;
%         close all         
%         figure;stackedplot(1:size(EEG_sub,2),EEG_sub') 
%         figure;stackedplot(1:size(EEG_sub,2),EEG_PW')
        
        % Save Part Specific information
        save([Session_names(ss).name ,'_',sprintf('CSD_ind_%d',CSD_ind), '_Delta_SDII_vConfusion_240minOv180min_CCWL20min_Xcrr_Pruned_vVII.mat'], 'EEG_PW', 'EEG_sub', 'EEG_mask', 'CSD_type', 'Events_type', 'Events', 'CSD_GT', 'srate', 'CSD_GT_total','CSD_Labels_total','T_start','T_end', '-v7.3');
        close all
        
        % Clear CSD and Event type variables for next Part
        clear CSD_type Events_type
        % End of 'Part'-specific loop
    end
    
    % Find the channel specific global mean across the Session and save
    Global_M = M_global_tot./M_global_Count;
    save([Session_names(ss).name ,'GlobalMean_Delta_SDII_vConfusion_240minOv180min_CCWL20min_Xcrr_Improved_Pruned_vVII.mat'],'Global_M')
    
    % Exit Session-specific loop
    cd ..
end



