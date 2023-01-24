clear all
close all
clc

%% Initialization:

Sub_ID = {'04-1167','04-1177','04-1206','04-1203','04-1201','04-1205','04-1229','04-1214','04-1215','04-1216','04-1219','04-1220'};

subj = 1;
DELTA_T = 2;


Detection_count_subj = {};
Detection_duration_subj = {};
Detection_vector_subj = {};
CSD_Count_subj = {};
Window_quality_subj = {};
Thr_tot_subj = {};
T_start_subj{subj} = {};
T_end_subj{subj} = {};

THR_1 = 0.300000000000000 ;
THR_2 = 0.600000000000000;
THR_3 = 0.689473684210526;
THR_4 = 2.0000;

parfor subj = 1:size(Sub_ID,2)
    cd(['D:\SD-II_POST_ICA\SD-II\',Sub_ID{subj},'\',Sub_ID{subj}])
    Session_names = dir('Patient*');
    Sessions_size = size(Session_names,1);
    
    Detection_count_windows = {};
    Detection_duration_windows = {};
    Detection_vector_windows = {};
    CSD_Count_windows = {};
    T_start_windows = {};
    T_end_windows = {};
    Window_quality_windows = {};
    
    
    for ss = 1:Sessions_size
        cd(Session_names(ss).name)
        Thr_tot = [];
        
        
        TPR_num_Temp = [];
        TPR_denum_Temp = [];
        TNR_num_Temp = [];
        TNR_denum_Temp = [];
        
        EEG_Imp = pop_loadset('filename',[Session_names(ss).name,'_Visualization_Delta_Imp.set'],'filepath',pwd);

        event_delete_ind = [];
        for event_ind=1:size(EEG_Imp.event,2)
            test_ind = strfind(string(EEG_Imp.event(event_ind).type),'Data Quality');
            test_ind = cat(1,test_ind,strfind(string(EEG_Imp.event(event_ind).type),'Detection'));
            test_ind = cat(1,test_ind,strfind(string(EEG_Imp.event(event_ind).type),'boundary'));
            test_ind = cat(1,test_ind,strfind(string(EEG_Imp.event(event_ind).type),'Electrode'));
            
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
                continue;
            end
            
        end
        
        
        EEG_Imp = pop_editeventvals(EEG_Imp,'delete',event_delete_ind);
        
        srate = EEG_Imp.srate;
        
        CSD_GT_total = zeros(1,size(EEG_Imp.data,2));
        Event_total = zeros(1,size(EEG_Imp.data,2));
        CSD_Labels_total = [];
        Event_Labels_total = [];
        for i=1:size(EEG_Imp.event,2)
            if(contains(EEG_Imp.event(i).type,'CSD') || contains(EEG_Imp.event(i).type,'CSD/ISD') || contains(EEG_Imp.event(i).type,'scCSD'))
                CSD_GT_total(round(max(EEG_Imp.event(i).latency,1)))=1;
                CSD_Labels_total = cat(1,CSD_Labels_total,{EEG_Imp.event(i).type});
            end
        end
        
        Block_size = 30*60*60*srate+1;
        Step_size = 1*60*60*srate;
        Min_size = 20*60*60*srate;
        EEG = EEG_Imp.data(1:19,:);
        
        
        if(size(EEG,2)<Min_size)
            cd ..
            continue;
        end
        
        if(size(EEG,2)<Block_size)
            Block_size = size(EEG,2);
        end
        
        Windows_Counter = 0;
        CSD_ind = 1-Step_size;
        while(CSD_ind<=(size(EEG,2)-Block_size+1))
            CSD_ind = CSD_ind + Step_size;
            Thr_tot = [];
            T_start = CSD_ind;
            T_end = min(CSD_ind+Block_size-1,size(EEG,2));
            
            
            
            
            if(T_end-T_start+1<Min_size)
                Step_size = 1*60*60*srate;
                continue;
                
            end
            
            
            CSD_GT = CSD_GT_total;
            
            
            EEG_mask = EEG(:,T_start:T_end);
 
            num_good_elec = (sum(EEG_mask,1)>=0.5*11);
            Window_quality = sum(sum(EEG_mask,1));
            
            if(sum(num_good_elec)<Min_size)
                Step_size = 1*60*60*srate;
                continue;
            end

            Step_size = 1*60*60*srate;
            
            Windows_Counter = Windows_Counter+1;
            CSD_GT = CSD_GT(1,max(T_start-60*60*srate,1):min(T_end+60*60*srate,size(EEG,2)));
            
            CSD_sub_ind = find(CSD_GT);
            
            CSD_Count = size(find(CSD_GT),2);
            for CSD_pruning_ind = 1:size(CSD_sub_ind,2)
                if(mean(sum(EEG_mask(:,max(CSD_sub_ind(CSD_pruning_ind)-5*60*srate+1,1)...
                        :min(CSD_sub_ind(CSD_pruning_ind)+5*60*srate-1,size(EEG_mask,2))),1))<0.1*11)
                    CSD_Count = CSD_Count-1;
                end
            end
            
            
            Detection_count_Temp = [];
            Detection_duration_Temp = [];
            Detection_vector_Temp = [];
            CSD_Count_Temp = [];
            Window_quality_temp = [];
            
            for Thr_1_ind = 1:size(THR_1,2)
                for Thr_2_ind = 1:size(THR_2,2)
                    for Thr_3_ind = 1:size(THR_3,2)
                        for Thr_4_ind = 1:size(THR_4,2)
                            Thr_1 = THR_1(Thr_1_ind)
                            Thr_2 = THR_2(Thr_2_ind)
                            Thr_3 = THR_3(Thr_3_ind)
                            Thr_4 = THR_4(Thr_4_ind)
                            
                            Thr_tot = cat(1,Thr_tot,[Thr_1,Thr_2,Thr_3,Thr_4]);
                            
                            
                            
                            EEG_vis_Event = load([Session_names(ss).name,'_',sprintf('Thr_1_%.2d',Thr_1)...
                                ,'_',sprintf('Thr_2_%.2d',Thr_2),'_',sprintf('Thr_3_%.2d',Thr_3),'_detection_SDIII_240minOv180min_CCWL20_Delta_LowResVelFixedvTEST_VelocitySave.mat']);
                            
                            Strt_tot_TOTAL_temp = EEG_vis_Event.Strt_tot_TOTAL;
                            Strt_tot_TOTAL_temp = Strt_tot_TOTAL_temp{1,Thr_4_ind};
                            Endd_tot_TOTAL_temp = EEG_vis_Event.Endd_tot_TOTAL;
                            Endd_tot_TOTAL_temp = Endd_tot_TOTAL_temp{1,Thr_4_ind};
                            
                            Detection_vector = zeros(1,60*60*60*srate+1);
                            Detection_count = 0;
                            Detection_duration = 0;
                            
                     
                            Strt = Strt_tot_TOTAL_temp*srate;
                            %                                     
                            Endd = Endd_tot_TOTAL_temp*srate;
                            
                            if(subj==1 && ss==2)%% Poor/No ECoG
                                Strt = [];
                                Endd = [];
                            end
                            %
                            if(subj==3 && ss==1 && CSD_ind<50000)%% Poor/No ECoG
                                Strt = [];
                                Endd = [];
                            end
                            
                            for ind_sub=1:size(Strt,1)
                                if(((Strt(ind_sub)>=T_start && Strt(ind_sub)<=T_end) ...
                                        || (Endd(ind_sub)>=T_start && Endd(ind_sub)<=T_end)) && mean(sum(EEG_mask(:,max(Strt(ind_sub)-T_start+1-Thr_4*60*srate+1,1)...
                                        :min(Endd(ind_sub)-T_start+1+Thr_4*60*srate-1,size(EEG_mask,2))),1))>=0*11)
                                    Detection_vector(1,max((Strt(ind_sub)-Thr_4*60*srate-T_start+1),1):...
                                        min((Endd(ind_sub)+Thr_4*60*srate-T_start+1),T_end-T_start+1)) = 1;
                                end
                                if(((Strt(ind_sub)>=T_start-60*60*srate && Strt(ind_sub)<=T_end+60*60*srate) ...
                                        || (Endd(ind_sub)>=T_start-60*60*srate && Endd(ind_sub)<=T_end+60*60*srate))...
                                        && mean(sum(EEG_mask(:,max(Strt(ind_sub)-T_start+1-Thr_4*60*srate+1,1)...
                                        :min(Endd(ind_sub)-T_start+1+Thr_4*60*srate-1,size(EEG_mask,2))),1))>=0.5*11)
                                    Detection_count = Detection_count + 1;
                                    Detection_duration = Detection_duration + ...
                                        (2*Thr_4) + (Endd(ind_sub) - Strt(ind_sub))/(60*srate);
                                end
                            end
                            
                            
                            Detection_count_Temp = cat(2,Detection_count_Temp,Detection_count);
                            Detection_duration_Temp = cat(2,Detection_duration_Temp,Detection_duration);
                            Detection_vector_Temp = cat(1,Detection_vector_Temp,Detection_vector);
                            CSD_Count_Temp = cat(2,CSD_Count_Temp,CSD_Count);
                            Window_quality_temp = cat(2,Window_quality_temp,Window_quality);
                            
                        end
                    end
                end
            end
            Detection_count_windows{ss,Windows_Counter} = Detection_count_Temp;
            Detection_duration_windows{ss,Windows_Counter} = Detection_duration_Temp;
            Detection_vector_windows{ss,Windows_Counter} = Detection_vector_Temp;
            CSD_Count_windows{ss,Windows_Counter} = CSD_Count_Temp;
            Window_quality_windows{ss,Windows_Counter} = Window_quality_temp;
            T_start_windows{ss,Windows_Counter} = T_start;
            T_end_windows{ss,Windows_Counter} = T_end;
        end
        
        cd ..
    end
    Detection_count_subj{subj} = Detection_count_windows;
    Detection_duration_subj{subj} = Detection_duration_windows;
    Detection_vector_subj{subj} = Detection_vector_windows;
    CSD_Count_subj{subj} = CSD_Count_windows;
    Window_quality_subj{subj} = Window_quality_windows;
    T_start_subj{subj} = T_start_windows;
    T_end_subj{subj} = T_end_windows;
    Thr_tot_subj{subj} = Thr_tot;
end


save('Regression_Results_30hrs_vF.mat', 'Detection_count_subj', ...
    'Detection_duration_subj', 'CSD_Count_subj', 'Thr_tot_subj', 'Window_quality_subj', 'Detection_vector_subj', 'T_start_subj', 'T_end_subj', '-v7.3');

%%
clear all
load Regression_Results_30hrs_vF.mat


Thr_ind = 1;

figure;
hold on
CSD_Count = [];
Window_quality = [];
Detection_count = [];
Detection_Duration = [];
Detection_vector = [];
Session_Count = [];
Subj_Count = [];
T_start_tot = [];
T_end_tot = [];




for subj_ind = 1:12
    subj = subj_ind;
    for i=1:size(CSD_Count_subj{1,subj},1)
        for j=1:size(CSD_Count_subj{1,subj},2)
            if(~isempty(CSD_Count_subj{1,subj}{i,j}))
                Session_Count = cat(1,Session_Count,i);
                Subj_Count = cat(1,Subj_Count,subj_ind);
                T_start_tot = cat(1,T_start_tot,T_start_subj{1,subj}{i,j}(Thr_ind));
                T_end_tot = cat(1,T_end_tot,T_end_subj{1,subj}{i,j}(Thr_ind));
                CSD_Count = cat(1,CSD_Count,CSD_Count_subj{1,subj}{i,j}(Thr_ind));
                Window_quality = cat(1,Window_quality,Window_quality_subj{1,subj}{i,j}(Thr_ind));
                Detection_count = cat(1,Detection_count,Detection_count_subj{1,subj}{i,j}(Thr_ind));
                Detection_Duration = cat(1,Detection_Duration,Detection_duration_subj{1,subj}{i,j}(Thr_ind));
                Detection_vector = cat(1,Detection_vector,Detection_vector_subj{1,subj}{i,j}(Thr_ind,:));
            end
        end
    end
end
plot3(CSD_Count,Detection_Duration,Window_quality,'*')


[CSD_Count,CSD_ind] = sort(CSD_Count);
Window_quality = Window_quality(CSD_ind);
Detection_vector = Detection_vector(CSD_ind,:);
Session_Count = Session_Count(CSD_ind,:);
Subj_Count = Subj_Count(CSD_ind,:);
T_start_tot = T_start_tot(CSD_ind,:);
T_end_tot = T_end_tot(CSD_ind,:);


After_dist_tot = {};
Befor_dist_tot = {};
Length_det_tot = {};

figure;
hold on

for i=1:size(CSD_Count,1)
        
    
    Detection_ConComp = bwconncomp(Detection_vector(i,:));
    
    After_dist = [];
    Befor_dist = [];
    Length_det = [];
    for j=2:size(Detection_ConComp.PixelIdxList,2)-1
        
        After_dist = cat(1,After_dist,Detection_ConComp.PixelIdxList{1,j+1}(1)-Detection_ConComp.PixelIdxList{1,j}(end));
        Befor_dist = cat(1,Befor_dist,Detection_ConComp.PixelIdxList{1,j}(1)-Detection_ConComp.PixelIdxList{1,j-1}(end));
        Length_det = cat(1,Length_det,size(Detection_ConComp.PixelIdxList{1,j},1));
        
        if(Detection_ConComp.PixelIdxList{1,j+1}(1)-Detection_ConComp.PixelIdxList{1,j}(end)>240*60*64 && ...
                Detection_ConComp.PixelIdxList{1,j}(1)-Detection_ConComp.PixelIdxList{1,j-1}(end)>240*60*64 ...
                && size(Detection_ConComp.PixelIdxList{1,j},1)<20*60*64)
            Detection_vector(i,Detection_ConComp.PixelIdxList{1,j}(1):...
                Detection_ConComp.PixelIdxList{1,j}(end)) = 0;
        end
    end
    
    Detection_ConComp = bwconncomp(Detection_vector(i,:));
    
    if(size(Detection_ConComp.PixelIdxList,2)>=2)
        
        After_dist = cat(1,After_dist,Detection_ConComp.PixelIdxList{1,2}(1)-Detection_ConComp.PixelIdxList{1,1}(end));
        Befor_dist = cat(1,Befor_dist,Detection_ConComp.PixelIdxList{1,2}(1)-Detection_ConComp.PixelIdxList{1,1}(end));
        Length_det = cat(1,Length_det,size(Detection_ConComp.PixelIdxList{1,1},1));
        
        
        if(Detection_ConComp.PixelIdxList{1,2}(1)-Detection_ConComp.PixelIdxList{1,1}(end)>240*60*64 ...
                && size(Detection_ConComp.PixelIdxList{1,1},1)<20*60*64)
            Detection_vector(i,Detection_ConComp.PixelIdxList{1,1}(1):...
                Detection_ConComp.PixelIdxList{1,1}(end)) = 0;
        end
        
        After_dist = cat(1,After_dist,Detection_ConComp.PixelIdxList{1,end}(1)-Detection_ConComp.PixelIdxList{1,end-1}(end));
        Befor_dist = cat(1,Befor_dist,Detection_ConComp.PixelIdxList{1,end}(1)-Detection_ConComp.PixelIdxList{1,end-1}(end));
        Length_det = cat(1,Length_det,size(Detection_ConComp.PixelIdxList{1,end},1));
        
        if(Detection_ConComp.PixelIdxList{1,end}(1)-Detection_ConComp.PixelIdxList{1,end-1}(end)>240*60*64 ...
                && size(Detection_ConComp.PixelIdxList{1,end},1)<20*60*64)
            Detection_vector(i,Detection_ConComp.PixelIdxList{1,end}(1):...
                Detection_ConComp.PixelIdxList{1,end}(end)) = 0;
        end
        plot3(min(Befor_dist,After_dist)/(60*64),Length_det/(60*64),repmat(CSD_Count(i),size(Length_det,1),1),'*')
        
    end
    
    
    
    After_dist_tot{i} = After_dist;
    Befor_dist_tot{i} = Befor_dist;
    Length_det_tot{i} = Length_det;
end

%%
Detection_vector_temp = Detection_vector;

for stitchingparam = 240
    
    Detection_count = [];
    Detection_temp = [];
    CSD_Count_temp = [];
    Subj_Count_temp = [];
    Session_Count_temp = [];
    for i=1:size(CSD_Count,1)
        
        Detection_ConComp = bwconncomp(Detection_vector_temp(i,:));
        
        
        for j=1:size(Detection_ConComp.PixelIdxList,2)-1
            if(Detection_ConComp.PixelIdxList{1,j+1}(1)-Detection_ConComp.PixelIdxList{1,j}(end)<stitchingparam*60*64)
                Detection_vector_temp(i,Detection_ConComp.PixelIdxList{1,j}(end):...
                    Detection_ConComp.PixelIdxList{1,j+1}(1)) = 1;
            end
        end
        
        Detection_ConComp = bwconncomp(Detection_vector_temp(i,:));
        for j=1:size(Detection_ConComp.PixelIdxList,2)
            if(size(Detection_ConComp.PixelIdxList{1,j},1)<1*60*64)
                Detection_vector_temp(i,Detection_ConComp.PixelIdxList{1,j}(1):...
                    Detection_ConComp.PixelIdxList{1,j}(end)) = 0;
            end
        end
        
        Detection_ConComp_F = bwconncomp(Detection_vector_temp(i,:));
        for j=1:size(Detection_ConComp_F.PixelIdxList,2)
            Detection_temp = cat(1,Detection_temp,size(Detection_ConComp_F.PixelIdxList{1,j},1));
            CSD_Count_temp = cat(1,CSD_Count_temp,CSD_Count(i));
            Subj_Count_temp = cat(1,Subj_Count_temp,Subj_Count(i));
            Session_Count_temp = cat(1,Session_Count_temp,Session_Count(i));
        end
        if(size(Detection_ConComp_F.PixelIdxList,2)==0)
            
            Detection_temp = cat(1,Detection_temp,0);
            CSD_Count_temp = cat(1,CSD_Count_temp,CSD_Count(i));
            Subj_Count_temp = cat(1,Subj_Count_temp,Subj_Count(i));
            Session_Count_temp = cat(1,Session_Count_temp,Session_Count(i));
            
        end
        
        
        Detection_count = cat(1,Detection_count,size(Detection_ConComp_F.PixelIdxList,2));
    end
    
    figure;plot(CSD_Count_temp,Detection_temp/(60*64),'*')
    figure;plot3(CSD_Count_temp,Detection_temp/(60*64),Subj_Count_temp,'*')
    
    
    Detection_Duration = sum(Detection_vector_temp,2)/(60*64);
    
    WL_Q = Window_quality;
    WL = (T_end_tot-T_start_tot)/(60*60*64);
    figure;plot3(CSD_Count,Detection_Duration,WL_Q,'*')
    figure;plot3(CSD_Count,Detection_Duration,WL,'*')
    figure;plot3(CSD_Count(WL>21 & WL_Q>6),Detection_Duration(WL>21 & WL_Q>6),WL_Q(WL>21 & WL_Q>6),'*')
    
    CSD_Count_normalize = CSD_Count./WL_Q;
    figure;plot3(CSD_Count_normalize,Detection_Duration,WL,'*')
    figure;plot3(CSD_Count(WL_Q>8.3),Detection_Duration(WL_Q>8.3),WL(WL_Q>8.3),'*')
    
end

%%
Detection_vector = Detection_vector_temp;

fun = @(x)sum((x(1)*CSD_Count.^2+x(2)*CSD_Count+x(3)-Detection_Duration).^2);

x0 = [-1,20,0];
A = [-150,-1,0];
b = 0;
x = fmincon(fun,x0,A,b)
figure;plot(CSD_Count,Detection_Duration,'*')
hold on
Detection_Duration_est = x(1)*CSD_Count.^2+x(2)*CSD_Count+x(3);
plot(CSD_Count,Detection_Duration_est,'-r')

r = Detection_Duration-Detection_Duration_est;
SSE = sum(r.^2);
SST = sum((Detection_Duration-mean(Detection_Duration)).^2);
R2 = 1 - SSE/SST
RMSE = sqrt(mean((Detection_Duration_est-Detection_Duration).^2))

%%

fun = @(x)sum((x(1)*sqrt(CSD_Count)+x(2)-Detection_Duration).^2);

x0 = [1,200];
x = fmincon(fun,x0)
figure;plot(CSD_Count,Detection_Duration,'*')
hold on
Detection_Duration_est = x(1)*sqrt(CSD_Count)+x(2);
plot(CSD_Count,Detection_Duration_est,'-r')

r = Detection_Duration-Detection_Duration_est;
SSE = sum(r.^2);
SST = sum((Detection_Duration-mean(Detection_Duration)).^2);
R2 = 1 - SSE/SST
RMSE = sqrt(mean((Detection_Duration_est-Detection_Duration).^2))

CSD_Count_est = ((Detection_Duration - x(2))/x(1)).^2;
RMSE_CSD_Count = sqrt(mean((CSD_Count_est-CSD_Count).^2))




%%
LR = {'R','L','R','R','R','R','L','R','R','R','R','R'};
%possible good low freq SDs: 51
for i=15:size(CSD_Count,1)
    
    cd(['D:\SD-II_POST_ICA\SD-II\',Sub_ID{Subj_Count(i)},'\',Sub_ID{Subj_Count(i)}])
    Session_names = dir('Patient*');
    cd(Session_names(Session_Count(i)).name)
    
    
    EEG_vis_PW = pop_loadset('filename',[Session_names(Session_Count(i)).name,'_Visualization_Delta_PW_xcrr.set'],'filepath',pwd);
    EEG_vis = pop_loadset('filename',[Session_names(Session_Count(i)).name,'_Visualization_Delta.set'],'filepath',pwd);
    EEG_Imp = pop_loadset('filename',[Session_names(Session_Count(i)).name,'_Visualization_Delta_Imp.set'],'filepath',pwd);
    EEG = EEG_Imp.data(1:19,:);
    event_delete_ind = [];
    for event_ind=1:size(EEG_vis.event,2)
        test_ind = strfind(string(EEG_vis.event(event_ind).type),'Data Quality Normal');
        %         test_ind = cat(1,test_ind,strfind(string(EEG_vis.event(event_ind).type),'Data Quality'));
        %         test_ind = cat(1,test_ind,strfind(string(EEG_vis.event(event_ind).type),'boundary'));
        %         test_ind = cat(1,test_ind,strfind(string(EEG_vis.event(event_ind).type),'Electrode'));
        
        target_detected = 0;
        for ss_ind = 1:size(test_ind,1)
            if(size(test_ind,1)>1)
                if(~isempty(test_ind{ss_ind,1}) && mean(sum(EEG(:,max(round(max(EEG_Imp.event(event_ind).latency,1))-5*60*srate+1,1)...
                        :min(round(max(EEG_Imp.event(event_ind).latency,1))+5*60*srate-1,size(EEG,2))),1))>=0.1*11)
                    
                    target_detected = 1;
                    
                end
            else
                if(~isempty(test_ind) && mean(sum(EEG(:,max(round(max(EEG_Imp.event(event_ind).latency,1))-5*60*srate+1,1)...
                        :min(round(max(EEG_Imp.event(event_ind).latency,1))+5*60*srate-1,size(EEG,2))),1))>=0.1*11)
                    target_detected = 1;
                end
            end
        end
        
        if(target_detected==0)
            event_delete_ind = cat(1,event_delete_ind,event_ind);
            %                         EEG_vis.event(event_ind) = [];
            continue;
        end
        
    end
    
    %     EEG_vis = pop_editeventvals(EEG_vis,'delete',event_delete_ind);
    %     EEG_vis_PW = pop_editeventvals(EEG_vis_PW,'delete',event_delete_ind);
    srate = EEG_vis_PW.srate;
    
    
    %     EEG_vis.data = EEG_vis.data(:,T_start_tot(i):T_end_tot(i));
    
    STRT = EEG_vis.xmin + (T_start_tot(i)/srate)
    %      EEG_vis_temp.xmin + (T_end_tot(i)/srate);
    
    
    Detection_ConComp_F = bwconncomp(Detection_vector(i,:));
    for j=1:size(Detection_ConComp_F.PixelIdxList,2)
        
        Strt = EEG_vis_PW.xmin + ((T_start_tot(i)+Detection_ConComp_F.PixelIdxList{1,j}(1))/srate);
        Endd = EEG_vis_PW.xmin + ((T_start_tot(i)+Detection_ConComp_F.PixelIdxList{1,j}(end))/srate);
        
        EEG_vis_PW = pop_editeventvals(EEG_vis_PW,'insert',{1,[],[],[]},'changefield',{1,'type','Start_Detection'},'changefield',{1,'latency',Strt},'changefield',{1,'duration',1});
        EEG_vis_PW = pop_editeventvals(EEG_vis_PW,'insert',{1,[],[],[]},'changefield',{1,'type','Stop_Detection'},'changefield',{1,'latency',Endd},'changefield',{1,'duration',1});
        
    end
    
    EEG_vis.event = EEG_vis_PW.event;
    EEG_vis_PW.event = EEG_vis_PW.event;
    
    
    crani_ind = 1:19;
    if(LR{Subj_Count(i)}=='L')
        LR_sgn = -1;
        crani_ind([1,3,4,6,8,10,11,13,15,17,18]) = [];
        Ch_id_removed = crani_ind;
        Ch_referential = [8,6,4,10,18,1,3,15,13,17,11,Ch_id_removed,20:25];
    else
        Ch_id_removed = crani_ind([1,4,6,8,11,13,15,18]);
        Ch_referential = [9,7,5,10,19,2,3,16,14,17,12,Ch_id_removed,20:25];
    end
    
    %     Ch_referential = [9,7,5,10,19,2,3,16,14,17,12,8,6,4,18,1,15,13,11,20:25];
    
    EEG_vis_PW.data = EEG_vis_PW.data(Ch_referential,:);
    EEG_vis_PW.chanlocs = EEG_vis_PW.chanlocs(Ch_referential);
    EEG_vis_PW.nbchan = size(EEG_vis_PW.data,1);
    
    EEG_vis.data = EEG_vis.data(Ch_referential,:);
    EEG_vis.chanlocs = EEG_vis.chanlocs(Ch_referential);
    EEG_vis.nbchan = size(EEG_vis.data,1);
    
    %     for ind_removal=1:size(Ch_id_removed,2)
    %         Ch_id_removed(ind_removal) = find(Ch_referential==Ch_id_removed(ind_removal));
    %     end
    
    EEG_vis_PW = pop_select( EEG_vis_PW,'nochannel',12:19);
    EEG_vis = pop_select( EEG_vis,'nochannel',12:19);
    
    
    pop_eegplot( EEG_vis, 1, 1, 1);
    pop_eegplot( EEG_vis_PW, 1, 1, 1);
    
    
    
end





% figure;stackedplot((1:size(Detection_vector,2))/(60*64),Detection_vector(1:42,:)')
%
%
% figure;stackedplot((1:size(Detection_vector,2))/(60*64),Detection_vector(48:end,:)')












