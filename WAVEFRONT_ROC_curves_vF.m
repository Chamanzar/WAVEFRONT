clear all
clc

%%
% ROC curve analysis and SD velocity estimation and corresponding plots
%
% See: README.txt and [1] for more info.

% [1] Alireza Chamanzar, Jonathan Elmer, Lori Shutter, Jed Hartings, Pulkit Grover,
%  "Noninvasive and reliable automated detection of spreading depolarization in severe traumatic brain injury using scalp EEG",
%   Submitted to Nature Comms Med, 2022.

% Author: Alireza 	Date: 2022/01/06 12:00:00 	Revision: 0.1
% Copyright: This project is licensed - see the LICENSE.md file for details


%% Initialization:

Sub_ID = {'04-1229','04-1206','04-1201','04-1167','04-1205','04-1203','04-1177','04-1214','04-1215','04-1216','04-1220','04-1219'};
subj = 1;
DELTA_T = 60;

TPR_num_TOT = [];
TPR_denum_TOT = [];
TNR_num_TOT = [];
TNR_denum_TOT = [];



THR_1 = 0.05:0.05:0.6;%mm
THR_3 = linspace(0.1,0.9,20);
THR_1 = THR_1([1,2,3,4,5,6,7,8,9]);%mm
THR_2 = 0.6:0.1:0.7;
THR_3 = THR_3(10:15);
THR_4 = 2:3;

%%

% Average_Vel_tot_TOT = {};

for subj = 1:size(Sub_ID,2)
    cd(['D:\SD-II_POST_ICA\SD-II\',Sub_ID{subj},'\',Sub_ID{subj}])
    Session_names = dir('Patient*');
    Sessions_size = size(Session_names,1);
    TPR_num_Temp_temp = zeros(1,size(THR_1,2)*size(THR_2,2)*size(THR_3,2)*size(THR_4,2));
    TPR_denum_Temp_temp = 0;
    TNR_num_Temp_temp = 0;
    TNR_denum_Temp_temp = 0;
    for ss = 1:Sessions_size
%         cd(Session_names(ss).name)
        
        Thr_tot = [];
        
        TPR_num_Temp = [];
        TPR_denum_Temp = [];
        TNR_num_Temp = [];
        TNR_denum_Temp = [];
        
        for Thr_1_ind = 1:size(THR_1,2)
            for Thr_2_ind = 1:size(THR_2,2)
                for Thr_3_ind = 1:size(THR_3,2)
                        Thr_1 = THR_1(Thr_1_ind)
                        Thr_2 = THR_2(Thr_2_ind)
                        Thr_3 = THR_3(Thr_3_ind)
                        
                                                        load([Session_names(ss).name,'\',Session_names(ss).name ,'_',sprintf('Thr_1_%.2d',Thr_1)...
                                ,'_',sprintf('Thr_2_%.2d',Thr_2),'_',sprintf('Thr_3_%.2d',Thr_3), '_detection_SDIII_240minOv180min_CCWL20_Delta_LowResVelFixedvTEST_VelocitySave.mat'])
                            
                        TPR_num = double(TPR_num_TOTAL_temp(1:size(THR_4,2)))';
                        TPR_denum = double(TPR_denum_TOTAL_temp(1:size(THR_4,2)))';
                        TNR_num = double(TNR_num_TOTAL_temp(1:size(THR_4,2)))';
                        TNR_denum = double(TNR_denum_TOTAL_temp(1:size(THR_4,2)))';
                        for Thr_4_ind = 1:size(THR_4,2)
                            Thr_4 = THR_4(Thr_4_ind)
                            Thr_tot = cat(1,Thr_tot,[Thr_1,Thr_2,Thr_3,Thr_4]);                        
                        end
                        
                        TPR_num_Temp = cat(2,TPR_num_Temp,TPR_num);
                        TPR_denum_Temp = cat(2,TPR_denum_Temp,TPR_denum);
                        TNR_num_Temp = cat(2,TNR_num_Temp,TNR_num);
                        TNR_denum_Temp = cat(2,TNR_denum_Temp,TNR_denum);
      
                end
            end
        end
        
        TPR_num_Temp_temp = TPR_num_Temp_temp + TPR_num_Temp;
        TPR_denum_Temp_temp = TPR_denum_Temp_temp + TPR_denum_Temp;
        TNR_num_Temp_temp = TNR_num_Temp_temp + TNR_num_Temp;
        TNR_denum_Temp_temp = TNR_denum_Temp_temp + TNR_denum_Temp;
        
    end
    TPR_num_TOT = cat(1,TPR_num_TOT,TPR_num_Temp_temp);
    TPR_denum_TOT = cat(1,TPR_denum_TOT,TPR_denum_Temp_temp);
    TNR_num_TOT = cat(1,TNR_num_TOT,TNR_num_Temp_temp);
    TNR_denum_TOT = cat(1,TNR_denum_TOT,TNR_denum_Temp_temp);
end

%%

Y_tot = [];
X_tot = [];


X_Tr_avg = [];
Y_Tr_avg = [];

X_Val_avg = [];
Y_Val_avg = [];


X_Val_num = [];
X_Val_denum = [];
Y_Val_num = [];
Y_Val_denum = [];


X_Tr_num = [];
X_Tr_denum = [];
Y_Tr_num = [];
Y_Tr_denum = [];


X_Val_TOT = [];
Y_Val_TOT = [];

X_Tr_TOT = [];
Y_Tr_TOT = [];

Thr_train_TOT = [];

ind_tempThr_TOT = 0;
for TPR_Thr = 0.50:0.00005:0.953
    ind_tempThr_TOT = ind_tempThr_TOT+1;
    X_Val_tot = [];
    Y_Val_tot = [];
    X_Tr_tot = [];
    Y_Tr_tot = [];
    Thr_train_tot = [];
    TNR_num_Val = [];
    TNR_denum_Val = [];
    TPR_num_Val = [];
    TPR_denum_Val = [];
    TNR_num_Tr = [];
    TNR_denum_Tr = [];
    TPR_num_Tr = [];
    TPR_denum_Tr = [];
    Sub_ID_tot = [];

    ind_tempThr = 0;
    Val_non_CSD_num_TOT = 0;
    Val_CSD_num_TOT = 0;
    Tr_non_CSD_num_TOT = 0;
    Tr_CSD_num_TOT = 0;
    for i=1:12
        for j=i+1:12
            Tr_ind = 1:12;
            Val_ind = [i,j];
            Tr_ind(Val_ind) = [];

            X_Tr = 1-(sum(TNR_num_TOT(Tr_ind,:),1)./sum(TNR_denum_TOT(Tr_ind,:),1));
            [Y_Tr,ind] = sort(sum(TPR_num_TOT(Tr_ind,:),1)./sum(TPR_denum_TOT(Tr_ind,:),1));
            X_Tr = X_Tr(ind);
            PPV_Tr = sum(TPR_num_TOT(Tr_ind,:),1)./(sum(TPR_num_TOT(Tr_ind,:),1)+...
                (sum(TNR_denum_TOT(Tr_ind,:),1)-sum(TNR_num_TOT(Tr_ind,:),1)));
            PPV_Tr = PPV_Tr(ind);

            Thr_train = Thr_tot(ind,:);

            [~,Trained_Param_ind_L] =  max(X_Tr(Y_Tr<=TPR_Thr));
            [~,Trained_Param_ind_U] =  min(X_Tr(Y_Tr>=TPR_Thr));
            ind_temp = find(Y_Tr>=TPR_Thr);
            Trained_Param_ind_U = ind_temp(Trained_Param_ind_U);
            ind_temp = find(Y_Tr<=TPR_Thr);
            Trained_Param_ind_L = ind_temp(Trained_Param_ind_L);
            if(Trained_Param_ind_L~=Trained_Param_ind_U)
                X_Tr_U = X_Tr(Trained_Param_ind_U);
                X_Tr_L = X_Tr(Trained_Param_ind_L);
                X_Tr_mu = (X_Tr_U+X_Tr_L)/2;

            end

            [~,Trained_Param_ind] =  min(X_Tr(Y_Tr>TPR_Thr));
            ind_temp = find(Y_Tr>TPR_Thr);
            Trained_Param_ind = ind_temp(Trained_Param_ind);
            Y_Tr = Y_Tr(Trained_Param_ind);
            X_Tr = X_Tr(Trained_Param_ind);

            Thr_train = Thr_train(Trained_Param_ind,:);
            
            Val_Param_ind = find((Thr_tot(:,1)==Thr_train(1)) & (Thr_tot(:,2)==Thr_train(2)) ...
                & (Thr_tot(:,3)==Thr_train(3)) & (Thr_tot(:,4)==Thr_train(4)));
            X_Val = 1-(sum(TNR_num_TOT(Val_ind,Val_Param_ind),1)./sum(TNR_denum_TOT(Val_ind,Val_Param_ind),1));
            Y_Val = sum(TPR_num_TOT(Val_ind,Val_Param_ind),1)./sum(TPR_denum_TOT(Val_ind,Val_Param_ind),1);

            TNR_num_Val = cat(1,TNR_num_Val,sum(TNR_num_TOT(Val_ind,Val_Param_ind),1));
            TNR_denum_Val = cat(1,TNR_denum_Val,sum(TNR_denum_TOT(Val_ind,Val_Param_ind),1));
            TPR_num_Val = cat(1,TPR_num_Val,sum(TPR_num_TOT(Val_ind,Val_Param_ind),1));
            TPR_denum_Val = cat(1,TPR_denum_Val,sum(TPR_denum_TOT(Val_ind,Val_Param_ind),1));

            TNR_num_Tr = cat(1,TNR_num_Tr,sum(TNR_num_TOT(Tr_ind,Val_Param_ind),1));
            TNR_denum_Tr = cat(1,TNR_denum_Tr,sum(TNR_denum_TOT(Tr_ind,Val_Param_ind),1));
            TPR_num_Tr = cat(1,TPR_num_Tr,sum(TPR_num_TOT(Tr_ind,Val_Param_ind),1));
            TPR_denum_Tr = cat(1,TPR_denum_Tr,sum(TPR_denum_TOT(Tr_ind,Val_Param_ind),1));


            X_Val_tot = cat(1,X_Val_tot,X_Val);
            Y_Val_tot = cat(1,Y_Val_tot,Y_Val);

            X_Tr_tot = cat(1,X_Tr_tot,X_Tr);
            Y_Tr_tot = cat(1,Y_Tr_tot,Y_Tr);
            
            Sub_ID(1,[i,j])
            Sub_ID_tot = cat(1,Sub_ID_tot,[i,j]);
            Tr_CSD_num = sum(TPR_denum_TOT(Tr_ind,Trained_Param_ind),1)
            Tr_non_CSD_num = sum(TNR_denum_TOT(Tr_ind,Trained_Param_ind),1)
            Y_Tr
            X_Tr
            Val_CSD_num = sum(TPR_denum_TOT(Val_ind,Val_Param_ind),1)
            Val_non_CSD_num = sum(TNR_denum_TOT(Val_ind,Val_Param_ind),1)
            Y_Val
            X_Val
            
            Val_non_CSD_num_TOT = Val_non_CSD_num_TOT + Val_non_CSD_num;
            Val_CSD_num_TOT = Val_CSD_num_TOT + Val_CSD_num;
            
            Tr_non_CSD_num_TOT = Tr_non_CSD_num_TOT + Tr_non_CSD_num;
            Tr_CSD_num_TOT = Tr_CSD_num_TOT + Tr_CSD_num;



            Thr_train_tot = cat(1,Thr_train_tot,Thr_train);

        end
    end

    X_Val_TOT = cat(2,X_Val_TOT,X_Val_tot);
    Y_Val_TOT = cat(2,Y_Val_TOT,Y_Val_tot);

    X_Tr_TOT = cat(2,X_Tr_TOT,X_Tr_tot);
    Y_Tr_TOT = cat(2,Y_Tr_TOT,Y_Tr_tot);

    X_Tr_avg = cat(1,X_Tr_avg,1-(sum(TNR_num_Tr,1)./sum(TNR_denum_Tr,1)));
    Y_Tr_avg = cat(1,Y_Tr_avg,sum(TPR_num_Tr,1)./sum(TPR_denum_Tr,1));

    X_Val_avg = cat(1,X_Val_avg,1-(sum(TNR_num_Val,1)./sum(TNR_denum_Val,1)));
    Y_Val_avg = cat(1,Y_Val_avg,sum(TPR_num_Val,1)./sum(TPR_denum_Val,1));
    
    X_Val_num = cat(2,X_Val_num, TNR_num_Val);
    X_Val_denum = cat(2,X_Val_denum, TNR_denum_Val);
    Y_Val_num = cat(2,Y_Val_num, TPR_num_Val);
    Y_Val_denum = cat(2,Y_Val_denum, TPR_denum_Val);
    
    
    X_Tr_num = cat(2,X_Tr_num, TNR_num_Tr);
    X_Tr_denum = cat(2,X_Tr_denum, TNR_denum_Tr);
    Y_Tr_num = cat(2,Y_Tr_num, TPR_num_Tr);
    Y_Tr_denum = cat(2,Y_Tr_denum, TPR_denum_Tr);
    
    
    Thr_train_TOT = cat(3,Thr_train_TOT,Thr_train_tot);

end


%%

Y_Val_TOT(isnan(Y_Val_TOT)) = 1;


%% Averaging ROC across thresholds: 
Y_Val_avg = zeros(size(Thr_tot,1),1);
X_Val_avg = zeros(size(Thr_tot,1),1);
Y_Val_size_TOT = zeros(size(Thr_tot,1),1);
X_Val_size_TOT = zeros(size(Thr_tot,1),1);
X_Val_TOT = {};
Y_Val_TOT = {};
X_Val_size = {};
Y_Val_size = {};

for j = 1:size(Thr_tot,1)
    X_Val_avg_num = 0;
    X_Val_avg_denum = 0;
    Y_Val_avg_num = 0;
    Y_Val_avg_denum = 0;
    Y_Val_avg_size = 0;
    X_Val_avg_size = 0;
        
    for i=1:size(Thr_train_TOT,1)
        Val_Param_ind = find((Thr_train_TOT(i,1,:)==Thr_tot(j,1)) & (Thr_train_TOT(i,2,:)==Thr_tot(j,2)) ...
            & (Thr_train_TOT(i,3,:)==Thr_tot(j,3)) & (Thr_train_TOT(i,4,:)==Thr_tot(j,4)));
       
        
        X_Val_TOT{i,j} = [];
        Y_Val_TOT{i,j} =[];
        X_Val_size{i,j} = [];
        Y_Val_size{i,j} = [];
        
        if(isempty(Val_Param_ind))
            continue
        end
        
        X_Val_TOT{i,j} = 1-(sum(X_Val_num(i,Val_Param_ind),2)./sum(X_Val_denum(i,Val_Param_ind),2));
        Y_Val_TOT{i,j} = sum(Y_Val_num(i,Val_Param_ind),2)./sum(Y_Val_denum(i,Val_Param_ind),2);
        X_Val_size{i,j} = sum(X_Val_denum(i,Val_Param_ind),2);
        Y_Val_size{i,j} = sum(Y_Val_denum(i,Val_Param_ind),2);
        
        X_Val_avg_num = X_Val_avg_num + sum(X_Val_num(i,Val_Param_ind),2);
        X_Val_avg_denum = X_Val_avg_denum + sum(X_Val_denum(i,Val_Param_ind),2);
        Y_Val_avg_num = Y_Val_avg_num + sum(Y_Val_num(i,Val_Param_ind),2);
        Y_Val_avg_denum = Y_Val_avg_denum + sum(Y_Val_denum(i,Val_Param_ind),2);
        
        X_Val_avg_size = X_Val_avg_size + (X_Val_denum(i,1));
        Y_Val_avg_size = Y_Val_avg_size + (Y_Val_denum(i,1));
        
    end
    
    if(X_Val_avg_denum>0)
        X_Val_avg(j) = 1-(X_Val_avg_num/X_Val_avg_denum);
    else
        X_Val_avg(j) = 1;
    end
    if(Y_Val_avg_denum>0)
        Y_Val_avg(j) = (Y_Val_avg_num/Y_Val_avg_denum);
    else
        Y_Val_avg(j) = 1;
    end
    Y_Val_size_TOT(j) = Y_Val_avg_size;
    X_Val_size_TOT(j) = X_Val_avg_size;
end



X_Val_TOT_bootstrped = [];
Y_Val_TOT_bootstrped = [];
CI_X_Val = [];
CI_Y_Val = [];
ts_TOT = [];
for j=1:size(X_Val_size,2)
    rng('shuffle') %For reproducibility
    temp = [];
    size_temp = [];
    for i=1:size(X_Val_size,1)
        if(isempty(X_Val_TOT{i,j}))
            continue;
        end
        temp = cat(2,temp,X_Val_TOT{i,j});
        size_temp = cat(2,size_temp,X_Val_size{i,j});
    end
    
    if(isempty(size_temp))
         ts_TOT = cat(2,ts_TOT,0);
         X_Val_TOT_bootstrped = cat(2,X_Val_TOT_bootstrped,zeros(100,1));
         Y_Val_TOT_bootstrped = cat(2,Y_Val_TOT_bootstrped,zeros(100,1));
         CI_X_Val = cat(1,CI_X_Val,0);
         CI_Y_Val = cat(1,CI_Y_Val,0);
        continue;
    end
    ts = tinv([0.025  0.975],max(1,size(size_temp,2)-1));      % T-Score
    CI_X_Val = cat(1,CI_X_Val,ts(1)*std(temp,0,2)/sqrt(size(size_temp,2)));
    if(size(size_temp,2)>1)
    X_Val_TOT_bootstrped = cat(2,X_Val_TOT_bootstrped,bootstrp(100,@mean,temp','Weights',size_temp'));
    else
        X_Val_TOT_bootstrped = cat(2,X_Val_TOT_bootstrped,temp*ones(100,1));
    end
    
    temp = [];
    size_temp = [];
    for i=1:size(Y_Val_size,1)
        if(isempty(Y_Val_TOT{i,j}))
            continue;
        end
        temp = cat(2,temp,Y_Val_TOT{i,j});
        temp(isnan(temp)) = 1;
        size_temp = cat(2,size_temp,Y_Val_size{i,j});
    end
    ts = tinv([0.025  0.975],max(1,size(size_temp,2)-1));      % T-Score
    ts_TOT = cat(2,ts_TOT,ts(1));
    CI_Y_Val = cat(1,CI_Y_Val,ts(1)*std(temp,0,2)/sqrt(size(size_temp,2)));
    if(size(size_temp,2)>1)
    Y_Val_TOT_bootstrped = cat(2,Y_Val_TOT_bootstrped,bootstrp(100,@mean,temp','Weights',size_temp'));
    else
        Y_Val_TOT_bootstrped = cat(2,Y_Val_TOT_bootstrped,temp*ones(100,1));
    end
end

SE_X_Val = std(X_Val_TOT_bootstrped,0,1)';%(std(X_Tr_TOT,0,1)/sqrt(size(X_Tr_TOT,1)))';
SE_Y_Val = std(Y_Val_TOT_bootstrped,0,1)';%(std(Y_Tr_TOT,0,1)/sqrt(size(Y_Tr_TOT,1)))';
CI_Y_Val = ts_TOT'.*SE_Y_Val;                      % Confidence Intervals
CI_X_Val = ts_TOT'.*SE_X_Val; 


Y_Tr_avg = zeros(size(Thr_tot,1),1);
X_Tr_avg = zeros(size(Thr_tot,1),1);
Y_Tr_size_TOT = zeros(size(Thr_tot,1),1);
X_Tr_size_TOT = zeros(size(Thr_tot,1),1);
X_Tr_TOT = {};
Y_Tr_TOT = {};
X_Tr_size = {};
Y_Tr_size = {};

for j = 1:size(Thr_tot,1)
    X_Tr_avg_num = 0;
    X_Tr_avg_denum = 0;
    Y_Tr_avg_num = 0;
    Y_Tr_avg_denum = 0;
        
    for i=1:size(Thr_train_TOT,1)
        Val_Param_ind = find((Thr_train_TOT(i,1,:)==Thr_tot(j,1)) & (Thr_train_TOT(i,2,:)==Thr_tot(j,2)) ...
            & (Thr_train_TOT(i,3,:)==Thr_tot(j,3)) & (Thr_train_TOT(i,4,:)==Thr_tot(j,4)));
       
        
        X_Tr_TOT{i,j} = [];
        Y_Tr_TOT{i,j} =[];
        X_Tr_size{i,j} = [];
        Y_Tr_size{i,j} = [];
        
        if(isempty(Val_Param_ind))
            continue
        end
        
        X_Tr_TOT{i,j} = 1-(sum(X_Tr_num(i,Val_Param_ind),2)./sum(X_Tr_denum(i,Val_Param_ind),2));
        Y_Tr_TOT{i,j} = sum(Y_Tr_num(i,Val_Param_ind),2)./sum(Y_Tr_denum(i,Val_Param_ind),2);
        X_Tr_size{i,j} = sum(X_Tr_denum(i,Val_Param_ind),2);
        Y_Tr_size{i,j} = sum(Y_Tr_denum(i,Val_Param_ind),2);
        
        X_Tr_avg_num = X_Tr_avg_num + sum(X_Tr_num(i,Val_Param_ind),2);
        X_Tr_avg_denum = X_Tr_avg_denum + sum(X_Tr_denum(i,Val_Param_ind),2);
        Y_Tr_avg_num = Y_Tr_avg_num + sum(Y_Tr_num(i,Val_Param_ind),2);
        Y_Tr_avg_denum = Y_Tr_avg_denum + sum(Y_Tr_denum(i,Val_Param_ind),2);
        
    end
    
    if(X_Tr_avg_denum>0)
        X_Tr_avg(j) = 1-(X_Tr_avg_num/X_Tr_avg_denum);
    else
        X_Tr_avg(j) = 1;
    end
    if(Y_Tr_avg_denum>0)
        Y_Tr_avg(j) = (Y_Tr_avg_num/Y_Tr_avg_denum);
    else
        Y_Tr_avg(j) = 1;
    end
    Y_Tr_size_TOT(j) = Y_Tr_avg_denum;
    X_Tr_size_TOT(j) = X_Tr_avg_denum;
end



X_Tr_TOT_bootstrped = [];
Y_Tr_TOT_bootstrped = [];
CI_X_Tr = [];
CI_Y_Tr = [];
ts_TOT = [];
for j=1:size(X_Tr_size,2)
    rng('shuffle') %For reproducibility
    temp = [];
    size_temp = [];
    for i=1:size(X_Tr_size,1)
        if(isempty(X_Tr_TOT{i,j}))
            continue;
        end
        temp = cat(2,temp,X_Tr_TOT{i,j});
        size_temp = cat(2,size_temp,X_Tr_size{i,j});
    end
    
    if(isempty(size_temp))
         ts_TOT = cat(2,ts_TOT,0);
         X_Tr_TOT_bootstrped = cat(2,X_Tr_TOT_bootstrped,zeros(100,1));
         Y_Tr_TOT_bootstrped = cat(2,Y_Tr_TOT_bootstrped,zeros(100,1));
         CI_X_Tr = cat(1,CI_X_Tr,0);
         CI_Y_Tr = cat(1,CI_Y_Tr,0);
        continue;
    end
    ts = tinv([0.025  0.975],max(1,size(size_temp,2)-1));      % T-Score
    CI_X_Tr = cat(1,CI_X_Tr,ts(1)*std(temp,0,2)/sqrt(size(size_temp,2)));
    if(size(size_temp,2)>1)
    X_Tr_TOT_bootstrped = cat(2,X_Tr_TOT_bootstrped,bootstrp(100,@mean,temp','Weights',size_temp'));
    else
        X_Tr_TOT_bootstrped = cat(2,X_Tr_TOT_bootstrped,temp*ones(100,1));
    end
    
    temp = [];
    size_temp = [];
    for i=1:size(Y_Tr_size,1)
        if(isempty(Y_Tr_TOT{i,j}))
            continue;
        end
        temp = cat(2,temp,Y_Tr_TOT{i,j});
        temp(isnan(temp)) = 1;
        size_temp = cat(2,size_temp,Y_Tr_size{i,j});
    end
    ts = tinv([0.025  0.975],max(1,size(size_temp,2)-1));      % T-Score
    ts_TOT = cat(2,ts_TOT,ts(1));
    CI_Y_Tr = cat(1,CI_Y_Tr,ts(1)*std(temp,0,2)/sqrt(size(size_temp,2)));
    if(size(size_temp,2)>1)
    Y_Tr_TOT_bootstrped = cat(2,Y_Tr_TOT_bootstrped,bootstrp(100,@mean,temp','Weights',size_temp'));
    else
        Y_Tr_TOT_bootstrped = cat(2,Y_Tr_TOT_bootstrped,temp*ones(100,1));
    end
end

SE_X_Tr = std(X_Tr_TOT_bootstrped,0,1)';%(std(X_Tr_TOT,0,1)/sqrt(size(X_Tr_TOT,1)))';
SE_Y_Tr = std(Y_Tr_TOT_bootstrped,0,1)';%(std(Y_Tr_TOT,0,1)/sqrt(size(Y_Tr_TOT,1)))';
CI_Y_Tr = ts_TOT'.*SE_Y_Tr;                      % Confidence Intervals
CI_X_Tr = ts_TOT'.*SE_X_Tr; 


ind_Tr = [];
FPR_min = min(X_Tr_avg);
FPR_max = max(X_Tr_avg);
for FPR_Thr = FPR_min:0.001:min(FPR_max,1-1e-4)
    [~,TPR_ind] = max(Y_Tr_avg(X_Tr_avg<=FPR_Thr));
    ind_temp = find(X_Tr_avg<=FPR_Thr);
    ind_temp = ind_temp(TPR_ind);
    ind_Tr = cat(1,ind_Tr,ind_temp);
end

ind_Tr = unique(ind_Tr);


ind_TOTAL = ind_Tr;%intersect(ind_Val,ind_Tr);

[X_Tr_avg,ind_srt] = sort(X_Tr_avg(ind_TOTAL));
Y_Tr_avg = Y_Tr_avg(ind_TOTAL);
Y_Tr_avg = Y_Tr_avg(ind_srt);
Y_Tr_size_TOT = Y_Tr_size_TOT(ind_TOTAL);
Y_Tr_size_TOT = Y_Tr_size_TOT(ind_srt);
X_Tr_size_TOT = X_Tr_size_TOT(ind_TOTAL);
X_Tr_size_TOT = X_Tr_size_TOT(ind_srt);
CI_X_Tr = CI_X_Tr(ind_TOTAL);
CI_X_Tr = CI_X_Tr(ind_srt);
CI_Y_Tr = CI_Y_Tr(ind_TOTAL);
CI_Y_Tr = CI_Y_Tr(ind_srt);

X_Val_avg = X_Val_avg(ind_TOTAL);
X_Val_avg = X_Val_avg(ind_srt);
Y_Val_avg = Y_Val_avg(ind_TOTAL);
Y_Val_avg = Y_Val_avg(ind_srt);
Y_Val_size_TOT = Y_Val_size_TOT(ind_TOTAL);
Y_Val_size_TOT = Y_Val_size_TOT(ind_srt);
X_Val_size_TOT = X_Val_size_TOT(ind_TOTAL);
X_Val_size_TOT = X_Val_size_TOT(ind_srt);
CI_X_Val = CI_X_Val(ind_TOTAL);
CI_X_Val = CI_X_Val(ind_srt);
CI_Y_Val = CI_Y_Val(ind_TOTAL);
CI_Y_Val = CI_Y_Val(ind_srt);


ind_Val = [];
FPR_min = min(X_Val_avg);
FPR_max = max(X_Val_avg);
for FPR_Thr = FPR_min:0.001:min(FPR_max,1-1e-4)
    [~,TPR_ind] = max(Y_Val_avg(X_Val_avg<=FPR_Thr));
    ind_temp = find(X_Val_avg<=FPR_Thr);
    ind_temp = ind_temp(TPR_ind);
    ind_Val = cat(1,ind_Val,ind_temp);
end

ind_Val = unique(ind_Val);


ind_TOTAL = ind_Val;%intersect(ind_Val,ind_Tr);

[X_Tr_avg,ind_srt] = sort(X_Tr_avg(ind_TOTAL));
Y_Tr_avg = Y_Tr_avg(ind_TOTAL);
Y_Tr_avg = Y_Tr_avg(ind_srt);
Y_Tr_size_TOT = Y_Tr_size_TOT(ind_TOTAL);
Y_Tr_size_TOT = Y_Tr_size_TOT(ind_srt);
X_Tr_size_TOT = X_Tr_size_TOT(ind_TOTAL);
X_Tr_size_TOT = X_Tr_size_TOT(ind_srt);
CI_X_Tr = CI_X_Tr(ind_TOTAL);
CI_X_Tr = CI_X_Tr(ind_srt);
CI_Y_Tr = CI_Y_Tr(ind_TOTAL);
CI_Y_Tr = CI_Y_Tr(ind_srt);

[X_Val_avg,ind_srt] = sort(X_Val_avg(ind_TOTAL));
Y_Val_avg = Y_Val_avg(ind_TOTAL);
Y_Val_avg = Y_Val_avg(ind_srt);
Y_Val_size_TOT = Y_Val_size_TOT(ind_TOTAL);
Y_Val_size_TOT = Y_Val_size_TOT(ind_srt);
X_Val_size_TOT = X_Val_size_TOT(ind_TOTAL);
X_Val_size_TOT = X_Val_size_TOT(ind_srt);
CI_X_Val = CI_X_Val(ind_TOTAL);
CI_X_Val = CI_X_Val(ind_srt);
CI_Y_Val = CI_Y_Val(ind_TOTAL);
CI_Y_Val = CI_Y_Val(ind_srt);

%% Better ROC plots:

figure;
hold on
p = fill([(X_Val_avg'), fliplr(X_Val_avg')], [(Y_Val_avg'+CI_Y_Val'), fliplr(Y_Val_avg'-CI_Y_Val')], 'r');   % Shaded Confidence Intervals
p.FaceColor = [1 0.8 0.8];      
p.EdgeColor = 'none'; 
p = fill([(X_Val_avg'+CI_X_Val'), fliplr(X_Val_avg'-CI_X_Val')], [(Y_Val_avg'), fliplr(Y_Val_avg')], 'r');   % Shaded Confidence Intervals
p.FaceColor = [1 0.8 0.8];      
p.EdgeColor = 'none';
hold on

p = fill([(X_Tr_avg'), fliplr(X_Tr_avg')], [(Y_Tr_avg'+CI_Y_Tr'), fliplr(Y_Tr_avg'-CI_Y_Tr')], 'b');   % Shaded Confidence Intervals
p.FaceColor = [1 0.8 0.8];      
p.EdgeColor = 'none'; 
p = fill([(X_Tr_avg'+CI_X_Tr'), fliplr(X_Tr_avg'-CI_X_Tr')], [(Y_Tr_avg'), fliplr(Y_Tr_avg')], 'b');   % Shaded Confidence Intervals
p.FaceColor = [0.8 0.8 1];      
p.EdgeColor = 'none';
hold on
plot(X_Tr_avg,Y_Tr_avg,'-b.','LineWidth',2)
hold on
plot(X_Val_avg,Y_Val_avg,'-r.','LineWidth',2)

ppv = 0.50;
FPR = 0:0.0001:0.05;
NOP = X_Val_size_TOT(6)/Y_Val_size_TOT(6);
TPR = NOP.*(ppv/(1-ppv))*FPR;
hold on
plot(FPR,TPR,'-g')

PPV_error = 1-((Y_Val_avg.*Y_Val_size_TOT)./(Y_Val_avg.*Y_Val_size_TOT+X_Val_avg.*X_Val_size_TOT));

dist = sqrt(PPV_error.^2+(1-Y_Val_avg).^2);


%% Getting the average velocity of the detections (false and true positives):
Average_Vel_tot_TOT = [];
i_TOT = [];
% figure;
% hold on
Vel_ind = ind_Tr(ind_Val(ind_srt));
Vel_ind = Vel_ind(5);%%%%%%%%%%%%%%%OPTIMAL OPERATION POINT%%%%%%%%%
for i=1:size(Thr_train_TOT,1)
    for subj = 1:size(Sub_ID_tot,2)
        cd(['E:\SD-II_POST_ICA\SD-II\',Sub_ID{Sub_ID_tot(i,subj)},'\',Sub_ID{Sub_ID_tot(i,subj)}])
        Session_names = dir('Patient*');
        Sessions_size = size(Session_names,1);
        for ss = 1:Sessions_size
%             cd(Session_names(ss).name)
           
                        Thr_1 = Thr_tot(Vel_ind,1)
                        Thr_2 = Thr_tot(Vel_ind,2)
                        Thr_3 = Thr_tot(Vel_ind,3)
                        
                        load([Session_names(ss).name,'\',Session_names(ss).name ,'_',sprintf('Thr_1_%.2d',Thr_1)...
                            ,'_',sprintf('Thr_2_%.2d',Thr_2),'_',sprintf('Thr_3_%.2d',Thr_3), '_detection_SDIII_240minOv180min_CCWL20_Theta_LowResVelFixedvTEST_VelocitySaveII.mat'])
                        Part_names = dir([Session_names(ss).name,'\',Session_names(ss).name,'_CSD_ind_*_Theta_SDII_vConfusion_240minOv180min_CCWL20min_Xcrr_Pruned_vVII.mat']);
                        Part_size = size(Part_names,1);
                        Part_name_char = string({Part_names.name});
                        part_ind = strfind(Part_name_char,'ind');
                        of_ind = strfind(Part_name_char,'Theta');
                        Part_num = [];
                        if(Part_size>1)
                            for j=1:Part_size
                                Part_num = cat(1,Part_num,str2double(Part_name_char{j}(part_ind{j}+4:of_ind{j}-2)));
                            end
                        else
                            Part_num = str2double(Part_name_char{1}(part_ind+4:of_ind-2));
                        end
                        [~,Part_ind] = sort(Part_num);
                        load_temp = load([Session_names(ss).name,'\',Part_names(Part_ind(1)).name])
                        CSD_GT_total = load_temp.CSD_GT_total;
                        Strt_tot_TOTAL_temp = Strt_tot_TOTAL;
                        Strt_tot = Strt_tot_TOTAL_temp{1,find(THR_4==Thr_tot(Vel_ind,4))};
                        Endd_tot_TOTAL_temp = Endd_tot_TOTAL;
                        Endd_tot = Endd_tot_TOTAL_temp{1,find(THR_4==Thr_tot(Vel_ind,4))};
                        Average_Vel_tot_TEMP = Average_Vel_tot_TOTAL{1,find(THR_4==Thr_tot(Vel_ind,4))};

                        Strt = Strt_tot;
                        Strt = Strt;
                        Endd = Endd_tot;
                        Endd = Endd;
                        for ind_sub=1:size(Strt,1)
                           Detection_ind = find(CSD_GT_total(max(1,Strt(ind_sub)*64-60*60*64):min(Endd(ind_sub)*64+60*60*64,size(CSD_GT_total,2))));
                           if(~isempty(Detection_ind))
                              Average_Vel_tot_TOT = cat(1,Average_Vel_tot_TOT,Average_Vel_tot_TEMP(ind_sub)); 
                              
                           end
                               
                        end
                        
        end
    end
end
Average_Vel_tot_TOT(Average_Vel_tot_TOT<0.5) = [];
figure;histogram(Average_Vel_tot_TOT,linspace(min(Average_Vel_tot_TOT),max(Average_Vel_tot_TOT),20));
min(Average_Vel_tot_TOT)
max(Average_Vel_tot_TOT)
mean(Average_Vel_tot_TOT)
CI = 1.96*std(Average_Vel_tot_TOT)/sqrt(size(Average_Vel_tot_TOT,1))
