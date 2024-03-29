function [] = Detection_parallelized_LowRes_vF(Sub_ID,ss_start,ss_End,LR)

%% Function inputs  
% Sub_ID: Patient ID
% ss_start, ss_end: Session Range for Processing
% LR: 'L' or 'R' representing craniectomy hemisphere 

% clear all
% close all
% clc

%%
% Main SD detection code for WAVEFRONT, along with performance calculations
%
% See: README.txt and [1] for more info.

% [1] Alireza Chamanzar, Jonathan Elmer, Lori Shutter, Jed Hartings, Pulkit Grover,
%  "Noninvasive and reliable automated detection of spreading depolarization in severe traumatic brain injury using scalp EEG",
%   Submitted to Nature Comms Med, 2022.

% Author: Alireza 	Date: 2022/01/05 12:00:00 	Revision: 0.1
% Copyright: This project is licensed - see the LICENSE.md file for details


% Sub_ID = '04-1201';
% LR='R';
% ss_start = 1;
% ss_End = 1;
% im_scale = 0.1;


%% Initialization of Parameters:
crani_ind = 1:19;
if(LR=='L')
    LR_sgn = -1;
    crani_ind = crani_ind([1,3,4,6,8,10,11,13,15,17,18]);
else
    LR_sgn = 1;
    crani_ind([1,4,6,8,11,13,15,18]) = [];
end

im_scale = 0.0945;
DELTA_T = 60;

% Grid search for set of parameters in WAVEFRONT:
% THR_1 = 0.05:0.05:0.8;%mm
% THR_3 = linspace(0.1,0.9,20);
% % THR_1 = THR_1([1,2,3,4,5,6,7,8,9]);%mm
% THR_2 = 0.6:0.1:0.8;
% % THR_3 = THR_3(10:15);
% THR_4 = 2:10;

% OPTIMAL set of parameters (for quick run through the code):
THR_1 = 0.300000000000000 ;
THR_2 = 0.600000000000000;
THR_3 = 0.689473684210526;
THR_4 = 2.0000;


cd(['D:\SD-II_POST_ICA\SD-II\',Sub_ID,'\',Sub_ID])

Session_names = dir('Patient*');
Sessions_size = size(Session_names,1);

TPR_num = 0;
TPR_denum = 0;
TNR_num = 0;
TNR_denum = 0;
Tot_detectable_CSD = 0;
Detected_CSD = 0;
CSD_event_quality_tot = [];

% Beginning of Session-Specific Loop
for ss = ss_start:ss_End
    cd(Session_names(ss).name)
    
    current_path = pwd;
    EEG_vis = pop_loadset('filename',[Session_names(ss).name,'_Visualization_Delta.set'],'filepath',current_path);
    
    % Filter unnecessary events  
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
            continue;
        end
        
    end
    
    EEG_vis = pop_editeventvals(EEG_vis,'delete',event_delete_ind);
    
    [ind_r,ind_c] = find(isnan(EEG_vis.data));
    EEG_vis.data(ind_r,ind_c) = 0;
    
    % Load Global Mean data (Session-wide average of EEG_PW)
    Global_M = load([Session_names(ss).name ,'GlobalMean_Delta_SDII_vConfusion_240minOv180min_CCWL20min_Xcrr_Improved_Pruned_vVII.mat']);
    Global_M = Global_M.Global_M;
    Part_names = dir([Session_names(ss).name,'_CSD_ind_*_Delta_SDII_vConfusion_240minOv180min_CCWL20min_Xcrr_Pruned_vVII.mat']);

    Part_size = size(Part_names,1);
    Part_name_char = string({Part_names.name});
    part_ind = strfind(Part_name_char,'ind');
    of_ind = strfind(Part_name_char,'Delta');
    Part_num = [];
    if(Part_size>1)
        for i=1:Part_size
            Part_num = cat(1,Part_num,str2double(Part_name_char{i}(part_ind{i}+4:of_ind{i}-2)));
        end
    else
        Part_num = str2double(Part_name_char{1}(part_ind+4:of_ind-2));
    end
    [~,Part_ind] = sort(Part_num);

    % Load EEG Electrode Spatial Locations 
    UPMC_TBI_EEG_Loc = load('D:\SilenceMap-v1.0\Chamanzar-SilenceMap-668f8ab\EEG_prep\UPMC_TBI_EEG_Loc_SDII_ordered.mat');
    UPMC_TBI_EEG_Loc = UPMC_TBI_EEG_Loc.UPMC_TBI_EEG_Loc;
    
    BCNT_tot = {};
    CNT_tot = {};
    CSD_GT_tot = {};
    T_start_tot = {};
    T_end_tot = {};
    Events_tot = {};
    CSD_type_tot = {};
    Events_type_tot = {};
    num_good_elec_tot = {};
    CSD_GT_total = {};
    CSD_Labels_total = {};
    Strt_tot = {};
    Endd_tot = {};
    T_start_tot_temp = {};
    T_end_tot_temp = {};
    BCNT_tot_temp = {};
    CNT_tot_temp = {};
    CSD_GT_tot_temp = {};
    num_good_elec_tot_temp = {};
    CSD_type_tot_tempp = {};
    Endd_tot_temp = {};
    Strt_tot_temp = {};
    Thr_counter_tot = {};
    
    %% Beginning of Part-Specific Loop
    for CSD_ind_orig = 1:Part_size
        
        
        T_start_tot_temp_temp = {};
        T_end_tot_temp_temp = {};
        BCNT_tot_temp_temp = {};
        CNT_tot_temp_temp = {};
        CSD_GT_tot_temp_temp = {};
        num_good_elec_tot_temp_temp = {};
        CSD_type_tot_temp_temp = {};
        Endd_tot_temp_temp = {};
        Average_Vel_tot_temp_temp = {};
        Strt_tot_temp_temp = {};
        Thr_counter_tot_temp = {};
        
        
        
        fprintf(['Starting ' sprintf('CSD_ind_%d',CSD_ind_orig) '\n'])
        Name = Session_names(ss).name;
        load_temp = load(Part_names(Part_ind(CSD_ind_orig)).name)
        CSD_GT = load_temp.CSD_GT;
        CSD_GT_total{CSD_ind_orig} = load_temp.CSD_GT_total;
        CSD_Labels_total{CSD_ind_orig} = load_temp.CSD_Labels_total;
        CSD_type = load_temp.CSD_type;
        EEG_PW = load_temp.EEG_PW;
        EEG_mask = load_temp.EEG_mask;
        Events = load_temp.Events;
        Events_type = load_temp.Events_type;
        T_end = load_temp.T_end
        T_start = load_temp.T_start
        srate = load_temp.srate;
        
        
        %% Spatial projection on a 2D plane:
        sensor_locs_n = UPMC_TBI_EEG_Loc(1:size(UPMC_TBI_EEG_Loc,1)-1,1:3)/100;
        
        %Reorder the electrodes based on SD-II dataset:
        [~,sensor_ind] = sort(UPMC_TBI_EEG_Loc(1:size(UPMC_TBI_EEG_Loc,1)-1,4),'ascend');
        sensor_locs_n = sensor_locs_n(sensor_ind,:);
        
        % Get Part specific event info including Event type and timestamp
        if(isempty(CSD_type))
            CSD_type = {'1_No CSD'};
            CSD_GT(1) = 1;
        end
        CSD_ind = Part_num(Part_ind(CSD_ind_orig));
        
        Sub_sampl = 30*srate;
        EEG_PW(EEG_PW<0) = 0;        
        EEG_PW = EEG_PW(:,1:Sub_sampl:end);
        EEG_mask = EEG_mask(:,1:Sub_sampl:end);
        CSD_GT_ind = find(CSD_GT);
        CSD_GT_ind = round(CSD_GT_ind/Sub_sampl);
        CSD_GT_ind(CSD_GT_ind==0) = 1;
        CSD_GT = zeros(1,size(EEG_PW,2));
        CSD_GT(CSD_GT_ind) = 1;
        CSD_GT_temp = CSD_GT;
                
        %%
        Normal_events_ind = (strcmp(Events_type,'Data Quality Normal'));
        Events_ind = find(Events);
        Events(Events_ind(Normal_events_ind)) = 0;
        Events_type((Normal_events_ind)) = [];
        
        Events_ind = find(Events);
        Events_ind = round(Events_ind/Sub_sampl);
        Events = zeros(1,size(EEG_PW,2));
        for event_id=1:size(Events_ind,2)
            if(Events_ind(event_id)==0)
                Events_ind(event_id)=1;
            end
            temp = Events_ind(event_id);
            while(Events(temp)==1 && temp<size(Events,2))
                temp = temp+1;
            end
            Events(temp) = 1;
        end
        
        Events_temp = Events;
        
        % Calculate spatial loc of electrodes in cart. and spherical coords

        Alphax = 0;%*(rand(1)-0.5)+n*0;
        Alphay = 0;%80*(rand(1)-0.5)+n*0;
        Alphaz = 90;
        
        Rx = [1 0 0; 0 cos((pi/180)*Alphax) -sin((pi/180)*Alphax); 0 sin((pi/180)*Alphax) cos((pi/180)*Alphax)];
        Ry = [cos((pi/180)*Alphay) 0 sin((pi/180)*Alphay); 0 1 0; -sin((pi/180)*Alphay) 0 cos((pi/180)*Alphay)];
        Rz = [cos((pi/180)*Alphaz) -sin((pi/180)*Alphaz) 0; sin((pi/180)*Alphaz) cos((pi/180)*Alphaz) 0; 0 0 1];
        
        
        Rtot = Rz * Ry * Rx;
        sensor_locs_n = (Rtot*sensor_locs_n')';
        
        
        [phi,theta,r]=cart2sph(sensor_locs_n(:, 1), sensor_locs_n(:, 2), sensor_locs_n(:, 3));
        phi((phi<0))=phi((phi<0))+2*pi;
        
        
        theta((theta>=0))=((pi/2)-theta((theta>=0)));
        theta((theta<0))=(-1)*(theta((theta<0))-(pi/2));
        
        theta=round((180/pi)*theta);
        phi=round((180/pi)*phi);
        phi((phi==0))= phi((phi==0))+1;
        phi((phi==360))= 1;
        theta((theta==0))= theta((theta==0))+1;
        
        % Project Data into 2D space using above spatial conversions 
        EEG_median = (median(EEG_PW(crani_ind,:),1));
        I_EEG = ones(180,359,size(EEG_PW,2));
        I_mask = zeros(180,359,size(EEG_mask,2));
        
        for k=1:size(EEG_PW,2)
            I_EEG(:,:,k) = I_EEG(:,:,k)*EEG_median(1,k);
            for i=1:size(EEG_PW,1)
                if(LR_sgn==1)
                    if(theta(i)>0 && theta(i)<=180 && phi(i)>0 && phi(i)<185)
                        
                        phi_temp = phi(i);
                        if(phi(i)==1 && theta(i)<20)
                            phi_temp = 90;
                        end
                        
                        I_EEG(theta(i),phi_temp,k) = EEG_PW(i,k);
                        I_mask(theta(i),phi_temp,k) = EEG_mask(i,k);
                    end
                else
                    if(theta(i)>0 && theta(i)<=180 && (phi(i)>175 || phi(i)<5) && phi(i)<360)
                        phi_temp = phi(i);
                        if(phi(i)==1 && theta(i)<20)
                            phi_temp = 90;
                        end
                        if(phi(i)>175)
                            phi_temp = 180 - (phi(i)-180);
                        end
                        
                        
                        I_EEG(theta(i),phi_temp,k) = EEG_PW(i,k);
                        I_mask(theta(i),phi_temp,k) = EEG_mask(i,k);
                    end
                end
            end
        end
                
        %% EEG spatial interpolation:
        
        % Smooth 'Sparse' spatial projections, creating 'I_smoothed'
        %Gaussian:
        figure;mesh(I_EEG(:,:,309))
        h = ones(30);
        for k=1:size(I_EEG,3)
            I_EEG(:,:,k) = imgaussfilt(I_EEG(:,:,k),20,'Padding',EEG_median(1,k));%resample(imgaussfilt(temp,10,'Padding',0),size(I_EEG,1),size(temp,1),5);
            I_mask(:,:,k) = imgaussfilt( I_mask(:,:,k),20,'Padding',0);%resample(imgaussfilt(temp,10,'Padding',0),size(I_mask,1),size(temp,1),5);
        end
        I_mask = double(I_mask>0);
        figure;mesh(I_EEG(:,:,309))
        I_EEG = I_EEG(1:100,:,:);

        %
        %% Temporal subsampling of frames II:
        Thr_counter = 0;
        Thr_ind = 0;

        % Begin Parameter-specific loops, starting with Thr_3
        for Thr_3_ind = 1:size(THR_3,2)
            Thr_3 = THR_3(Thr_3_ind)
            
            II = I_EEG(end:-1:1,:,:);
            num_good_elec = sum(EEG_mask(crani_ind,:),1);
            II_median = (median(reshape(II,[],1,size(II,3)),1));
            
            % Thresholding the frames to obtain binary images (I_BW):
            Thr = Thr_3;
            hh = max(max(II));
            ll = min(min(II));
            temp = II;
            II(temp<=II_median+Thr*(hh-II_median)) = 0;
            II(temp>II_median+Thr*(hh-II_median)) = 255;
            II(:,:,num_good_elec<5)=0;
            CC_BW = squeeze(sum(sum(II,1),2))';
            II(:,:,CC_BW>=0.2*255*size(II,1)*size(II,2))=0;
            
            
            II_resize = [];
            for t=1:size(II,3)
                II_resize = cat(3,II_resize,imresize(squeeze(II(:,:,t)),im_scale));
            end

            II = II_resize;

            value = II;
            activity = (II(end:-1:1,:,:));

            %% Calculating Optical flows:
            
            Ttotal = size(activity,3);
            %Orange color
            C = [1 .5 0];
            
            Stp = 1;
            
            disp('Computing Optical Flow');
            % Compute optical flow 
            opticalFlow = opticalFlowHS;
            of = estimateFlow(opticalFlow, value(:,:,1));
            % the variable 'of' below stores optical flows at each 
            % x,y as real and imaginary values.
            bb = {};
            cc = {};
            for t = 1:Ttotal-Stp
                of(t)=estimateFlow(opticalFlow, value(:,:,t+Stp));
                label = bwlabel(activity(:,:,t));
                st = regionprops( label,'centroid', 'Area', 'BoundingBox' );
                bb{t} = round(cat(1,st.BoundingBox)); % Bounding box for connected components
                cc{t} = round(cat(1,st.Centroid)); % Centroid for connected components
            end
            
            
            %% Calculating the histogram of orientations for all BBoxes:
            vel_x = zeros(size(II,1),size(II,2),size(of,2));
            vel_y = zeros(size(II,1),size(II,2),size(of,2));
            Theta_scale = 100/(size(II,1)-1);
            Phi_scale = 359/(size(II,2)-1);
            for i=1:size(of,2)
                vel_x(:,:,i) = of(i).Vx(end:-1:1,:).*cos((pi/180)*repmat(linspace(90,-10,size(II,1))',1,size(II,2)))*75*Phi_scale*(pi/180);
                vel_y(:,:,i) = of(i).Vy(end:-1:1,:)*75*Theta_scale*(pi/180);
            end
            
            
            
            % Old code for plotting optical flows
%             Fr=[];
%             
%             for r=309
%                 figure;
%                 %                                 r = 225;%t=15min;
%                 set(gcf, 'Position', get(0, 'Screensize'));
%                 imshow(value(end:-1:1,:,r),'InitialMagnification','fit');
%                 hold on
%                 %ploting BBoxs:
%                 for i = 1:size(bb{r},1)
%                     rectangle('Position',bb{r}(i,:),'EdgeColor','r','LineWidth',4)
%                 end
%                 quiver(of(r).Vx(end:-1:1,:).*cos((pi/180)*repmat(linspace(90,-10,size(II,1))',1,size(II,2))),(-of(r).Vy(end:-1:1,:)),1,'color',C,'LineWidth',2)
%                 Fr = cat(1,Fr,getframe(gcf));
%                 hold off
%             end
            
            % Calculate magnitude and orientation of optical flows 
            velocity = {};
            for t =1:size(vel_x,3)
                if isempty(bb{t})
                    velocity{t,1} = [];
                    continue;
                end
                for i = 1:size(bb{t},1)
                    if((bb{t}(i,3)* bb{t}(i,4))<1) % Checking area of bounding box is >4
                        continue;
                    end
                    
                    orient =  atan2d(vel_y(bb{t}(i,2):bb{t}(i,2)+bb{t}(i,4)-1,bb{t}(i,1):bb{t}(i,1)+bb{t}(i,3)-1,t),...
                        vel_x(bb{t}(i,2):bb{t}(i,2)+bb{t}(i,4)-1,bb{t}(i,1):bb{t}(i,1)+bb{t}(i,3)-1,t));
                    mag = (vel_y(bb{t}(i,2):bb{t}(i,2)+bb{t}(i,4)-1,bb{t}(i,1):bb{t}(i,1)+bb{t}(i,3)-1,t).^2+...
                        vel_x(bb{t}(i,2):bb{t}(i,2)+bb{t}(i,4)-1,bb{t}(i,1):bb{t}(i,1)+bb{t}(i,3)-1,t).^2).^0.5;
                    orient(orient<0) = 360 + orient(orient<0);
                    orient(mag==0) = -1;
                    
                    [Hof,idx] = histcounts(reshape(orient,1,[]),[linspace(0,360,9)-22.5,360+22.5]);
                    Mag = mag;
                    
                    velocity{t,i} = struct(...
                        'bbox',bb{t}(i,:),...
                        'cc',cc{t}(i,:),...
                        'mag',Mag, ...
                        'hof',Hof,...
                        'idx',idx,...
                        'angle',orient);
                end
            end
            
            
            %% Histogram @ 10min(86):
%             % Old code for plotting histogram of BBoxs
%                         r=309;% hidtogram of BBox @ 10 min
%                         figure;polarhistogram('BinEdges',(pi/180)*[linspace(0,360,9)-22.5,360+22.5],'BinCounts',velocity{r, 1}.hof)
%             
            %% Quantization of orientations and creating OBBox:
            velocity_tot = velocity;
            % Begin parameter 'Thr_1' specific-loop
            for Thr_1_ind = 1:size(THR_1,2)
                
                Thr_1 = THR_1(Thr_1_ind)
                velocity = velocity_tot;
                
                Fr = [];
                VX = zeros(size(II,1),size(II,2),size(velocity,1));
                VY = zeros(size(II,1),size(II,2),size(velocity,1));
                SCR_tot = zeros(size(II,1),size(II,2),size(velocity,1));
                ANgle = zeros(size(II,1),size(II,2),size(velocity,1));
                Bin_Cntr = linspace(0,360,9);
                Bin_Cntr(end) = [];
                
                bbb = {};
                ccc = {};
                BB = {};
                for t=1:size(velocity,1)
                    Vx = zeros(size(vel_x,1),size(vel_x,2),size(velocity,2));
                    Vy = zeros(size(vel_y,1),size(vel_y,2),size(velocity,2));
                    scr_tot = zeros(size(vel_x,1),size(vel_x,2),size(velocity,2));
                    Angls = -1*ones(size(Vx,1),size(Vx,2),size(velocity,2));
                    for i=1:size(velocity,2)
                        if isempty(velocity{t,i})
                            BB{t,1} = [];
                            continue;
                        end
                        
                        vel = velocity{t,i}.angle;
                        vel(vel==-1) = [];
                        
                        mag = velocity{t,i}.mag;
                        mag(mag==0) = [];
                        
                        % Score BBox's based upon maag and orietnation 
                        % compare to Thr_1 to remove pop-ups and fades 
                        scr = sqrt(mean(mean(mag.*cos(vel*pi/180))).^2 + mean(mean(mag.*sin(vel*pi/180))).^2);
                        
                        if(isempty(vel))
                            scr = 0;
                        end
                        
                        if(scr<Thr_1)
                            velocity{t,i}.angle = -1 *  ones(size(velocity{t,i}.angle));
                        end
                        
                        Hof=idx(velocity{t,i}.hof>=0.5*max(velocity{t,i}.hof));
                        velocity{t,i}.hof=Hof;
                        orient=velocity{t,i}.angle;
                        Orient = -1*ones(size(orient));
                        for j=1:size(Hof,2)
                            Orient(orient>=Hof(j) & orient<(Hof(j)+45))=Hof(j)+22.5;%%%%Central point of the bin
                        end
                        Orient(Orient==360) = 0;
                        
                        velocity{t,i}.angle=Orient;
                        mag = velocity{t,i}.mag;
                        mag(Orient==-1)=0;
                        velocity{t,i}.mag = mag;
                        
                        temp = zeros(size(vel_x,1),size(vel_x,2));
                        temp(velocity{t,i}.bbox(2):velocity{t,i}.bbox(2)+velocity{t,i}.bbox(4)-1,velocity{t,i}.bbox(1):velocity{t,i}.bbox(1)+velocity{t,i}.bbox(3)-1)=(mag.*cos((pi/180)*Orient));
                        Vx(:,:,i) = temp;
                        temp = zeros(size(vel_y,1),size(vel_y,2));
                        temp(velocity{t,i}.bbox(2):velocity{t,i}.bbox(2)+velocity{t,i}.bbox(4)-1,velocity{t,i}.bbox(1):velocity{t,i}.bbox(1)+velocity{t,i}.bbox(3)-1)=(mag.*sin((pi/180)*Orient));
                        Vy(:,:,i) = temp;
                        temp_Angle = -1*ones(size(vel_x,1),size(vel_x,2));
                        temp_Angle(velocity{t,i}.bbox(2):velocity{t,i}.bbox(2)+velocity{t,i}.bbox(4)-1,velocity{t,i}.bbox(1):velocity{t,i}.bbox(1)+velocity{t,i}.bbox(3)-1)=(Orient);
                        Angls(:,:,i) = temp_Angle;
                        scr_tot(velocity{t,i}.bbox(2):velocity{t,i}.bbox(2)+velocity{t,i}.bbox(4)-1,velocity{t,i}.bbox(1):velocity{t,i}.bbox(1)+velocity{t,i}.bbox(3)-1,i) = scr;
                        
                    end
                    
                    for i=1:size(velocity,2)
                        [i_r,i_c] = find(Angls(:,:,i)~=-1);
                        Angls(i_r,i_c,1) = Angls(i_r,i_c,i);
                    end
                    Angls = Angls(:,:,1);
                    Vx = sum(Vx,3);
                    Vy = sum(Vy,3);
                    scr_tot = sum(scr_tot,3);
                    
                    %Creating the OBBoxes (Oriented Bounding Boxes):
                    ANgle(:,:,t) = Angls;
                    VX(:,:,t) = Vx;
                    VY(:,:,t) = Vy;
                    SCR_tot(:,:,t) = scr_tot;
                    C=unique(Angls);
                    C(C==-1)=[];
                    if (size(C,2)>0)
                        for q=1:size(C,1)
                            temp=Angls;
                            temp(temp~=C(q))=0;
                            temp(temp~=0)=255;
                            label(:,:,t)=bwlabel(temp);
                            st = regionprops( label(:,:,t),'centroid', 'BoundingBox' );
                            bbb{t} = round(cat(1,st.BoundingBox)); % Bounding box for connected components of Orient
                            ccc{t} = round(cat(1,st.Centroid)); % Bounding box for connected components of Orient
                            
                            if (isempty(bbb{t})==0)
                                p=0;
                                for o = 1:size(bbb{t},1)
                                    if(bbb{t}(o,3)*bbb{t}(o,4)>1)
                                        p=p+1;
                                        BB{t,find(Bin_Cntr==C(q))}(p,:) = cat(2,bbb{t}(o,:),C(q),ccc{t}(o,:));
                                    end
                                end
                            end
                        end
                    end
                end
                
                if(size(BB,1)<size(velocity,1))
                    for t=size(BB,1)+1:size(velocity,1)
                        for dum=1:size(BB,2)
                            BB{t,dum}=[];
                        end
                    end
                end
                
                %% Scoring the OBBox based on consistency in orientation & speed:
                
                % Begin 'Thr_2' specific-loop
                for Thr_2_ind = 1:size(THR_2,2)

                    % Begin 'Thr_4' specific-loop           
                    for Thr_4_ind = 1:size(THR_4,2)
                        
                        Thr_2 = THR_2(Thr_2_ind)
                        
                        Thr_4 = THR_4(Thr_4_ind)
                        
                        Threshld_Frame=Thr_2;%0.1:0.2:0.9
                        SCR = BB;
                        Vel_SCR = BB;
                        L_wind = floor(Thr_4*(60*srate/(Sub_sampl)));%temporal neighborhood 2min
                        Gap_q = floor(0*(60*srate/(Sub_sampl)));%temporal neighborhood 2min
                        
                        for t=L_wind+1:size(BB,1)-L_wind
                            for i=1:size(BB,2)
                                SCR{t,i} = zeros(size(SCR{t,i},1),1);
                                Vel_SCR{t,i} = zeros(size(Vel_SCR{t,i},1),1);
                                if (isempty(BB{t,i}))
                                    continue;
                                end
                                for j=1:size(BB{t,i},1)
                                    VV = mean(mean(sqrt((VX(BB{t,i}(j,2):BB{t,i}(j,2)+BB{t,i}(j,4)-1,BB{t,i}(j,1):BB{t,i}(j,1)+BB{t,i}(j,3)-1,t)).^2 ...
                                        +(VY(BB{t,i}(j,2):BB{t,i}(j,2)+BB{t,i}(j,4)-1,BB{t,i}(j,1):BB{t,i}(j,1)+BB{t,i}(j,3)-1,t)).^2)));
                                    if(((60*srate/(Sub_sampl))*VV)>8 || ((60*srate/(Sub_sampl))*VV)<0.5)%speed of propagation between 0.5 to 8mm/min -> CSD candidate
                                        SCR{t,i}(j)=0;
                                        Vel_SCR{t,i}(j)=-1;
                                        continue;
                                    end
                                    Vel_SCR{t,i}(j) = Vel_SCR{t,i}(j) + ((60*srate/(Sub_sampl))*VV);
                                    q_tot = [-L_wind:L_wind];
                                    f_SCR = zeros(size(q_tot,2),1);
                                    for q_ind=1:size(q_tot,2)
                                        q = q_tot(q_ind);
                                        ii = i;
                                        if (isempty(BB{t+q,ii})==0)
                                            for w=1:size(BB{t+q,ii},1)
                                                
                                                %avg speed of OBBox:
                                                VV = mean(mean(sqrt(VX(BB{t+q,ii}(w,2):BB{t+q,ii}(w,2)+BB{t+q,ii}(w,4)-1,BB{t+q,ii}(w,1):BB{t+q,ii}(w,1)+BB{t+q,ii}(w,3)-1,t)).^2+ ...
                                                    (VY(BB{t+q,ii}(w,2):BB{t+q,ii}(w,2)+BB{t+q,ii}(w,4)-1,BB{t+q,ii}(w,1):BB{t+q,ii}(w,1)+BB{t+q,ii}(w,3)-1,t)).^2));
                                                SCR_BB = mean(mean(SCR_tot(BB{t+q,ii}(w,2):BB{t+q,ii}(w,2)+BB{t+q,ii}(w,4)-1,BB{t+q,ii}(w,1):BB{t+q,ii}(w,1)+BB{t+q,ii}(w,3)-1,t)));
                                                
                                                %Consistency check in speed and orientation:
                                                Theta_degree = linspace(90,-10,size(II,1));
                                                
                                                if(abs(BB{t+q,ii}(w,6)-BB{t,i}(j,6))<2)
                                                    Delta_BB = (abs(BB{t+q,ii}(w,6)-BB{t,i}(j,6))*75*Phi_scale*(pi/180)*cos((pi/180)*Theta_degree(min(BB{t+q,ii}(w,7),BB{t,i}(j,7)))))^2+...
                                                    (75*Theta_scale*(pi/180)*abs(BB{t+q,ii}(w,7)-BB{t,i}(j,7)))^2;
                                                else
                                                    Delta_BB = (75*Phi_scale*(pi/180)*sum(cos((pi/180)*Theta_degree(round(linspace(min(BB{t+q,ii}(w,7),BB{t,i}(j,7)),max(BB{t+q,ii}(w,7),BB{t,i}(j,7)),abs(BB{t+q,ii}(w,6)-BB{t,i}(j,6)))))),2))^2+...
                                                    (75*Theta_scale*(pi/180)*abs(BB{t+q,ii}(w,7)-BB{t,i}(j,7)))^2;
                                                end
                                                
                                                if(BB{t+q,ii}(w,5)==BB{t,i}(j,5) && ...
                                                        (Delta_BB<(ceil(70*Sub_sampl/(srate*60))^2)) && ...
                                                        (Delta_BB>(ceil(0*Sub_sampl/(srate*60))^2)) && ...
                                                        ((60*srate/(Sub_sampl))*VV)<=8 && (((60*srate/(Sub_sampl))*VV)>=0.5))% && ((LR_sgn*BB{t+q,ii}(w,1))<(LR_sgn*181)))
                                                    
                                                    %OBBox score:
                                                    SCR{t,i}(j)=SCR{t,i}(j)+1;%SCR_BB;%(BB{t+q,ii}(w,4)*BB{t+q,ii}(w,3));
                                                    Vel_SCR{t,i}(j) = Vel_SCR{t,i}(j) + ((60*srate/(Sub_sampl))*VV);
                                                    %frame score (number of frames contributing
                                                    %in SCR):
                                                    f_SCR(q_ind) = 1;
                                                end
                                            end
                                        end
                                    end
                                    if (size(find(f_SCR),1)<Threshld_Frame*size(f_SCR,1))% check if the presence of matching BBoxes is consistent over 2*L_wind frames
                                        SCR{t,i}(j) = 0;
                                        Vel_SCR{t,i}(j) = -1;
                                    else
                                        Vel_SCR{t,i}(j) = Vel_SCR{t,i}(j)/(SCR{t,i}(j)+1);
                                    end
                                end
                            end
                        end
                        
                        %Some adjustment at boundaries:
                        for t=size(SCR,1)-L_wind+1:size(SCR,1)
                            for i=1:size(SCR,2)
                                SCR{t,i} = [];
                                Vel_SCR{t,i} = [];
                            end
                        end
                        
                        for t=1:L_wind
                            for i=1:size(SCR,2)
                                SCR{t,i} = [];
                                Vel_SCR{t,i} = [];
                            end
                        end
                        
                                                
                        Threshld = 0.01;%:0.05:0.99
                        
                        BSCR = SCR; %binary score of OBBox
                        CNT = zeros(size(SCR,1),1);%total score of frames
                        CNT_Vel = zeros(size(SCR,1),1);%total score of frames
                        
                        for t=1:size(SCR,1)
                            temp=[];
                            for ii=1:size(SCR,2)
                                if (isempty(SCR{t,ii}))
                                    continue;
                                end
                                temp=cat(1,temp,max(SCR{t,ii},[],1));
                            end
                            if (isempty(temp))
                                continue;
                            end
                            SCMax=max(temp);
                            for i=1:size(SCR,2)
                                if (isempty(SCR{t,i}))
                                    continue;
                                end
                                for j=1:size(SCR{t,i},1)
                                    if(SCR{t,i}(j)<Threshld*SCMax)% reject OBBoxes with very small score
                                        BSCR{t,i}(j)=0;
                                    else
                                        BSCR{t,i}(j)=1;
                                        CNT(t) = CNT(t)+SCR{t,i}(j);
                                        CNT_Vel(t) = CNT_Vel(t)+(SCR{t,i}(j)*Vel_SCR{t,i}(j));
                                    end
                                end
                            end
                        end
                        CNT_Vel(CNT~=0) = CNT_Vel(CNT~=0)./CNT(CNT~=0);
                        
                        %% Temporal stitching:
                        % reject frames with very small total score
                        CNT(CNT<0.05*median(CNT)) = 0;
                        CNT_Vel(CNT==0) = -1;
                        CNT_T = bwconncomp(CNT);
                        CNT_TT = [];
                        for i=1:size(CNT_T.PixelIdxList,2)
                            CNT_TT = cat(1,CNT_TT,CNT_T.PixelIdxList{1,i});
                        end
                                        
                        BCNT = zeros(size(CNT)); % binary score of frames
                        WL = 1; % time window length threshold = 30 sec
                        
                        %stitching time window = 2 min
                        WL_stitching = floor(2*60*srate/(Sub_sampl)); 
                        CNT_e = CNT;
                        for t=WL:size(CNT_e,1)-WL+1
                            temp = find(CNT_e(t-WL+1:t+WL-1));
                            if(size(temp,1)>=WL)
                                BCNT(t-WL+1) = 1;
                            end
                        end
                        
                        % Reject small time windows & consider overlapping 
                        % connected components (time samps) at boundaries:
                        BCNT_T = bwconncomp(BCNT);
                        
                        for j=1:size(BCNT_T.PixelIdxList,2)
                            if(BCNT_T.PixelIdxList{1,j}(end,1)-BCNT_T.PixelIdxList{1,j}(1,1)+1<WL)
                                BCNT(BCNT_T.PixelIdxList{1,j}) = 0;
                                continue;
                            end
                            Indx = 0;
                            Strt = BCNT_T.PixelIdxList{1,j}(1);
                            for i=size(CNT_T.PixelIdxList,2):-1:1
                                while(size(find(CNT_T.PixelIdxList{1,i}(:,1)==(Strt-WL_stitching:Strt)))>0)
                                    Indx = i;
                                    Strt = CNT_T.PixelIdxList{1,Indx}(1,1)-1;
                                end
                            end
                            
                            if(Indx~=0)
                                BCNT(CNT_T.PixelIdxList{1,Indx}(1,1):BCNT_T.PixelIdxList{1,j}(1)) = 1;
                            end
                            
                            
                            Indx = 0;
                            Endd = BCNT_T.PixelIdxList{1,j}(end);
                            for i=1:size(CNT_T.PixelIdxList,2)
                                while(size(find(CNT_T.PixelIdxList{1,i}(:,1)==(Endd:Endd+WL_stitching)))>0)
                                    Indx = i;
                                    Endd = CNT_T.PixelIdxList{1,Indx}(end,1)+1;
                                end
                            end
                            
                            if(Indx~=0)
                                BCNT(BCNT_T.PixelIdxList{1,j}(end):CNT_T.PixelIdxList{1,Indx}(end,1)) = 1;
                            end
                            
                        end
                        BCNT_T = bwconncomp(BCNT);                        
                        
                        %% Find connected comps of BCNT 
 
                        CNT_T = zeros(size(CNT));
                        CNT_T(CNT_TT) = 1;
                        
                        BCNT_Vis = bwconncomp(BCNT);
                        Strt = [];
                        Endd = [];
                        Average_Vel = [];
                        for i=1:size(BCNT_Vis.PixelIdxList,2)
                            Vel_temp = CNT_Vel(BCNT_Vis.PixelIdxList{1,i}(1,1):BCNT_Vis.PixelIdxList{1,i}(end,1));
                            Average_Vel = cat(1,Average_Vel,mean(Vel_temp(Vel_temp>0)));
                            Strt = cat(1,Strt,uint32((T_start+Sub_sampl*BCNT_Vis.PixelIdxList{1,i}(1,1))/srate));
                            Endd = cat(1,Endd,uint32((T_start+Sub_sampl*BCNT_Vis.PixelIdxList{1,i}(end,1))/srate));
                        end
                        Thr_ind = Thr_ind+1;
                        Average_Vel_tot_temp_temp{Thr_ind} = Average_Vel;
                        Strt_tot_temp_temp{Thr_ind} = Strt;
                        Endd_tot_temp_temp{Thr_ind} = Endd;                       

                        %% Performance calculation:
                        % Find all detection predictions for given 
                        % ID/Session/Part using parameters Thr1-4 
                        CSD_GT = CSD_GT_temp;
                        Events = Events_temp;
                        CSD_GT(end) = [];
                        Events(end) = [];
                        EEG_mask_tot = EEG_mask(:,1:size(EEG_mask,2)-1);
                        num_good_elec_craniectomy = (sum(EEG_mask_tot(crani_ind,:),1))';
                        CC_BW = squeeze(sum(sum(value,1),2))';
                        
                        num_good_elec = num_good_elec_craniectomy;
                        
                        Trial_boundary_indicators = zeros(size(BCNT));
                        delta_trial = 30*(srate*1)/Sub_sampl;
                        L_trial = 2*(srate*60)/Sub_sampl;
                        Trial_boundary_indicators(delta_trial:delta_trial:(size(BCNT,1)-L_trial)) = 1;
                        Trial_boundary_indicators(1) = 1;
                        Trial_boundary_ind = find(Trial_boundary_indicators);
                        Trial_boundary_ind(1) = 0;
                        T_start_tot_sub = zeros(1,size(Trial_boundary_ind,1));
                        T_end_tot_sub = zeros(1,size(Trial_boundary_ind,1));
                        BCNT_tot_sub = zeros(L_trial,(size(Trial_boundary_ind,1)));
                        CNT_tot_sub = zeros(L_trial,(size(Trial_boundary_ind,1)));
                        CSD_GT_tot_sub = zeros(L_trial,(size(Trial_boundary_ind,1)));
                        num_good_elec_tot_sub = zeros(L_trial,(size(Trial_boundary_ind,1)));
                        CSD_type_tot_sub = cell(1,size(Trial_boundary_ind,1));
                        for i_sub = 1:size(Trial_boundary_ind,1)
                            Strt = Trial_boundary_ind(i_sub)+1;
                            Endd = (Trial_boundary_ind(i_sub)+L_trial);
                            
                            T_start_sub = T_start+((i_sub-1)*delta_trial*Sub_sampl);
                            T_end_sub = (T_start_sub+(L_trial*Sub_sampl))+10*srate;
                            
                            T_start_tot_sub(i_sub) = T_start_sub;
                            T_end_tot_sub(i_sub) = T_end_sub;
                            BCNT_tot_sub(:,i_sub) = BCNT(Strt:Endd);
                            CNT_tot_sub(:,i_sub) = CNT_Vel(Strt:Endd);
                            CSD_GT_tot_sub(:,i_sub) = CSD_GT(Strt:Endd);
                            CSD_type_counter = find(CSD_GT);
                            CSD_type_ind_strt = min(find(CSD_type_counter>=Strt));
                            CSD_type_ind_Endd = max(find(CSD_type_counter<=Endd));
                            
                            if(~isempty(CSD_type_ind_strt) && ~isempty(CSD_type_ind_Endd))
                                if(CSD_type_ind_strt<=CSD_type_ind_Endd)
                                    CSD_type_tot_sub(i_sub) = {CSD_type(CSD_type_ind_strt:CSD_type_ind_Endd)};
                                end
                            end

                            num_good_elec_tot_sub(:,i_sub) = num_good_elec(Strt:Endd);
                        end

                        % Save all threshold-combination results
                        T_start_tot_temp_temp{Thr_ind} = T_start_tot_sub;
                        T_end_tot_temp_temp{Thr_ind} = T_end_tot_sub;
                        BCNT_tot_temp_temp{Thr_ind} = reshape(BCNT_tot_sub,[],1);
                        CNT_tot_temp_temp{Thr_ind} = reshape(CNT_tot_sub,[],1);
                        CSD_GT_tot_temp_temp{Thr_ind} = reshape(CSD_GT_tot_sub,[],1);
                        num_good_elec_tot_temp_temp{Thr_ind} = reshape(num_good_elec_tot_sub,[],1);
                        CSD_type_tot_temp_temp{Thr_ind} = CSD_type_tot_sub;

                    end
                end
            end
        end
        % Save each threshold-combination specific result for each Part
        Average_Vel_tot_temp{CSD_ind_orig} = Average_Vel_tot_temp_temp;
        Strt_tot_temp{CSD_ind_orig} = Strt_tot_temp_temp;
        Endd_tot_temp{CSD_ind_orig} = Endd_tot_temp_temp;
        T_start_tot_temp{CSD_ind_orig} = T_start_tot_temp_temp;
        T_end_tot_temp{CSD_ind_orig} = T_end_tot_temp_temp;
        BCNT_tot_temp{CSD_ind_orig} = BCNT_tot_temp_temp;
        CNT_tot_temp{CSD_ind_orig} = CNT_tot_temp_temp;
        CSD_GT_tot_temp{CSD_ind_orig} = CSD_GT_tot_temp_temp;
        num_good_elec_tot_temp{CSD_ind_orig} = num_good_elec_tot_temp_temp;
        CSD_type_tot_tempp{CSD_ind_orig} = CSD_type_tot_temp_temp;
        
    end

    % End of Part-specific loop

    currentFolder = pwd;
    cd ..
    %% Performance Metric Calculations
    TPR_num_TOTAL = {};
    TPR_denum_TOTAL = {};
    TNR_num_TOTAL = {};
    TNR_denum_TOTAL = {};
    CSD_event_quality_tot_TOTAL = {};

    Thr_ind = 0;
    CSD_GT_total_temp = CSD_GT_total{1};
    CSD_Labels_total_temp = CSD_Labels_total{1};

    % Loop back through thres-combinations for TPR/FNR calculations at each
    % Session of given Patient-ID
    for Thr_3_ind = 1:size(THR_3,2)
        for Thr_1_ind = 1:size(THR_1,2)
            for Thr_2_ind = 1:size(THR_2,2)
                TPR_num_TOTAL_temp = [];
                TPR_denum_TOTAL_temp = [];
                TNR_num_TOTAL_temp = [];
                TNR_denum_TOTAL_temp = [];
                CSD_event_quality_tot_TOTAL_temp = [];
                Endd_tot_TOTAL = {};
                Strt_tot_TOTAL = {};
                
                Thr_1 = THR_1(Thr_1_ind)
                Thr_2 = THR_2(Thr_2_ind)
                Thr_3 = THR_3(Thr_3_ind)
                
                for Thr_4_ind = 1:size(THR_4,2)
                    Endd_tot_TOTAL_temp = [];
                    Strt_tot_TOTAL_temp = [];
                    Average_Vel_tot_TOTAL_temp = [];

                    Thr_4 = THR_4(Thr_4_ind)
                    CSD_GT_total = CSD_GT_total_temp;
                    CSD_Labels_total = CSD_Labels_total_temp;
                    
                    BCNT_tot = {};
                    CNT_tot = {};
                    CSD_GT_tot = {};
                    T_start_tot = {};
                    T_end_tot = {};
                    CSD_type_tot = {};
                    num_good_elec_tot = {};
                    Strt_tot = {};
                    Endd_tot = {};
                    Average_Vel_tot = {};
                                
                    Thr_ind = Thr_ind+1;
                    for CSD_ind_orig = 1:Part_size
                        T_start_tot{CSD_ind_orig}  = T_start_tot_temp{CSD_ind_orig}{Thr_ind};
                        T_end_tot{CSD_ind_orig} = T_end_tot_temp{CSD_ind_orig}{Thr_ind};
                        BCNT_tot{CSD_ind_orig} = BCNT_tot_temp{CSD_ind_orig}{Thr_ind};
                        CNT_tot{CSD_ind_orig} = CNT_tot_temp{CSD_ind_orig}{Thr_ind};
                        CSD_GT_tot{CSD_ind_orig} = CSD_GT_tot_temp{CSD_ind_orig}{Thr_ind};
                        num_good_elec_tot{CSD_ind_orig} = num_good_elec_tot_temp{CSD_ind_orig}{Thr_ind};
                        CSD_type_tot{CSD_ind_orig} = CSD_type_tot_tempp{CSD_ind_orig}{Thr_ind};
                        Average_Vel_tot{CSD_ind_orig} = Average_Vel_tot_temp{CSD_ind_orig}{Thr_ind};
                        Strt_tot{CSD_ind_orig} = Strt_tot_temp{CSD_ind_orig}{Thr_ind};
                        Endd_tot{CSD_ind_orig} = Endd_tot_temp{CSD_ind_orig}{Thr_ind};
                    end
                    
                    TPR_num = 0;
                    TPR_denum = 0;
                    TNR_num = 0;
                    TNR_denum = 0;
                    Tot_detectable_CSD = 0;
                    Detected_CSD = 0;
                    CSD_event_quality_tot = [];              
                    
                    CSD_type_tot_temp = [];
                    for type_i = 1:size(CSD_type_tot,2)
                        tmp = CSD_type_tot{1,type_i};
                        for type_j = 1:size(tmp,2)
                            tmp_temp = tmp{1,type_j};
                            for type_l = 1:size(tmp_temp,1)
                                CSD_type_tot_temp = cat(2,CSD_type_tot_temp,tmp_temp(type_l));
                            end
                        end
                    end
                    CSD_type_tot = CSD_type_tot_temp;
                    
                    T_start_TOT = [];
                    T_end_TOT = [];
                    BCNT_TOT = [];
                    CNT_TOT = [];
                    CSD_GT_TOT = [];
                    num_good_elec_TOT = [];
                    for ind_CSD = 1:Part_size
                        T_start_TOT = cat(2,T_start_TOT,T_start_tot{ind_CSD});
                        T_end_TOT = cat(2,T_end_TOT,T_end_tot{ind_CSD});
                        BCNT_TOT = cat(1,BCNT_TOT,BCNT_tot{ind_CSD});
                        CNT_TOT = cat(1,CNT_TOT,CNT_tot{ind_CSD});
                        CSD_GT_TOT = cat(1,CSD_GT_TOT,CSD_GT_tot{ind_CSD});
                        num_good_elec_TOT = cat(1,num_good_elec_TOT,num_good_elec_tot{ind_CSD});
                    end
                    
                    T_start_tot = T_start_TOT;
                    T_end_tot = T_end_TOT;
                    BCNT_tot = BCNT_TOT;
                    CNT_tot = CNT_TOT;
                    CSD_GT_tot = CSD_GT_TOT;
                    num_good_elec_tot = num_good_elec_TOT;
                    srate = 64;
                    Sub_sampl = 30*srate;
                    
                    %% Calculation of correct alarm and detection rate:
                    
                    Trial_boundary_indicators = zeros(size(BCNT_tot));
                    delta_trial = size(Trial_boundary_indicators,1)/size(T_start_tot,2);
                    Trial_boundary_indicators(delta_trial:delta_trial:end) = 1;
                    Trial_boundary_indicators(1) = 1;
                    Trial_boundary_ind = find(Trial_boundary_indicators);
                    Trial_boundary_ind(1) = 0;
                    
                    
                    CSD_GT_tot_ind = find(CSD_GT_tot);
                    NO_CSD_events_ind = (contains(CSD_type_tot,'No CSD'));
                    CSD_GT_tot(CSD_GT_tot_ind(NO_CSD_events_ind)) = 0;
                    CSD_type_tot(NO_CSD_events_ind) = [];
                    
                    ISD_events_ind = (contains(CSD_type_tot,'_ISD'));
                    CSD_GT_tot_ind = find(CSD_GT_tot);
                    CSD_GT_tot_with_ISD = CSD_GT_tot;
                    CSD_GT_tot(CSD_GT_tot_ind(ISD_events_ind)) = 0;
                    CSD_type_tot(ISD_events_ind) = [];

                    num_good_elec_tot_Avg = num_good_elec_tot;
                    num_good_elec_tot(CSD_GT_tot_ind(ISD_events_ind)) = 0;
                    
                    
                    CSD_type_tot_unique =  unique(CSD_type_tot);
                    CSD_GT_tot_ind = find(CSD_GT_tot);
                    
                    if(~isempty(CSD_Labels_total))
                        Actual_ISD_events_ind = (contains(CSD_Labels_total,'_ISD'));
                        Actual_CSD_GT_ind_with_ISD = find(CSD_GT_total);
                        CSD_GT_total(Actual_CSD_GT_ind_with_ISD(Actual_ISD_events_ind)) = 0;
                    end
                    Actual_CSD_GT_ind = find(CSD_GT_total);
                    
                    BCNT_tot_ind = find(BCNT_tot);
                    CSD_blocks_size = size(Trial_boundary_ind,1)-1-size(NO_CSD_events_ind,2);
                    det_rate = zeros(size(Trial_boundary_ind,1)-1,1);
                    correct_alarms = zeros(size(Trial_boundary_ind,1)-1,1);
                    CNT_Vel_temp = -1*ones(size(Trial_boundary_ind,1)-1,1);
                    bad_event_ind = zeros(size(Trial_boundary_ind,1)-1,1);
                    true_negative = zeros(size(Trial_boundary_ind,1)-1,1);
                    CSD_event_quality = zeros(size(Trial_boundary_ind,1)-1,1);
                    WL = DELTA_T*(60*srate/Sub_sampl);%20min
                    WL_orig = DELTA_T*(60*srate);%20min
                    CSD_intervals = [];
                    CSD_Det_Neg_ind = zeros(size(Trial_boundary_ind,1)-1,1);
                    
                    
                    CSD_block_count = 0;
                    Det_block_count = 0;
                    correct_alarm = 0;
                    
                    % Compare detection times with GT CSD times 
                    parfor i=1:size(Trial_boundary_ind,1)-1
                        Strt = Trial_boundary_ind(i)+1;
                        Endd = Trial_boundary_ind(i+1);
                        det_temp_bad = find(mean(num_good_elec_tot(Strt:Endd))<0.5*(size(crani_ind,2)));
                        if(~isempty(det_temp_bad))
                            bad_event_ind(i) = 1;
                        end
                        CSD_event_quality(i) = mean(num_good_elec_tot(Strt:Endd));
                        CSD_block = ~isempty(find(CSD_GT_tot(Strt:Endd)));
                        Det_block = ~isempty(find(BCNT_tot(Strt:Endd)));
                        if(CSD_block)
                            CSD_block_count = CSD_block_count+1;
                            CSD_Det_Neg_ind(i) = 1;
                            First_CSD_ind = Actual_CSD_GT_ind(min(find(Actual_CSD_GT_ind>=(T_start_tot(i)-Sub_sampl))));
                            Last_CSD_ind = Actual_CSD_GT_ind(max(find(Actual_CSD_GT_ind<=(T_end_tot(i)+Sub_sampl))));
                            if(isempty(Last_CSD_ind))
                                Last_CSD_ind = First_CSD_ind;
                            end
                            if(isempty(First_CSD_ind))
                                First_CSD_ind = Last_CSD_ind;
                            end
                            i_prev = max(find(First_CSD_ind-WL_orig>T_end_tot))+1;
                            i_next = min(find(Last_CSD_ind+WL_orig<T_start_tot))-1;
                            if(~isempty(i_prev))
                                CSD_start = Trial_boundary_ind(i_prev)+1 + round((First_CSD_ind-WL_orig-T_start_tot(i_prev))/Sub_sampl);
                            else
                                CSD_start = Strt;
                            end
                            if(~isempty(i_next))
                                CSD_end = Trial_boundary_ind(i_next+1) - round((T_end_tot(i_next)-Last_CSD_ind-WL_orig)/Sub_sampl);
                            else
                                CSD_end = Endd;
                            end
                            det_temp = find(BCNT_tot(CSD_start:CSD_end));
                            if(~isempty(det_temp))
                                det_rate(i) = 1;
                            end
                        elseif(Det_block)
                            Det_block_count = Det_block_count+1;
                            CSD_Det_Neg_ind(i) = 2;
                            First_BCNT_ind = T_start_tot(i)+(BCNT_tot_ind(min(find(BCNT_tot_ind>=Strt)))-Strt)*Sub_sampl;
                            Last_BCNT_ind = T_end_tot(i)-(Endd-BCNT_tot_ind(max(find(BCNT_tot_ind<=Endd))))*Sub_sampl;
                            i_prev = find(First_BCNT_ind-WL_orig>T_end_tot, 1, 'last' )+1;
                            i_next = min(find(Last_BCNT_ind+WL_orig<T_start_tot))-1;
                            if(~isempty(i_prev))
                                BCNT_start = Trial_boundary_ind(i_prev)+1 + round((First_BCNT_ind-WL_orig-T_start_tot(i_prev))/Sub_sampl);
                            else
                                BCNT_start = Strt;
                            end
                            if(~isempty(i_next))
                                BCNT_end = Trial_boundary_ind(i_next+1) - round((T_end_tot(i_next)-Last_BCNT_ind-WL_orig)/Sub_sampl);
                            else
                                BCNT_end = Endd;
                            end
                            
                            correct_alarm_temp = find(CSD_GT_tot_with_ISD(BCNT_start:BCNT_end));
                            det_temp_bad = find(mean(num_good_elec_tot(Strt:Endd))<0.5*(size(crani_ind,2)));
                            if(~isempty(correct_alarm_temp))
                                correct_alarms(i) = 1;
                                Vel_temp = CNT_tot(Strt:Endd);
                                CNT_Vel_temp(i) = mean(Vel_temp(Vel_temp~=-1)); 
                            end
                            if(~isempty(det_temp_bad))
                                bad_event_ind(i) = 1;
                            end
                        else
                            CSD_Det_Neg_ind(i) = 3;
                            i_prev = max(find(T_start_tot(i)-WL_orig>T_end_tot))+1;
                            i_next = min(find(T_end_tot(i)+WL_orig<T_start_tot))-1;
                            if(~isempty(i_prev))
                                Check_start = Trial_boundary_ind(i_prev)+1;
                            else
                                Check_start = Strt;
                            end
                            if(~isempty(i_next))
                                Check_end = Trial_boundary_ind(i_next+1);
                            else
                                Check_end = Endd;
                            end
                            
                            false_negative_temp = find(CSD_GT_tot_with_ISD(Check_start:Check_end));
                            det_temp_bad = find(mean(num_good_elec_tot(Strt:Endd))<0.5*(size(crani_ind,2)));
                            if(isempty(false_negative_temp))
                                true_negative(i) = 1;
                            end
                            if(~isempty(det_temp_bad))
                                bad_event_ind(i) = 1;
                            end
                        end
                    end
                    
                    % Calculate summary stats of given Session+Thres-combn. 
                    bad_Ps = size(find(CSD_Det_Neg_ind==1 & bad_event_ind==1),1);
                    Total_Ns = (size(find(true_negative),1)+size(find(CSD_Det_Neg_ind==2),1)-size(find(correct_alarms),1));
                    CSD_Det_Neg_ind(bad_event_ind==1) = [];
                    det_rate(bad_event_ind==1) = [];
                    correct_alarms(bad_event_ind==1) = [];
                    CNT_Vel_temp(bad_event_ind==1) = [];
                    true_negative(bad_event_ind==1) = [];
                    CSD_event_quality(bad_event_ind==1) = [];
                    CSD_event_quality = [CSD_event_quality(CSD_Det_Neg_ind==1);CSD_event_quality(find(correct_alarms))];
                    
                    TPR = (size(find(det_rate),1))...
                        /(size(find(CSD_Det_Neg_ind==1),1));
                    TNR = (size(find(true_negative),1))...
                        /(size(find(true_negative),1)+size(find(CSD_Det_Neg_ind==2),1)-size(find(correct_alarms),1));
                    FPR = 1-TNR;
                    FNR = 1-TPR;
                    
                    TPR_num = TPR_num+(size(find(det_rate),1));
                    TPR_denum = TPR_denum+(size(find(CSD_Det_Neg_ind==1),1));
                    TNR_num = TNR_num+(size(find(true_negative),1));
                    TNR_denum = TNR_denum+(size(find(true_negative),1)+size(find(CSD_Det_Neg_ind==2),1)-size(find(correct_alarms),1));
                    bad_Ns = Total_Ns-(size(find(true_negative),1)+size(find(CSD_Det_Neg_ind==2),1)-size(find(correct_alarms),1));
                                       
                    CSD_event_quality_tot = [size(bad_event_ind,1),bad_Ns,bad_Ps];%mean(CSD_event_quality_tot);
                    
                    TPR_num_TOTAL_temp = cat(1,TPR_num_TOTAL_temp,uint32(TPR_num));
                    TPR_denum_TOTAL_temp = cat(1,TPR_denum_TOTAL_temp, uint32(TPR_denum));
                    TNR_num_TOTAL_temp = cat(1, TNR_num_TOTAL_temp, uint32(TNR_num));
                    TNR_denum_TOTAL_temp = cat(1,TNR_denum_TOTAL_temp, uint32(TNR_denum));
                    CSD_event_quality_tot_TOTAL_temp = cat(1, CSD_event_quality_tot_TOTAL_temp, uint32(CSD_event_quality_tot));
                    
                    CNT_Vel_temp(CNT_Vel_temp==-1) = []; 
                    
                    for Strt_ind=1:size(Strt_tot,2)
                        if(~isempty(Strt_tot(Strt_ind)))
                            Average_Vel = Average_Vel_tot(Strt_ind);
                            Strt = Strt_tot(Strt_ind);
                            Endd = Endd_tot(Strt_ind);
                            Average_Vel_tot_TOTAL_temp = cat(1,Average_Vel_tot_TOTAL_temp,Average_Vel{1,1});
                            Strt_tot_TOTAL_temp = cat(1,Strt_tot_TOTAL_temp,Strt{1,1});
                            Endd_tot_TOTAL_temp = cat(1,Endd_tot_TOTAL_temp,Endd{1,1});
                        end
                    end
                    Average_Vel_tot_TOTAL{Thr_4_ind} = CNT_Vel_temp;%Average_Vel_tot_TOTAL_temp;
                    Strt_tot_TOTAL{Thr_4_ind} = Strt_tot_TOTAL_temp;
                    Endd_tot_TOTAL{Thr_4_ind} = Endd_tot_TOTAL_temp;
                end
                
                % Save Session + Thres_1_2_3 combination specific files
                % (Multiple Thr_4 results are saved to same file)
                save([currentFolder,'\',Session_names(ss).name ,'_',sprintf('Thr_1_%.2d',Thr_1)...
                    ,'_',sprintf('Thr_2_%.2d',Thr_2),'_',sprintf('Thr_3_%.2d',Thr_3), '_detection_SDIII_240minOv180min_CCWL20_Delta_LowResVelFixedvTEST_VelocitySaveVf.mat'],...
                    'TPR_num_TOTAL_temp','TPR_denum_TOTAL_temp','TNR_num_TOTAL_temp','TNR_denum_TOTAL_temp','CSD_event_quality_tot_TOTAL_temp','Strt_tot_TOTAL','Endd_tot_TOTAL','Average_Vel_tot_TOTAL')
            end
        end
    end
end
