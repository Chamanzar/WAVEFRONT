close all
clear all
clc
%%
% This code loads the raw dataset (ECG, PLETH, ECoG, and EEG signals), recorded using MOBERG Amp, converted using CNS Envision, loading/adding the SD and ICU events/annotations, and band-pass filtering:
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
current_path = pwd;
cd(eegpath)
eeglab
cd(current_path)
Patient_ID = current_path(end-6:end);
Session_names = dir('Patient*');
Sessions_size = size(Session_names,1);

for ss = 1:Sessions_size
    cd(Session_names(ss).name)
    current_path = pwd;
    EEG_Imp = dir('EEG,Composite,Impedance,*.mat')
    load(EEG_Imp(1).name)
    
    comp_elements_Imp = comp_elements;
    measurement_data_Imp = measurement_data;
    start_date_time_Imp = start_date_time;
    time_vector_Imp = time_vector;
    EEG_names = dir('EEG,Composite,SampleSeries*.mat');
    Part_size = size(EEG_names,1);
    for pp = 1: Part_size
        Name = [Patient_ID,'_',Session_names(ss).name,'_',EEG_names(pp).name(46:end-4)];
        load(EEG_names(pp).name)
        data_qual_time_EEG = data_qual_time;
        data_qual_str_EEG = data_qual_str;
        for cc = 1:size(comp_elements,2)
            Q_events = zeros(size(time_vector));
            Q_events(ismember(time_vector,data_qual_time))=1000;
            fid =  fopen('Event.txt');
            Event_time_tot = [];
            Event_label_tot = [];
            Events = textscan(fid,'%s','Delimiter',{'*'});
            fclose(fid);
            start_date_time = datetime(start_date_time,'Format','dd MMMM yyyy, HH:mm:ss.SSS');
            for ee = 1:2:size(Events{1,1},1)
                Event_label_tot = cat(1,Event_label_tot,Events{1,1}(ee));
                Events_stamps_datetime = datetime(Events{1,1}(ee+1),'Format','yyyy-MM-dd HH:mm:ss.SSS');
                Event_time = (etime(datevec(Events_stamps_datetime),datevec(start_date_time)));
                Event_time_tot = cat(1,Event_time_tot,Event_time);
                T_events = zeros(size(time_vector_Imp));
                [~,ind_min] = min(abs(time_vector_Imp-Event_time));
                T_events(ind_min)=1000;
                               hold on
                               plot(time_vector_Imp,T_events,'y')
                               text(time_vector_Imp(ind_min),1000,Events{1,1}(ee))
            end
        end
        hold off
        
        
        close all
        %% Read in bdf file, average data to mastoids (electrodes 129 130), bandpass filter 1-100Hz, and assign channel names
        % read in edf
        EEG = pop_biosig('H:\SD-II\04-1167\04-1167\LabchartECG.edf', 'ref',[] ,'refoptions',{'keepref' 'on'});
        EEG.data = measurement_data;
        EEG.times = time_vector*1000;%time in ms
        EEG.nbchan = size(measurement_data,1);
        EEG.comments = [];
        EEG.srate = 256;%1/(time_vector(2)-time_vector(1));
        EEG.xmin = time_vector(1);
        EEG.xmax = time_vector(end);
        EEG.pnts = size(time_vector,2);
        EEG.etc.T0 = datevec(start_date_time);
        for i=1:size(EEG.data,1)
            EEG.chanlocs(i).labels = char(comp_elements(i));
            %    EEG.chanlocs(i).urchan = i;
        end
        EEG = eeg_checkset( EEG );
        EEG=pop_chanedit(EEG, 'lookup',[eegpath '/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp']); 
        
        ECoG = pop_select( EEG,'channel',{'ECoG1' 'ECoG2' 'ECoG3' 'ECoG4' 'ECoG5' 'ECoG6'});
        EEG = pop_select( EEG,'nochannel',{'ECoG1' 'ECoG2' 'ECoG3' 'ECoG4' 'ECoG5' 'ECoG6' 'ECoG7' 'ECoG8'});
        
        
        % filter
        EEG = pop_eegfiltnew(EEG, 0.01,50,85000,0,[],1);
        ECoG = pop_eegfiltnew(ECoG, 0.01,50,85000,0,[],1);
    
        ECoG = eeg_checkset( ECoG );
        EEG = eeg_checkset( EEG );
        %% visually check electrodes and interpolate if necessary. If all are fine, then run...
        
        EEG.setname= Name;
        EEG = eeg_checkset( EEG );
        ECoG.setname= Name;
        ECoG = eeg_checkset( ECoG );
   
        EEG = pop_saveset( EEG, 'filename',[Name , '_EEG_filtered_withDC.set'],'filepath',current_path);
        ECoG = pop_saveset( ECoG, 'filename',[Name , '_ECoG_filtered_withDC.set'],'filepath',current_path);
        
        % STOP
        %% Loading ECG, PLETH, Resp, adding events (SDs, quality, etc.):
        
        ECG_names = dir('ECG,II,SampleSeries*.mat');
        if(~isempty(ECG_names))
            ECG_size = size(ECG_names,1);
            measurement_data_ECG = [];
            start_date_time_ECG = [];
            time_vector_ECG = [];
            data_qual_time_ECG = [];
            data_qual_str_ECG = [];
            for i = 1: ECG_size
                load(ECG_names(i).name)
                start_date_time_ECG = cat(1,start_date_time_ECG,start_date_time);
                Delta_step = time_vector(2)-time_vector(1);
                if(i>1)
                    Delta_tp = (time_vector(1)-time_vector_ECG(end));
                    extended_time_vector = (time_vector_ECG(end)...
                        +Delta_step:Delta_step:time_vector(1)-Delta_step);
                else
                    Delta_tp = (Delta_step);
                    extended_time_vector = [];
                end
                if(Delta_tp<0)
                    overlap_start = find(ismember(time_vector_ECG(end),time_vector));
                    time_vector_ECG = cat(2,time_vector_ECG,time_vector(overlap_start+1:end));
                    measurement_data_ECG =  cat(1,measurement_data_ECG,measurement_data(overlap_start+1:end,:));
                else
                    time_vector_ECG = cat(2,time_vector_ECG,cat(2,extended_time_vector,time_vector));
                    measurement_data_ECG =  cat(1,zeros(size(extended_time_vector,2),1),measurement_data_ECG,measurement_data);
                end
                
                data_qual_time_ECG = cat(2,data_qual_time_ECG,data_qual_time);
                data_qual_str_ECG = cat(2,data_qual_str_ECG,data_qual_str);
            end
            
            [time_vector_ECG,ind] = sort(time_vector_ECG,2);
            measurement_data_ECG = measurement_data_ECG(ind);
            
            time_vector_ECG_shifted = [time_vector_ECG(1)-Delta_step,time_vector_ECG(1:end-1)];
            ind = find((time_vector_ECG - time_vector_ECG_shifted)>2*Delta_step);
            for t=1:size(ind,2)
                time_vector_extension = time_vector_ECG(ind(t)-1):Delta_step:time_vector_ECG(ind(t));
                time_vector_ECG = cat(2,time_vector_ECG(1:ind(t)-2),time_vector_extension,time_vector_ECG(ind(t)+1:end));
                measurement_data_ECG =  cat(1,measurement_data_ECG(1:ind(t)-2),zeros(size(time_vector_extension,2),1),measurement_data_ECG(ind(t)+1:end));
                ind = ind+size(time_vector_extension,2)-2;
            end
        else
            measurement_data_ECG = EEG.data(1,:)' * 0;
            time_vector_ECG = EEG.times/1000;
        end
        
        PLETH_names = dir('PLETH,na,SampleSeries,Integer,IntelliVue,data*.mat');
        if(~isempty(PLETH_names))
            PLETH_size = size(PLETH_names,1);
            measurement_data_PLETH = [];
            start_date_time_PLETH = [];
            time_vector_PLETH = [];
            data_qual_time_PLETH = [];
            data_qual_str_PLETH = [];
            for i = 1: PLETH_size
                load(PLETH_names(i).name)
                start_date_time_PLETH = cat(1,start_date_time_PLETH,start_date_time);
                Delta_step = time_vector(2)-time_vector(1);
                if(i>1)
                    Delta_tp = (time_vector(1)-time_vector_PLETH(end));
                    extended_time_vector = (time_vector_PLETH(end)...
                        +Delta_step:Delta_step:time_vector(1)-Delta_step);
                else
                    Delta_tp = (Delta_step);
                    extended_time_vector = [];
                end
                if(Delta_tp<0)
                    overlap_start = find(ismember(time_vector_PLETH(end),time_vector));
                    time_vector_PLETH = cat(2,time_vector_PLETH,time_vector(overlap_start+1:end));
                    measurement_data_PLETH =  cat(1,measurement_data_PLETH,measurement_data(overlap_start+1:end,:));
                else
                    time_vector_PLETH = cat(2,time_vector_PLETH,cat(2,extended_time_vector,time_vector));
                    measurement_data_PLETH =  cat(1,zeros(size(extended_time_vector,2),1),measurement_data_PLETH,measurement_data);
                end
                
                data_qual_time_PLETH = cat(2,data_qual_time_PLETH,data_qual_time);
                data_qual_str_PLETH = cat(2,data_qual_str_PLETH,data_qual_str);
            end
            
            [time_vector_PLETH,ind] = sort(time_vector_PLETH,2);
            measurement_data_PLETH = measurement_data_PLETH(ind);
            
            time_vector_PLETH_shifted = [time_vector_PLETH(1)-Delta_step,time_vector_PLETH(1:end-1)];
            ind = find((time_vector_PLETH - time_vector_PLETH_shifted)>2*Delta_step);
            for t=1:size(ind,2)
                time_vector_extension = time_vector_PLETH(ind(t)-1):Delta_step:time_vector_PLETH(ind(t));
                time_vector_PLETH = cat(2,time_vector_PLETH(1:ind(t)-2),time_vector_extension,time_vector_PLETH(ind(t)+1:end));
                measurement_data_PLETH =  cat(1,measurement_data_PLETH(1:ind(t)-2),zeros(size(time_vector_extension,2),1),measurement_data_PLETH(ind(t)+1:end));
                ind = ind+size(time_vector_extension,2)-2;
            end
        else
            measurement_data_PLETH = EEG.data(1,:)' * 0;
            time_vector_PLETH = EEG.times/1000;
        end
        
        
        
        RESP_names = dir('RESP,na,SampleSeries,Integer,IntelliVue,data*.mat');
        if(~isempty(RESP_names))
            RESP_size = size(RESP_names,1);
            measurement_data_RESP = [];
            start_date_time_RESP = [];
            time_vector_RESP = [];
            data_qual_time_RESP = [];
            data_qual_str_RESP = [];
            for i = 1: RESP_size
                load(RESP_names(i).name)
                start_date_time_RESP = cat(1,start_date_time_RESP,start_date_time);
                Delta_step = time_vector(2)-time_vector(1);
                if(i>1)
                    Delta_tp = (time_vector(1)-time_vector_RESP(end));
                    extended_time_vector = (time_vector_RESP(end)...
                        +Delta_step:Delta_step:time_vector(1)-Delta_step);
                else
                    Delta_tp = (Delta_step);
                    extended_time_vector = [];
                end
                if(Delta_tp<0)
                    overlap_start = find(ismember(time_vector_RESP(end),time_vector));
                    time_vector_RESP = cat(2,time_vector_RESP,time_vector(overlap_start+1:end));
                    measurement_data_RESP =  cat(1,measurement_data_RESP,measurement_data(overlap_start+1:end,:));
                else
                    time_vector_RESP = cat(2,time_vector_RESP,cat(2,extended_time_vector,time_vector));
                    measurement_data_RESP =  cat(1,zeros(size(extended_time_vector,2),1),measurement_data_RESP,measurement_data);
                end
                
                data_qual_time_RESP = cat(2,data_qual_time_RESP,data_qual_time);
                data_qual_str_RESP = cat(2,data_qual_str_RESP,data_qual_str);
            end
            
            
            [time_vector_RESP,ind] = sort(time_vector_RESP,2);
            measurement_data_RESP = measurement_data_RESP(ind);
            
            time_vector_RESP_shifted = [time_vector_RESP(1)-Delta_step,time_vector_RESP(1:end-1)];
            ind = find((time_vector_RESP - time_vector_RESP_shifted)>2*Delta_step);
            for t=1:size(ind,2)
                time_vector_extension = time_vector_RESP(ind(t)-1):Delta_step:time_vector_RESP(ind(t));
                time_vector_RESP = cat(2,time_vector_RESP(1:ind(t)-2),time_vector_extension,time_vector_RESP(ind(t)+1:end));
                measurement_data_RESP =  cat(1,measurement_data_RESP(1:ind(t)-2),zeros(size(time_vector_extension,2),1),measurement_data_RESP(ind(t)+1:end));
                ind = ind+size(time_vector_extension,2)-2;
            end
        else
            measurement_data_RESP = EEG.data(1,:)' * 0;
            time_vector_RESP = EEG.times/1000;
        end

        if(~isempty(ECG_names))
            measurement_data_ECG = resample(measurement_data_ECG,256,500);
        end
        if(~isempty(PLETH_names))
            measurement_data_PLETH = resample(measurement_data_PLETH,256,125);
        end
        if(~isempty(RESP_names))
            measurement_data_RESP = resample(measurement_data_RESP,2560,625);
        end
        
        Delta_T = 1/256;
        Delta_To = ((EEG.times(1)/1000)-time_vector_ECG(1));
        Delta_To = round(Delta_To / Delta_T);
        if(Delta_To>=0)
            measurement_data_ECG = measurement_data_ECG(Delta_To+1:end);
        else
            measurement_data_ECG = cat(1,zeros(abs(Delta_To),1),measurement_data_ECG);
        end
        
        if(size(measurement_data_ECG,1)>=EEG.pnts)
            measurement_data_ECG = measurement_data_ECG(1:EEG.pnts);
        else
            measurement_data_ECG = cat(1,measurement_data_ECG,zeros((EEG.pnts-size(measurement_data_ECG,1)),1));
        end
        
        Delta_To = ((EEG.times(1)/1000)-time_vector_PLETH(1));
        Delta_To = round(Delta_To / Delta_T);
        if(Delta_To>=0)
            measurement_data_PLETH = measurement_data_PLETH(Delta_To+1:end);
        else
            measurement_data_PLETH = cat(1,zeros(abs(Delta_To),1),measurement_data_PLETH);
        end
        if(size(measurement_data_PLETH,1)>=EEG.pnts)
            measurement_data_PLETH = measurement_data_PLETH(1:EEG.pnts);
        else
            measurement_data_PLETH = cat(1,measurement_data_PLETH,zeros((EEG.pnts-size(measurement_data_PLETH,1)),1));
        end
        
        Delta_To = ((EEG.times(1)/1000)-time_vector_RESP(1));
        Delta_To = round(Delta_To / Delta_T);
        if(Delta_To>=0)
            measurement_data_RESP = measurement_data_RESP(Delta_To+1:end);
        else
            measurement_data_RESP = cat(1,zeros(abs(Delta_To),1),measurement_data_RESP);
        end
        if(size(measurement_data_RESP,1)>=EEG.pnts)
            measurement_data_RESP = measurement_data_RESP(1:EEG.pnts);
        else
            measurement_data_RESP = cat(1,measurement_data_RESP,zeros((EEG.pnts-size(measurement_data_RESP,1)),1));
        end
        
        % Add the ECG, PLETH, and RESP channels to EEG:
        EEG.data(end+1,:) = measurement_data_ECG;
        EEG.data(end+1,:) = measurement_data_PLETH;
        EEG.data(end+1,:) = measurement_data_RESP;
        EEG.nbchan = size(EEG.data,1);
        if ~isempty(EEG.chanlocs)
            EEG.chanlocs(end+1).labels = 'ECG';
            EEG.chanlocs(end+1).labels = 'PLETH';
            EEG.chanlocs(end+1).labels = 'RESP';
        end
        
        EEG_pruned = pop_select(EEG,'point',[2,size(EEG.data,2)]);
        EEG.event = EEG_pruned.event;
        EEG.event.type = 'Dummy';
        EEG.event.latency = 0;
        EEG.event.duration = 0;
        start_date_time = datetime(start_date_time,'Format','dd MMMM yyyy, HH:mm:ss.SSS');
        
        %Insert quality and other events:
        for j=1%:size(data_qual_str_EEG,1)
            for i=1:size(data_qual_str_EEG(j,:),2)
                if(data_qual_time_EEG(i)>=EEG.xmin)
                    EEG = pop_editeventvals(EEG,'insert',{1,[],[],[]},'changefield',{1,'type',char(data_qual_str_EEG(j,i))},'changefield',{1,'latency',data_qual_time_EEG(i)},'changefield',{1,'duration',1});
                end
            end
        end
        
        for i=1:size(Event_time_tot,1)
            if(Event_time_tot(i)>=EEG.xmin)
                EEG = pop_editeventvals(EEG,'insert',{1,[],[],[]},'changefield',{1,'type',char(Event_label_tot(i))},'changefield',{1,'latency',Event_time_tot(i)},'changefield',{1,'duration',1});
            end
        end
        
        
        %Insert CSD events (15min duration assumed for CSD episode):
        T = readtable('Pitt ECoG SD-II.xls');
        T = T(9:end,:);
        CSD_ind = find(start_date_time<T.Date_Time);
        T = T(CSD_ind,:);
        CSD_stamps = char(between(start_date_time,T.Date_Time));
        for i=1:size(T,1)
            Delta_To = CSD_stamps(i,:);
            m_ind = strfind(Delta_To,'m');
            h_uind = strfind(Delta_To,'h');
            d_ind = strfind(Delta_To,'d');
            space_ind = strfind(Delta_To,' ');
            h_lind = space_ind(find(h_uind>space_ind, 1, 'last' ));
            if(isempty(d_ind))
                
                dt_i = str2double(Delta_To(m_ind+2:end-1)) + str2double(Delta_To(h_uind+2:m_ind-1))*60 + ...
                    str2double(Delta_To(h_lind+1:h_uind-1))*60*60;
            else
                dt_i = str2double(Delta_To(m_ind+2:end-1)) + str2double(Delta_To(h_uind+2:m_ind-1))*60 + ...
                    str2double(Delta_To(h_lind+1:h_uind-1))*60*60 + str2double(Delta_To(d_ind-1))*24*60*60;
            end
            if(dt_i>=EEG.xmin)
                EEG = pop_editeventvals(EEG,'insert',{1,[],[],[]},'changefield',{1,'type',char(T.Comment(i))},'changefield',{1,'latency',dt_i},'changefield',{1,'duration',900});
            end
            
        end
        
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        EEG_Only = pop_select( EEG,'nochannel',{'ECG' 'RESP' 'PLETH'});
        pop_eegplot( EEG_Only, 1, 1, 1);
        
        All_Physio = pop_select( EEG,'channel',{'ECG' 'RESP' 'PLETH'});
        
        All_Physio = pop_eegfiltnew(All_Physio, 58,62,1000,1,[],1);
        EEG.data(end-2:end,:) = All_Physio.data;
        
        EEG.data(end+1,:) = mean(EEG.data(8:9,:),1);
        EEG.data(end+1,:) = EEG.data(8,:) - EEG.data(9,:);
        EEG.nbchan = size(EEG.data,1);
        if ~isempty(EEG.chanlocs)
            EEG.chanlocs(end+1).labels = 'EOGv';
            EEG.chanlocs(end+1).labels = 'EOGh';
        end
        
        EOG = pop_select( EEG,'channel',{'EOGh' 'EOGv'});
        EOG = pop_eegfiltnew(EOG, 0.1,20,10000,0,[],1);
        EEG.data(end-1:end,:) = EOG.data;
        
        EEG = eeg_checkset( EEG );
        EEG = pop_saveset( EEG, 'filename',[Name , '_preica_withDC.set'],'filepath',current_path);
        close all
        
        
    end
    
    cd ..
end





