close all
clear all
clc
%%
% This code loads the raw dataset (ECG, PLETH, ECoG, and EEG signals), recorded using MOBERG Amp, converted using CNS Envision, loading/adding the SD and ICU events/annotations, and band-pass filtering:
%
% See: README.txt and [1] for more info.

% [1] Alireza Chamanzar, Jonathan Elmer, Lori Shutter, Jed Hartings, Pulkit Grover,
%  "Noninvasive and reliable automated detection of spreading depolarization in severe traumatic brain injury using scalp EEG",
%   Submitted to Nature Comms Med, 2022.

% Author: Alireza 	Date: 2021/12/19 12:00:48 	Revision: 0.1
% Copyright: This project is licensed - see the LICENSE.md file for details

% Documented/Commented by John McNamee
%% Preprocessing steps for EEG data:
eegpath = 'C:\Users\Alireza\Downloads\eeglab_current\eeglab2020_0'; % Initialize eeglab
current_path = pwd;
cd(eegpath)  % Jump to eeg path for starting eeglab 
eeglab  % Start eeglab 
cd(current_path) % Run code in Patient ID specific folder 
Patient_ID = current_path(end-6:end); % Extract patient ID from path
Session_names = dir('Patient*');  % Find list of Session in Patient ID Folder
Sessions_size = size(Session_names,1); 

for ss = 1:Sessions_size  % Loop through all found sessions for pateint ID
    cd(Session_names(ss).name) % Enter session specific folder
    current_path = pwd;  % Update current path to include Session extension
    EEG_Imp = dir('EEG,Composite,Impedance,*.mat')
    load(EEG_Imp(1).name) % Load Session Specific Impdence file (saved within session folder)
    % Assign all four loaded variables as 'Impedence variables'
    comp_elements_Imp = comp_elements;
    measurement_data_Imp = measurement_data;
    start_date_time_Imp = start_date_time;
    time_vector_Imp = time_vector;
    EEG_names = dir('EEG,Composite,SampleSeries*.mat'); % Locate all 'parts' of the recorded session
    Part_size = size(EEG_names,1);
    for pp = 1: Part_size  % Loop through all 'part' files which contain raw EEG data
        Name = [Patient_ID,'_',Session_names(ss).name,'_',EEG_names(pp).name(46:end-4)];
        load(EEG_names(pp).name)  % For each part 'pp' load the data
        data_qual_time_EEG = data_qual_time;  % Assign loaded varaibles as 'EEG Variables'
        data_qual_str_EEG = data_qual_str;
        for cc = 1:size(comp_elements,2)
            Q_events = zeros(size(time_vector));
            Q_events(ismember(time_vector,data_qual_time))=1000;
            fid =  fopen('Event.txt');  % Load 'Event.txt' which contains CSD Ground Truth information
            Event_time_tot = [];  % Initialize vectors for storing event data
            Event_label_tot = [];
            Events = textscan(fid,'%s','Delimiter',{'*'});  % Read events properly as inidividual events
            fclose(fid);
            % Newer eeglab versions may require 'MMM' rather than 'MMMM' as provided below
            %start_date_time = datetime(start_date_time,'Format','dd MMM yyyy, HH:mm:ss.SSS');  % Format start datatime into desired format
            start_date_time = datetime(start_date_time,'Format','dd MMMM yyyy, HH:mm:ss.SSS');  % Format start datatime into desired format
            for ee = 1:2:size(Events{1,1},1)  % This loop plots events which are read from previous file
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
        % Read in EDF for initalizaing 'EEG' data structure with desired fields
        EEG = pop_biosig('H:\SD-II\04-1167\04-1167\LabchartECG.edf', 'ref',[] ,'refoptions',{'keepref' 'on'});
        EEG.data = measurement_data; % Assign data loaded in line 41 to EEG structure
        EEG.times = time_vector*1000; % time in ms
        EEG.nbchan = size(measurement_data,1); % Measure number of channels
        EEG.comments = []; % Initialize comments vector
        EEG.srate = 256; % calculated by 1/(time_vector(2)-time_vector(1))
        EEG.xmin = time_vector(1);  % Find min and max time points
        EEG.xmax = time_vector(end);
        EEG.pnts = size(time_vector,2);  % Find total number of data points
        EEG.etc.T0 = datevec(start_date_time);  % Record start time 
        for i=1:size(EEG.data,1)  % Label channels names via info from comp elements
            EEG.chanlocs(i).labels = char(comp_elements(i));
            %    EEG.chanlocs(i).urchan = i;
        end
        EEG = eeg_checkset( EEG ); % eeglab call to ensure data is formatted correctly 
        % Load channel location information based on standard 10-5 configuration 
        EEG=pop_chanedit(EEG, 'lookup',[eegpath '/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp']); 
        
        % Seperate EEG and ECoG data into seperate variable structs
        ECoG = pop_select( EEG,'channel',{'ECoG1' 'ECoG2' 'ECoG3' 'ECoG4' 'ECoG5' 'ECoG6'});
        EEG = pop_select( EEG,'nochannel',{'ECoG1' 'ECoG2' 'ECoG3' 'ECoG4' 'ECoG5' 'ECoG6' 'ECoG7' 'ECoG8'});
        
        % Bandpass filter raw data between [0.001 - 50] Hz
        % NOTE: newer eeglabs may return error for below params
        % if given error modify to EEG = pop_eegfiltnew(EEG, 0.01,50,85000,0,[],1)
        EEG = pop_eegfiltnew(EEG, 0.001,50,85000,0,[],1);
        ECoG = pop_eegfiltnew(ECoG, 0.001,50,85000,0,[],1);
        % eeglab call for ensuring data is formatted correctly 
        ECoG = eeg_checkset( ECoG );
        EEG = eeg_checkset( EEG );
        %% visually check electrodes and interpolate if necessary. If all are fine, then run...
        EEG.setname= Name; % Assign appropriate names to both EEG and ECoG structs 
        EEG = eeg_checkset( EEG ); % eeglab format check 
        ECoG.setname= Name;
        ECoG = eeg_checkset( ECoG );
        % Save bandpass filtered EEG and ECoG data as seperate 'set' files
        EEG = pop_saveset( EEG, 'filename',[Name , '_EEG_filtered_withDC.set'],'filepath',current_path);
        ECoG = pop_saveset( ECoG, 'filename',[Name , '_ECoG_filtered_withDC.set'],'filepath',current_path);
        
        % STOP
        %% Loading ECG, PLETH, Resp, adding events (SDs, quality, etc.):
        % This large code block works to process additional Physio data
        ECG_names = dir('ECG,II,SampleSeries*.mat'); % Begin w ECG files
        if(~isempty(ECG_names)) % If any ECG files are present enter processing loop
            ECG_size = size(ECG_names,1); % Measure number of ECG files present 
            measurement_data_ECG = []; % Initialize vectors
            start_date_time_ECG = [];
            time_vector_ECG = [];
            data_qual_time_ECG = [];
            data_qual_str_ECG = [];
            for i = 1: ECG_size % Loop throug i number of ECG files
                load(ECG_names(i).name) % Load ECG file 'i'
                start_date_time_ECG = cat(1,start_date_time_ECG,start_date_time); % Save start data of file
                Delta_step = time_vector(2)-time_vector(1); % Calculate time step of file i 
                if(i>1)     % If this not our first file
                    % Find time diff across ECG files (looping for time overlap betweeen files)
                    Delta_tp = (time_vector(1)-time_vector_ECG(end)); 
                    % If large gap between new files and exisiting data fill in gap with extended time vec
                    extended_time_vector = (time_vector_ECG(end)... 
                        +Delta_step:Delta_step:time_vector(1)-Delta_step);
                else         % If loading first ECG file 
                    Delta_tp = (Delta_step); % Set delta tp as delta step
                    extended_time_vector = []; % Initialize extended time vector
                end
                if(Delta_tp<0) % If time overlap is detected only append new information
                    overlap_start = find(ismember(time_vector_ECG(end),time_vector));
                    % Find where time overlap ends, then add all vals beyond that index for time and data
                    time_vector_ECG = cat(2,time_vector_ECG,time_vector(overlap_start+1:end));
                    measurement_data_ECG =  cat(1,measurement_data_ECG,measurement_data(overlap_start+1:end,:));
                else % In cases where there delta_tp>0, simply append loaded time vector and measurement data
                    time_vector_ECG = cat(2,time_vector_ECG,cat(2,extended_time_vector,time_vector));
                    measurement_data_ECG =  cat(1,zeros(size(extended_time_vector,2),1),measurement_data_ECG,measurement_data);
                end
                % Append all data quality times and strings to Session-wide data variable
                data_qual_time_ECG = cat(2,data_qual_time_ECG,data_qual_time);
                data_qual_str_ECG = cat(2,data_qual_str_ECG,data_qual_str);
            end
            % Ensure data is properly time sorted by soritng times and then data in same order 
            [time_vector_ECG,ind] = sort(time_vector_ECG,2);
            measurement_data_ECG = measurement_data_ECG(ind);
            % Shift time vector by delta step
            time_vector_ECG_shifted = [time_vector_ECG(1)-Delta_step,time_vector_ECG(1:end-1)];
            % Check for data beyond expected ECG time vector
            ind = find((time_vector_ECG - time_vector_ECG_shifted)>2*Delta_step);
            for t=1:size(ind,2) % If 'out of range' data found
                time_vector_extension = time_vector_ECG(ind(t)-1):Delta_step:time_vector_ECG(ind(t));
                % Calculate how much to extend the data and time vectors by 
                time_vector_ECG = cat(2,time_vector_ECG(1:ind(t)-2),time_vector_extension,time_vector_ECG(ind(t)+1:end));
                measurement_data_ECG =  cat(1,measurement_data_ECG(1:ind(t)-2),zeros(size(time_vector_extension,2),1),measurement_data_ECG(ind(t)+1:end));
                ind = ind+size(time_vector_extension,2)-2; % Update value of ind
            end
        else % If no ECG files are avaliable for given Patient-ID/session 
             %set dummy vals for measurement data and time vectors
            measurement_data_ECG = EEG.data(1,:)' * 0;
            time_vector_ECG = EEG.times/1000;
        end
        % End of ECG processing loop
        
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






