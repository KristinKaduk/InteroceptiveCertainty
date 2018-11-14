function Data2trials_ECG()
%% Kristin start:03072017

%To DO: Weich-Kodieren - z.B. Name PARAMETER.source_file_s(5:end-4)]
clear all;
clc;
%% PARAMETER
addToExistingTable = 1;
Computer           = 'PC308';  %PC on which you are working
ExtractionPeriod   = 'HBC_Period';  %Baseline_Period  %in which TimePeriod the HPTs are detected
NameStudy1          = '3ndDataCollection_30trials_CertaintyScale';  %2ndDataCollection_30trials %1stDataCollection_5trials
NameStudy2          = 'HPT_certaintyscale'; %HPT_certaintyscale  2nd_HPT_30Trials
%% PATH
% CHANGE PATH according to the location where you  stored the data
switch Computer
    case 'PC308'
       path_HPT_eeg = ['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Interoceptive_Certainty\Data\ECG\HPT\raw\NewData'];
    case 'PC202' %Anna's computer            
        path_HPT_eeg = ['D:\Students\Anna\HeartbeatCountingTask\Data\', NameStudy1, '\ECG\', NameStudy2, '\raw\NewData'];
end

cd(path_HPT_eeg) % "Current Folder" ändert sich oben links automatisch zu dem angegebenen Pfad
Files = dir('*.eeg'); %Current Folder - zeige alle .eeg Daten
AllData_Combined = [];% initiate a variable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% PREPROCESSING: EACH SUBJECT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% USING Fieldtrip %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subject: ID as indSubject = subject you want to analyse
    indSubject =1; %:size(Files,1)
    cd(path_HPT_eeg)
    PARAMETER.stim_not_clear = 'yes'; %'no'- trigger protocol available
    PARAMETER.source_file_s  = Files(indSubject).name; %output is e.g. HBT_RIAS03.eeg
    PARAMETER.filter_s       = 'yes';
    PARAMETER.detrend_s      = 'yes';


    % 1. define the triggers for cutting the data into trials;
    %reading and segmenting the data such that the experimental conditions are represented as trials in a
    %data strucure. % to read: http://www.fieldtriptoolbox.org/tutorial/preprocessing
    cfg = [];
    cfg.dataset = PARAMETER.source_file_s;
    cfg.trialdef.eventtype = 'Stimulus'; %Trigger
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Extract the HBC_period OR extract the BaselinePeriod    %%  %%  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    switch ExtractionPeriod
        case 'HBC_Period'
           % extract the HBC_period
            cfg.trialdef.eventvalue_Start = 'S 41';
            cfg.trialdef.eventvalue_End   = 'S 42';
            cfg.trialfun = 'ft_trialfun_extractBetweenTwoTriggers'; % searching between two triggers
        case 'Baseline_Period'
            % extract the BaselinePeriod
            cfg.trialfun = 'KK_ft_trialfun_extractBetweenTwoTriggers_OneTrial'; %
            cfg.trialdef.eventvalue_Start = 'S 2';
            cfg.trialdef.eventvalue_End   = 'S 3';
        case 'NoTask'
             PARAMETER.stim_not_clear = 'no'; %'no'-> trigger protocol available
             cfg                    = [];
             cfg.dataset            = PARAMETER.source_file_s;
             cfg.trialdef.eventtype = '?';  
        end
    %% process
    cfg = ft_definetrial(cfg);
  %OUTPUT to control:found 660 events; created 30 trials  
    %% 2. Read the data and separate them into trials around the trigger of interest
    % http://www.fieldtriptoolbox.org/reference/ft_preprocessing?s[]=cfg&s[]=bsfreq
    cfg.continuous = 'yes';
    cfg.demean     = 'no'; %whether to apply baseline correction (default = 'no')
    %cfg.padding    = 10;   %defines the duration to which the data in the trial will be padded:
                            %not only the segment of data that corresponds to the trial is read from the file,
                            %but also some extra data around it. After filtering this,
                            %padding is removed from the data and preprocessing only returns the segment of interest
    cfg.bsfilter   = PARAMETER.filter_s;
    cfg.bsfreq     = [48 52]; %  bandstop frequency range, specified as [low high] in Hz
    cfg.detrend    = PARAMETER.detrend_s;
    cfg.lpfilter   = 'yes';% or 'no'  lowpass filter (default = 'no')
    cfg.lpfreq     = 150;
    %cfg.hpfreq     = 0.8;
    %cfg.absdiff    = 'yes';
    cfg.channel    = {'ECG'}; %{ 'ECG' 'Ear_L' 'Eye_L' }
    cfg.feedback   = 'no';
    data_tr        = ft_preprocessing(cfg);

    %% 3. check data
           da = [];
           [~] = ft_databrowser(da, data_tr);
    
    %% 4. save data
    
    switch Computer
        case 'PC308'
            path_HPT_Segmented = ['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Interoceptive_Certainty\Data\ECG\HPT\segmented'];
        case 'PC202' %Anna's computer
            path_HPT_Segmented = ['D:\Students\Anna\HeartbeatCountingTask\Data\', NameStudy1, '\ECG\', NameStudy2, '\segmented'];
    end
    
    
    cd(path_HPT_Segmented);
    switch ExtractionPeriod
        case 'Baseline_Period'
            name = [PARAMETER.source_file_s(1:3), PARAMETER.source_file_s(4:end-4), '_Segm_Baseline_Period'];
            name2 = ['Data_',name(5:end)]; 

        case 'HBC_Period'
            name = [PARAMETER.source_file_s(1:3), PARAMETER.source_file_s(4:end-4), '_Segm'];
            name2 = ['Data_',name(5:end)]; 
    end
    
    save(name,'data_tr');       %save preprocessed and segmented file -> will be changed for the artefact detection
    save(name2,'data_tr');      %save preprocessed and segmented file 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% DETECT R-Peak (Heartbeat) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = data_tr;
Da   = [];
Data = [];
%% detect R-peaks AUTOMATICALLY for each trial
switch ExtractionPeriod
    case 'Baseline_Period'
        data.sampleinfo = data_tr.sampleinfo;
        data.trial = data_tr.trial;
        data.time = data_tr.time;
        cfg = [];
        cfg.trl  = data_tr.cfg.trl;
        cfg.continuous = 'yes';
        cfg.artfctdef.ecg.method  = 'zvalue'; %peak-detection method
        cfg.artfctdef.ecg.cutoff  = 3;% peak-threshold
        cfg.artfctdef.ecg.pretim  = 0.05; %pre-artifact rejection-interval in seconds
        cfg.artfctdef.ecg.psttim  = 0.3; %post-artifact rejection-interval in seconds
        cfg.artfctdef.ecg.inspect  = 'ECG';% peak-threshold
        [cfg, artifact] = KK_Baseline_ft_artifact_ecg(cfg,data);
        Data.automaticDetected_beats_BaselinePeriod(1)          = size(artifact, 1);
        Data.beats_BaselinePeriod(1)                            = size(artifact, 1);
        Data.duration(1)       = length(data_tr.time{1});
        Data.eventStart(1)     = data_tr.sampleinfo(1);
        Data.eventFinsish(1)   = data_tr.sampleinfo(end);
        
    case 'HBC_Period'
        for indTrials =1: size(data_tr.trial,2)
            data.sampleinfo = data_tr.sampleinfo(indTrials,:);
            data.trial = data_tr.trial{1,indTrials};
            data.time = data_tr.time(indTrials);
            % medfilt1(signal, N) -> N needs to be found experimentally
            data.trial = data.trial - medfilt1(data.trial,data.fsample); %FILTER!!!
           % figure(2); plot(data.trial);figure(2);
            data.trial = [{data.trial}]; 

            cfg = [];
            cfg.trl  = data_tr.cfg.trl(indTrials,:);
            cfg.continuous = 'no';
            cfg.artfctdef.ecg.method  = 'zvalue'; %peak-detection method
            cfg.artfctdef.ecg.cutoff  = 3;% peak-threshold
            cfg.artfctdef.ecg.pretim  = 0.05; %pre-artifact rejection-interval in seconds
            cfg.artfctdef.ecg.psttim  = 0.3; %post-artifact rejection-interval in seconds
            cfg.artfctdef.ecg.inspect  = 'ECG';% peak-threshold
           %% SPECIAL CASE: In two dataset, one trial is not working!!! -> a work around it
%             if strmatch(PARAMETER.source_file_s(4:end-4), 'MIYE17') && indTrials == 18
%             Da.automaticDetected_beats(indTrials,1)          = 0;
%             Da.beats(indTrials,1)          = 0;
%             Da.duration(indTrials,1)       = data_tr.time{1,indTrials}(end) - data_tr.time{1,indTrials}(1);
%             Da.eventStart(indTrials,1)     = data_tr.time{1,indTrials}(1);
%             Da.eventFinsish(indTrials,1)   = data_tr.time{1,indTrials}(end);      
%                 
%             else
               % !!! automatic heartbeat detection
            [cfg, artifact] = KK_ft_artifact_ecg_InputRemovedLowFreq(cfg,data);
            Da.automaticDetected_beats(indTrials,1)          = size(artifact, 1);
            Da.beats(indTrials,1)          = size(artifact, 1);
            Da.duration(indTrials,1)       = data_tr.time{1,indTrials}(end) - data_tr.time{1,indTrials}(1);
            Da.eventStart(indTrials,1)     = data_tr.time{1,indTrials}(1);
            Da.eventFinsish(indTrials,1)   = data_tr.time{1,indTrials}(end);
%              end
        end
        
end


%CLOSE all graph: otherwise screw ups later on while running the pipeline

%% MANUAL CHECK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Manual check of the automatic Heartbeatdetection %%%%%%%%%%%%

%only for HBC_Period nicht baseline: go trial by trial through the data
%!! change the number of trials !!

indTrials =15;   %: size(data_tr.trial,2)
data.sampleinfo = data_tr.sampleinfo(indTrials,:);
data.trial = data_tr.trial{1,indTrials};
data.trial = data.trial - medfilt1(data.trial,data.fsample); %FILTER!!!
data.trial = [{data.trial}]; 
data.time = data_tr.time(indTrials);

cfg = [];
cfg.trl  = data_tr.cfg.trl(indTrials,:);
cfg.continuous = 'no';
cfg.artfctdef.ecg.method  = 'zvalue';   %peak-detection method
cfg.artfctdef.ecg.cutoff  = 3;          %peak-threshold
cfg.artfctdef.ecg.pretim  = 0.05;       %pre-artifact rejection-interval in seconds
cfg.artfctdef.ecg.psttim  = 0.3;        %post-artifact rejection-interval in seconds
cfg.artfctdef.ecg.inspect  = 'ECG';     %peak-threshold

[cfg, artifact] = KK_ft_artifact_ecg(cfg,data);

%%     %% check manually the detected R-peaks -> each peak should be marked in RED
% If an R-peak is not automatically detected, add it manually:
% Without any tool selected -> click right side on ECG ->  mark the area of
% the peak and click inside the area

cfg = [];
cfg.datafile   = PARAMETER.source_file_s;
cfg.artfctdef.ECG.artifact= artifact;
cfg.headerfile = PARAMETER.source_file_s;
cfg.continuous = 'no';
cfg.channel    = {'all'};
cfg.viewmode   = 'vertical';

%% !! Schau ob alle Herzschlaege automatisch erkannt wurden
cfg_manualArtefact = ft_databrowser(cfg,data);

%% !! HIER ueberspeicherst du die automatisch erfassten Herzschlaege
Da.beats(indTrials,1)          = size(cfg_manualArtefact.artfctdef.ECG.artifact, 1); %fuer HBC_Period
save(['DataStrc_HBCperiod', PARAMETER.source_file_s(4:end-4)],'Da');


%nach allen Trials einer Person hier weiter machen%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% transform the structur with the information about heartbeat into behavior data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Table_HPT = struct2table(Da);
Table_HPT %for direct visualization

%Table_BaselinePeriod = struct2table(Data);
%Table_CT = struct2table(CT);

%load( [path_HPT_Segmented, '\DataStrc_HBCperiod_', PARAMETER.source_file_s(5:end-4) , '.mat']);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Behavior Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch Computer
        case 'PC308'
            path_HPT_Segmented = ['C:\Users\kkaduk\Desktop\Kristin\Projects\Interception_Perception_Metacognition\Pilotstudy_HeartrateCountingTask\Data\3ndDataCollection_30trials_CertaintyScale\EMG\HPT_certaintyscale\segmented'];
        case 'PC202' %Anna's computer
            path_HPT_Segmented = ['D:\Students\Anna\HeartbeatCountingTask\Data\', NameStudy1, '\ECG\', NameStudy2, '\segmented'];
    end
cd(path_HPT_Segmented)
load( [path_HPT_Segmented, '\DataStrc_HBCperiod_', PARAMETER.source_file_s(5:end-4) , '.mat']);
Table_HPT = struct2table(Da);
switch Computer
    case 'PC308'
        path = ['C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Interoceptive_Certainty\Data\behavior\HPT'];
        cd(path);
        
    case 'PC202' %Anna's computer   
        path = ['D:\Students\Anna\HeartbeatCountingTask\Data\', NameStudy1, '\behavior\', NameStudy2];
        cd(path);
end


%% Load all Data
Files_beh_Subject = dir([PARAMETER.source_file_s(5:end-4), '*']);

%% Find the behavior data file matching the HBT-file
Filename = [];
indFiles = 1;

while indFiles <= size(Files_beh_Subject,1)
    if strmatch(PARAMETER.source_file_s(12:end-4), Files_beh_Subject(indFiles).name(8:9)) 
       Filename = Files_beh_Subject(indFiles).name;
       break
    else
             Filename = [];
    end
    indFiles = indFiles +1;
end

% Manual check by looking at
% PARAMETER.source_file_s(5:end-4)

Filename
%% loop & concatentate (verketten) files of one day in a table

DataOneFile  = [];
trial  = [];
load(Filename)
Tab_Trial = struct2table(trial);
Tab_Trial = Tab_Trial(2:Tab_Trial.trial(end),:);
Tab_Trial = Tab_Trial(:,1:end-7); %some not helpful variables

Subject   = Filename(1:6);
Session   = Filename(9:18);
Task      = Filename(20:22);
Tab_Trial.trial     = (1:Tab_Trial.trial(end)-1)';

if strmatch(class(Tab_Trial.counted_heartbeats), 'cell')
    Tab_Trial.counted_heartbeats = str2double(Tab_Trial.counted_heartbeats);
end

if strmatch(class(Tab_Trial.wager_choosen_post), 'char')
    Tab_Trial.wager_choosen_post =  str2num(Tab_Trial.wager_choosen_post);
end
if strmatch(class(Tab_Trial.wager_choosen_pre), 'char')
    Tab_Trial.wager_choosen_pre =  str2num(Tab_Trial.wager_choosen_pre);
end


T                   = table(repmat(Subject,size(Tab_Trial(:,1),1),1), repmat(Session,size(Tab_Trial(:,1),1),1),repmat(Task,size(Tab_Trial(:,1),1),1) ,'VariableNames', {'Subject' 'Session' 'Task'} );
Table_HPT.duration  = round(Table_HPT.duration);
Table_HPT = Table_HPT(:,2:end); %remove the automatic detected HB
DataOneFile         = [T,Table_HPT, Tab_Trial]; % finale Version für ein Subject

%% combine the new subject with the previous Data to a huge table
if addToExistingTable == 1
    switch Computer
        case 'PC308' 
            path_save = 'C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Interoceptive_Certainty\Data\combined\';
        case 'PC202' %Anna's computer
               path_save = 'D:\Students\Anna\HeartbeatCountingTask\Data\3ndDataCollection_30trials_CertaintyScale\combined\';
    end   
    
  load( [path_save, 'HPT_CertaintyScale_Concatenated.mat']);
  AllData_Combined         = [DataOneFile; AllData_Combined];
  save([path_save, 'HPT_CertaintyScale_Concatenated' ],'AllData_Combined');

end

%% save the combined behavior & ECG-data to: 
%
    path_save = 'C:\Users\kkaduk\Desktop\Kristin\Projects\Metacognition_Interoception_Human\Interoceptive_Certainty\Data\combined\';
    save([path_save, PARAMETER.source_file_s(5:end-4),'_Table_HBCperiod_behavior' ],'DataOneFile');
    
%% Finally, transform the mat-file to a txt-file   
    FileName = 'HPT_CertaintyScale_Concatenated.txt';

    writetable(AllData_Combined, [path_save, FileName], 'Delimiter', ' ')
    writetable(AllData_Combined, ['C:\Users\kkaduk\Dropbox\promotion\Projects\Metacogition_Interoception_Human\HPT_Certainty\Data\' FileName], 'Delimiter', ' ')

end
