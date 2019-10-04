%% -------------------------------------------------------------------- %%
%           project file Beesearch _ template
%------------------------------------------------------------------------%
% berryv.dberg@gmail.com
%------------------------------------------------------------------------%

set(0,'DefaultFigureColormap',feval('jet'));
close all;clear all;clc;
restoredefaultpath;
set(0,'defaultuicontrolfontname','Arial');
set(0,'defaultaxesfontname','Arial');
set(0,'fixedwidthfontname','Arial');

global RUN

if ispc
    RUN.rootpath = 'C:\Users\Berry\Dropbox\projects\projects_current\beesearch_workshop\';
    RUN.analysisDir = fullfile(RUN.rootpath,'analysis');cd(RUN.analysisDir)
    RUN.dataPath = fullfile(RUN.rootpath,'data');
    RUN.matlabAppsPath = RUN.rootpath;
end

RUN.fieldtrip = fullfile(RUN.matlabAppsPath,'fieldtrip-20170918');
RUN.eeglab = fullfile(RUN.matlabAppsPath,'eeglab14_1_2b');
addpath(fullfile(RUN.matlabAppsPath,'useful_functions'));% contains many functions such as
%   detection of latency, decoding, topoplots, etc.

%---- set preproc params ----%
RUN.preproc.remIC = 'no';
RUN.preproc.cue.epochs = {'6' '7' '8' '9'};
RUN.preproc.target.epochs =  cellfun(@num2str,num2cell(10:100),'un',0);
RUN.preproc.epochLength = [-1 4]; % choose long epochs, useful for time frequency analysis

%% Define subject information here
RUN.subjectID = {'1','2','3','4','6','8','9','10','11','12','13','14','15','16'};
RUN.filename = {'1_beesearch_2019-05-07_14-31-25.cnt','2_beesearch_2019-05-09_10-48-45.cnt', ...
    '3_beesearch_2019-05-09_14-47-57.cnt','4_beesearch_2019-05-13_10-38-54.cnt', ...
    '6_beesearch_2019-05-14_14-12-14.cnt','8_beesearch_2019-05-27_10-32-15.cnt', ...
    '9_beesearch_2019-05-27_14-30-59.cnt','10_beesearch_2019-05-28_10-33-38.cnt', ...
    '11_beesearch_2019-05-28_14-36-21.cnt','12_beesearch_2019-05-29_10-57-15.cnt', ...
    '13_beesearch_2019-05-29_14-33-21.cnt','14_beesearch_2019-06-04_14-50-15.cnt', ...
    '15_beesearch_2019-06-05_10-40-36.cnt','16_beesearch_2019-06-05_14-45-16.cnt', ...
    };


%% Initialize EEGLab
addpath(RUN.eeglab)
eeglab


%% RUN ICA
for iSub = 1:length(RUN.subjectID)
    file = fullfile(RUN.dataPath,'raw',RUN.subjectID{iSub},RUN.filename{iSub}); % make this into a loop

    EEG = pop_loadeep_v4(file); % load data
    EEG = pop_chanedit(EEG, 'lookup',fullfile(fileparts(which('eeglab.m')), ...
        '/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp')); % add channel locations
    EEG2 = EEG; % this is just an extra structure that we can use to compare raw to filtered data. 
    EEG = pop_eegfiltnew(EEG,0.5,0); 
    EEG = pop_eegfiltnew(EEG,0,30); %
    EEG = pop_epoch( EEG, RUN.preproc.cue.epochs, RUN.preproc.epochLength);
    EEG = pop_rmbase( EEG, [-200 0]);
    
    % do yourself RUN ICA here after removing major artefacts  - reduce the duration of epoch to speed up
    EEG = pop_eegthresh(EEG,1,find(~strcmp({EEG.chanlocs.labels},'EOG')) ,-1*140,1250,-2.5,2.5,2,0); % this trick exludes the eyeblinks
    tmp = pop_rejepoch(EEG,EEG.reject.rejthresh,0); 
    tmp = pop_runica(tmp, 'icatype', 'runica');
    
    % copy the important stuff
    ICrun.icaweights    = tmp.icaweights;
    ICrun.icasphere     = tmp.icasphere;
    ICrun.icawinv       = tmp.icawinv;
    ICrun.icachansind   = tmp.icachansind;
    save(fullfile(RUN.dataPath,'preproc', [RUN.subjectID{iSub} '_ICrun.mat']),'ICrun')
end

%% PREPROCESSING
for iSub = 1:length(RUN.subjectID)
    file = fullfile(RUN.dataPath,'raw',RUN.subjectID{iSub},RUN.filename{iSub}); % make this into a loop
    EEG = pop_loadeep_v4(file);
    EEG = pop_chanedit(EEG, 'lookup',fullfile(fileparts(which('eeglab.m')), ...
        '/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp'));
    eventList = EEG.event; % copy event information
    tmp = num2cell(1:size({EEG.event.latency},2));
    [EEG.event.indx] = tmp{:};
    EEG.urevent = EEG.event;
    EEG = pop_eegfiltnew(EEG,0.1,0); % FILTERS! (please look into the difference of causal and non-causal fitlers) -> this should be done as early as possible [kernell is quite big]

    EEG = pop_eegfiltnew(EEG,0,30); % low pass filter [you can do this much later on the average data as well] -> kernell is much smaller
    EEG = pop_epoch( EEG, RUN.preproc.cue.epochs, RUN.preproc.epochLength);
    EEG = pop_rmbase( EEG, [-200 0]);
    
    
    load(fullfile(RUN.dataPath,'preproc', [RUN.subjectID{iSub} '_ICrun.mat'])) % What I do is the following:
    % 1: RUN ICA on data that is filtered at 0.5Hz or 1Hz -> ICA doesnt do
    % well with slow wave drifts
    % 2: copy the ICA weights to a dataset that is filtered at 0.1 or lower
    % (there is brain activity in slow wave stuff [CDA and CNV])
    % 3: remove those components that reflect eyeblinks [max 1/2 per subject]
    
    EEG.icaweights = ICrun.icaweights;
    EEG.icasphere = ICrun.icasphere;
    EEG.icawinv =ICrun.icawinv;
    EEG.icachansind =ICrun.icachansind;
    
    
    pop_topoplot(EEG,0,1:10)
    pop_eegplot(EEG,0)
    pop_eegplot(EEG,1)
    EEG = pop_subcomp(EEG,RUN.ics{iSub}, 1);
    keyboard;
    
    % the next part is key - management of trial information - I call this a trial structure
    % port codes
    % cue: 6 7 8 9
    % target: 11:88
    % response 1 2 3 4
    
    trialStruct.cueType = zeros(length(EEG.epoch),1);
    trialStruct.cueLatency = zeros(length(EEG.epoch),1);
    trialStruct.targetType = zeros(length(EEG.epoch),1);
    trialStruct.targetLatency = zeros(length(EEG.epoch),1);
    trialStruct.responseType = zeros(length(EEG.epoch),1);
    trialStruct.responseLatency = zeros(length(EEG.epoch),1);
    
    
    for i = 1:length(EEG.epoch)
        % get the urevent index
        urIndx = cell2mat(EEG.epoch(i).eventindx(cell2mat(EEG.epoch(i).eventlatency) == 0));
        
        % okay, this is the cue
        trialStruct.cueType(i) =  str2double(eventList(urIndx).type);
        trialStruct.cueLatency(i) =  eventList(urIndx).latency;
        
        % lets add the target (which happens somewhere around 1200ms post cue)
        % how do we do this?
        % 1) look forward which event occured around the expected time:
        tmpIdx = cell2mat({eventList.latency}) > eventList(urIndx).latency &  cell2mat({eventList.latency}) < eventList(urIndx).latency + 1300;
        tmpIdx = find( tmpIdx & ismember({eventList.type},RUN.preproc.target.epochs));
        
        
        % 1) find the latency of the cue
        if ~isempty(tmpIdx)
            trialStruct.targetType(i) =  str2double(eventList(tmpIdx).type);
            trialStruct.targetLatency(i) = eventList(tmpIdx).latency;
        end
        % 2) Assignment: find the latency of the response
        
        tmpIdx = cell2mat({eventList.latency}) > eventList(urIndx).latency &  cell2mat({eventList.latency}) < eventList(urIndx).latency + 1900;
        tmpIdx = find( tmpIdx & ismember({eventList.type},{'1' '2' '3' '4'}));
        if ~isempty(tmpIdx)
            trialStruct.responseType(i) =  str2double(eventList(tmpIdx).type);
            trialStruct.responseLatency(i) = eventList(tmpIdx).latency;
        end
    end
    refchannel = ismember({EEG.chanlocs.labels}, {'M1', 'M2'});
    EEG = pop_reref( EEG, refchannel, 'keepref','on');
    
    trialStruct.targetLatency - trialStruct.cueLatency; % huh! something is wrong! What is it?
    data = eeglab2fieldtrip(EEG,'preprocessing'); % I use fieldtrip from this moment onwards
    data.trialStruct = trialStruct;
    save(fullfile(RUN.dataPath,'preproc', [RUN.subjectID{iSub} '_data.mat']),'data')
end

%% Artefact detection & data export. 
for iSub = 1:length(RUN.subjectID)
    load(fullfile(RUN.dataPath,'preproc', [RUN.subjectID{iSub} '_data.mat']))
    
    threshold = 150;
    stepThresh = 50;
    timeIntv = [-0.5 1.5];
    trialStruct = data.trialStruct;
    
    % threshold
    trialStruct.artThresh = zeros(length(data.trial),1);
    for iTrl = 1:length(data.trial);
        intvIdx = data.time{iTrl}>timeIntv(1) & data.time{iTrl} <timeIntv(2);
        if(any(max(abs(data.trial{iTrl}(:,intvIdx)))> threshold));
            trialStruct.artThresh(iTrl) = 1;
        end
    end
    disp(['threshold of ' num2str(threshold) ' detected '  num2str(sum(trialStruct.artThresh)) ' artefacts'] )
    
    % step function
    %------ horizontal eye movements -----%
    trialStruct.artheog = ft_artstep(data, [-200 1200], stepThresh, 100, 50, [1 2 3]);
    disp(['stepfunc of ' num2str(50) ' detected '  num2str(sum(trialStruct.artheog)) ' artefacts'] )
    trialStruct.artefactEEG = any([...
        trialStruct.artThresh],2);
    
    trialStruct.subjectID = ones(length(trialStruct.artheog),1) *  str2double(RUN.subjectID{iSub});
    trialStruct = dataset(trialStruct);
    
    cfg = [];
    cfg.keeptrials = 'yes';
    data = ft_timelockanalysis(cfg, data);
    timeIDX =  data.time>0.7 & data.time<1.1;
    chanIDX = ismember(data.label,'Cz');
    trialStruct.CNV = squeeze(mean(mean(data.trial(:,chanIDX,timeIDX),2),3));
    
    % here I export the data to R for multilevel analysis
    export(trialStruct, 'File',fullfile(RUN.dataPath,'dataExport', [RUN.subjectID{iSub} '_trialStruct.csv']),'Delimiter',',');
    save(fullfile(RUN.dataPath,'preproc', [RUN.subjectID{iSub} '_trialStruct.mat']),'trialStruct')
end



%% FROM HERE we have the choice to go to freq space/ connectivity space/ linear models/ etc.... 
% for now, lets focus on the timelocked average data. 


%% AVERAGE
for iSub = 1:length(RUN.subjectID)
    load(fullfile(RUN.dataPath,'preproc', [RUN.subjectID{iSub} '_data.mat']))
    load(fullfile(RUN.dataPath,'preproc', [RUN.subjectID{iSub} '_trialStruct.mat']))
    
    
    booleans.control = ismember(trialStruct.cueType,6) & ~ trialStruct.artThresh ;
    booleans.neither = ismember(trialStruct.cueType,7) & ismember(trialStruct.responseType,[3 4]) & ~ trialStruct.artThresh ;
    booleans.charity = ismember(trialStruct.cueType,8) & ~ trialStruct.artThresh ;
    booleans.yourself = ismember(trialStruct.cueType,9) & ~ trialStruct.artThresh;
    booleans.targetLeft = ismember(trialStruct.targetType,50:100) & ~ trialStruct.artThresh ;
    booleans.targetRight = ismember(trialStruct.targetType,10:50) & ~ trialStruct.artThresh;
    
    
    
    addpath(RUN.fieldtrip)
    % add more conditions here
    fn = fieldnames(booleans);
    
    for iConditions = 1 : length(fn)
        cfg=[];
        cfg.trials = booleans.(fn{iConditions});
        if sum(cfg.trials) > 10
            avgData = ft_timelockanalysis (cfg, data);
            timeLock.(fn{iConditions}){iSub}=avgData;
        end
    end
end

%% GRANDAVERAGE
fn = fieldnames(timeLock);

for i = 1:length(fn)
    cfg = [];
    cfg.keepindividual = 'yes';
    gaER.(fn{i}) =  ft_timelockgrandaverage(cfg,timeLock.(fn{i}){:});
end

%% PREP LAYOUT
layout.elec = timeLock.charity{1}.elec;
layout.elec.pnt = layout.elec.chanpos;

cfg = [];
cfg.elec = layout.elec;
cfg.rotate = 90;
cfg.layout = ft_prepare_layout(cfg);
layout = ft_prepare_layout(cfg);
layout.elec.pnt = layout.elec.chanpos;
ft_layoutplot(cfg);

%% SINGLEPLOT
figure;
cfg = [];
cfg.xlim = [-0.2 1.5];
cfg.layout = layout;
cfg.channel = {'Cz' 'P1', 'P2'};
ft_singleplotER(cfg,gaER.control,gaER.neither,gaER.charity,gaER.yourself);
legend({'control','neither','charity','yourself'})

%% CONTRAST
contrastER.yourselfMcontrol = gaER.control;
contrastER.yourselfMcontrol.individual = gaER.control.individual - gaER.yourself.individual;

%% TOPOMAPS
cfg= [];
cfg.latencyInc = 0.1;
cfg.latency = 0:cfg.latencyInc:1.2;cfg.latency = [cfg.latency' cfg.latency'+cfg.latencyInc];
cfg.perspective = {'top' 'back' 'left' 'right'};
cfg.ncontours = 10;
cfg.granularity = 0.1;
cfg.clim = [-2 2];
cfg.plotchannellab = 'no';
cfg.elecsize = 10;
cfg.layout = layout.elec;
jp_topoplotFT(cfg,gaER.control);
colormap jet



