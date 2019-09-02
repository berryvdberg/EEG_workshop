%% -------------------------------------------------------------------- %%
%           project file Beesearch
%------------------------------------------------------------------------%
% berryv.dberg@gmail.com
%------------------------------------------------------------------------%

set(0,'DefaultFigureColormap',feval('jet'));
close all;clear all;clc;
restoredefaultpath;
set(0,'defaultuicontrolfontname','Arial');
set(0,'defaultaxesfontname','Arial');
set(0,'fixedwidthfontname','Arial');

global RUN;
RUN.debug = 'no';
if isunix
    RUN.roothpath = '/home/berry/Dropbox/projects/projects_current/beesearch_workshop/';
    RUN.analysisDir = fullfile(RUN.roothpath,'analysis');cd(RUN.analysisDir)
    RUN.dataPath = fullfile(RUN.roothpath,'data');
    RUN.matlabAppsPath = RUN.roothpath;
end

if ispc
    RUN.roothpath = 'C:\Users\Berry\Dropbox\projects\projects_current\beesearch_workshop\';
    RUN.analysisDir = fullfile(RUN.roothpath,'analysis');cd(RUN.analysisDir)
    RUN.dataPath = fullfile(RUN.roothpath,'data');
    RUN.matlabAppsPath =RUN.roothpath;
end



RUN.fieldtrip = fullfile(RUN.matlabAppsPath,'fieldtrip-20170918');
RUN.eeglab = fullfile(RUN.matlabAppsPath,'eeglab14_1_2b');

addpath(fullfile(RUN.matlabAppsPath,'useful_functions'))

%---- set preproc params ----%
RUN.preproc.remIC = 'no';
RUN.preproc.cue.epochs = {'6' '7' '8' '9'};
RUN.preproc.target.epochs =  cellfun(@num2str,num2cell(10:100),'un',0);
RUN.preproc.epochLength = [-2.5 2.5];


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


%% load data and do some basic filtering
for iSub = 1:length(RUN.subjectID)
    file = fullfile(RUN.dataPath,'raw',RUN.subjectID{iSub},RUN.filename{iSub}); % make this into a loop
    
    
    EEG = pop_loadeep_v4(file);
    
    EEG = pop_chanedit(EEG, 'lookup',fullfile(fileparts(which('eeglab.m')), ...
        '/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp'));
    
    eventList = EEG.event; % copy event information
    tmp = num2cell(1:size({EEG.event.latency},2));
    [EEG.event.indx] = tmp{:};
    EEG.urevent = EEG.event;
    EEG2 = EEG;
    EEG = pop_eegfiltnew(EEG,0.5,0); % filters! lets discuss this a bit
    EEG = pop_eegfiltnew(EEG,0,30);
    EEG = pop_epoch( EEG, RUN.preproc.cue.epochs, RUN.preproc.epochLength);
    EEG = pop_rmbase( EEG, [-200 0]);
    
    % do yourself RUN ICA here after removing major artefacts  - reduce the duration of epoch to speed up
    EEG = pop_eegthresh(EEG,1,find(~strcmp({EEG.chanlocs.labels},'EOG')) ,-1*140,1250,-2.5,2.5,2,0);
    tmp = pop_rejepoch(EEG,EEG.reject.rejthresh,0);
    tmp = pop_runica(tmp, 'icatype', 'runica');
    
    % copy the important stuff
    ICrun.icaweights    =  tmp.icaweights;
    ICrun.icasphere     =   tmp.icasphere;
    ICrun.icawinv       =     tmp.icawinv;
    ICrun.icachansind   = tmp.icachansind;
    save(fullfile(RUN.dataPath,'preproc', [RUN.subjectID{iSub} '_ICrun.mat']),'ICrun')
end

%%

for iSub = 1:length(RUN.subjectID)
    file = fullfile(RUN.dataPath,'raw',RUN.subjectID{iSub},RUN.filename{iSub}); % make this into a loop
    EEG = pop_loadeep_v4(file);
    
    EEG = pop_chanedit(EEG, 'lookup',fullfile(fileparts(which('eeglab.m')), ...
        '/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp'));
    
    
    eventList = EEG.event; % copy event information
    tmp = num2cell(1:size({EEG.event.latency},2));
    [EEG.event.indx] = tmp{:};
    EEG.urevent = EEG.event;
    EEG2 = EEG;
    EEG = pop_eegfiltnew(EEG,0.1,0); % filters! lets discuss this a bit
    EEG = pop_eegfiltnew(EEG,0,30);
    EEG = pop_epoch( EEG, RUN.preproc.cue.epochs, RUN.preproc.epochLength);
    EEG = pop_rmbase( EEG, [-200 0]);
    
    
    % EEG.icaweights = ICrun.icaweights;
    % EEG.icasphere = ICrun.icasphere;
    % EEG.icawinv =ICrun.icawinv;
    % EEG.icachansind =ICrun.icachansind;
    %
    %
    % pop_topoplot(EEG,0,1:10)
    % pop_eegplot(EEG,0)
    % pop_eegplot(EEG,1)
    % EEG = pop_subcomp(EEG,[5], 1);
    
    
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
    
    trialStruct.targetLatency - trialStruct.cueLatency; % huh! something is wrong! What is it?
    data = eeglab2fieldtrip(EEG,'preprocessing'); % I use fieldtrip from this moment onwards
    data.trialStruct = trialStruct;
    save(fullfile(RUN.dataPath,'preproc', [RUN.subjectID{iSub} '_data.mat']),'data')
end

%% Artefact detection
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
    
    
    export(trialStruct, 'File',fullfile(RUN.dataPath,'dataExport', [RUN.subjectID{iSub} '_trialStruct.csv']),'Delimiter',',');
    save(fullfile(RUN.dataPath,'preproc', [RUN.subjectID{iSub} '_trialStruct.mat']),'trialStruct')
end

%% DO YOURSELF! implement a loop dat loads the data and create the average per subject
% hint: timeLock.control{iSub}


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


%% CREATE THE GRANDAVERAGE 
fn = fieldnames(timeLock);

for i = 1:length(fn)
   cfg = [];
   cfg.keepindividual = 'yes';
   gaER.(fn{i}) =  ft_timelockgrandaverage(cfg,timeLock.(fn{i}){:});
end


%%
layout.elec = timeLock.charity{1}.elec;
layout.elec.pnt = layout.elec.chanpos;

cfg = [];
cfg.elec = layout.elec;
cfg.rotate = 90;
cfg.layout = ft_prepare_layout(cfg);
layout = ft_prepare_layout(cfg);
ft_layoutplot(cfg);


%%

figure;
cfg = [];
cfg.xlim = [-0.2 1.5];
cfg.layout = layout;
cfg.channel = {'Cz' 'P1', 'P2'};
ft_singleplotER(cfg,gaER.control,gaER.neither,gaER.charity,gaER.yourself);
legend({'control','neither','charity','yourself'})

%% 



%% plotting
cfg.xlim = [1.2 1.8];
cfg.ylim = [-20 10]
cfg.channel = {'PO3' 'O1' 'PO7' 'PO5'};
figure;subplot(121);ft_singleplotER(cfg,timeLock.targetLeft,timeLock.targetRight);legend('left', 'right')
cfg.channel = {'PO6' 'O2' 'PO8' 'PO4'};
subplot(122);ft_singleplotER(cfg,timeLock.targetLeft,timeLock.targetRight);legend('left', 'right')


%%
% create difference wave
timeLock.leftMINUSRight  = timeLock.targetRight;
timeLock.leftMINUSRight.avg  = timeLock.targetLeft.avg -  timeLock.targetRight.avg;

cfg = [];
cfg.baselinewindow = [1 1.2];
cfg.demean  = 'yes';
ft_preprocessing(cfg,timeLock.leftMINUSRight)

timeLock.leftMINUSRight.individual = timeLock.leftMINUSRight.avg;
timeLock.leftMINUSRight.time = 1000* timeLock.leftMINUSRight.time;
layout.elec = timeLock.charity.elec;
layout.elec.pnt = layout.elec.chanpos;

cfg = [];
cfg.elec = layout.elec;
cfg.rotate = 90;
cfg.layout = ft_prepare_layout(cfg);
ft_layoutplot(cfg);
%%
cfg= [];
cfg.latencyInc = 100;
cfg.latency = 1200:cfg.latencyInc:1800;cfg.latency = [cfg.latency' cfg.latency'+cfg.latencyInc];
cfg.perspective = {'top' 'back' 'left' 'right'};
cfg.ncontours = 10;
cfg.granularity = 0.1;
cfg.clim = [-2 2];
cfg.plotchannellab = 'no';
cfg.elecsize = 10;
cfg.layout = layout.elec;
jp_topoplotFT(cfg,timeLock.leftMINUSRight,timeLock.leftMINUSRight,timeLock.leftMINUSRight,timeLock.leftMINUSRight);
colormap jet


%% Plot filtering

unfiltered = timeLock.control;
unfiltered.avg = zeros(size(unfiltered.avg));
unfiltered.avg(:,600:660) = 1;


cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 4;
cfg.hpfiltord = 4;
cfg.hpfiltdir = 'onepass';
filtered = ft_preprocessing(cfg,unfiltered);
cfg = [];
cfg.ylim = [-1 2];
cfg.channel = {'Oz'};
figure;subplot(121);ft_singleplotER(cfg,filtered,unfiltered);


cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 2;
cfg.hpfiltord = 4;
cfg.hpfiltdir = 'twopass';
filtered = ft_preprocessing(cfg,unfiltered);

cfg = [];
cfg.ylim = [-1 2];
cfg.channel = {'Oz'};
subplot(122);ft_singleplotER(cfg,filtered,unfiltered);

%% TimeFreq decomp
for iSub = 1:length(RUN.subjectID)
    load(fullfile(RUN.dataPath,'preproc', [RUN.subjectID{iSub} '_data.mat']))
    load(fullfile(RUN.dataPath,'preproc', [RUN.subjectID{iSub} '_trialStruct.mat']))
    
    
    booleans.control = ismember(trialStruct.cueType,6) & ~ trialStruct.artThresh ;
    booleans.neither = ismember(trialStruct.cueType,7) & ismember(trialStruct.responseType,[3 4]) & ~ trialStruct.artThresh ;
    booleans.charity = ismember(trialStruct.cueType,8) & ~ trialStruct.artThresh ;
    booleans.yourself = ismember(trialStruct.cueType,9) & ~ trialStruct.artThresh;
    booleans.targetLeft = ismember(trialStruct.targetType,50:100) & ~ trialStruct.artThresh ;
    booleans.targetRight = ismember(trialStruct.targetType,10:50) & ~ trialStruct.artThresh;
    
    cfg  = [];
    cfg.keeptrials = 'yes';
    cfg.output     = 'pow';
    cfg.method     = 'mtmconvol';
    cfg.taper = 'hanning';
    cfg.foi =  1:40;
    cfg.toi        = -0.5:0.2:1.5;
    cfg.t_ftimwin =   7./cfg.foi;
    cfg.keeptapers = 'no';
    cfg.pad        = 'maxperlen';
    TFR.all = ft_freqanalysis(cfg, data);
    TFR.all.powspctrm = log10(TFR.all.powspctrm);
    
    
    fn = fieldnames(booleans);
    for iConditions = 1 : length(fn)
        cfg = [];
        cfg.trials = booleans.(fn{iConditions})';
        TFR.(fn{iConditions}){iSub} = ft_freqdescriptives(cfg,TFR.all);
    end
    
    
end


%%
TFR = rmfield(TFR,'all');
%%
fn = fieldnames(TFR);



for i = 1:length(fn)
   cfg = [];
   cfg.keepindividual = 'yes';
   
   gaTFR.(fn{i}) =  ft_freqgrandaverage(cfg,TFR.(fn{i}){:});
   
   cfg = [];
   cfg.baseline = [-0.5 -0.2];
   gaTFR.(fn{i}) =  ft_freqbaseline(cfg,gaTFR.(fn{i}));

end






%%
gaTFR.yourselfMINUScontrol = gaTFR.yourself;
gaTFR.yourselfMINUScontrol.powspctrm = gaTFR.yourself.powspctrm - gaTFR.control.powspctrm;
gaTFR.yourselfMINUSneither= gaTFR.yourself;
gaTFR.yourselfMINUSneither.powspctrm = gaTFR.yourself.powspctrm - gaTFR.neither.powspctrm;
gaTFR.yourselfMINUScharity= gaTFR.yourself;
gaTFR.yourselfMINUScharity.powspctrm = gaTFR.yourself.powspctrm - gaTFR.charity.powspctrm;
gaTFR.charityMINUSneither= gaTFR.yourself;
gaTFR.charityMINUSneither.powspctrm = gaTFR.charity.powspctrm - gaTFR.neither.powspctrm;




cfg = [];
cfg.channel = {'Cz'};
cfg.zlim = [-0.1 0.1];
cfg.layout = layout;
figure;subplot(1,3,1);ft_singleplotTFR(cfg,gaTFR.yourselfMINUScharity)
subplot(1,3,2);ft_singleplotTFR(cfg,gaTFR.yourselfMINUSneither)
subplot(1,3,3);ft_singleplotTFR(cfg,gaTFR.charityMINUSneither)

