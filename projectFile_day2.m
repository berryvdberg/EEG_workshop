global RUN
RUN.rootpath = 'C:\Users\berry\Desktop\beesearch_workshop\day1';
RUN.eeglab = 'eeglab14_1_2b';
RUN.datapath = fullfile(RUN.rootpath,'data');

cd(RUN.rootpath);
b=2;
RUN.b = b;
% hey marlon, hieronder ben ik nog wat onzeker over XYZ

a =1;
plosOne(a)

%---- set plotting params ----%
%RUN.linetype = 
%RUN.linecolor = 


%---- set preproc params ----%
RUN.preproc.remIC = 'no';
RUN.preproc.cue.epochs = {'6' '7' '8' '9'};
RUN.preproc.target.epochs =  cellfun(@num2str,num2cell(10:100),'un',0);
RUN.preproc.epochLength = [-1 4]; % choose long epochs, useful for time frequency analysis

%% Define subject information here
RUN.subjectID = {'1'};
RUN.filename = {'1_beesearch_2019-05-07_14-31-25.cnt' };

%% you only need EEGlab for preprocessing, after that we go on a fieltrip

addpath(RUN.eeglab)
eeglab

%%
figure;imagesc(EEG.data(:,1:1000))


%% Preprocessing
iSub = 1;
file = fullfile(RUN.datapath,RUN.subjectID{iSub},RUN.filename{iSub}); % make this into a loop later
EEG = pop_loadeep_v4(file);
EEG = pop_chanedit(EEG, 'lookup',fullfile(fileparts(which('eeglab.m')), ...
    '/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp'));

eventList = EEG.event; % copy event information
% tmp = num2cell(1:size({EEG.event.latency},2));
% [EEG.event.indx] = tmp{:};
% EEG.urevent = EEG.event;

% make pivot table here


%%

EEG = pop_eegfiltnew(EEG,0.1,0); % FILTERS! (please look into the difference of causal and non-causal fitlers) -> this should be done as early as possible [kernell is quite big]

EEG = pop_eegfiltnew(EEG,0,30); % low pass filter [you can do this much later on the average data as well] -> kernell is much smaller
EEG = pop_epoch( EEG, RUN.preproc.cue.epochs, RUN.preproc.epochLength);
EEG = pop_rmbase( EEG, [-200 0]);


%load(fullfile(RUN.dataPath,'preproc', [RUN.subjectID{iSub} '_ICrun.mat'])) % What I do is the following:
% 1: RUN ICA on data that is filtered at 0.5Hz or 1Hz -> ICA doesnt do
% well with slow wave drifts
% 2: copy the ICA weights to a dataset that is filtered at 0.1 or lower
% (there is brain activity in slow wave stuff [CDA and CNV])
% 3: remove those components that reflect eyeblinks [max 1/2 per subject]

% EEG.icaweights = ICrun.icaweights;
% EEG.icasphere = ICrun.icasphere;
% EEG.icawinv =ICrun.icawinv;
% EEG.icachansind =ICrun.icachansind;


% pop_topoplot(EEG,0,1:10)
% pop_eegplot(EEG,0)
% pop_eegplot(EEG,1)
% EEG = pop_subcomp(EEG,RUN.ics{iSub}, 1);
% keyboard;

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
% data = eeglab2fieldtrip(EEG,'preprocessing'); % I use fieldtrip from this moment onwards
% data.trialStruct = trialStruct;
% save(fullfile(RUN.dataPath,'preproc', [RUN.subjectID{iSub} '_data.mat']),'data')
% 


%% why regresssion coefficients are just as valid to construct ERPs as mean
% and the quickest way to get the coefs
% https://nl.mathworks.com/help/matlab/data_analysis/linear-regression.html
load accidents
x = hwydata(:,14); %Population of states
y = hwydata(:,4);
X = [ones(length(x),1) x];
b = X\y

% do the above for a single channel.... 

figure;plot(x,y)




