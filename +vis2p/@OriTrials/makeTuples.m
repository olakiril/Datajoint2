function makeTuples( obj, key)

import vis2p.*

% get stim file
stim = fetch1(VisStims(key),'stim_file');

% define starting event 
if stim.params.constants.stimulusTime>stim.params.constants.stimFrames*1000/60
    startEvent = 'showSubStimulus';
else
    startEvent = 'showStimulus';
end

% Get stim data
trial_conditions = horzcat(stim.params.trials.conditions);
conditions = horzcat(stim.params.conditions);
startIndx = find(strcmp(vertcat(stim.eventTypes),startEvent));

% get the selected trials
selectedConditions = 1:length(conditions);

trials = cell(1,length(selectedConditions));
for icond = 1:length(selectedConditions) 
    trials{icond} = find(trial_conditions == selectedConditions(icond));
end
trials = cell2mat(trials);
times = horzcat(stim.events.syncedTimes);
idx = ([stim.events.types]==startIndx);
times = times(idx);
stimuli = [stim.params.conditions.orientation];

% loop through repeats
for iRepeat = 1:length(trials)
    key.repeat_num = iRepeat;
    key.ori_times = times(trials(iRepeat));
    key.ori_num = stimuli(trial_conditions(trials(iRepeat)));

    % insert data
    insert( obj, key );
end
