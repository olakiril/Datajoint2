function makeTuples( obj, key)

import vis2p.*

% get stim file
stim = fetch1(VisStims(key),'stim_file');

% Get stim data
trial_conditions = horzcat(stim.params.trials.conditions);
conditions = horzcat(stim.params.conditions);
startIndx = find(strcmp(vertcat(stim.eventTypes),'showStimulus'));
stopIndx = find(strcmp(vertcat(stim.eventTypes),'endStimulus'));

% get the selected trials
selectedConditions = 1:length(conditions);

trials = cell(1,length(selectedConditions));
for icond = 1:length(selectedConditions) 
    trials{icond} = find(trial_conditions == selectedConditions(icond));
end
trials = cell2mat(trials);

% loop through repeats
for iRepeat = 1:length(trials)
    key.repeat_num = iRepeat;
    key.movie_times = [stim.events(trials(iRepeat)).syncedTimes(stim.events(trials(iRepeat)).types == startIndx) ...
            stim.events(trials(iRepeat)).syncedTimes(stim.events(trials(iRepeat)).types == stopIndx)];
    key.ori_in = stim.params.conditions(stim.params.trials(trials(iRepeat)).conditions).orientationIn;
    key.ori_out = stim.params.conditions(stim.params.trials(trials(iRepeat)).conditions).orientationOut;

    % insert data
    insert( obj, key );
end
