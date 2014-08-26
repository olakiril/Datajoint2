function makeTuples( obj, key)

import vis2p.*

% get stim file
stim = fetch1(VisStims(key),'stim_file');

% Get stim data
trial_conditions = horzcat(stim.params.trials.conditions);
conditions = horzcat(stim.params.conditions);

% get the selected trials
selectedConditions = 1:length(conditions);

trials = cell(1,length(selectedConditions));
for icond = 1:length(selectedConditions)
    trials{icond} = find(trial_conditions == selectedConditions(icond));
end
trials = cell2mat(trials);
stimuli = [stim.params.conditions.orientation];

% loop through repeats
old_key = key;
key.total_trials = length(trials);
key.uni_ori = length(unique(stimuli));

% insert data
insert( obj, key );
makeTuples(OriTrials,old_key);

end
