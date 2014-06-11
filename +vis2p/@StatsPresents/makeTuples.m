function makeTuples( obj, key, stim, movieStat,trialTrigger )

% Get stim data
trial_conditions = horzcat(stim.params.trials.condition);
conditions = horzcat(stim.params.conditions);
startIndx = find(strcmp(vertcat(stim.eventTypes),'showStimulus'));

% get the selected trials
if isfield(conditions,'movieStat') && isfield(conditions,'trialTrigger')
     selectedConditions = find([conditions.movieStat] == movieStat & [conditions.trialTrigger] == trialTrigger);
elseif isfield(conditions,'movieStat')
    selectedConditions = find([conditions.movieStat] == movieStat);
else
    selectedConditions = 1:length(conditions);
end
trials = cell(1,length(selectedConditions));
for icond = 1:length(selectedConditions) 
    trials{icond} = find(trial_conditions == selectedConditions(icond));
end
trials = cell2mat(trials);

% loop through repeats
for iRepeat = 1:length(trials)
    key.repeat_num = iRepeat;
    if isfield(conditions,'movieNumber')
        key.movie_num = conditions(trial_conditions(trials(iRepeat))).movieNumber;
    else % handle cases where there was only one movie shown
        key.movie_num = stim.params.constants.movieNumber;
    end
    key.movie_times = stim.params.trials(trials(iRepeat)).presentedStimTimes *1000 + ...
            stim.events(trials(iRepeat)).syncedTimes(stim.events(trials(iRepeat)).types == startIndx);

    % insert data
    insert( obj, key );
end
