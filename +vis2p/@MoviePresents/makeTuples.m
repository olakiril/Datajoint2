function makeTuples( obj, key, stim)

% Get stim data
trials = horzcat(stim.params.trials.condition);
conditions = horzcat(stim.params.conditions);
startIndx = find(strcmp(vertcat(stim.eventTypes),'showStimulus'));

% loop through repeats
for iRepeat = 1:length(trials)
    key.repeat_num = iRepeat;
    key.movie_start_time = conditions(trials(iRepeat)).clipStartTimes;
    key.movie_times = (stim.params.trials(iRepeat).presentedStimTimes - key.movie_start_time)*1000  + ...
            stim.events(iRepeat).syncedTimes(stim.events(iRepeat).types == startIndx);

    % insert data
    insert( obj, key );
end

