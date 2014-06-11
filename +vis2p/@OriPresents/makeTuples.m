function makeTuples( obj, key)

% get stim file
stim = fetch1(VisStims(key),'stim_file');

% get the different stimulus conditions
uniStims = stim.params.constants.orientation;

% loop through unique dotSizes and stimFrames
for iOri = 1:length(uniStims)
    
    % Get stim data
    if strcmp(fetch1(VisStims(key),'exp_type'),'MultDimExperiment');
        trial_conditions = horzcat(stim.params.trials.conditions);
    else
        trial_conditions = horzcat(stim.params.trials.condition);
    end
    
    startIndx = find(strcmp(vertcat(stim.eventTypes),'showStimulus'));
    endIndx = find(strcmp(vertcat(stim.eventTypes),'endStimulus'));
    
    % get the selected trials
    trials = find(trial_conditions == iOri);
    
    % loop through repeats
    for iRepeat = 1:length(trials)
        key.repeat_num = iRepeat;
        key.ori_num = iOri;
        key.ori_times = [stim.events(trials(iRepeat)).syncedTimes(stim.events(iRepeat).types == startIndx) ...
            stim.events(trials(iRepeat)).syncedTimes(stim.events(iRepeat).types == endIndx)];
        
        % insert data
        insert( obj, key );
    end
end


