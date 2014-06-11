function makeTuples( obj, key, stim, indx )

% Get stim data
startIndx = find(strcmp(vertcat(stim.eventTypes),'showStimulus'));
dotTimes = cell(1,length(stim.params.trials));
for iTrial = 1:length(stim.params.trials)
    syncedTimes = stim.events(iTrial).syncedTimes;
    startTime = syncedTimes(startIndx == stim.events(iTrial).types);
    stimsL = length(stim.params.trials(iTrial).dotColors);
    timeWindow = stim.params.trials(1).delayTime / stimsL;
    dotTimes{iTrial} = startTime + (0:stimsL-1) * timeWindow;
end
dotTimes = cell2mat(dotTimes);
dotLocations = horzcat(stim.params.trials.dotLocations);
dotColors = horzcat(stim.params.trials.dotColors);

% loop through repeats
for iRepeat = 1:length(indx)
    key.repeat_num = iRepeat;
    key.dot_time = round(dotTimes(indx(iRepeat)));
    key.dot_loc_x = dotLocations{indx(iRepeat)}(1,1);
    key.dot_loc_y = dotLocations{indx(iRepeat)}(2,1);
    key.dot_color = dotColors{indx(iRepeat)};
    
    if key.dot_loc_x < 0 || key.dot_loc_y < 0 
        continue
    end
    % insert data
    insert( obj, key );
end





