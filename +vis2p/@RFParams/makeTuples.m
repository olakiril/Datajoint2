function makeTuples( obj, key )

% get stim file
stim = fetch1(vis2p.VisStims(key),'stim_file');

% initialize
dotSize = cell(1,length(stim.params.trials));
stimFrames = cell(1,length(stim.params.trials));

% put all the trials in one big cell
for i = 1:length(stim.params.trials)
    l = length(stim.params.trials(i).dotColors);
    if  datenum(key.exp_date, 'yyyy-mm-dd')> datenum('2012-05-01','yyyy-mm-dd')
        dotSize{i} = repmat(stim.params.constants.dotSize,[1 l]);
        stimFrames{i} = repmat(stim.params.constants.stimFrames,[1 l]);
    else
        dotSize{i} = repmat(stim.params.trials(i).dotSize,[1 l]);
        stimFrames{i} = repmat(stim.params.trials(i).stimFrames,[1 l]);
    end
end

% find the unique stimuli
stims = [cell2mat(dotSize) ;cell2mat(stimFrames)]';
[uniStims,foo,ii] = unique(stims,'rows'); %#ok<ASGLU>

% loop through unique dotSizes and stimFrames
for i = 1:size(uniStims,1)
    indx = find(ii == i);
    key.dot_size = uniStims(i,1);
    key.stim_frames = uniStims(i,2);
    
    % make sure that we pass only the common key
    tuple = key;
    tuple.total_repetitions = length(indx);
    tuple.loc_total_num = length(unique(cell2mat(horzcat(stim.params.trials.dotLocations))','rows'));
    insert( obj, tuple );
    
    % fill the dependant tables
    makeTuples(vis2p.RFPresents,key,stim,indx);
end



