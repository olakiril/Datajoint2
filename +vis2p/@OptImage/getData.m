function  [data,times,trials] = getData(obj,key)

% getData(obj)
%
% gets the Optical data
%
% MF 2015-06

import vis2p.*
    
[name, date] = fetch1( Scans(key),'file_name','exp_date' );
if ~isempty(strfind(name,'.h5'));fend = [];else fend = '.h5';end
path = ['M:/IntrinsicImaging/' datestr(date, 'yymmdd')];
filename = getLocalPath([path '/' name fend]);
[data,Fs] = getOpticalData(filename);
times = 1/Fs:1/Fs:size(data,1)/Fs;

if nargout>2
    % get stim
    stim = fetch1(VisStims(key),'stim_file');

    % get the different stimulus conditions
    uniStims = stim.params.constants.location;

    % Get stim data
    trial_conditions = horzcat(stim.params.trials.conditions);
    allconditions = horzcat(stim.params.conditions.location);
    startIndx = find(strcmp(vertcat(stim.eventTypes),'showStimulus'));
    endIndx = find(strcmp(vertcat(stim.eventTypes),'endStimulus'));

    % get the selected trials info
    conditions = nan(size(uniStims));
    for icond = 1:size(uniStims,2)
        conditions(:,icond) = find(allconditions(1,:) == uniStims(1,icond) & ...
            allconditions(2,:) == uniStims(2,icond));
    end
    trials = struct;
    for itrial = 1:length(trial_conditions)
        trials(itrial).start =  stim.events(itrial).syncedTimes(stim.events(itrial).types == startIndx)/1000;
        trials(itrial).end =  stim.events(itrial).syncedTimes(stim.events(itrial).types == endIndx)/1000;
        trials(itrial).condIdx = find(sum(trial_conditions(itrial) == conditions));
    end
end
