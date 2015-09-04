function makeTuples( obj, key )

% params 
resp_delay = 0;
resp_period = 750;

%
import vis2p.*

% get stim file
stim = fetch1(VisStims(key),'stim_file');

% define starting event
if stim.params.constants.stimulusTime>stim.params.constants.stimFrames*1000/60
    startEvent = 'showSubStimulus';
else
    startEvent = 'showStimulus';
end
endEvent = 'endStimulus';

% Get stim data
startIndx = find(strcmp(vertcat(stim.eventTypes),startEvent));
endIndx = find(strcmp(vertcat(stim.eventTypes),endEvent));

times = horzcat(stim.events.syncedTimes);

startTimes = times([stim.events.types]==startIndx);
endTimes = times([stim.events.types]==endIndx);

% get trace
traces = fetch1(Traces(key),'trace');

% Load times of the trace
times = fetch1(VisStims(key) & OriTrials,'frame_timestamps');

% find trace segments
OnTraces = nan(1,length(startTimes));
for iTimes = 1:length(startTimes)
    OnTraces(iTimes) = mean(traces(times>startTimes(iTimes)+ resp_delay & ...
        times<(startTimes(iTimes)+resp_delay+resp_period),:));
end

OffTraces = nan(1,length(startTimes));
for iTimes = 1:length(startTimes)
    OffTraces(iTimes) = mean(traces(times>endTimes(iTimes)+ resp_delay & ...
        times<(endTimes(iTimes)+resp_delay+resp_period),:));
end

[~,pResp] = ttest2(OnTraces,OffTraces);

key.pResp = pResp;
key.OnTraces = single(OnTraces);
key.OffTraces = single(OffTraces);

% insert data
insert( obj, key );