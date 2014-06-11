function makeTuples( obj, key )

% get Traces
traces = fetchn(Traces('trace_opt = 6 and masknum > -1').*Scans(key),'trace');
fps = fetch1(Movies(key),'fps');
traces = [traces{:}];

% Load times of the traces
times = fetch1(VisStims(key),'frame_timestamps');

% stim times and types of movies shown
stimTimes = fetchn(StatsPresents(key),'movie_times');

% find trace segments for each stimulus and remove 0.5 from each
% side
tracesC = cell(1,length(stimTimes));
for iTimes = 1:length(stimTimes)
    tracesC{iTimes} = traces(times > (stimTimes{iTimes}(1)) & ...
        times < (stimTimes{iTimes}(end)),:);
end
traces = cell2mat(tracesC');
clear tracesC
clear times

% binarize
s = std(traces);
traces = bsxfun(@(t,s) t > key.thr_factor * s,traces,s);

% bin
d = max(1,round(key.binsize/1000*fps));
k = ones(d,1)/d;
traces = conv2(double(traces),k,'valid');
traces = traces(1:d:end,:);
traces = traces>0;
tracesS = traces;
for i = 1:size(traces,2)
    tracesS(:,i) = traces(randperm(size(traces,1)),i);
end

% calculate simultaneously active neurons
patterns = sum(traces,2);
patternsS = sum(tracesS,2);

% calculate probabilities
key.pp = nan(size(traces,2),1);
for iNeurons = 1:size(traces,2)
    key.pp(iNeurons) = mean(patterns == iNeurons);
end

key.ppS = nan(size(traces,2),1);
for iNeurons = 1:size(traces,2)
    key.ppS(iNeurons) = mean(patternsS == iNeurons);
end

% insert data
insert( obj, key );