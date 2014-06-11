function data = getTraceStruct(obj, key )

% check if input is already a key
if nargin<2
    key = fetch(obj);
end

% get traces of all neurons of a site and fps
T = Traces(key);

% select cells
k.exp_date= key.exp_date;
k.scan_idx= key.scan_idx;
if ~isempty(MaskTracesQuality(k))
    T = T.*MaskTracesQuality(k,['ca_event_snr>' num2str(key.trace_qual)]);
end
k.trace_opt= key.trace_opt;
if ~isempty(RFFit(k))
    T = T.*RFFit(k,['snr >' num2str(key.rf_snr_thr)]);
end
if ~isempty(RFStats(k))
    T = T.*RFStats(k,['onpoff_p <' num2str(key.rf_p_thr)]);
end
tracesR = fetchn(T,'trace');

% get RFs 
[on off] = fetchn(RFMap('rf_opt_num = 3').*T,'on_rf','off_rf');

% get
tracesR = [tracesR{:}];
fps    = fetch1( Movies(key), 'fps' );

% Load times of the traces
times = fetch1(VisStims(key).*StatsPresents,'frame_timestamps');

% equalize traces/times %%%%%% CHECK WHY
times = times(1,1:min([size(times,2) size(tracesR,1)]));
tracesR = tracesR(1:min([size(times,2) size(tracesR,1)]),:);

% stim times and types of movies shown
stimTimes = fetchn(StatsPresents(key),'movie_times');
movieNums = fetchn(StatsPresents(key),'movie_num');

% find trace segments for each stimulus and remove 0.5 from each
% side
traces = cell(1,length(stimTimes));
ttimes = cell(1,length(stimTimes));
for iTimes = 1:length(stimTimes)
    traces{iTimes} = tracesR(times > (stimTimes{iTimes}(1)) & ...
        times < (stimTimes{iTimes}(end)),:);
    ttimes{iTimes} = times(times > (stimTimes{iTimes}(1)) & ...
        times < (stimTimes{iTimes}(end)))';
end

% remove incomplete trials
L = cell2mat(cellfun(@size,traces,'UniformOutput',0)');
indx = L(:,1) >= mean(L(:,1))*9/10;
traces = traces(indx);
ttimes = ttimes(indx);
stimTimes = stimTimes(indx);
movieNums = movieNums(indx);
L = L(indx);

% equalize
ttimes = cellfun(@(x) x(1:min(L),:),ttimes,'UniformOutput',0);
traces = cellfun(@(x) x(1:min(L),:),traces,'UniformOutput',0);

% rebin to approximate binsize (downsampling by an integer factor)
traces = cat(3,traces{:});
d = max(1,round(key.binsize/1000*fps));
k = ones(d,1)/d;
traces = convn(traces,k,'valid');
traces = double(traces(1:d:end,:,:));
ttimes = cat(2,ttimes{:});
ttimes = convn(ttimes,k,'valid');
ttimes = double(ttimes(1:d:end,:));

data.rfmap.on = on;
data.rfmap.off = off;
data.actual_binsize = d*1000/fps;
data.exp_date = key.exp_date;
data.scan_idx = key.scan_idx;
data.trace_opt = key.trace_opt;
data.movie_path = '/lab/users/philipp/stimuli/MouseMovie/Mac';
data.movie_type = key.movie_type;
data.trials = [];
for iTrial = 1:size(traces,3)
    data.trials(iTrial).times = ttimes(:,iTrial)';
    data.trials(iTrial).traces = traces(:,:,iTrial)';
    data.trials(iTrial).movie_num = movieNums(iTrial);
    data.trials(iTrial).stim_times = stimTimes{iTrial};
end

