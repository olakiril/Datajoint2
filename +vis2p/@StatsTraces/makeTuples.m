function makeTuples( obj, key )


% get trace
trace = fetch1(Traces(key),'trace');
fps    = fetch1( Movies(key), 'fps' );

% Load times of the trace
times = fetch1(VisStims(key),'frame_timestamps');
times = times(1:length(trace));

%stim times and types
stimTimes = fetchn(StatsPresents(key),'movie_times');
movieTypes = fetchn(StatsPresents(key),'movie_num');
movieTypes = bsxfun(@eq,movieTypes,unique(movieTypes)');

% find trace segments
traces = cell(1,length(stimTimes));
for iTimes = 1:length(stimTimes)
    traces{iTimes} = trace(times>stimTimes{iTimes}(1) & times<stimTimes{iTimes}(end));
end

% remove incomplete trials
L = cell2mat(cellfun(@size,traces,'UniformOutput',0)');
indx = L(:,1) >= mean(L(:,1))*9/10;
traces = traces(indx);
L = L(indx);
movieTypes = movieTypes(indx,:);

% equalize
traces = cellfun(@(x) x(1:min(L),:),traces,'UniformOutput',0);

% collapse segments
traces = cell2mat(traces);

% % make sure there are no edge effects by removing 0.5 sec from the start and
% % 0.5 sec from the end of the trials
% traces = traces(round(fps/2):end - round(fps/2),:);

% if size(traces,2) <= 5
%     return
% end

tuple = key;

%bin
d = max(1,round(key.binsize/1000*fps));
k = ones(d,1)/d;
trace = conv2(traces,k,'valid');

if key.undersample
    trace = trace(1:d:end,:);
end

% get the zscore
traceZ = zscore(trace);

% calculate all combinations of correlations
mtrace = cell(size(movieTypes,2),1);
cmbIn = cell(size(movieTypes,2),1);
cmbOut = cell(size(movieTypes,2),1);
for iMovie = 1:size(movieTypes,2)
    % get the index of the same movie trials
    trialIn = find(movieTypes(:,iMovie));
    trialOut = find(~movieTypes(:,iMovie));
    
    % compute permuted correlations
    cmbIn{iMovie} = combnk (trialIn,2);
    if ~isempty(trialOut)
        cmbOut{iMovie} = setxor(setxor(combnk([trialIn;trialOut],2), ...
            cmbIn{iMovie},'rows'),combnk(trialOut,2),'rows');
        randOut = randperm(size(cmbOut{iMovie},1));
        cmbOut{iMovie} = cmbOut{iMovie}(randOut(1:size(cmbIn{iMovie},1)),:);
    end
    % mean trials
    mtrace{iMovie} = mean(trace(:,movieTypes(:,iMovie)),2);
end
cmbIn = cell2mat(cmbIn);
cmbOut = cell2mat(cmbOut);

% Correlate and calculate significance
[tuple.inCorr tuple.inCorrP] = corr(reshape(traceZ(:,cmbIn(:,1)),[],1),...
    reshape(traceZ(:,cmbIn(:,2)),[],1));
if ~isempty(trialOut)
    [tuple.outCorr tuple.outCorrP] = corr(reshape(traceZ(:,cmbOut(:,1)),[],1),...
        reshape(traceZ(:,cmbOut(:,2)),[],1));
else
    tuple.outCorr = 0;
    tuple.outCorrP =0;
end

% reshape trace
mtrace = cell2mat(mtrace);

% calculate Lifetime sparseness
L = length(mtrace);
tuple.life_sparse_avg = (1 - ((sum(mtrace)/L)^2/sum(mtrace.^2/L)))/(1 - (1/L));

L = length(trace(:));
tuple.life_sparse = (1 - ((sum(trace(:))/L)^2/sum(trace(:).^2/L)))/(1 - (1/L));

% calculate autocorrelation
[auto_corr.corr, auto_corr.lag, auto_corr.confidence] = autocorr(mtrace);
tuple.auto_corr = auto_corr;

% calculate kurtosis
tuple.kurtosis_avg = kurtosis(mtrace);
tuple.kurtosis = nanmean(kurtosis(trace));
tuple.kurtosis(isinf(tuple.kurtosis)) = 0;

% calculate mean
tuple.mean = mean(trace,2);

% calculate variance
tuple.variance = var(trace,0,2);

% Proportion of frames with 0 response
% binarize
s = std(trace(:));
tuple.pzero = sum(bsxfun(@(t,s) t < key.event_thr * s,trace(:),s)) / numel(trace);

s = std(mtrace(:));
tuple.pzero_avg = sum(bsxfun(@(t,s) t < key.event_thr * s,mtrace,s)) / numel(mtrace);

% insert data
insert( obj, tuple );





