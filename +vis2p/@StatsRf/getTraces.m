function traces = getTraces(obj, key )

% check if input is already a key
if nargin<2
    key = fetch(obj);
end

% get traces of all neurons of a site and fps
T = Traces(key);

% select cells
k.mouse_id = key.mouse_id;
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

% get
tracesR = [tracesR{:}];
fps    = fetch1( Movies(key), 'fps' );

% Load times of the traces
times = fetch1(VisStims(key),'frame_timestamps');

% equalize traces/times %%%%%% CHECK WHY
times = times(1,1:min([size(times,2) size(tracesR,1)]));
tracesR = tracesR(1:min([size(times,2) size(tracesR,1)]),:);

% stim times and types of movies shown
stimTimes = fetchn(StatsPresents(key),'movie_times');
movieNums = fetchn(StatsPresents(key),'movie_num');
movieTypes = bsxfun(@eq,movieNums,unique(movieNums)');

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

if isempty(traces) || isempty(fetch(StatsImages(key)))
    display('Nothing to compute!')
    return
end

% remove incomplete trials
L = cell2mat(cellfun(@size,traces,'UniformOutput',0)');
indx = L(:,1) >= mean(L(:,1))*9/10;
traces = traces(indx);
ttimes = ttimes(indx);
stimTimes = stimTimes(indx);
movieNums = movieNums(indx);
L = L(indx);
movieTypes = movieTypes(indx,:);

% equalize
ttimes = cellfun(@(x) x(1:min(L),:),ttimes,'UniformOutput',0);
traces = cellfun(@(x) x(1:min(L),:),traces,'UniformOutput',0);

% rebin to approximate binsize (downsampling by an integer factor)
traces = cat(3,traces{:});
d = max(1,round(key.binsize/1000*fps));
key.actual_binsize = d*1000/fps;
k = ones(d,1)/d;
traces = convn(traces,k,'valid');
traces = double(traces(1:d:end,:,:));
ttimes = cat(2,ttimes{:});
ttimes = convn(ttimes,k,'valid');
ttimes = double(ttimes(1:d:end,:));

% correct for calcium onset delay 
ttimes = ttimes - key.delay;

% image statistics
[ik im is ic ib mn] = fetchn(StatsImages(key),'im_kurtosis','im_mean',...
    'im_std','im_pwz_corr','im_bispectrum','movie_num');
im_k = nan(size(ttimes,1),length(stimTimes));im_m = im_k;im_s = im_k;im_c = im_k;im_b = im_k;
for iTrial = 1:length(stimTimes)
    for iStim = 1:size(ttimes,1)
        if iStim == size(ttimes,1)
            step = stimTimes{iTrial} >= ttimes(iStim,iTrial);
        else
            step = stimTimes{iTrial} < ttimes(iStim+1,iTrial);
        end
        indx = stimTimes{iTrial} >= ttimes(iStim,iTrial) & step;
        im_k(iStim,iTrial) = nanmean(ik{mn == movieNums(iTrial)} (indx));
        im_m(iStim,iTrial) = nanmean(im{mn == movieNums(iTrial)} (indx));
        im_s(iStim,iTrial) = nanmean(is{mn == movieNums(iTrial)} (indx));
        im_c(iStim,iTrial) = nanmean(nansum((ic{mn == movieNums(iTrial)}(indx,1:50))));
        im_b(iStim,iTrial) = nanmean(ib{mn == movieNums(iTrial)} (indx));
    end
end

% remove trials with no stimulus
indx = ~isnan(sum(im_k,2));
traces = permute(traces,[1 3 2]);
traces = permute(traces(indx,:,:),[1 3 2]);

% mean across same stimulus of different repetitions and collapse segments
t = cell(2,1);k = t;m = t;s = t;c = t;b = t;
for iMovie = 1:size(movieTypes,2)
    t{iMovie} = mean(traces(:,:,movieTypes(:,iMovie)),3);
    k{iMovie} = mean(im_k(:,movieTypes(:,iMovie)),2);
    m{iMovie} = mean(im_m(:,movieTypes(:,iMovie)),2);
    s{iMovie} = mean(im_s(:,movieTypes(:,iMovie)),2);
    c{iMovie} = mean(im_c(:,movieTypes(:,iMovie)),2);
    b{iMovie} = mean(im_b(:,movieTypes(:,iMovie)),2);
end
% traces = [stim cells]
traces = cell2mat(t);
