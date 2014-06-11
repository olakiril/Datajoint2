function makeTuples( obj, key )

% image statistics
[ik im is ic ib mn] = fetchn(StatsImages(key,'statparam = 1'),'im_kurtosis','im_mean',...
    'im_std','im_pwz_corr','im_bispectrum','movie_num');
sim_traces = fetchn(StatsSim(key),'sim_traces');
try
    ori_traces = fetch1(StatsOri(key,'ori_param = 1'),'R');
catch
    ori_traces = 0;
end

% make sure all the data are already computed
if length(mn) ~= length(sim_traces)
      disp('SimTraces table is not computed!')
    return
elseif length(mn) ~= size(ori_traces,1)
    disp('OriTraces table is not computed!')
    return
end

% get traces of all neurons of a site and fps
T = Traces(key);
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

im_k = nan(size(ttimes,1),length(stimTimes));im_m = im_k;im_s = im_k;im_c = im_k;im_b = im_k;
str = nan(size(ttimes,1),length(stimTimes),size(sim_traces{1},1));
otr = nan(size(ttimes,1),length(stimTimes),size(ori_traces,2));
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
        im_b(iStim,iTrial) = nanmean(ib{mn == movieNums(iTrial)} (indx));
        im_c(iStim,iTrial) = nanmean(nanmean(ic{mn == movieNums(iTrial)}(indx,1:min([50 size(ic{1},2)])))); % was  im_c(iStim,iTrial) = nanmean(nansum((ic{mn == movieNums(iTrial)}(indx,1:50))));
        str(iStim,iTrial,:) = nanmean(sim_traces{mn == movieNums(iTrial)} (:,indx),2);
        otr(iStim,iTrial,:) = nanmean(ori_traces(mn == movieNums(iTrial),:,indx),3);
    end
end

% remove trials with no stimulus
indx = ~isnan(sum(im_k,2));
im_k = im_k(indx,:);
im_m = im_m(indx,:);
im_s = im_s(indx,:);
im_c = im_c(indx,:);
im_b = im_b(indx,:);
traces = permute(traces(indx,:,:),[1 3 2]);
str = str(indx,:,:);
otr = otr(indx,:,:);

% collapse segments of different repetitions
t = cell(2,1);k = t;m = t;s = t;c = t;st = t;ot = t;b = t;
for iMovie = 1:size(movieTypes,2)
    t{iMovie} = traces(:,movieTypes(:,iMovie),:);
    k{iMovie} = im_k(:,movieTypes(:,iMovie));
    m{iMovie} = im_m(:,movieTypes(:,iMovie));
    s{iMovie} = im_s(:,movieTypes(:,iMovie));
    c{iMovie} = im_c(:,movieTypes(:,iMovie));
    b{iMovie} = im_b(:,movieTypes(:,iMovie));
    st{iMovie} = str(:,movieTypes(:,iMovie),:);
    ot{iMovie} = otr(:,movieTypes(:,iMovie),:);
end
minin = @(t) cellfun(@(x) x(:,1:min(cellfun(@(x) size(x,2),t)),:),t,'Uniformoutput',0);
t = minin(t);
k = minin(k);m = minin(m);s = minin(s);c = minin(c);st = minin(st);ot = minin(ot);b = minin(b);
key.traces = permute(cat(4,t{:}),[1 2 4 3]);
key.sim_traces = permute(cat(4,st{:}),[1 2 4 3]);
key.ori_traces = permute(cat(4,ot{:}),[1 2 4 3]);
key.im_k = cat(3,k{:});
key.im_m = cat(3,m{:});
key.im_s = cat(3,s{:});
key.im_c = cat(3,c{:});
key.im_b = cat(3,b{:});

% insert data
insert( obj, key );

