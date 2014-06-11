function makeTuples( obj, key )

% get traces of all neurons of a site and fps
T = Traces(key);
if ~isempty(RFFit(key))
    T = T.*RFFit(key,['snr >' num2str(key.rf_snr_thr)]);
end
if ~isempty(RFStats(key))
    T = T.*RFStats(key,['onpoff_p <' num2str(key.rf_p_thr)]);
end
if ~isempty(MaskTracesQuality(key))
    T = T.*MaskTracesQuality(key,['ca_event_snr>' num2str(key.trace_qual)]);
end
tracesR = fetchn(T,'trace');

tracesR = [tracesR{:}];
fps    = fetch1( Movies(key), 'fps' );

% Load times of the traces
times = fetch1(VisStims(key),'frame_timestamps');

% equalize traces/times %%%%%% CHECK WHY
times = times(1,1:min([size(times,2) size(tracesR,1)]));
tracesR = tracesR(1:min([size(times,2) size(tracesR,1)]),:);

% stim times and types of movies shown
stimTimes = fetchn(StatsPresents(key),'movie_times');
movieTypes = fetchn(StatsPresents(key),'movie_num');
movieTypes = bsxfun(@eq,movieTypes,unique(movieTypes)');

% find trace segments for each stimulus and remove 0.5' from each
% side
traces = cell(1,length(stimTimes));
for iTimes = 1:length(stimTimes)
    traces{iTimes} = tracesR(times > (stimTimes{iTimes}(1) + 500) & ...
        times < (stimTimes{iTimes}(end) - 500 ),:);
end

if isempty(traces)
    display('Nothing to compute!')
    return
end

% remove incomplete trials
L = cell2mat(cellfun(@size,traces,'UniformOutput',0)');
indx = L(:,1) >= mean(L(:,1))*9/10;
traces = traces(indx);
L = L(indx);
movieTypes = movieTypes(indx,:);

% equalize trial length
traces = cellfun(@(x) x(1:min(L),:),traces,'UniformOutput',0);
traces = cat(3,traces{:});

% rebin to approximate binsize (downsampling by an integer factor)
d = max(1,round(key.binsize/1000*fps));
k = ones(d,1)/d;
traces = convn(traces,k,'valid');
traces = double(traces(1:d:end,:,:));

% equilize trials
mintr = min(sum(movieTypes,1));
trialIndx = nan(size(movieTypes,1),1);
for iMovie = 1:size(movieTypes,2)
    trialIndx(movieTypes(:,iMovie)) = 1:sum(movieTypes(:,iMovie));
end
indx = trialIndx <= mintr;
trialIndx = trialIndx(indx);
traces = traces(:,:,indx);

% organize in [classes cells trials]
data = nan(size(traces,1),size(traces,2),mintr,size(traces,3)/mintr);
for iTrial = 1:mintr
    data(:,:,iTrial,1: size(traces,3)/mintr) = traces(:,:,trialIndx == iTrial);
end
traces = permute(data,[2 3 1 4]);
traces = reshape(traces,size(traces,1),size(traces,2),[]);
traces = permute(traces,[1 3 2]);

% measure mean and variance of the responses
m = mean(traces,3);
v = var(traces,[],3);

% randomize traces
randtraces = traces(:);
rindx = 1:numel(traces);
for i = 1:1000
    rindx = rindx(randperm(numel(traces)));
end
randtraces = reshape(randtraces(rindx),size(traces,1),size(traces,2),size(traces,3));

% get rid of covariance structure by randomizing trials
straces = randdim(traces,3);

% randomize variance across mean
pm = repmat(m,[1 1 size(traces,3)]);
vr = traces(:) - pm(:);
vtraces = pm + reshape(vr(rindx),size(traces,1),size(traces,2),size(traces,3));

% replace variance with opposite stimulus
k = key;
if strcmp(k.movie_type,'phase')
    k.movie_type = 'natural';
elseif strcmp(k.movie_type,'natural')
    k.movie_type = 'phase';
end
traces2 = fetchTraces(StatsSites(k));
traces2 = traces2(1:size(traces,1),1:size(traces,2),1:size(traces,3));
[t2 idx2] = sort(mean(traces2,3),2);
pm2 = repmat(t2,[1 1 size(traces,3)]);
[x,~,z] = ind2sub(size(t2),1:numel(t2));
idx = sub2ind(size(t2),x,idx2(:)',z)';
mat = repmat(idx,1,size(traces2,3)) + numel(idx)*(meshgrid(1:size(traces2,3),1:size(idx,1))-1);
vr = reshape(traces2(mat(:)) - pm2(:),size(traces));

[mt idx] = sort(m,2);
pm = repmat(mt,[1 1 size(traces,3)]);
vstraces = pm + vr;
[x,~] = ind2sub(size(t2),1:numel(t2));
i(sub2ind(size(t2),x,idx(:)')) = 1:numel(t2);
mat = repmat(i',1,size(traces2,3)) + numel(i)*(meshgrid(1:size(traces2,3),1:size(i',1))-1);
vstraces = reshape(vstraces(mat(:)),size(traces));


% do it
D = nan(size(traces,1),(size(traces,2)^2 - size(traces,2))/2);
for iCell = 1:size(traces,1);
    d = nan((size(traces,2)^2 - size(traces,2))/2,key.repetitions);
    for iRep = 1:key.repetitions
        cellindx = randperm(size(traces,1));
        data = traces(cellindx(1:iCell),:,:);
        d(:,iRep) = pdist(mean(data,3)');
    end
    D(iCell,:) = mean(d,2); 
end

key.class_mean = m;
key.class_var = v;
key.performance = classloop(traces,key.repetitions);
key.sperformance = classloop(straces,key.repetitions);
key.rperformance = classloop(randtraces,key.repetitions);
key.vperformance = classloop(vtraces,key.repetitions);
key.vsperformance = classloop(vstraces,key.repetitions);
key.distances = D;

% insert data
insert( obj, key );

function J = classloop(traces,rep)
J = nan(size(traces,1),1);
for iCell = 1:size(traces,1);
    p = nan(rep,size(traces,3),size(traces,2));
    for iRep = 1:rep
        cellindx = randperm(size(traces,1));
        data = traces(cellindx(1:iCell),:,:);
        for iTrial = 1:size(traces,3)
            ind = true(size(traces,3),1);
            ind(iTrial) = false;
            r = mean(data(:,:,ind),3);
            for iClass = 1:size(traces,2)
                dist = pdist2(r',data(:,iClass,iTrial)');
                [foo indx] = min(dist); %#ok<ASGLU>
                p(iRep,iTrial,iClass) = indx == iClass;
           
            end
        end
    end
    J(iCell) = mean(p(:));
end

