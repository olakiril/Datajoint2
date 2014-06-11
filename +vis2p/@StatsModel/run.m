function [J rJ] = run( obj, classes,classreps )

key = fetch(obj);

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
[stimTimes movieTypes] = fetchn(StatsPresents(key),'movie_times','movie_num');
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

% organize in [cells classes trials]
data = nan(size(traces,1),size(traces,2),mintr,size(traces,3)/mintr);
for iTrial = 1:mintr
    data(:,:,iTrial,1: size(traces,3)/mintr) = traces(:,:,trialIndx == iTrial);
end
traces = permute(data,[2 3 1 4]);
traces = reshape(traces,size(traces,1),size(traces,2),[]);
traces = permute(traces,[1 3 2]);

% randomize traces
randtraces = traces(:);
rindx = 1:numel(traces);
for i = 1:1000
    rindx = rindx(randperm(numel(traces)));
end
randtraces = reshape(randtraces(rindx),size(traces,1),size(traces,2),size(traces,3));

% select classes
if nargin<2
    classes = size(traces,2);
    classreps = 1;
elseif nargin<3
    classreps = 30;
end

%% do it
J = nan(size(traces,1),key.repetitions,classreps);
rJ = nan(size(traces,1),key.repetitions,classreps);

for iCell = 1:size(traces,1);
    p = nan(key.repetitions,classreps,size(traces,3),classes);
    rp = nan(key.repetitions,classreps,size(traces,3),classes);
    for iClassRep = 1:classreps
        classindx = randperm(size(traces,2));
        Data = traces(:,classindx(1:classes),:);
        Randdata = randtraces(:,classindx(1:classes),:);
        for iRep = 1:key.repetitions
            cellindx = randperm(size(traces,1));
            data = Data(cellindx(1:iCell),:,:);
            randdata = Randdata(cellindx(1:iCell),:,:);
            for iTrial = 1:size(traces,3)
                ind = true(size(traces,3),1);
                ind(iTrial) = false;
                r = mean(data(:,:,ind),3);
                rr = mean(randdata(:,:,ind),3);
                for iClass = 1:classes
                    dist = pdist2(r',data(:,iClass,iTrial)');
                    [foo indx] = min(dist); %#ok<ASGLU>
                    dist = pdist2(rr',randdata(:,iClass,iTrial)'); % randdata was data 2011/08/20
                    [foo rindx] = min(dist); %#ok<ASGLU>
                    p(iRep,iClassRep,iTrial,iClass) = indx == iClass;
                    rp(iRep,iClassRep,iTrial,iClass) = rindx == iClass;
                    
                end
            end
        end
    end
    J(iCell,:,:) = mean(mean(p,4),3);
    rJ(iCell,:,:) = mean(mean(rp,4),3);
end