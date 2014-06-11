function makeTuples( obj, key )

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
uniMov = unique(movieNums);
movieTypes = bsxfun(@eq,movieNums,unique(movieNums)');

% find trace segments for each stimulus and remove 0.5 from each
% side
traces = cell(1,length(stimTimes));
ttimes = cell(1,length(stimTimes));
for iTimes = 1:length(stimTimes)
    traces{iTimes} = tracesR(times >= (stimTimes{iTimes}(1)) & ...
        times <= (stimTimes{iTimes}(end)),:);
    ttimes{iTimes} = times(times >= (stimTimes{iTimes}(1)) & ...
        times <= (stimTimes{iTimes}(end)))';
end

if isempty(traces)
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
movieTypes = movieTypes(indx,:);
L = L(indx);

% get movies
path = getLocalPath('/lab/users/philipp/stimuli/MouseMovie/Mac');
type = 'nat';if strcmp(key.movie_type,'phase');type = 'phs';end;
movies = cell(2,1);
for iMov = 1:length(uniMov)
    movR = mmreader([path '/mov' num2str(uniMov(iMov)) '_' type '.avi']); %#ok<TNMLP>
    frame = imresize(mean(read(movR, 1),3),0.1);
    movie = nan(size(frame,1),size(frame,2),max(cellfun(@length,stimTimes)));
    for iFrame = 1:size(movie,3)
        movie(:,:,iFrame) = imresize(mean(read(movR, iFrame),3),0.1);
    end
    movie = permute(movie,[3 1 2]);
    d = max(1,round(key.binsize/1000*movR.FrameRate));
    k = ones(d,1)/d;
    movie = convn(movie,k,'valid');
    movies{iMov} = movie(1:d:end,:,:); 
end
stimTimes = cellfun(@(x) x(1:d:end),stimTimes,'Uniformoutput',0);

% equalize
ttimes = cellfun(@(x) x(1:min(L),:),ttimes,'UniformOutput',0);
traces = cellfun(@(x) x(1:min(L),:),traces,'UniformOutput',0);

% rebin to approximate binsize (downsampling by an integer factor)
trace = cat(3,traces{:});
ttimes = cat(2,ttimes{:});
traces = nan(size(movies{1},1),size(trace,2),size(trace,3));
for iTrial = 1:length(stimTimes)
    traces(:,:,iTrial) = interp1(ttimes(:,iTrial),trace(:,:,iTrial),stimTimes{iTrial},'cubic');
end

% mean across same stimulus of different repetitions and collapse segments
t = cell(2,1);
for iMovie = 1:size(movieTypes,2)
    t{iMovie} = mean(traces(:,:,movieTypes(:,iMovie)),3);
end
tra = cell2mat(t);
clear t

% do it
statsRFmap = zeros(size(traces,2),size(movies{1},2),size(movies{1},3));
for iCell = 1:size(traces,2)
    for iTrial = 1:size(traces,3)
        [c p] = corr(traces(4:end,iCell,iTrial),...
            movies{uniMov == movieNums(iTrial)}(1:end-3,:));
        c(p>0.001) = 0;
        rfcorr(iTrial,:) = c;
    end
    statsRFmap(iCell,:,:) = reshape(rfcorr,size(movies{1},2),size(movies{1},3));
end


% insert data
insert( obj, key );

