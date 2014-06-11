function [T binsize] = getSpikes(obj,varargin)

params.binsize = 100; % msec
ds = 150;

params = getParams(params,varargin);

%get key
keys = fetch(obj);

T = cell(1,length(keys));
for k = 1:length(keys)
    
    key = keys(k);
    
    [path,name] = fetch1( Experiments*Scans(key), 'directory','file_name' );
    filename = getLocalPath([path '/' name '.Htt']);
    tt = ah_readTetData(filename);
    spikeTimes = tt.t;
    
    % get trace
    fps    = fetch1( Movies(key), 'fps' );
    
    % Load times of the trace
    times = fetch1(VisStims(key),'frame_timestamps');
    times = linspace(times(1),times(end),round(length(times)/ds))';
    tracesM = zeros(size(times));
    for it = 1:length(spikeTimes)
        indx = find(spikeTimes(it)>times,1,'last');
        tracesM(indx) = tracesM(indx) + 1;
    end
    
    if isempty(params.movies)
        Sobj = StatsPresents(key);
    else
        strmov = ['movie_num = ' num2str(params.movies(1))];
        for iMov = 2:length(params.movies)
            nstr = [' or movie_num = ' num2str(params.movies(iMov))];
            strmov(end+1:end+length(nstr)) = nstr;
        end
         Sobj = StatsPresents(key,strmov);
    end
    
    movies = fetchn(Sobj,'movie_num');
    uniMovies = unique(movies');
    
    %stim times and types
    stimTimes = fetchn(Sobj,'movie_times');
    
    % find trace segments
    traces = cell(1,length(stimTimes));
    for iTimes = 1:length(stimTimes)
        traces{iTimes} = tracesM(times>stimTimes{iTimes}(1) & times<stimTimes{iTimes}(end),:);
    end
    
    % equalize
    tracessize = cell2mat(cellfun(@size,traces,'UniformOutput',0)');
    indx = tracessize(:,1) > mean(tracessize(:,1)) -  std(tracessize(:,1));
    traces = traces(indx);
    tracessize = tracessize(indx,1);
    traces = cellfun(@(x) x(1:min(tracessize(:,1)),:),traces,'UniformOutput',0);
    Trace = cell(1,length(uniMovies));
    for iMovie = 1:length(uniMovies)
        t = traces(movies(indx,1) == uniMovies(iMovie));
        Trace{iMovie} = cat(3,t{:});
    end
    c = cellfun(@size,Trace,'UniformOutput',0);
    if  mod(length([c{:}]),2) || sum(cell2mat(c)) == 0;continue;end
    Trace = cellfun(@squeeze,Trace,'UniformOutput',0);
    tracessize = cell2mat(cellfun(@size,Trace,'UniformOutput',0)');
    traces = cellfun(@(x) x(:,1:min(tracessize(:,2))),Trace,'UniformOutput',0);
    traces = cell2mat(traces');
    
    % downsample
    d = round(params.binsize/1000/ds*fps);
    traces = cumsum(traces);
    
    traces = traces(1:d:end,:);
    traces = traces(2:end,:) - traces(1:end-1,:);
    
    T{k} = traces;
    
end

if length(keys) == 1
    T = cell2mat(T);
end

