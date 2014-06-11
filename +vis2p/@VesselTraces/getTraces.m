function [T, binsize,Tr,uniMovies] = getTraces(obj,varargin)

% function [T binsize] = getTraces(obj,varargin)
%
% gets the traces for the StatsSites experiment
% [cells time trials]
%
% MF 2012-02

params.bin = []; % msec
params.thr = []; % sd's
params.method = 'conv'; % method of trace resampling
params.key = 'vmasknum>0';
params.compute = 0; % overide obj when computing StatsSites
params.compute2 = 0;
params.minLength = 8; % min seconds of movie for concatenation
params.minTr = 2; % minimum trials for concatenation
params.collapse = 0; % collapse trials of same site recorgings
params.unimovie = 0; % separate traces for different movies
params.half = 0; % get the second half of the movie

params = getParams(params,varargin);

%get key
keys = fetch( Movies.*obj);

T = cell(1,length(keys));
binsize = T;
for ikey = 1:length(keys); key = keys(ikey);
    
    % get Traces
   Tr = VesselTraces(key,params.key);

    
    % get trace
    traces = fetchn(Tr,'radii_trace');
    tracesM = cell2mat(traces');
    fps    = fetch1( Movies(key), 'fps' )/10;
    d = max(1,round(4*fps));
    kk = ones(d,1)/d;
    tracesd = convn(tracesM,kk,'same');
    tracesM(2:end-2,:) = tracesd(2:end-2,:);

    % Load times of the trace
    times = fetchn(VesselTraces(key,params.key),'timestamps');
    times = times{1};
    
    % equalize traces 2 times if length difference is 1
    if abs(length(tracesM) - length(times)) == 1
        ml = min([size(tracesM,1)  length(times)]);
        tracesM = tracesM(1:ml,:);
        times = times(1:ml);
    end
        
    %stim times and types
    key.movie_type = params.movie_type;
    [stimTimes, movies] = fetchn(StatsPresents(key),'movie_times','movie_num');
    uniMovies = unique(movies');
    
    % find trace segments
    traces = cell(1,length(stimTimes));
    for iTimes = 1:length(stimTimes)
        traces{iTimes} = tracesM(times>stimTimes{iTimes}(1) & ...
            times<stimTimes{iTimes}(end),:);
    end
    
    % remove incomplete trials
    tracessize = cell2mat(cellfun(@size,traces,'UniformOutput',0)');
    indx = tracessize(:,1) >= mean(tracessize(:,1))*9/10;% 2012-09-24 >=mean(tracessize(:,1)) -  2*std(tracessize(:,1));
    traces = traces(indx);
    tracessize = tracessize(indx,1);
    
    % equalize trial length
    traces = cellfun(@(x) x(1:min(tracessize(:,1)),:),traces,'UniformOutput',0);
    traces = cat(3,traces{:});
    
    % organize classes
    Trace = cell(1,length(uniMovies));
    for iMovie = 1:length(uniMovies)
        Trace{iMovie} = traces(:,:,movies(indx,1) == uniMovies(iMovie));
        if params.half
             Trace{iMovie} = Trace{iMovie}(:,round(size(Trace{iMovie},2)/2)+1:end,:);
        end
    end
    
    
    % equalize trials
    tracessize = [];
    tracessize(:,1) = cell2mat(cellfun(@(x) size(x,1),Trace,'UniformOutput',0)');
    tracessize(:,2) = cell2mat(cellfun(@(x) size(x,2),Trace,'UniformOutput',0)');
    tracessize(:,3) = cell2mat(cellfun(@(x) size(x,3),Trace,'UniformOutput',0)');
    if sum(tracessize(:)) == 0;continue;end % skip if traces are empty
    traces = cellfun(@(x) x(:,:,1:min(tracessize(:,3))),Trace,'UniformOutput',0);
    if ~params.unimovie
        traces = permute(cell2mat(traces'),[2 1 3]);
    else
         traces = cellfun(@(x) permute(x,[2 1 3]),traces,'UniformOutput',0);
    end
    
    T{ikey} = traces;
end

% collapse same sites
if length(keys)>1 && params.collapse
    keys = arrayfun(@(x) structfun(@(x) num2str(x),x,'uniformoutput',0),keys);
    if params.compute || params.compute2;
        keys = struct2cell(arrayfun(@(x) rmfield(x,'stim_idx'),keys));
    end
    ind = [];for i = 1:size(keys,1);ind(i) = length(unique(keys(i,:)))==1;end
    if ind 
        sT = cellfun(@(x) size(x,3),T);
        sC = cellfun(@(x) size(x,2),T);
        ind = sT>=params.minTr &  bin*sC/1000/length(uniMovies)>params.minLength;
        T = T(ind);sT = sT(ind);sC = sC(ind);
        traces = T{1}(:,1:min(sC),:);
        for iT = 2:length(T)
            traces(:,:,end+1:end+sT(iT)) = T{iT}(:,1:min(sC),:);
        end
        T = traces;
    end
end

% convert to array
if iscell(T) && length(T) == 1
    T = T{1};
end




