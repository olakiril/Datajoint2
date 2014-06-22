function [T binsize] = getExtraTracesJake(obj,type,varargin)

% function [T binsize] = getExtraTraces(obj,varargin)
%
% gets the traces for the StatsSites experiment
% [cells time trials]
%
% MF 2012-02

params.bin = 100; % msec
params.thr = 0; % sd's
params.method = 'conv'; % method of trace resampling
params.key = 'mouse_id>0 and masknum>0';
params.function = [];
params.compute = 0; % overide obj when computing StatsSites
params.compute2 = 0;
params.minLength = 8; % min seconds of movie for concatenation
params.minTr = 2; % minimum trials for concatenation
params.collapse = 0; % collapse trials of same site recorgings
params.unimovie = 0; % separate traces for different movies

params = getParams(params,varargin);

%get key
keys = fetch(obj);

import vis2p.*

T = cell(1,length(keys));
binsize = T;
for k = 1:length(keys)
    
    key = keys(k);
    
    % get trace
    traces = fetchn(ExtraTracesJake(key),type);
    if isempty(traces);continue;end
    tracesM = cell2mat(traces');
   
    % get framerate    
    fps    = fetch1( Movies(key), 'fps' );
    
    % Load times of the trace
    times = fetch1(VisStims(key),'frame_timestamps');
    
    % apply function if supplied
    if ~isempty(params.function)
        tracesM = params.function(tracesM);
    end
    
    % equalize traces 2 times if length difference is 1
    if abs(length(tracesM) - length(times)) == 1
       ml = min([size(tracesM,1)  length(times)]);
       tracesM = tracesM(1:ml,:);
       times = times(1:ml);
    end

    % rebin to approximate binsize (downsampling by an integer factor)
    [tracesM, d] = trresize(tracesM,fps,params.bin,params.method);
    times = trresize(times',fps,params.bin,params.method);
    binsize{k} = d*1000/fps;
    
    %stim times and types
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
    end
        
    % equilize trials
    tracessize = [];
    tracessize(:,1) = cell2mat(cellfun(@(x) size(x,1),Trace,'UniformOutput',0)');
    tracessize(:,2) = cell2mat(cellfun(@(x) size(x,2),Trace,'UniformOutput',0)');
    tracessize(:,3) = cell2mat(cellfun(@(x) size(x,3),Trace,'UniformOutput',0)');
    if sum(tracessize(:)) == 0;continue;end % skip if traces are empty
    traces = cellfun(@(x) x(:,:,1:min(tracessize(:,3))),Trace,'UniformOutput',0);
    traces = permute(cell2mat(traces'),[2 1 3]);

    % threshold traces
    if params.thr
        sz = size(traces);
        traces = reshape(traces,size(traces,1),[]);
        sd = std(traces,[],2);
        traces = reshape(bsxfun(@gt,traces,sd),sz);
    end
    
    T{k} = traces;
end

% collapse same sites
if length(keys)>1 && params.collapse
    keys = arrayfun(@(x) structfun(@(x) num2str(x),x,'uniformoutput',0),keys);
    if params.compute || params.compute2;
        keys = struct2cell(arrayfun(@(x) rmfield(x,'stim_idx'),keys));
    end
    ind = [];for i = 1:size(keys,1);ind(i) = length(unique(keys(i,:)))==1;end
    if ind 
        if iscell(T{1})
            T2 = cell(length(T{1}),1);
            for imov = 1:length(T{1})
                T3 = cellfun(@(x) x{imov},T,'UniformOutput',0);
                sT = cellfun(@(x) size(x,3),T3);
                sC = cellfun(@(x) size(x,2),T3);
                traces = T{1}{imov}(:,1:min(sC),:);
                for iT = 2:length(T)
                    traces(:,:,end+1:end+sT(iT)) = T{iT}{imov}(:,1:min(sC),:);
                end
                T2{imov} = traces;
            end
            T = T2;
        else
            sT = cellfun(@(x) size(x,3),T);
            sC = cellfun(@(x) size(x,2),T);
            ind = sT>=params.minTr &  params.bin*sC/1000/length(uniMovies)>params.minLength;
            T = T(ind);sT = sT(ind);sC = sC(ind);
            traces = T{1}(:,1:min(sC),:);
            for iT = 2:length(T)
                traces(:,:,end+1:end+sT(iT)) = T{iT}(:,1:min(sC),:);
            end
            T = traces;
        end
    end
end


if length(keys) == 1
    T = cell2mat(T);
end

