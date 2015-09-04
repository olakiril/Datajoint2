function [T,oris] = getTraces(obj,varargin)

% function [T ,{oriOut} {oriIn}] = getTraces(obj,varargin)
%
% gets the traces for the CenterSurround experiment
% [oriOut oriIn cells trials]
%
% MF 2014-06

params.thr = []; % sd's
params.method = [] ; % method of trace resampling
params.key = 'masknum>0';
params.compute = 0; % overide obj when computing StatsSites
params.compute2 = 0;
params.minLength = 8; % min seconds of movie for concatenation
params.minTr = 2; % minimum trials for concatenation
params.collapse = 0; % collapse trials of same site recorgings
params.off = 0;
params.resp_time = 1000;

params = getParams(params,varargin);

import vis2p.*

%get key
if params.compute
    keys = fetch(Scans(params.key)*CenterSurParams(params.key));
    for ikey = 1:length(keys);
        keys(ikey).trace_opt = params.key.trace_opt;
        keys(ikey).center_sur_opt = params.key.center_sur_opt;
        if isfield(params.key,'masknum');
            keys(ikey).masknum = params.key.masknum;
        end
    end
    params.key = 'masknum>0';
else  keys = fetch(obj);
end

T = cell(1,length(keys));
for k = 1:length(keys); key = keys(k);
    
    % get extra params
    [rsThr, Qual, rpThr, eThr,method]= fetch1(CenterSurParams(key),...
        'rf_snr_thr','trace_qual','rf_p_thr','event_thr','samp_method');
    if isempty(params.thr);params.thr = eThr;end
    if isempty(params.method);params.method = method;end
   
    % select cells on trace reconstruction quality
    if ~isnan(Qual); 
        Tr = Traces(key,[params.key ' and quality>' num2str(Qual)]);
    else
        Tr = Traces(key,params.key);
    end
    
    % select cells on rf quality
    if ~isnan(rsThr);Tr = Tr & RFFit(key,['snr >' num2str(rsThr)]);end
    if ~isnan(rpThr);Tr=Tr & RFStats(key,['onpoff_p<' num2str(rpThr)]);end
    
    % get trace
    traces = fetchn(Tr,'trace');
    if isempty(traces);disp noTraces!;T{k} = [];oris = [];continue;end
    tracesM = cell2mat(traces');
    fps    = fetch1( Movies(key), 'fps' );
    
    % Load times of the trace
    times = fetch1(VisStims(key) & CenterSurTrials,'frame_timestamps');
    
    % equalize traces 2 times if length difference is 1
    if abs(length(tracesM) - length(times)) == 1
        ml = min([size(tracesM,1)  length(times)]);
        tracesM = tracesM(1:ml,:);
        times = times(1:ml);
    end
    
    % stim times and types
    [stimTimes,oriIn,oriOut] = fetchn(CenterSurTrials(key),'movie_times','ori_in','ori_out');
    movies = [oriIn oriOut];
    uniMovies =unique(movies,'rows');
    
    % switch to off response
    if params.off
       [st,sidx] = sort(cell2mat(stimTimes));
       s = reshape(st(2:end-1),size(st,1)-1,2);
       s(end+1,:) = [st(end) st(end)+abs(mean(diff(s')))];
       stimTimes = mat2cell(s,ones(size(s,1),1),size(s,2));
       oriIn = oriIn(sidx(:,1));
       oriOut = oriOut(sidx(:,1));       
    end
    
    % find trace segments
    traces = cell(1,length(stimTimes));
    for iTimes = 1:length(stimTimes)
        traces{iTimes} = mean(tracesM(times>stimTimes{iTimes}(1) & ...
            times<(stimTimes{iTimes}(1)+params.resp_time),:));
    end
    
    % remove incomplete trials
    tracessize = cell2mat(cellfun(@size,traces,'UniformOutput',0)');
    indx = tracessize(:,1) >= mean(tracessize(:,1))*9/10;
    traces = traces(indx);
    tracessize = tracessize(indx,1);
    
    % equalize trial length
    traces = cellfun(@(x) x(1:min(tracessize(:,1)),:),traces,'UniformOutput',0);
    traces = cat(3,traces{:});
    
    % organize classes
    Trace = cell(1,length(uniMovies));
    for iMovie = 1:length(uniMovies)
        Trace{iMovie} = traces(:,:,all(bsxfun(@eq,movies(indx,:),uniMovies(iMovie,:))'));
    end
    
    % equalize trials
    tracessize = [];
    tracessize(:,1) = cell2mat(cellfun(@(x) size(x,1),Trace,'UniformOutput',0)');
    tracessize(:,2) = cell2mat(cellfun(@(x) size(x,2),Trace,'UniformOutput',0)');
    tracessize(:,3) = cell2mat(cellfun(@(x) size(x,3),Trace,'UniformOutput',0)');
    if sum(tracessize(:)) == 0;continue;end % skip if traces are empty
    traces = cellfun(@(x) x(:,:,1:min(tracessize(:,3))),Trace,'UniformOutput',0);

    traces = reshape(cat(1,traces{:}),length(unique(oriIn)),[],size(traces{1},2),size(traces{1},3));
    oris{k,2} = unique(oriIn);
    oris{k,1} = unique(oriOut);
    
    % threshold traces
    if ~isnan(params.thr)
        sz = size(traces);
        traces = reshape(traces,size(traces,1),[]);
        sd = std(traces,[],2);
        traces = reshape(bsxfun(@gt,traces,sd),sz);
    end
    T{k} = traces;
end

% convert to array
if iscell(T) && length(T) == 1
    T = T{1};
end




