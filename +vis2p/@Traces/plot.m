function plot(obj,opts)

if nargin<2
    opts = [3 17 22];
end

% get traces and frames
[traces(:,:,1),traces(:,:,2),traces(:,:,3)] = ...
    fetchn(MaskTraces('masknum>0').*obj,...
    'calcium_trace','red_trace','annulus_trace');
traceTypes = size(traces,3);
traces = ([traces{:}]);  % traces are in columns
traces = reshape(traces,size(traces,1),[],traceTypes);

fps = fetch1(Movies.*obj,'fps'); %#ok<NASGU>

key = fetch(Scans.*obj);
TRACES = cell(length(opts),1);
for iTopt = 1:length(opts)
    key.trace_opt = opts(iTopt);
    tracesP = fetchn(Traces(key,'masknum>0'),'trace');
    if isempty(tracesP)
        traceOpt = fetch(TracesOpt(key),'*'); %#ok<NASGU>
        options = fetch1(TracesOpt(key),'trace_computation');
        options = strread(options, '%s','delimiter',',');
        
        % filter traces
        tracesP = traces;
        for iopt = 1:length(options)
            [tracesP, qual, traceOpt]= eval([options{iopt} '(tracesP,fps,traceOpt)']);
        end
    else
        tracesP = cell2mat(tracesP');
    end
    TRACES{iTopt} = tracesP;
end

figure
s(1) = subplot(211);
plot(bsxfun(@plus,TRACES{1},(1:size(TRACES{1},2))/2),'color',[0.2 0.2 0.2])
hold on
plot(bsxfun(@plus,TRACES{2},(1:size(TRACES{2},2))/2),'color',[1 0 0])

s(2) = subplot(212);
plot(bsxfun(@plus,TRACES{1},(1:size(TRACES{2},2))/2),'color',[0.2 0.2 0.2])
hold on
plot(bsxfun(@plus,TRACES{3},(1:size(TRACES{2},2))/2),'color',[1 0 0])

% s(3) = subplot(223);
% plot(bsxfun(@plus,TRACES{1},(1:size(TRACES{2},2))/2),'color',[0.2 0.2 0.2])
% hold on
% plot(bsxfun(@plus,TRACES{4},(1:size(TRACES{2},2))/2),'color',[1 0 0])
% linkaxes(s)
%
% s(4) = subplot(224);
% plot(bsxfun(@plus,TRACES{1},(1:size(TRACES{2},2))/2),'color',[0.2 0.2 0.2])
% hold on
% plot(bsxfun(@plus,TRACES{5},(1:size(TRACES{2},2))/2),'color',[1 0 0])
linkaxes(s)

