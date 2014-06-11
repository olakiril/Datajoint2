keys = fetch(MaskGroup.*(StatArea.*StatsSites('exp_date>"2013-08-05"')));


for ikey = 1:length(keys);
    
    
    key = keys(ikey);
    clear traces
    % get traces and frames
    [masknums,traces(:,:,1),traces(:,:,2),traces(:,:,3)] = ...
        fetchn(MaskTraces(key,'masknum <> -1 and masknum <> -3 and masknum <> 0'),...
        'masknum','calcium_trace','red_trace','annulus_trace');
    fps = fetch1(Movies(key),'fps');
    opts = [3 17 22 25 24];
    TRACES = cell(length(opts),1);
    for iTopt = 1:length(opts)
        key.trace_opt = opts(iTopt);
        traceOpt = fetch(TracesOpt(key),'*'); %#ok<NASGU>
        options = fetch1(TracesOpt(key),'trace_computation');
        options = strread(options, '%s','delimiter',',');
        
        % filter traces
        traceTypes = size(traces,3);
        tracesP = double([traces{:}]);  % traces are in columns
        tracesP = reshape(tracesP,size(tracesP,1),[],traceTypes);
        for iopt = 1:length(options)
            [tracesP, qual]= eval([options{iopt} '(tracesP,fps,traceOpt)']);
        end
        TRACES{iTopt} = tracesP;
    end
end

%%
close all
figure
s(1) = subplot(221);
plot(bsxfun(@plus,TRACES{1},(1:size(TRACES{1},2))/2),'color',[0.2 0.2 0.2])
hold on
plot(bsxfun(@plus,TRACES{2},(1:size(TRACES{2},2))/2),'color',[1 0 0])

s(2) = subplot(222);
plot(bsxfun(@plus,TRACES{1},(1:size(TRACES{2},2))/2),'color',[0.2 0.2 0.2])
hold on
plot(bsxfun(@plus,TRACES{3},(1:size(TRACES{2},2))/2),'color',[1 0 0])

s(3) = subplot(223);
plot(bsxfun(@plus,TRACES{1},(1:size(TRACES{2},2))/2),'color',[0.2 0.2 0.2])
hold on
plot(bsxfun(@plus,TRACES{4},(1:size(TRACES{2},2))/2),'color',[1 0 0])
linkaxes(s)

s(4) = subplot(224);
plot(bsxfun(@plus,TRACES{1},(1:size(TRACES{2},2))/2),'color',[0.2 0.2 0.2])
hold on
plot(bsxfun(@plus,TRACES{5},(1:size(TRACES{2},2))/2),'color',[1 0 0])
linkaxes(s)

