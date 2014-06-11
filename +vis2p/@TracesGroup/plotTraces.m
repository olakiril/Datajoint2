function plotTraces(obj,key,trace_opt,colors)

if nargin<4
    if nargin<3
        trace_opt = [3 17];
    end
    colors = lines(length(trace_opt));
end

% get trace
fps = fetchn( Movies(key), 'fps');
clf
hold on


for iTrace = 1:length(trace_opt)
    traces = fetchn(Traces(key,['trace_opt =' num2str(trace_opt(iTrace))]),'trace');

    traces = cell2mat(traces');
traces = traces*2;
    plot(bsxfun(@plus,traces*2,1:size(traces,2)),'color',colors(iTrace,:));
end

set(gca,'XTick',0:fps*100:size(traces,1),'XTickLabel',0:100:size(traces,1)/fps)

xlabel('time(sec)')
ylabel('cell #')
set(gcf,'Color',[1 1 1])

