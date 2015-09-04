function plotTraces(obj,key,trace_opt,colors)

import vis2p.*

if nargin<4
    if nargin<3
        trace_opt = [3 17];
    end
    colors = lines(length(trace_opt));
end

% get trace
fps = fetchn( Movies(key).*TracesGroup, 'fps');
% clf
% subplot(4,1,1:3)
hold on


for iTrace = 1:length(trace_opt)
    traces = fetchn(Traces(key,['trace_opt =' num2str(trace_opt(iTrace))]),'trace');

    traces = cell2mat(traces');
        plot((0:1:(size(traces,1)-1))/fps,bsxfun(@plus,traces*1,1:size(traces,2)),'color',colors(iTrace,:));
end

set(gca,'XTick',0:fps*100:size(traces,1),'XTickLabel',0:100:size(traces,1)/fps)

xlabel('time(sec)')
ylabel('cell #')
set(gcf,'Color',[1 1 1])

