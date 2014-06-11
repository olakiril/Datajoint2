function plotSite(obj,tracetype)

if length(nargin)<2
    tracetype = 'calcium_trace';
end

%get key
keys = fetch(obj);
k = 1;
while k <length(keys)
    
    key = keys(k);
    
    % get trace
    traces = fetchn(MaskTraces(key),tracetype);
    [fps x y]   = fetchn( Movies(key), 'fps','alignXshifts','alignYshifts');
    
    traces = cell2mat(traces');
    scale = 1;
    numpoints = 1:0.5:ceil(size(traces,2)/2)+1;
    plot(bsxfun(@plus,traces*scale,numpoints(1:size(traces,2))));
    set(gca,'XTick',0:fps*100:size(traces,1),'XTickLabel',0:100:size(traces,1)/fps,...
        'YTick',2.5:2:size(traces,2)/2,'YTickLabel',4:4:size(traces,2))
    hold on
    plot(x{1}-2)
     plot(y{1}-4,'r')
    xlabel('time(sec)')
    ylabel('cell #')
    set(gcf,'Color',[1 1 1])
    key
    reply = input('','s');
    if isempty(reply)
        k = k+1;
    elseif k~=1 && ~isempty(reply)
        k = k - 1;
    end
    clf(1)
end