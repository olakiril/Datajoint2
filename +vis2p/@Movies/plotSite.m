function plotSite(obj,tracetype)

if nargin<2 || isempty(tracetype)
    tracetype = 'calcium_trace';
end

%get key
keys = fetch(obj);
if length(keys)~=1
    k = 1;
else
    k = 0;
end
while k < length(keys)
    if length(keys)==1
        k = 1;
    end
    key = keys(k);
    
    if strcmp(tracetype,'norm')
        gtraces = cell2mat(fetchn(MaskTraces(key,'masknum > 0'), 'calcium_trace')');
        rtraces = cell2mat(fetchn(MaskTraces(key,'masknum > 0'),'annulus_trace')');
        traces = gtraces- rtraces;
    else
        % get trace
        traces = fetchn(MaskTraces(key,'masknum > 0'),tracetype);
traces2 = fetchn(Traces(key,'trace_opt = 17 and masknum > 0'),'trace');
        traces = cell2mat(traces');
           traces2 = cell2mat(traces2');

    end
    
    [fps, x, y]   = fetchn( Movies(key), 'fps','alignXshifts','alignYshifts');
    try
    ball = fetch1(ExtraTraces(key),'speed_trace');
    catch
        ball = [];
    end
    traces = calcDFoF(double(traces),fps);
    scale = 1;
    numpoints = 1:0.5:ceil(size(traces,2)/2)+1;
     hold on
      plot(bsxfun(@plus,traces2*scale,numpoints(1:size(traces2,2))),'k');
    plot(bsxfun(@plus,traces*scale,numpoints(1:size(traces,2))));
    set(gca,'XTick',0:fps*100:size(traces,1),'XTickLabel',0:100:size(traces,1)/fps,...
        'YTick',2.5:2:size(traces,2)/2,'YTickLabel',4:4:size(traces,2))
   
    plot(normalize(x{1})-2)
    plot(normalize(y{1})-4,'r')
    plot(ball-6,'g')
    xlabel('time(sec)')
    ylabel('cell #')
    set(gcf,'Color',[1 1 1])
    key
    if length(keys)~= 1
        reply = input('','s');
        if isempty(reply)
            k = k+1;
        elseif k~=1 && ~isempty(reply)
            k = k - 1;
        end
        clf(1)
    end
end