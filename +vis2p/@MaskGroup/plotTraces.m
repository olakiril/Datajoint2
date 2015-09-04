function [TRACE, FPS] = plotTraces(obj,keys,varargin)

params.factor = 2;
params.highPass = 0.1;
params.lowPass = 10;
params.hz = 20;
params.pca = [];
try
    if  strcmp(fetch1(Scans(keys),'scan_prog'),'AOD') && isempty(params.pca);
        params.pca = 1;
    else
        params.pca = false;
    end
end

params = getParams(params,varargin);

import vis2p.*
keys = fetch(MaskGroup(keys));
% figure

for ikey = 1:length(keys)
    key = keys(ikey)%#ok<NOPRT>
    
    
    fps = fetchn( Movies(key), 'fps');
    traces = fetchn(MaskTraces(key),'calcium_trace');
    traces = double(cat(2,traces{:}));
    
    if strcmp(fetch1(Scans(key),'scan_prog'),'AOD');
        traces = traces+495424;
    end
    
    if exists(Experiments('dyes="Twitch"') & Scans(key))
        traces2 = fetchn(MaskTraces(key),'red_trace');
        traces2 = double(cat(2,traces2{:}));
        traces2 = traces2+495424;
        traces = traces2./traces;
%         traces = traces2;

    end
    
    k = hamming(round(fps/params.highPass)*2+1);
    k = k/sum(k);
    traces = traces./convmirr(traces,k)-1;  %  dF/F where F is low pass
    
    if params.pca
        [c, p] = princomp(traces);
        traces = p(:,params.pca+1:end)*c(:,params.pca+1:end)';
    end
    
    if params.lowPass
        k = hamming(round(fps/params.lowPass)*2+1);
        k = k/sum(k);
        traces = convmirr(traces,k);
    end
        
    if params.hz
        traces = trresize(traces,fps,1/params.hz*1000);
        fps = params.hz;
    end
    
    plot((0:1:(size(traces,1)-1))/fps,bsxfun(@plus,traces*params.factor,1:size(traces,2)),'color',[0.1 0.1 0.1]);
%     set(gca,'XTick',0:fps*10:size(traces,1),'XTickLabel',0:10:size(traces,1)/fps)
    
    xlabel('time(sec)')
    ylabel('cell #')
    set(gcf,'Color',[1 1 1])
    
    if nargout
        TRACE = traces;
        FPS = fps;
    end
end
