function [traces, keys] = getCaTraces(obj,varargin)

params.trace_opt = 17;
params.raw = 0;
params.rawpc = 0;

params = getParams(params,varargin);

import vis2p.*

masknum = fetchn(obj,'masknum');
traces = cell(length(masknum),1);
keys = fetch(obj);

if masknum==0
    disp('oh')
    keys(1)
end
for i = 1:length(masknum);
        keys(i).masknum = masknum(i);
    key = keys(i);
    
    if params.raw
        try traces{i} = fetch1(MaskTraces(key),'calcium_trace');end
    elseif params.rawpc
        fps = fetch1(Movies(key),'fps');
        k = rmfield(key,'masknum');
        [trac, masknums] = fetchn(MaskTraces(k,'masknum>0'),'calcium_trace','masknum');
        trac = cell2mat(trac');
%         trac = double(trac);
%         k = hamming(round(fps/0.05)*2+1);    k = k/sum(k);
        trac = trac + abs(min(trac(:)))+eps;
%         trac = trac./convmirr(trac,k)-1;  %  dF/F where F is low pass
%         trac(isnan(trac)) = 0;
        [u,~] = eigs( cov(double(trac)),1 );
        trac = (trac - trac*u*u');
        traces{i} = trac(:,masknums==masknum(i));
    else
        key.trace_opt = params.trace_opt;
        traces{i} = fetch1(Traces(key),'trace');
    end
end
