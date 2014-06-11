function [traces,out] = gnormAll(traces,fps,traceOpts) %#ok<INUSD,*INUSD>

traces = traces(:,:,1)./repmat(mean(traces(:,:,3),2),1,size(traces,2));


if nargout>1
    out = [];
end