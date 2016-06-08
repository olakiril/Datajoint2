function [traces,out] = gnorm(traces,fps,traceOpts) %#ok<INUSD,*INUSD>

traces = traces(:,:,1)./traces(:,:,3);


if nargout>1
    out = [];
end